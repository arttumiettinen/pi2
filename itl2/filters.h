#pragma once

#include <algorithm>

#include "image.h"
#include "neighbourhood.h"
#include "math/mathutils.h"
#include "pointprocess.h"
#include "fft.h"
#include "utilities.h"
#include "fastmaxminfilters.h"
#include "median.h"

namespace itl2
{
	

	/**
	Filter input image and place result to output image.
	@param pixel_t Pixel data type in input image.
	@param out_t Pixel data type in output image.
	@param processNeighbourhood Processing function that produces output value given image of a neighbourhood and corresponding mask.
	@param img Input image.
	@param nbType Type of neighbourhood to use.
	@param nbRadius Radius of the neighbourhood.
	@param out Output image.
	*/
	template<typename pixel_t, typename out_t, typename NumberUtils<pixel_t>::FloatType processNeighbourhood(const Image<pixel_t>& nb, const Image<pixel_t>& mask)>
	void filter(const Image<pixel_t>& img, Image<out_t>& out, Vec3c nbRadius, NeighbourhoodType nbType, BoundaryCondition bc)
	{
		out.mustNotBe(img);
		out.ensureSize(img);

		// Zero radius in those dimensions that are not in use
		for (size_t n = img.dimensionality(); n < nbRadius.size(); n++)
			nbRadius[n] = 0;

		Image<pixel_t> mask;
		createNeighbourhoodMask(nbType, nbRadius, mask);

		size_t totalProcessed = 0;
		#pragma omp parallel if(!omp_in_parallel() && img.pixelCount() > PARALLELIZATION_THRESHOLD)
		{
			
			Image<pixel_t> nb(mask.dimensions());

			if (img.dimensionality() >= 3)
			{
				#pragma omp for
				for (coord_t z = 0; z < img.depth(); z++)
				{
					for (coord_t y = 0; y < img.height(); y++)
					{
						for (coord_t x = 0; x < img.width(); x++)
						{
							getNeighbourhood(img, Vec3c(x, y, z), nbRadius, nb, bc);

							out(x, y, z) = pixelRound<out_t>(processNeighbourhood(nb, mask));
						}
					}

					showThreadProgress(totalProcessed, img.depth());
				}
			}
			else
			{
				#pragma omp for
				for (coord_t y = 0; y < img.height(); y++)
				{
					for (coord_t x = 0; x < img.width(); x++)
					{
						getNeighbourhood(img, Vec3c(x, y, 0), nbRadius, nb, bc);

						out(x, y, 0) = pixelRound<out_t>(processNeighbourhood(nb, mask));
					}
				}
			}
		}
	}

	/**
	Filter input image and place result to output image.
	@param pixel_t Pixel data type in input image.
	@param out_t Pixel data type in output image.
	@param param_t Type of parameter.
	@param processNeighbourhood Processing function that produces output value given image of a neighbourhood, corresponding mask, and parameter value.
	@param img Input image.
	@param nbType Type of neighbourhood to use.
	@param nbRadius Radius of the neighbourhood.
	@param out Output image.
	@param parameter Parameter that is given directly to processNeighbourhood function.
	*/
	template<typename pixel_t, typename out_t, typename param_t, typename NumberUtils<pixel_t>::FloatType processNeighbourhood(const Image<pixel_t>& nb, const Image<pixel_t>& mask, param_t parameter)>
	void filter(const Image<pixel_t>& img, Image<out_t>& out, Vec3c nbRadius, const param_t parameter, NeighbourhoodType nbType, BoundaryCondition bc)
	{
		out.mustNotBe(img);
		out.ensureSize(img);

		// Zero radius in those dimensions that are not in use
		for (size_t n = img.dimensionality(); n < nbRadius.size(); n++)
			nbRadius[n] = 0;

		Image<pixel_t> mask;
		createNeighbourhoodMask(nbType, nbRadius, mask);

		size_t totalProcessed = 0;
		#pragma omp parallel if(!omp_in_parallel() && img.pixelCount() > PARALLELIZATION_THRESHOLD)
		{

			Image<pixel_t> nb(mask.dimensions());

			if (img.dimensionality() >= 3)
			{
				#pragma omp for
				for (coord_t z = 0; z < img.depth(); z++)
				{
					for (coord_t y = 0; y < img.height(); y++)
					{
						for (coord_t x = 0; x < img.width(); x++)
						{
							getNeighbourhood(img, Vec3c(x, y, z), nbRadius, nb, bc);
							out(x, y, z) = pixelRound<out_t>(processNeighbourhood(nb, mask, parameter));
						}
					}

					showThreadProgress(totalProcessed, img.depth());
				}
			}
			else
			{
				#pragma omp for
				for (coord_t y = 0; y < img.height(); y++)
				{
					for (coord_t x = 0; x < img.width(); x++)
					{
						getNeighbourhood(img, Vec3c(x, y, 0), nbRadius, nb, bc);
						out(x, y, 0) = pixelRound<out_t>(processNeighbourhood(nb, mask, parameter));
					}
				}
			}
		}
	}

	namespace internals
	{
		/**
		Separable filtering helper (for one image+parameter).
		*/
		template<typename pixel_t, typename param_t, typename NumberUtils<pixel_t>::FloatType processNeighbourhood(const Image<pixel_t>& nb, const Image<pixel_t>& mask, param_t param)>
		void sepFilterOneDimension(Image<pixel_t>& img, coord_t r, size_t dim, param_t param, BoundaryCondition bc, bool showProgressInfo = true)
		{
			coord_t N = 2 * r + 1;
			Image<pixel_t> mask(N);
			setValue(mask, (pixel_t)1);

			size_t counter = 0;

			if (dim == 0)
			{
				#pragma omp parallel if(!omp_in_parallel() && img.pixelCount() > PARALLELIZATION_THRESHOLD)
				{
					Image<pixel_t> buffer(N);
					#pragma omp for
					for (coord_t z = 0; z < img.depth(); z++)
					{
						for (coord_t y = 0; y < img.height(); y++)
						{
							// Init for this line
							pixel_t edgeVal = bc == BoundaryCondition::Zero ? pixel_t() : img(0, y, z);
							setValue(buffer, edgeVal);
							for (coord_t x = 0; x < std::min(r + 1, img.width()); x++)
							{
								buffer(r + x) = img(x, y, z);
							}

							edgeVal = bc == BoundaryCondition::Zero ? pixel_t() : img(img.width() - 1, y, z);
							for (coord_t x = std::min(r + 1, img.width()); x < N - r; x++)
							{
								buffer(r + x) = edgeVal;
							}

							// Process
							for (coord_t x = 0; x < img.width() - r - 1; x++)
							{
								img(x, y, z) = pixelRound<pixel_t>(processNeighbourhood(buffer, mask, param));

								for (coord_t i = 0; i < N - 1; i++)
									buffer(i) = buffer(i + 1);
								buffer(N - 1) = img(x + r + 1, y, z);
							}

							// Process end of line
							for (coord_t x = std::max((coord_t)0, img.width() - r - 1); x < img.width(); x++)
							{
								img(x, y, z) = pixelRound<pixel_t>(processNeighbourhood(buffer, mask, param));

								for (coord_t i = 0; i < N - 1; i++)
									buffer(i) = buffer(i + 1);
								buffer(N - 1) = edgeVal;
							}

						}

						showThreadProgress(counter, img.depth(), showProgressInfo);
					}
				}
			}
			else if (dim == 1)
			{
				#pragma omp parallel if(!omp_in_parallel() && img.pixelCount() > PARALLELIZATION_THRESHOLD)
				{
					Image<pixel_t> buffer(N);
					#pragma omp for
					for (coord_t z = 0; z < img.depth(); z++)
					{
						for (coord_t x = 0; x < img.width(); x++)
						{
							// Init for this line
							pixel_t edgeVal = bc == BoundaryCondition::Zero ? pixel_t() : img(x, 0, z);
							setValue(buffer, edgeVal);
							for (coord_t y = 0; y < std::min(r + 1, img.height()); y++)
							{
								buffer(r + y) = img(x, y, z);
							}

							edgeVal = bc == BoundaryCondition::Zero ? pixel_t() : img(x, img.height() - 1, z);
							for (coord_t y = std::min(r + 1, img.height()); y < N - r; y++)
							{
								buffer(r + y) = edgeVal;
							}

							// Process
							for (coord_t y = 0; y < img.height() - r - 1; y++)
							{
								img(x, y, z) = pixelRound<pixel_t>(processNeighbourhood(buffer, mask, param));

								for (coord_t i = 0; i < N - 1; i++)
									buffer(i) = buffer(i + 1);
								buffer(N - 1) = img(x, y + r + 1, z);
							}

							// Process end of line
							for (coord_t y = std::max((coord_t)0, img.height() - r - 1); y < img.height(); y++)
							{
								img(x, y, z) = pixelRound<pixel_t>(processNeighbourhood(buffer, mask, param));

								for (coord_t i = 0; i < N - 1; i++)
									buffer(i) = buffer(i + 1);
								buffer(N - 1) = edgeVal;
							}

						}

						showThreadProgress(counter, img.depth(), showProgressInfo);
					}
				}
			}
			else if (dim == 2)
			{
				#pragma omp parallel if(!omp_in_parallel() && img.pixelCount() > PARALLELIZATION_THRESHOLD)
				{
					Image<pixel_t> buffer(N);
					#pragma omp for
					for (coord_t y = 0; y < img.height(); y++)
					{
						for (coord_t x = 0; x < img.width(); x++)
						{
							// Init for this line
							pixel_t edgeVal = bc == BoundaryCondition::Zero ? pixel_t() : img(x, y, 0);
							setValue(buffer, edgeVal);
							for (coord_t z = 0; z < std::min(r + 1, img.depth()); z++)
							{
								buffer(r + z) = img(x, y, z);
							}

							edgeVal = bc == BoundaryCondition::Zero ? pixel_t() : img(x, y, img.depth() - 1);
							for (coord_t z = std::min(r + 1, img.depth()); z < N - r; z++)
							{
								buffer(r + z) = edgeVal;
							}

							// Process
							for (coord_t z = 0; z < img.depth() - r - 1; z++)
							{
								img(x, y, z) = pixelRound<pixel_t>(processNeighbourhood(buffer, mask, param));

								for (coord_t i = 0; i < N - 1; i++)
									buffer(i) = buffer(i + 1);
								buffer(N - 1) = img(x, y, z + r + 1);
							}

							// Process end of line
							for (coord_t z = std::max((coord_t)0, img.depth() - r - 1); z < img.depth(); z++)
							{
								img(x, y, z) = pixelRound<pixel_t>(processNeighbourhood(buffer, mask, param));

								for (coord_t i = 0; i < N - 1; i++)
									buffer(i) = buffer(i + 1);
								buffer(N - 1) = edgeVal;
							}

						}

						showThreadProgress(counter, img.height(), showProgressInfo);
					}
				}
			}
			else
			{
				throw ITLException("Invalid dimension.");
			}
		}

		template<typename pixel_t, typename NumberUtils<pixel_t>::FloatType processNeighbourhood(const Image<pixel_t>& nb, const Image<pixel_t>& mask)> typename NumberUtils<pixel_t>::FloatType paramRemover(const Image<pixel_t>& nb, const Image<pixel_t>& mask, int param)
		{
			return processNeighbourhood(nb, mask);
		}

		template<typename pixel_t, typename NumberUtils<pixel_t>::FloatType processNeighbourhood(const Image<pixel_t>& nb, const Image<pixel_t>& mask)>
		void sepFilterOneDimension(Image<pixel_t>& img, coord_t r, size_t dim, BoundaryCondition bc, bool showProgressInfo = true)
		{
			internals::sepFilterOneDimension<pixel_t, int, internals::paramRemover<pixel_t, processNeighbourhood> >(img, r, dim, 0, bc, showProgressInfo);
		}


		/**
		Separable filtering (for image and parameter image).
		NOTE: This is mostly repetition of other versions. Consider techniques to have only single function.
		@return Count of pixels whose value changed.
		*/
		template<typename pixel1_t, typename pixel2_t, typename NumberUtils<pixel1_t>::FloatType processNeighbourhood(const Image<pixel1_t>& nb1, const Image<pixel2_t>& nb2)>
		size_t sepFilterOneDimension2Images(Image<pixel1_t>& img1, const Image<pixel2_t>& img2, coord_t r, size_t dim, BoundaryCondition bc)
		{
			img1.checkSize(img2);
			coord_t N = 2 * r + 1;

			size_t counter = 0;
			size_t changed = 0;

			if (dim == 0)
			{
				#pragma omp parallel if(!omp_in_parallel() && img1.pixelCount() > PARALLELIZATION_THRESHOLD)
				{
					Image<pixel1_t> buffer1(N);
					Image<pixel2_t> buffer2(N);
					#pragma omp for reduction(+:changed)
					for (coord_t z = 0; z < img1.depth(); z++)
					{
						for (coord_t y = 0; y < img1.height(); y++)
						{
							// Init for this line
							pixel1_t edgeVal1 = bc == BoundaryCondition::Zero ? pixel1_t() : img1(0, y, z);
							pixel2_t edgeVal2 = bc == BoundaryCondition::Zero ? pixel2_t() : img2(0, y, z);
							setValue(buffer1, edgeVal1);
							setValue(buffer2, edgeVal2);
							for (coord_t x = 0; x < std::min(r + 1, img1.width()); x++)
							{
								buffer1(r + x) = img1(x, y, z);
								buffer2(r + x) = img2(x, y, z);
							}

							edgeVal1 = bc == BoundaryCondition::Zero ? pixel1_t() : img1(img1.width() - 1, y, z);
							edgeVal2 = bc == BoundaryCondition::Zero ? pixel2_t() : img2(img2.width() - 1, y, z);
							for (coord_t x = std::min(r + 1, img1.width()); x < N - r; x++)
							{
								buffer1(r + x) = edgeVal1;
								buffer2(r + x) = edgeVal2;
							}

							// Process
							for (coord_t x = 0; x < img1.width() - r - 1; x++)
							{
								pixel1_t np = pixelRound<pixel1_t>(processNeighbourhood(buffer1, buffer2));
								if (np != img1(x, y, z))
									changed++;
								img1(x, y, z) = np;

								for (coord_t i = 0; i < N - 1; i++)
								{
									buffer1(i) = buffer1(i + 1);
									buffer2(i) = buffer2(i + 1);
								}
								buffer1(N - 1) = img1(x + r + 1, y, z);
								buffer2(N - 1) = img2(x + r + 1, y, z);
							}

							// Process end of line
							for (coord_t x = std::max((coord_t)0, img1.width() - r - 1); x < img1.width(); x++)
							{
								pixel1_t np = pixelRound<pixel1_t>(processNeighbourhood(buffer1, buffer2));
								if (np != img1(x, y, z))
									changed++;
								img1(x, y, z) = np;

								for (coord_t i = 0; i < N - 1; i++)
								{
									buffer1(i) = buffer1(i + 1);
									buffer2(i) = buffer2(i + 1);
								}
								buffer1(N - 1) = edgeVal1;
								buffer2(N - 1) = edgeVal2;
							}

						}

						showThreadProgress(counter, img1.depth());
					}
				}
			}
			else if (dim == 1)
			{
				#pragma omp parallel if(!omp_in_parallel() && img1.pixelCount() > PARALLELIZATION_THRESHOLD)
				{
					Image<pixel1_t> buffer1(N);
					Image<pixel2_t> buffer2(N);
					#pragma omp for reduction(+:changed)
					for (coord_t z = 0; z < img1.depth(); z++)
					{
						for (coord_t x = 0; x < img1.width(); x++)
						{
							// Init for this line
							pixel1_t edgeVal1 = bc == BoundaryCondition::Zero ? pixel1_t() : img1(x, 0, z);
							pixel2_t edgeVal2 = bc == BoundaryCondition::Zero ? pixel2_t() : img2(x, 0, z);
							setValue(buffer1, edgeVal1);
							setValue(buffer2, edgeVal2);
							for (coord_t y = 0; y < std::min(r + 1, img1.height()); y++)
							{
								buffer1(r + y) = img1(x, y, z);
								buffer2(r + y) = img2(x, y, z);
							}

							edgeVal1 = bc == BoundaryCondition::Zero ? pixel1_t() : img1(x, img1.height() - 1, z);
							edgeVal2 = bc == BoundaryCondition::Zero ? pixel2_t() : img2(x, img2.height() - 1, z);
							for (coord_t y = std::min(r + 1, img1.height()); y < N - r; y++)
							{
								buffer1(r + y) = edgeVal1;
								buffer2(r + y) = edgeVal2;
							}

							// Process
							for (coord_t y = 0; y < img1.height() - r - 1; y++)
							{
								pixel1_t np = pixelRound<pixel1_t>(processNeighbourhood(buffer1, buffer2));
								if (np != img1(x, y, z))
									changed++;
								img1(x, y, z) = np;

								for (coord_t i = 0; i < N - 1; i++)
								{
									buffer1(i) = buffer1(i + 1);
									buffer2(i) = buffer2(i + 1);
								}
								buffer1(N - 1) = img1(x, y + r + 1, z);
								buffer2(N - 1) = img2(x, y + r + 1, z);
							}

							// Process end of line
							for (coord_t y = img1.height() - r - 1; y < img1.height(); y++)
							{
								pixel1_t np = pixelRound<pixel1_t>(processNeighbourhood(buffer1, buffer2));
								if (np != img1(x, y, z))
									changed++;
								img1(x, y, z) = np;

								for (coord_t i = 0; i < N - 1; i++)
								{
									buffer1(i) = buffer1(i + 1);
									buffer2(i) = buffer2(i + 1);
								}
								buffer1(N - 1) = edgeVal1;
								buffer2(N - 1) = edgeVal2;
							}

						}

						showThreadProgress(counter, img1.depth());
					}
				}
			}
			else if (dim == 2)
			{
				#pragma omp parallel if(!omp_in_parallel() && img1.pixelCount() > PARALLELIZATION_THRESHOLD)
				{
					Image<pixel1_t> buffer1(N);
					Image<pixel2_t> buffer2(N);
					#pragma omp for reduction(+:changed)
					for (coord_t y = 0; y < img1.height(); y++)
					{
						for (coord_t x = 0; x < img1.width(); x++)
						{
							// Init for this line
							pixel1_t edgeVal1 = bc == BoundaryCondition::Zero ? pixel1_t() : img1(x, y, 0);
							pixel2_t edgeVal2 = bc == BoundaryCondition::Zero ? pixel2_t() : img2(x, y, 0);
							setValue(buffer1, edgeVal1);
							setValue(buffer2, edgeVal2);
							for (coord_t z = 0; z < std::min(r + 1, img1.depth()); z++)
							{
								buffer1(r + z) = img1(x, y, z);
								buffer2(r + z) = img2(x, y, z);
							}

							edgeVal1 = bc == BoundaryCondition::Zero ? pixel1_t() : img1(x, y, img1.depth() - 1);
							edgeVal2 = bc == BoundaryCondition::Zero ? pixel2_t() : img2(x, y, img2.depth() - 1);
							for (coord_t z = std::min(r + 1, img1.depth()); z < N - r; z++)
							{
								buffer1(r + z) = edgeVal1;
								buffer2(r + z) = edgeVal2;
							}

							// Process
							for (coord_t z = 0; z < img1.depth() - r - 1; z++)
							{
								pixel1_t np = pixelRound<pixel1_t>(processNeighbourhood(buffer1, buffer2));
								if (np != img1(x, y, z))
									changed++;
								img1(x, y, z) = np;

								for (coord_t i = 0; i < N - 1; i++)
								{
									buffer1(i) = buffer1(i + 1);
									buffer2(i) = buffer2(i + 1);
								}
								buffer1(N - 1) = img1(x, y, z + r + 1);
								buffer2(N - 1) = img2(x, y, z + r + 1);
							}

							// Process end of line
							for (coord_t z = img1.depth() - r - 1; z < img1.depth(); z++)
							{
								pixel1_t np = pixelRound<pixel1_t>(processNeighbourhood(buffer1, buffer2));
								if (np != img1(x, y, z))
									changed++;
								img1(x, y, z) = np;

								for (coord_t i = 0; i < N - 1; i++)
								{
									buffer1(i) = buffer1(i + 1);
									buffer2(i) = buffer2(i + 1);
								}
								buffer1(N - 1) = edgeVal1;
								buffer2(N - 1) = edgeVal2;
							}

						}

						showThreadProgress(counter, img1.height());
					}
				}
			}
			else
			{
				throw ITLException("Invalid dimension.");
			}

			return changed;
		}
	}

	/**
	Filter each dimension of the input image using the same processing function (separable filtering).
	Neighbourhood type is always Rectangular.
	@param pixel_t Pixel data type in the image. If the processing function produces real numbers, make sure that this data type can store them with adequate accuracy.
	@param processNeighbourhood Processing function that produces output value image of a neighbourhood and corresponding mask.
	@param img Image to filter. The image is filtered in-place.
	@param nbRadius Radius of the neighbourhood.
	*/
	template<typename pixel_t, typename NumberUtils<pixel_t>::FloatType processNeighbourhood(const Image<pixel_t>& nb, const Image<pixel_t>& mask)>
	void sepFilter(Image<pixel_t>& img, const Vec3c& nbRadius, BoundaryCondition bc, bool showProgressInfo = true)
	{
		for (size_t n = 0; n < std::max<size_t>(1, img.dimensionality()); n++)
		{
			internals::sepFilterOneDimension<pixel_t, processNeighbourhood>(img, nbRadius[n], n, bc, showProgressInfo);
		}
	}

	/**
	Filter each dimension of the input image using the same processing function (separable filtering).
	Neighbourhood type is always Rectangular.
	@param pixel_t Pixel data type in the image. If the processing function produces real numbers, make sure that this data type can store them with adequate accuracy.
	@param processNeighbourhood Processing function that produces output value image of a neighbourhood and corresponding mask.
	@param img Image to filter. The image is filtered in-place.
	@param nbRadius Radius of the neighbourhood.
	@param param Parameter for each dimension.
	*/
	template<typename pixel_t, typename param_t, typename NumberUtils<pixel_t>::FloatType processNeighbourhood(const Image<pixel_t>& nb, const Image<pixel_t>& mask, param_t param)>
	void sepFilter(Image<pixel_t>& img, const Vec3c& nbRadius, const Vec3<param_t>& params, BoundaryCondition bc, bool showProgressInfo = true)
	{
		// NOTE: (see also sepgauss)
		// If taking derivative with Nearest boundary condition, we must filter up to dimension where the derivative is being taken, otherwise we just return (2D) filtered version of the original.
		//	but if the derivative dimension > image dimensionality, the result will be zero.
		// If takind derivative with Zero boundary condition, the situation is the same except the result won't be zero.
		// If blurring 2D image with Nearest boundary condition, the third dimension can be processed or not.
		// If blurring 2D image with Zero boundary condition, then the third dimension MUST NOT be processed or the output becomes darker than the original.

		for (size_t n = 0; n < std::max<size_t>(1, img.dimensionality()); n++)
		{
			internals::sepFilterOneDimension<pixel_t, param_t, processNeighbourhood>(img, nbRadius[n], n, params[n], bc, showProgressInfo);
		}
	}


	/**
	Filter each dimension of the input image using the same processing function (separable filtering).
	Neighbourhood type is always Rectangular.
	@param pixel_t Pixel data type in the image. If the processing function produces real numbers, make sure that this data type can store them with adequate accuracy.
	@param processNeighbourhood Processing function that produces output value image of a neighbourhood and corresponding mask.
	@param img Image to filter. The image is filtered in-place.
	@param nbRadius Radius of the neighbourhood.
	@return Count of pixels whose value changed.
	*/
	template<typename pixel1_t, typename pixel2_t, typename NumberUtils<pixel1_t>::FloatType processNeighbourhood(const Image<pixel1_t>& nb1, const Image<pixel2_t>& nb2)>
	size_t sepFilterImageImage(Image<pixel1_t>& img1, const Image<pixel2_t>& img2, const Vec3c& nbRadius, BoundaryCondition bc)
	{
		img1.checkSize(img2);
		size_t changed = 0;
		for (size_t n = 0; n < std::max<size_t>(1, img1.dimensionality()); n++)
		{
			changed += internals::sepFilterOneDimension2Images<pixel1_t, pixel2_t, processNeighbourhood>(img1, img2, nbRadius[n], n, bc);
		}
		return changed;
	}

	namespace internals
	{
		// TODO: Filtering can be made faster for rectangular neighbourhoods where mask lookup is not required.
		// In practice the speed improvement is quite minor when compared to more complicated processes.

		template<typename pixel_t> typename NumberUtils<pixel_t>::FloatType meanOp(const Image<pixel_t>& nb, const Image<pixel_t>& mask)
		{
			typename NumberUtils<pixel_t>::FloatType sum = 0;
			typename NumberUtils<pixel_t>::RealFloatType count = 0;
			for (coord_t n = 0; n < nb.pixelCount(); n++)
			{
				typename NumberUtils<pixel_t>::FloatType val = (typename NumberUtils<pixel_t>::FloatType)nb(n);
				typename NumberUtils<pixel_t>::RealFloatType w = (typename NumberUtils<pixel_t>::RealFloatType)mask(n);

				sum += val * w;
				count += w;
			}
			return sum / count;
		}

		template<typename pixel_t> typename NumberUtils<pixel_t>::FloatType varianceOp(const Image<pixel_t>& nb, const Image<pixel_t>& mask)
		{
			/*
			// This algorithm works, but is not as stable as Welford's algorithm below.
			typename NumberUtils<pixel_t>::RealFloatType count = 0;
			typename NumberUtils<pixel_t>::FloatType total = 0;
			typename NumberUtils<pixel_t>::FloatType total2 = 0;

			// Use one value in nb as an estimate of mean to improve numerical stability.
			// See https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
			typename NumberUtils<pixel_t>::FloatType K =  nb(nb.pixelCount() / 2);

			for (coord_t n = 0; n < nb.pixelCount(); n++)
			{
				typename NumberUtils<pixel_t>::FloatType val = (typename NumberUtils<pixel_t>::FloatType)nb(n) - K;
				typename NumberUtils<pixel_t>::RealFloatType w = (typename NumberUtils<pixel_t>::RealFloatType)mask(n);
				count += w; // Note: w is either 0 or 1
				val *= w;
				total += val;
				total2 += val * val;
			}

			typename NumberUtils<pixel_t>::FloatType var = (typename NumberUtils<pixel_t>::FloatType)((total2 - (total * total / count)) / (count - 1));

			return var;
			*/

			// Welford's algorithm
			typename NumberUtils<pixel_t>::RealFloatType count = 0;
			typename NumberUtils<pixel_t>::FloatType mean = 0;
			typename NumberUtils<pixel_t>::FloatType M2 = 0;

			for (coord_t n = 0; n < nb.pixelCount(); n++)
			{
				typename NumberUtils<pixel_t>::FloatType x = (typename NumberUtils<pixel_t>::FloatType)nb(n);
				typename NumberUtils<pixel_t>::RealFloatType w = (typename NumberUtils<pixel_t>::RealFloatType)mask(n);
				if (w > 0)
				{
					count += 1;
					typename NumberUtils<pixel_t>::FloatType delta = x - mean;
					mean += delta / count;
					typename NumberUtils<pixel_t>::FloatType delta2 = x - mean;
					M2 += delta * delta2;
				}
			}

			typename NumberUtils<pixel_t>::FloatType var = M2 / (count - 1);

			return var;
		}

		/*
		Specific std dev operation avoids rounding to pixel data type before square root.
		*/
		template<typename pixel_t> typename NumberUtils<pixel_t>::FloatType stddevOp(const Image<pixel_t>& nb, const Image<pixel_t>& mask)
		{
			return std::sqrt(varianceOp(nb, mask));
		}

		template<typename pixel_t> typename NumberUtils<pixel_t>::FloatType medianOp(const Image<pixel_t>& nb, const Image<pixel_t>& mask)
		{
			std::vector<pixel_t> values;
			values.reserve(nb.pixelCount());

			for (coord_t n = 0; n < nb.pixelCount(); n++)
			{
				if (mask(n) != 0)
					values.push_back(nb(n));
			}

			return calcMedian(values);
		}

		inline float32_t nanMedianOp(const Image<float32_t>& nb, const Image<float32_t>& mask)
		{
			std::vector<float32_t> values;
			values.reserve(nb.pixelCount());

			for (coord_t n = 0; n < nb.pixelCount(); n++)
			{
				if (mask(n) != 0 && !NumberUtils<float32_t>::isnan(nb(n)))
					values.push_back(nb(n));
			}

			return calcMedian(values);
		}

		template<typename pixel_t> typename NumberUtils<pixel_t>::FloatType maskedMedianOp(const Image<pixel_t>& nb, const Image<pixel_t>& mask, pixel_t badValue)
		{
			std::vector<pixel_t> values;
			values.reserve(nb.pixelCount());

			for (coord_t n = 0; n < nb.pixelCount(); n++)
			{
				if (mask(n) != 0 && nb(n) != badValue)
					values.push_back(nb(n));
			}

			return calcMedian(values);
		}

		template<typename pixel_t> typename NumberUtils<pixel_t>::FloatType minOp(const Image<pixel_t>& nb, const Image<pixel_t>& mask)
		{
			typename NumberUtils<pixel_t>::FloatType res = std::numeric_limits<typename NumberUtils<pixel_t>::FloatType>::max();
			for (coord_t n = 0; n < nb.pixelCount(); n++)
			{
				typename NumberUtils<pixel_t>::FloatType val = (typename NumberUtils<pixel_t>::FloatType)nb(n);
				if (mask(n) != 0 && val < res)
					res = val;

			}
			return res;
		}

		template<typename pixel_t> typename NumberUtils<pixel_t>::FloatType maxOp(const Image<pixel_t>& nb, const Image<pixel_t>& mask)
		{
			typename NumberUtils<pixel_t>::FloatType res = std::numeric_limits<typename NumberUtils<pixel_t>::FloatType>::lowest();
			for (coord_t n = 0; n < nb.pixelCount(); n++)
			{
				typename NumberUtils<pixel_t>::FloatType val = (typename NumberUtils<pixel_t>::FloatType)nb(n);
				if (mask(n) != 0 && val > res)
					res = val;

			}
			return res;
		}


		template<typename seed_t, typename mask_t> typename NumberUtils<seed_t>::FloatType morphoRecOp(const Image<seed_t>& seeds, const Image<mask_t>& mask)
		{
			// If mask is zero or seed is nonzero, don't do anything
			Vec3c center((mask.dimensions() - Vec3c(1, 1, 1)) / 2);
			mask_t maskVal = mask(center);
			seed_t seedVal = seeds(center);
			if (maskVal == 0 || seedVal != 0)
				return (typename NumberUtils<seed_t>::FloatType)seeds(center);

			// Otherwise, take max of neighbourhood
			typename NumberUtils<seed_t>::FloatType res = (typename NumberUtils<seed_t>::FloatType)seeds(center);
			for (coord_t n = 0; n < seeds.pixelCount(); n++)
			{
				typename NumberUtils<seed_t>::FloatType val = (typename NumberUtils<seed_t>::FloatType)seeds(n);
				if (mask(n) != 0 && val > res)
					res = val;

			}
			return res;
		}

		/**
		Variance weighted mean filter.
		*/
		template<typename pixel_t> typename NumberUtils<pixel_t>::FloatType vaweOp(const Image<pixel_t>& nb, const Image<pixel_t>& mask, double noiseStdDev)
		{
			typename NumberUtils<pixel_t>::FloatType sum = 0.0;
			typename NumberUtils<pixel_t>::FloatType sum_sqr = 0.0;
			typename NumberUtils<pixel_t>::RealFloatType count = 0;
			for (coord_t n = 0; n < nb.pixelCount(); n++)
			{
				if (mask(n) != 0)
				{
					typename NumberUtils<pixel_t>::FloatType pix = (typename NumberUtils<pixel_t>::FloatType)nb(n);
					sum += pix;
					sum_sqr += pix * pix;
					count++;
				}
			}

			typename NumberUtils<pixel_t>::FloatType mean = sum / count;
			typename NumberUtils<pixel_t>::FloatType lvar = (sum_sqr - count * mean * mean) / (count - 1);
			typename NumberUtils<pixel_t>::FloatType gvar = (typename NumberUtils<pixel_t>::RealFloatType)noiseStdDev * (typename NumberUtils<pixel_t>::RealFloatType)noiseStdDev;

			typename NumberUtils<pixel_t>::FloatType multip;
			if (gvar <= lvar)
				multip = gvar / lvar;
			else
				multip = 1;

			typename NumberUtils<pixel_t>::FloatType pixel = (typename NumberUtils<pixel_t>::FloatType)nb(nb.pixelCount() / 2);
			return pixel - (multip * (pixel - mean));
		}


		/**
		Bilateral filter
		Mask is not used.
		*/
		template<typename pixel_t> typename NumberUtils<pixel_t>::FloatType bilateralOp(const Image<pixel_t>& nb, const Image<pixel_t>& mask, Vec2d spatialandrangesigma)
		{
			typename NumberUtils<pixel_t>::RealFloatType sigmas = (typename NumberUtils<pixel_t>::RealFloatType)spatialandrangesigma[0];
			typename NumberUtils<pixel_t>::RealFloatType sigmat = (typename NumberUtils<pixel_t>::RealFloatType)spatialandrangesigma[1];


			typename NumberUtils<pixel_t>::FloatType sum = 0.0;
			typename NumberUtils<pixel_t>::RealFloatType wsum = 0.0;

			Vec3c center((nb.dimensions() - Vec3c(1, 1, 1)) / 2);
			typename NumberUtils<pixel_t>::FloatType centerVal = (typename NumberUtils<pixel_t>::FloatType)nb(center);

			for (coord_t z = 0; z < nb.depth(); z++)
			{
				for (coord_t y = 0; y < nb.height(); y++)
				{
					for (coord_t x = 0; x < nb.width(); x++)
					{
						typename NumberUtils<pixel_t>::FloatType pix = (typename NumberUtils<pixel_t>::FloatType)nb(x, y, z);
						
						typename NumberUtils<pixel_t>::RealFloatType r2 = (Vec3c(x, y, z) - center).normSquared<typename NumberUtils<pixel_t>::RealFloatType>();
						//double r2 = (x - cx) * (x - cx) + (y - cy) * (y - cy) + (z - cz) * (z - cz);
						typename NumberUtils<pixel_t>::FloatType c = pix - centerVal;

						typename NumberUtils<pixel_t>::RealFloatType w = (typename NumberUtils<pixel_t>::RealFloatType)(
							(1 / (sigmas * sqrt(2 * PI)) * ::exp(-(r2) / (2 * sigmas * sigmas))) *	    // Spatial
							(1 / (sigmat * sqrt(2 * PI)) * ::exp(-(c * c) / (2 * sigmat * sigmat)))     // Range
							);
						// This approximation seems to create decent quality image but it is not really faster.
						//double w = (1 / (2 * PI)) * normPdfApprox(sqrt(r2), 0, sigmas)    // Spatial
						//	* normPdfApprox(c, 0, sigmat);								// Range


						sum += pix * w;
						wsum += w;
					}
				}
			}

			return sum / wsum;
		}

		/*
		Convolution for 1D case (for use in separable filtering).
		Mask is not used, so this works only with rectangular neighbourhood.
		Kernel size must match neighbourhood size.
		*/
		template<typename pixel_t> typename NumberUtils<pixel_t>::FloatType convolution1DOp(const Image<pixel_t>& nb, const Image<pixel_t>& mask, const Image<float32_t>* kernel)
		{
			typename NumberUtils<pixel_t>::FloatType sum = 0;
			coord_t N = kernel->pixelCount() - 1;
			for(coord_t n = 0; n <= N; n++)
			{
				sum += (typename NumberUtils<pixel_t>::FloatType)nb(n) * (typename NumberUtils<pixel_t>::RealFloatType)(*kernel)(N - n);
			}
			return sum;
		}

		/*
		Convolution.
		Mask is not used, so this works only with rectangular neighbourhood.
		Kernel size must match neighbourhood size.
		*/
		template<typename pixel_t> typename NumberUtils<pixel_t>::FloatType convolution3DOp(const Image<pixel_t>& nb, const Image<pixel_t>& mask, const Image<float32_t>* kernel)
		{
			typename NumberUtils<pixel_t>::FloatType sum = 0;
			Vec3c N = kernel->dimensions() - Vec3c(1, 1, 1);
			for (coord_t z = 0; z < nb.depth(); z++)
			{
				for (coord_t y = 0; y < nb.height(); y++)
				{
					for (coord_t x = 0; x < nb.width(); x++)
					{
						typename NumberUtils<pixel_t>::FloatType pix = (typename NumberUtils<pixel_t>::FloatType)nb(x, y, z);
						typename NumberUtils<pixel_t>::RealFloatType weight = (typename NumberUtils<pixel_t>::RealFloatType)(*kernel)(N.x - x, N.y - y, N.z - z);
						sum += pix * weight;
					}
				}
			}
			return sum;
		}
	}

	/**
	Helpers for Gaussian filtering.
	*/
	namespace internals
	{
		/**
		* Create 1-dimensional Gaussian kernel or its derivative.
		* @param sigma Standard deviation for each dimension.
		* @param nbRadius Radius of the kernel will be placed here.
		* @param kernel The kernel will be placed in this image.
		* @param derivative Set to true to return derivative of Gaussian.
		*/
		inline void gaussianKernel1D(double sigma, coord_t& nbRadius, Image<float32_t>& kernel, size_t derivativeOrder)
		{
			if (derivativeOrder > 2)
				throw ITLException("Invalid derivative order.");

			if (sigma > 0)
			{

				nbRadius = itl2::ceil(3.0 * sigma);

				kernel.ensureSize(2 * nbRadius + 1);

				float32_t kernelSum = 0.0;
				for (coord_t x = 0; x < kernel.width(); x++)
				{
					coord_t X = x - nbRadius;

					float32_t value = ::exp(-(X * X) / (float32_t)(2 * sigma * sigma));
					kernelSum += value;

					if (derivativeOrder == 1)
						value *= (float32_t)(-X / (sigma * sigma));
					else if (derivativeOrder == 2)
						value *= (float32_t)((X * X) / (sigma * sigma * sigma * sigma) - 1 / (sigma * sigma));

					kernel(x) = value;
				}

				// Normalize kernel
				for (coord_t n = 0; n < kernel.pixelCount(); n++)
				{
					kernel(n) /= kernelSum;
				}
			}
			else
			{
				// Do-nothing kernel
				kernel.ensureSize(1);
				kernel(0) = 1;
			}
		}

		/**
		* Create 3-dimensional Gaussian kernel or Gaussian derivative kernel.
		* @param sigma Standard deviation for each dimension.
		* @param nbRadius Radius of the kernel will be placed here.
		* @param kernel The kernel will be placed in this image.
		* @param derivativeDimension Dimension where derivative should be taken.
		* @param derivativeOrder Order of derivative of Gaussian. Set to zero value to disable derivative.
		*/
		inline void gaussianKernel3D(const Vec3d& sigma, Vec3c& nbRadius, Image<float32_t>& kernel, coord_t derivativeDimension1, coord_t derivativeDimension2, size_t dimensionality)
		{
			if (derivativeDimension1 > 2 || derivativeDimension2 > 2 || (derivativeDimension1 < 0 && derivativeDimension2 >= 0))
				throw ITLException("Invalid derivative dimension.");

			nbRadius = ceil(3.0 * sigma);

			// Zero radius in those dimensions that are not in use
			for (size_t n = dimensionality; n < nbRadius.size(); n++)
				nbRadius[n] = 0;

			Vec3c kernelSize = 2 * nbRadius + Vec3c(1, 1, 1);
			kernel.ensureSize(kernelSize);

			float32_t kernelSum = 0.0;
			for (coord_t z = 0; z < kernel.depth(); z++)
			{
				for (coord_t y = 0; y < kernel.height(); y++)
				{
					for (coord_t x = 0; x < kernel.width(); x++)
					{
						Vec3c X(x, y, z);
						X -= nbRadius;

						float32_t sum = 0.0;

						// Add contribution from each dimension.
						for (size_t i = 0; i < dimensionality; i++)
						{
							sum += (X[i] * X[i]) / (float32_t)(2 * sigma[i] * sigma[i]);
						}
						float32_t value = ::exp(-sum);
						kernelSum += value;

						if (derivativeDimension1 >= 0)
						{
							if (derivativeDimension2 == derivativeDimension1)
							{
								// dI^2 / dx_i^2 derivative
								value *= (float32_t)((X[derivativeDimension1] * X[derivativeDimension1]) / (sigma[derivativeDimension1] * sigma[derivativeDimension1] * sigma[derivativeDimension1] * sigma[derivativeDimension1]) - 1 / (sigma[derivativeDimension1] * sigma[derivativeDimension1]));
							}
							else
							{
								// dI/dx_i or dI^2 / dx_i dx_j derivative
								value *= (float32_t)(-X[derivativeDimension1] / (sigma[derivativeDimension1] * sigma[derivativeDimension1]));
								if(derivativeDimension2 >= 0)
									value *= (float32_t)(-X[derivativeDimension2] / (sigma[derivativeDimension2] * sigma[derivativeDimension2]));
							}
						}	
						
						kernel(x, y, z) = value;
					}
				}
			}

			// Normalize kernel
			for (coord_t n = 0; n < kernel.pixelCount(); n++)
			{
				kernel(n) /= kernelSum;
			}
		}

		inline void checkDerivativeDimension(coord_t dim)
		{
			if (dim < -1 || dim > 2)
				throw ITLException("Invalid derivative dimension.");
		}

		/**
		Non-separable Gaussian convolution.
		*/
		template<typename pixel_t> void gauss(const Image<pixel_t>& in, Image<pixel_t>& out, const Vec3d& sigma, coord_t derivativeDimension1, coord_t derivativeDimension2, BoundaryCondition bc)
		{
			checkDerivativeDimension(derivativeDimension1);
			checkDerivativeDimension(derivativeDimension2);

			Vec3c nbRadius;
			Image<float32_t> kernel;

			gaussianKernel3D(sigma, nbRadius, kernel, derivativeDimension1, derivativeDimension2, in.dimensionality());

			// Filter
			filter<pixel_t, pixel_t, const Image<float32_t>*, internals::convolution3DOp<pixel_t> >(in, out, nbRadius, &kernel, NeighbourhoodType::Rectangular, bc);
		}

		/**
		Separable Gaussian filtering in-place.
		*/
		template<typename pixel_t> void sepgauss(Image<pixel_t>& img, const Vec3d& sigma, coord_t derivativeDimension1, coord_t derivativeDimension2, BoundaryCondition bc, bool showProgressInfo = true)
		{
			checkDerivativeDimension(derivativeDimension1);
			checkDerivativeDimension(derivativeDimension2);

			if (derivativeDimension1 >= (coord_t)img.dimensionality() || derivativeDimension2 >= (coord_t)img.dimensionality())
			{
				// All derivatives in dimension more than image dimensionality are zero.
				setValue(img, 0);
				return;
			}

			Vec3<const Image<float32_t>* > kernels;
			Vec3c nbRadius;
			Image<float32_t> kernelImages[3];

			// Generate kernel for each direction
			for (coord_t n = 0; n < 3; n++)
			{
				coord_t derOrder = 0;
				if(n == derivativeDimension1 || n == derivativeDimension2)
					derOrder = derivativeDimension1 != derivativeDimension2 ? 1 : 2;

				gaussianKernel1D(sigma[n], nbRadius[n], kernelImages[n], derOrder);
				kernels[n] = &kernelImages[n];
			}

			sepFilter<pixel_t, const Image<float32_t>*, internals::convolution1DOp<pixel_t> >(img, nbRadius, kernels, bc, showProgressInfo);
		}

		/**
		Separable Gaussian filtering.
		Use only if data type has good enough accuracy.
		*/
		template<typename input_t, typename output_t> void sepgauss(const Image<input_t>& in, Image<output_t>& out, const Vec3d& sigma, coord_t derivativeDimension1, coord_t derivativeDimension2, BoundaryCondition bc, bool showProgressInfo = true)
		{
			setValue(out, in);
			sepgauss(out, sigma, derivativeDimension1, derivativeDimension2, bc, showProgressInfo);
		}
	}

	/**
	Gaussian filtering.
	@param in Input image (not modified).
	@param out Output image.
	@param sigma Standard deviation of the Gaussian kernel.
	@param derivativeDimension Set to negative value to calculate Gaussian filtering, and positive value less than image dimensionality to calculate image derivative in that direction using Gaussian convolution.
	@param allowOpt Allow separable filtering for 16-bit images and FFT filtering for floating point images.
	@param bc Boundary condition.
	*/
	template<typename pixel_t>
	typename std::enable_if<std::is_integral<pixel_t>::value && (sizeof(pixel_t) < 2)>::type
	gaussFilter(const Image<pixel_t>& in, Image<pixel_t>& out, const Vec3d& sigma, bool allowOpt = true, BoundaryCondition bc = BoundaryCondition::Nearest)
	{
		// Perform generic non-separable filtering
		internals::gauss(in, out, sigma, -1, -1, bc);
	}

	/**
	Gaussian filtering.
	@param in Input image (not modified).
	@param out Output image.
	@param sigma Standard deviation of the Gaussian kernel.
	@param derivativeDimension Set to negative value to calculate Gaussian filtering, and positive value less than image dimensionality to calculate image derivative in that direction using Gaussian convolution.
	@param allowOpt Allow separable optimization for 16-bit images and FFT optimization for floating point images.
	@param bc Boundary condition.
	*/
	template<typename pixel_t> void gaussFilter(const Image<pixel_t>& in, Image<pixel_t>& out, double sigma, bool allowOpt = true, BoundaryCondition bc = BoundaryCondition::Nearest)
	{
		gaussFilter(in, out, Vec3d(sigma, sigma, sigma), allowOpt, bc);
	}

	template<typename pixel_t>
	typename std::enable_if<std::is_integral<pixel_t>::value && (sizeof(pixel_t) >= 2)>::type
	gaussFilter(const Image<pixel_t>& in, Image<pixel_t>& out, const Vec3d& sigma, bool allowOpt, BoundaryCondition bc)
	{
		// Perform separable filtering if allowOpt is true
		// Perform non-separable filtering if allowOpt is false
		if (allowOpt)
		{
			internals::sepgauss(in, out, sigma, -1, -1, bc);
		}
		else
		{
			internals::gauss(in, out, sigma, -1, -1, bc);
		}
	}

	inline void gaussFilter(const Image<float32_t>& in, Image<float32_t>& out, const Vec3d& sigma, bool allowOpt, BoundaryCondition bc)
	{
		// Perform FFT filtering if allowOpt is true
		// Perform separable filtering if allowOpt is false
		if (allowOpt && bc == BoundaryCondition::Zero)
		{
			setValue(out, in);
			gaussFilterFFT(out, sigma);
		}
		else
		{
			internals::sepgauss(in, out, sigma, -1, -1, bc);
		}
	}

	/**
	Gaussian filtering in-place using separable algorithm.
	@param out Image to filter.
	@param sigma Standard deviation of the Gaussian kernel.
	@param derivativeDimension Set to negative value to calculate Gaussian filtering, and positive value less than image dimensionality to calculate image derivative in that direction using Gaussian convolution.
	@param bc Boundary condition.
	*/
	inline void gaussFilter(Image<uint16_t>& img, const Vec3d& sigma, BoundaryCondition bc)
	{
		internals::sepgauss(img, sigma, -1, -1, bc);
	}

	/**
	Gaussian filtering in-place using separable algorithm.
	@param out Image to filter.
	@param sigma Standard deviation of the Gaussian kernel.
	@param bc Boundary condition.
	*/
	inline void gaussFilter(Image<uint16_t>& img, double sigma, BoundaryCondition bc)
	{
		gaussFilter(img, Vec3d(sigma, sigma, sigma), bc);
	}

	/**
	Gaussian filtering in-place using separable algorithm.
	@param out Image to filter.
	@param sigma Standard deviation of the Gaussian kernel.
	@param bc Boundary condition.
	*/
	inline void gaussFilter(Image<uint32_t>& img, const Vec3d& sigma, BoundaryCondition bc)
	{
		internals::sepgauss(img, sigma, -1, -1, bc);
	}

	/**
	Gaussian filtering in-place using separable algorithm.
	@param out Image to filter.
	@param sigma Standard deviation of the Gaussian kernel.
	@param bc Boundary condition.
	*/
	inline void gaussFilter(Image<uint32_t>& img, double sigma, BoundaryCondition bc)
	{
		gaussFilter(img, Vec3d(sigma, sigma, sigma), bc);
	}

	/**
	Gaussian filtering in-place using separable algorithm.
	@param out Image to filter.
	@param sigma Standard deviation of the Gaussian kernel.
	@param derivativeDimension Set to negative value to calculate Gaussian filtering, and positive value less than image dimensionality to calculate image derivative in that direction using Gaussian convolution.
	@param bc Boundary condition.
	*/
	inline void gaussFilter(Image<uint64_t>& img, const Vec3d& sigma, BoundaryCondition bc)
	{
		internals::sepgauss(img, sigma, -1, -1, bc);
	}

	/**
	Gaussian filtering in-place using separable algorithm.
	@param out Image to filter.
	@param sigma Standard deviation of the Gaussian kernel.
	@param derivativeDimension Set to negative value to calculate Gaussian filtering, and positive value less than image dimensionality to calculate image derivative in that direction using Gaussian convolution.
	@param bc Boundary condition.
	*/
	inline void gaussFilter(Image<uint64_t>& img, double sigma, BoundaryCondition bc)
	{
		gaussFilter(img, Vec3d(sigma, sigma, sigma), bc);
	}

	/**
	Gaussian filtering in-place using separable algorithm.
	@param out Image to filter.
	@param sigma Standard deviation of the Gaussian kernel.
	@param derivativeDimension Set to negative value to calculate Gaussian filtering, and positive value less than image dimensionality to calculate image derivative in that direction using Gaussian convolution.
	@param bc Boundary condition.
	*/
	template<typename pixel_t>
	typename std::enable_if<std::is_floating_point<pixel_t>::value>::type
	gaussFilter(Image<pixel_t>& img, const Vec3d& sigma, BoundaryCondition bc)
	{
		internals::sepgauss(img, sigma, -1, -1, bc);
	}

	/**
	Gaussian filtering in-place using separable algorithm.
	@param out Image to filter.
	@param sigma Standard deviation of the Gaussian kernel.
	@param derivativeDimension Set to negative value to calculate Gaussian filtering, and positive value less than image dimensionality to calculate image derivative in that direction using Gaussian convolution.
	@param bc Boundary condition.
	*/
	template<typename pixel_t>
	typename std::enable_if<std::is_floating_point<pixel_t>::value>::type
	gaussFilter(Image<pixel_t>& img, double sigma, BoundaryCondition bc)
	{
		gaussFilter(img, Vec3d(sigma, sigma, sigma), bc);
	}





	/**
	Gaussian partial derivative.
	Calculates
	out = dI / dx_i, where I = in and i = derivativeDimension1, or
	out = dI^2 / (dx_i dx_j), where I = in, i = derivativeDimension1, and j = derivativeDimension2.
	If derivativeDimension2 < 0, only first derivative is calculated.
	@param in Input image (not modified).
	@param out Output image.
	@param sigma Standard deviation of the Gaussian kernel.
	@param derivativeDimension1 Dimension where first derivative should be calculated.
	@param derivativeDimension2 Dimension where second derivative should be calculated. Set to negative value to calculate first derivative only.
	@param bc Boundary condition.
	*/
	template<typename input_t, typename output_t>
	typename std::enable_if<std::is_signed<output_t>::value>::type
		gaussDerivative(const Image<input_t>& in, Image<output_t>& out, const Vec3d& sigma, coord_t derivativeDimension1, coord_t derivativeDimension2, BoundaryCondition bc = BoundaryCondition::Nearest, bool showProgressInfo = true)
	{
		internals::sepgauss(in, out, sigma, derivativeDimension1, derivativeDimension2, bc, showProgressInfo);
	}

	/**
	Gaussian partial derivative.
	Calculates
	out = dI / dx_i, where I = in and i = derivativeDimension1, or
	out = dI^2 / (dx_i dx_j), where I = in, i = derivativeDimension1, and j = derivativeDimension2.
	If derivativeDimension2 < 0, only first derivative is calculated.
	@param in Input image (not modified).
	@param out Output image.
	@param sigma Standard deviation of the Gaussian kernel.
	@param derivativeDimension1 Dimension where first derivative should be calculated.
	@param derivativeDimension2 Dimension where second derivative should be calculated. Set to negative value to calculate first derivative only.
	@param bc Boundary condition.
	*/
	template<typename input_t, typename output_t>
	typename std::enable_if<std::is_signed<output_t>::value>::type
		gaussDerivative(const Image<input_t>& in, Image<output_t>& out, double sigma, coord_t derivativeDimension1, coord_t derivativeDimension2, BoundaryCondition bc = BoundaryCondition::Nearest, bool showProgressInfo = true)
	{
		gaussDerivative(in, out, Vec3d(sigma, sigma, sigma), derivativeDimension1, derivativeDimension2, bc, showProgressInfo);
	}

	/**
	Gaussian partial derivative in-place.
	Calculates
	out = dI / dx_i, where I = in and i = derivativeDimension1, or
	out = dI^2 / (dx_i dx_j), where I = in, i = derivativeDimension1, and j = derivativeDimension2.
	If derivativeDimension2 < 0, only first derivative is calculated.
	@param in Input image (not modified).
	@param out Output image.
	@param sigma Standard deviation of the Gaussian kernel.
	@param derivativeDimension1 Dimension where first derivative should be calculated.
	@param derivativeDimension2 Dimension where second derivative should be calculated. Set to negative value to calculate first derivative only.
	@param bc Boundary condition.
	*/
	template<typename pixel_t>
	typename std::enable_if<std::is_signed<pixel_t>::value>::type
		gaussDerivative(Image<pixel_t>& img, const Vec3d& sigma, coord_t derivativeDimension1, coord_t derivativeDimension2, BoundaryCondition bc = BoundaryCondition::Nearest, bool showProgressInfo = true)
	{
		internals::sepgauss(img, sigma, derivativeDimension1, derivativeDimension2, bc, showProgressInfo);
	}

	/**
	Gaussian partial derivative in-place.
	Calculates
	out = dI / dx_i, where I = in and i = derivativeDimension1, or
	out = dI^2 / (dx_i dx_j), where I = in, i = derivativeDimension1, and j = derivativeDimension2.
	If derivativeDimension2 < 0, only first derivative is calculated.
	@param in Input image (not modified).
	@param out Output image.
	@param sigma Standard deviation of the Gaussian kernel.
	@param derivativeDimension1 Dimension where first derivative should be calculated.
	@param derivativeDimension2 Dimension where second derivative should be calculated. Set to negative value to calculate first derivative only.
	@param bc Boundary condition.
	*/
	template<typename pixel_t>
	typename std::enable_if<std::is_signed<pixel_t>::value>::type
		gaussDerivative(Image<pixel_t>& img, double sigma, coord_t derivativeDimension1, coord_t derivativeDimension2, BoundaryCondition bc = BoundaryCondition::Nearest, bool showProgressInfo = true)
	{
		gaussDerivative(img, Vec3d(sigma, sigma, sigma), derivativeDimension1, derivativeDimension2, bc, showProgressInfo);
	}




	/**
	High-pass filtering.
	@param in Input image (not modified).
	@param out Output image.
	@param sigma Standard deviation of the Gaussian kernel.
	@param shift Constant added to all pixel values. Use to shift filtered pixel values from zero to desired value. Useful especially for bandpass filtering unsigned images.
	@param allowOpt Allow separable filtering for 16-bit images and FFT filtering for floating point images.
	@param bc Boundary condition.
	*/
	template<typename pixel_t> void highpassFilter(const Image<pixel_t>& in, Image<pixel_t>& out, const Vec3d& sigma, pixel_t shift = 0, bool allowOpt = true, BoundaryCondition bc = BoundaryCondition::Nearest)
	{
		gaussFilter(in, out, sigma, allowOpt, bc);
		if (shift == 0)
		{
			// Calculate out = in - out
			invsubtract(out, in);
		}
		else
		{
			invsubtractAdd(out, in, shift);
		}
	}


	/**
	Median filter that does not consider nan values.
	@param in Input image.
	@param out Output image.
	@param nbRadius Radius of filtering neighbourhood.
	@param nbType Neighbourhood type.
	@param bc Boundary condition.
	*/
	inline void nanMedianFilter(const Image<float32_t> & in, Image<float32_t> & out, const Vec3c & nbRadius, NeighbourhoodType nbType = NeighbourhoodType::Ellipsoidal, BoundaryCondition bc = BoundaryCondition::Nearest)
	{
		filter<float32_t, float32_t, internals::nanMedianOp >(in, out, nbRadius, nbType, bc);
	}
	
	/**
	Median filter that does not consider nan values.
	@param in Input image.
	@param out Output image.
	@param nbRadius Radius of filtering neighbourhood.
	@param nbType Neighbourhood type.
	@param bc Boundary condition.
	*/
	inline void nanMedianFilter(const Image<float32_t>& in, Image<float32_t>& out, coord_t nbRadius, NeighbourhoodType nbType = NeighbourhoodType::Ellipsoidal, BoundaryCondition bc = BoundaryCondition::Nearest)
	{
		nanMedianFilter(in, out, Vec3c(nbRadius, nbRadius, nbRadius), nbType, bc);
	}


	/*
	Mean, variance, etc. filters
	*/
	// First define macro that creates two shorthand methods, first where neighbourhood radius is vector, and second where neighbourhood
	// radius is the same for all coordinate directions.

	// This version is for minimum and maximum filtering and includes separable filtering optimization for all data types.
#define DEFINE_FILTER_MINMAX(name, help) \
/** \
help \
\
In-place filtering supports only rectangular neighbourhood. \
@param in Image to filter. \
@param nbRadius Radius of filtering neighbourhood. \
@param bc Boundary condition. \
*/ \
template<typename pixel_t> void name##Filter(Image<pixel_t>& img, const Vec3c& nbRadius, BoundaryCondition bc = BoundaryCondition::Nearest) \
{ \
	sepFilter<pixel_t, internals::name##Op<pixel_t> >(img, nbRadius, bc); \
} \
/** \
help \
\
Separable filtering is used for all pixel data types for rectangular neighbourhoods. \
If allowOpt is true, spherical structuring elements larger in radius than 5 are approximated using periodic lines and van Herk algorithm. \
@param in Input image. \
@param out Output image. \
@param nbRadius Radius of filtering neighbourhood. \
@param nbType Neighbourhood type. \
@param bc Boundary condition. \
*/ \
template<typename pixel_t, typename out_t> void name##Filter(const Image<pixel_t>& in, Image<out_t>& out, const Vec3c& nbRadius, NeighbourhoodType nbType = NeighbourhoodType::Ellipsoidal, BoundaryCondition bc = BoundaryCondition::Nearest, bool allowOpt = true) \
{ \
	if(nbType == NeighbourhoodType::Rectangular) \
	{ \
		out.ensureSize(in); \
		setValue<out_t, pixel_t>(out, in); \
		name##Filter<out_t>(out, nbRadius, bc); \
	} \
	else if(allowOpt && nbType == NeighbourhoodType::Ellipsoidal && nbRadius.x == nbRadius.y && nbRadius.x == nbRadius.z && nbRadius.x >= 5) \
	{ \
		out.ensureSize(in); \
		setValue<out_t, pixel_t>(out, in); \
		name##FilterSphereApprox<out_t>(out, nbRadius.x, bc); \
	} \
	else \
	{ \
		filter<pixel_t, out_t, internals::name##Op<pixel_t> >(in, out, nbRadius, nbType, bc); \
	} \
} \
 \
/** \
help \
\
Separable filtering is used for all pixel data types for rectangular neighbourhoods. \
If allowOpt is true, spherical structuring elements larger in radius than 5 are approximated using periodic lines and van Herk algorithm. \
@param in Input image. \
@param out Output image. \
@param nbRadius Radius of filtering neighbourhood. \
@param nbType Neighbourhood type. \
@param bc Boundary condition. \
*/ \
template<typename pixel_t, typename out_t> void name##Filter(const Image<pixel_t>& in, Image<out_t>& out, coord_t nbRadius, NeighbourhoodType nbType = NeighbourhoodType::Ellipsoidal, BoundaryCondition bc = BoundaryCondition::Nearest, bool allowOpt = true) \
{ \
	name##Filter<pixel_t, out_t>(in, out, Vec3c(nbRadius, nbRadius, nbRadius), nbType, bc, allowOpt); \
}




// This version includes separable filtering optimization for float types only.
#define DEFINE_FILTER_SEP_FLOAT(name, help) \
/** \
help \
\
Separable filtering is used for floating point pixel data types for rectangular neighbourhoods. \
@param in Input image. \
@param out Output image. \
@param nbRadius Radius of filtering neighbourhood. \
@param nbType Neighbourhood type. \
@param bc Boundary condition. \
*/ \
template<typename pixel_t, typename out_t> void name##Filter(const Image<pixel_t>& in, Image<out_t>& out, const Vec3c& nbRadius, NeighbourhoodType nbType = NeighbourhoodType::Ellipsoidal, BoundaryCondition bc = BoundaryCondition::Nearest) \
{ \
	filter<pixel_t, out_t, internals::name##Op<pixel_t> >(in, out, nbRadius, nbType, bc); \
} \
 \
/** \
help \
\
Separable filtering is used for floating point pixel data types for rectangular neighbourhoods. \
@param in Input image. \
@param out Output image. \
@param nbRadius Radius of filtering neighbourhood. \
@param nbType Neighbourhood type. \
@param bc Boundary condition. \
*/ \
template<typename pixel_t, typename out_t> void name##Filter(const Image<pixel_t>& in, Image<out_t>& out, coord_t nbRadius, NeighbourhoodType nbType = NeighbourhoodType::Ellipsoidal, BoundaryCondition bc = BoundaryCondition::Nearest) \
{ \
	name##Filter<pixel_t, out_t>(in, out, Vec3c(nbRadius, nbRadius, nbRadius), nbType, bc); \
} \
/** \
help \
\
In-place filtering supports only rectangular neighbourhood. \
@param img Image to process. \
@param nbRadius Radius of filtering neighbourhood. \
@param nbType Neighbourhood type. \
@param bc Boundary condition. \
*/ \
inline void name##Filter(Image<float32_t>& img, const Vec3c& nbRadius, BoundaryCondition bc = BoundaryCondition::Nearest) \
{ \
	sepFilter<float32_t, internals::name##Op<float32_t> >(img, nbRadius, bc); \
} \
/** \
help \
\
In-place filtering supports only rectangular neighbourhood. \
@param img Image to process. \
@param nbRadius Radius of filtering neighbourhood. \
@param nbType Neighbourhood type. \
@param bc Boundary condition. \
*/ \
inline void name##Filter(Image<float32_t>& img, coord_t nbRadius, BoundaryCondition bc = BoundaryCondition::Nearest) \
{ \
	sepFilter<float32_t, internals::name##Op<float32_t> >(img, Vec3c(nbRadius, nbRadius, nbRadius), bc); \
} \
/** \
help \
\
Separable filtering is used for floating point pixel data types for rectangular neighbourhoods. \
@param in Input image. \
@param out Output image. \
@param nbRadius Radius of filtering neighbourhood. \
@param nbType Neighbourhood type. \
@param bc Boundary condition. \
*/ \
template<> inline void name##Filter(const Image<float32_t>& in, Image<float32_t>& out, const Vec3c& nbRadius, NeighbourhoodType nbType, BoundaryCondition bc) \
{ \
	if(nbType == NeighbourhoodType::Rectangular) \
	{ \
		out.ensureSize(in); \
		setValue<float32_t>(out, in); \
		name##Filter(out, nbRadius, bc); \
	} \
	else \
	{ \
		filter<float32_t, float32_t, internals::name##Op<float32_t> >(in, out, nbRadius, nbType, bc); \
	} \
}

	


// Version with no parameters and no separable optimization
#define DEFINE_FILTER(name, help) \
/** \
help \
@param in Input image. \
@param out Output image. \
@param nbRadius Radius of filtering neighbourhood. \
@param nbType Neighbourhood type. \
@param bc Boundary condition. \
*/ \
template<typename pixel_t, typename out_t> void name##Filter(const Image<pixel_t>& in, Image<out_t>& out, const Vec3c& nbRadius, NeighbourhoodType nbType = NeighbourhoodType::Ellipsoidal, BoundaryCondition bc = BoundaryCondition::Nearest) \
{ \
	filter<pixel_t, out_t, internals::name##Op<pixel_t> >(in, out, nbRadius, nbType, bc); \
} \
 \
/** \
help \
@param in Input image. \
@param out Output image. \
@param nbRadius Radius of filtering neighbourhood. \
@param nbType Neighbourhood type. \
@param bc Boundary condition. \
*/ \
template<typename pixel_t, typename out_t> void name##Filter(const Image<pixel_t>& in, Image<out_t>& out, coord_t nbRadius, NeighbourhoodType nbType = NeighbourhoodType::Ellipsoidal, BoundaryCondition bc = BoundaryCondition::Nearest) \
{ \
	name##Filter<pixel_t, out_t>(in, out, Vec3c(nbRadius, nbRadius, nbRadius), nbType, bc); \
}


// Version with one parameter and no separable optimization
#define DEFINE_FILTER_1PARAM(name, paramtype, help, paramhelp) \
/** \
help \
@param in Input image. \
@param out Output image. \
@param nbRadius Radius of filtering neighbourhood. \
@param parameter paramhelp \
@param nbType Neighbourhood type. \
@param bc Boundary condition. \
*/ \
template<typename pixel_t, typename out_t> void name##Filter(const Image<pixel_t>& in, Image<out_t>& out, const Vec3c& nbRadius, paramtype parameter, NeighbourhoodType nbType = NeighbourhoodType::Ellipsoidal, BoundaryCondition bc = BoundaryCondition::Nearest) \
{ \
	filter<pixel_t, out_t, paramtype, internals::name##Op<pixel_t> >(in, out, nbRadius, parameter, nbType, bc); \
} \
 \
/** \
help \
@param in Input image. \
@param out Output image. \
@param nbRadius Radius of filtering neighbourhood. \
@param parameter paramhelp \
@param nbType Neighbourhood type. \
@param bc Boundary condition. \
*/ \
template<typename pixel_t, typename out_t> void name##Filter(const Image<pixel_t>& in, Image<out_t>& out, coord_t nbRadius, paramtype parameter, NeighbourhoodType nbType = NeighbourhoodType::Ellipsoidal, BoundaryCondition bc = BoundaryCondition::Nearest) \
{ \
	name##Filter<pixel_t, out_t>(in, out, Vec3c(nbRadius, nbRadius, nbRadius), parameter, nbType, bc); \
}



	// Now define the filtering operations
	DEFINE_FILTER_MINMAX(min, Calculates minimum filtering.)
	DEFINE_FILTER_MINMAX(max, Calculates maximum filtering.)
	DEFINE_FILTER_SEP_FLOAT(mean, Calculates mean filtering.)
	DEFINE_FILTER(median, Calculates median filtering.)
	DEFINE_FILTER_1PARAM(maskedMedian, pixel_t, Calculates masked median filtering., Image value that should not be considered when calculating median.)

	#define COMMA ,
	DEFINE_FILTER_1PARAM(vawe, double, Calculates variance weighted mean filtering., Standard deviation of noise. For a rough order of magnitude estimateCOMMA measure standard deviation from a region that does not contain any features.)


	// Variance requires special handling
	/**
	Calculates variance filtering.

	Separable optimization is used for rectangular neighbourhoods.
	@param in Input image.
	@param out Output image.
	@param nbRadius Radius of filtering neighbourhood.
	@param nbType Neighbourhood type.
	@param bc Boundary condition (BoundaryCondition::Zero or Nearest).
	*/
	template<typename pixel_t, typename out_t> void varianceFilter(const Image<pixel_t>& in, Image<out_t>& out, const Vec3c& nbRadius, NeighbourhoodType nbType = NeighbourhoodType::Ellipsoidal, BoundaryCondition bc = BoundaryCondition::Nearest)
	{
		filter<pixel_t, out_t, internals::varianceOp<pixel_t> >(in, out, nbRadius, nbType, bc);
	}
	
	/**
	Calculates variance filtering.

	Separable optimization is used for rectangular neighbourhoods.
	@param in Input image.
	@param out Output image.
	@param nbRadius Radius of filtering neighbourhood.
	@param nbType Neighbourhood type.
	@param bc Boundary condition (BoundaryCondition::Zero or Nearest).
	*/
	template<typename pixel_t, typename out_t> void varianceFilter(const Image<pixel_t>& in, Image<out_t>& out, coord_t nbRadius, NeighbourhoodType nbType = NeighbourhoodType::Ellipsoidal, BoundaryCondition bc = BoundaryCondition::Nearest)
	{
		varianceFilter<pixel_t, out_t>(in, out, Vec3c(nbRadius, nbRadius, nbRadius), nbType, bc);
	}

	namespace internals
	{
		inline float32_t rectVarianceFinalization(float32_t out, float32_t tmp)
		{
			return (float32_t)tmp - (float32_t)out * (float32_t)out;
		}
	}

	template<> inline void varianceFilter(const Image<float32_t>& in, Image<float32_t>& out, const Vec3c& nbRadius, NeighbourhoodType nbType, BoundaryCondition bc)
	{
		if(nbType == NeighbourhoodType::Rectangular)
		{
			out.ensureSize(in);

			// Calculate mean filtering of in. For now on, out = mean(in)
			setValue(out, in);
			sepFilter<float32_t, internals::meanOp<float32_t> >(out, nbRadius, bc);

			// Calculate in^2. For now on, out contains in^2.
			Image<float32_t> tmp;
			tmp.ensureSize(in);
			setValue<float32_t>(tmp, in);
			multiply(tmp, tmp);

			// Calculate mean of orig^2. For now on, tmp = mean(in^2).
			sepFilter<float32_t, internals::meanOp<float32_t> >(tmp, nbRadius, bc);

			// Calculate variance with mean(in^2) - mean(in)^2 = tmp - out * out
			pointProcessImageImage<float32_t, float32_t, float32_t, internals::rectVarianceFinalization>(out, tmp, false);
		}
		else
		{
			filter<float32_t, float32_t, internals::varianceOp<float32_t> >(in, out, nbRadius, nbType, bc);
		}
	}


	/**
	Calculates standard deviation filtering.

	Separable optimization is used for rectangular neighbourhoods.
	@param in Input image.
	@param out Output image.
	@param nbRadius Radius of filtering neighbourhood.
	@param nbType Neighbourhood type.
	@param bc Boundary condition (BoundaryCondition::Zero or Nearest).
	*/
	template<typename pixel_t, typename out_t> void stddevFilter(const Image<pixel_t>& in, Image<out_t>& out, const Vec3c& nbRadius, NeighbourhoodType nbType = NeighbourhoodType::Ellipsoidal, BoundaryCondition bc = BoundaryCondition::Nearest)
	{
		filter<pixel_t, out_t, internals::stddevOp<pixel_t> >(in, out, nbRadius, nbType, bc);
	}

	/**
	Calculates standard deviation filtering.

	Separable optimization is used for rectangular neighbourhoods.
	@param in Input image.
	@param out Output image.
	@param nbRadius Radius of filtering neighbourhood.
	@param nbType Neighbourhood type.
	@param bc Boundary condition (BoundaryCondition::Zero or Nearest).
	*/
	template<typename pixel_t, typename out_t> void stddevFilter(const Image<pixel_t>& in, Image<out_t>& out, coord_t nbRadius, NeighbourhoodType nbType = NeighbourhoodType::Ellipsoidal, BoundaryCondition bc = BoundaryCondition::Nearest)
	{
		stddevFilter<pixel_t, out_t>(in, out, Vec3c(nbRadius, nbRadius, nbRadius), nbType, bc);
	}


	/**
	Opening operation.
	@param img Image that is to be processed.
	@param tmp Temporary image having the same size as the input image. If the size is wrong, it is changed by this method.
	@param nbType Neighbourhood type.
	@param nbRadius Neighbourhood radius.
	*/
	template<typename pixel_t> void openingFilter(Image<pixel_t>& img, Image<pixel_t>& tmp, const Vec3c& nbRadius, NeighbourhoodType nbType = NeighbourhoodType::Ellipsoidal, BoundaryCondition bc = BoundaryCondition::Nearest, bool allowOpt = true)
	{
		tmp.ensureSize(img);
		minFilter<pixel_t, pixel_t>(img, tmp, nbRadius, nbType, bc, allowOpt);
		maxFilter<pixel_t, pixel_t>(tmp, img, nbRadius, nbType, bc, allowOpt);
	}

	/**
	Opening operation.
	@param img Image that is to be processed.
	@param tmp Temporary image having the same size as the input image. If the size is wrong, it is changed by this method.
	@param nbType Neighbourhood type.
	@param nbRadius Neighbourhood radius.
	*/
	template<typename pixel_t> void openingFilter(Image<pixel_t>& img, Image<pixel_t>& tmp, coord_t nbRadius, NeighbourhoodType nbType = NeighbourhoodType::Ellipsoidal, BoundaryCondition bc = BoundaryCondition::Nearest, bool allowOpt = true)
	{
		openingFilter(img, tmp, Vec3c(nbRadius, nbRadius, nbRadius), nbType, bc, allowOpt);
	}

	/**
	Closing operation.
	@param img Image that is to be processed.
	@param tmp Temporary image having the same size as the input image. If the size is wrong, it is changed by this method.
	@param nbType Neighbourhood type.
	@param nbRadius Neighbourhood radius.
	*/
	template<typename pixel_t> void closingFilter(Image<pixel_t>& img, Image<pixel_t>& tmp, const Vec3c& nbRadius, NeighbourhoodType nbType = NeighbourhoodType::Ellipsoidal, BoundaryCondition bc = BoundaryCondition::Nearest, bool allowOpt = true)
	{
		tmp.ensureSize(img);
		maxFilter<pixel_t, pixel_t>(img, tmp, nbRadius, nbType, bc, allowOpt);
		minFilter<pixel_t, pixel_t>(tmp, img, nbRadius, nbType, bc, allowOpt);
	}

	/**
	Closing operation.
	@param img Image that is to be processed.
	@param tmp Temporary image having the same size as the input image. If the size is wrong, it is changed by this method.
	@param nbType Neighbourhood type.
	@param nbRadius Neighbourhood radius.
	*/
	template<typename pixel_t> void closingFilter(Image<pixel_t>& img, Image<pixel_t>& tmp, coord_t nbRadius, NeighbourhoodType nbType = NeighbourhoodType::Ellipsoidal, BoundaryCondition bc = BoundaryCondition::Nearest, bool allowOpt = true)
	{
		closingFilter(img, tmp, Vec3c(nbRadius, nbRadius, nbRadius), nbType, bc, allowOpt);
	}

	/**
	Calculates bilateral filtering.
	@param in Input image.
	@param out Output image.
	@param spatialSigma Standard deviation of Gaussian kernel used for spatial smoothing.
	@param rangeSigma Standard deviation of Gaussian kernel used to avoid smoothing edges of features. Order of magnitude must be similar to difference between gray levels of background and objects.
	@param bc Boundary condition (BoundaryCondition::Zero or Nearest).
	*/
	template<typename pixel_t, typename out_t> void bilateralFilter(const Image<pixel_t>& in, Image<out_t>& out, double spatialSigma, double rangeSigma, BoundaryCondition bc = BoundaryCondition::Nearest)
	{
		coord_t nbRadius = (coord_t)ceil(3 * spatialSigma);
		filter<pixel_t, out_t, Vec2d, internals::bilateralOp<pixel_t> >(in, out, Vec3c(nbRadius, nbRadius, nbRadius), Vec2d(spatialSigma, rangeSigma), NeighbourhoodType::Rectangular, bc);
	}



	/**
	Performs morphological reconstruction for binary, label, or grayscale image.
	Uses simple constrained dilation algorithm. Faster algorithms are be available in particular for binary and label images e.g. in 
	Vincent - Morphological Grayscale Reconstruction in Image Analysis: Applications and Efficient Algorithms.
	@param img Image that contains seed points. The result of the reconstruction is placed into this image.
	@param mask Image that contains the mask. The reconstruction is constrained to non-zero pixels of this image.
	*/
	template<typename pixel_t, typename mask_t> size_t morphoRec(Image<pixel_t>& img, const Image<mask_t>& mask)
	{
		size_t totalChanged = 0;
		size_t changed = 0;
		size_t round = 0;
		do
		{
			changed = sepFilterImageImage<pixel_t, mask_t, internals::morphoRecOp<pixel_t, mask_t> >(img, mask, Vec3c(1, 1, 1), BoundaryCondition::Zero);
			round++;
			std::cout << "Round " << round << ", " << changed << " changes." << std::endl;
			totalChanged += changed;
		}
		while (changed > 0);

		return totalChanged;
	}

	
	namespace tests
	{
		void filters();
		void stddevuint16();
		void separableOptimization();
		void gaussFilters();
		void bilateral();
	}
}
