#pragma once

#include <algorithm>

#include "image.h"
#include "neighbourhood.h"
#include "math/mathutils.h"
#include "pointprocess.h"
#include "fft.h"
#include "utilities.h"

using math::Vec3c;
using math::Vec2d;

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
	template<typename pixel_t, typename out_t, double processNeighbourhood(const Image<pixel_t>& nb, const Image<pixel_t>& mask)>
	void filter(const Image<pixel_t>& img, Image<out_t>& out, const math::Vec3c& nbRadius, NeighbourhoodType nbType, BoundaryCondition bc)
	{
		out.mustNotBe(img);
		out.ensureSize(img);

		Image<pixel_t> mask;
		createNeighbourhoodMask(nbType, nbRadius, mask);

		size_t totalProcessed = 0;
		#pragma omp parallel if(!omp_in_parallel() && img.pixelCount() > PARALLELIZATION_THRESHOLD)
		{
			
			Image<pixel_t> nb(mask.dimensions());

			#pragma omp for
			for (coord_t z = 0; z < img.depth(); z++)
			{
				for (coord_t y = 0; y < img.height(); y++)
				{
					for (coord_t x = 0; x < img.width(); x++)
					{
						getNeighbourhood(img, math::Vec3c(x, y, z), nbRadius, nb, bc);

						out(x, y, z) = math::pixelRound<out_t>(processNeighbourhood(nb, mask));
					}
				}

				showThreadProgress(totalProcessed, img.depth());
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
	template<typename pixel_t, typename out_t, typename param_t, double processNeighbourhood(const Image<pixel_t>& nb, const Image<pixel_t>& mask, param_t parameter)>
	void filter(const Image<pixel_t>& img, Image<out_t>& out, const math::Vec3c& nbRadius, const param_t parameter, NeighbourhoodType nbType, BoundaryCondition bc)
	{
		out.mustNotBe(img);
		out.ensureSize(img);

		Image<pixel_t> mask;
		createNeighbourhoodMask(nbType, nbRadius, mask);

		size_t totalProcessed = 0;
		#pragma omp parallel if(!omp_in_parallel() && img.pixelCount() > PARALLELIZATION_THRESHOLD)
		{

			Image<pixel_t> nb(mask.dimensions());

			#pragma omp for
			for (coord_t z = 0; z < img.depth(); z++)
			{
				for (coord_t y = 0; y < img.height(); y++)
				{
					for (coord_t x = 0; x < img.width(); x++)
					{
						getNeighbourhood(img, math::Vec3c(x, y, z), nbRadius, nb, bc);
						out(x, y, z) = math::pixelRound<out_t>(processNeighbourhood(nb, mask, parameter));
					}
				}

				showThreadProgress(totalProcessed, img.depth());
			}
		}
	}

	namespace internals
	{
		template<typename pixel_t, typename param_t, double processNeighbourhood(const Image<pixel_t>& nb, const Image<pixel_t>& mask, param_t param)>
		void sepFilterOneDimension(Image<pixel_t>& img, coord_t r, size_t dim, param_t param, BoundaryCondition bc)
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
							pixel_t edgeVal = bc == Zero ? pixel_t() : img(0, y, z);
							setValue(buffer, edgeVal);
							for (coord_t x = 0; x < math::min(r + 1, img.width()); x++)
							{
								buffer(r + x) = img(x, y, z);
							}

							edgeVal = bc == Zero ? pixel_t() : img(img.width() - 1, y, z);
							for (coord_t x = math::min(r + 1, img.width()); x < N - r; x++)
							{
								buffer(r + x) = edgeVal;
							}

							// Process
							for (coord_t x = 0; x < img.width() - r - 1; x++)
							{
								img(x, y, z) = math::pixelRound<pixel_t>(processNeighbourhood(buffer, mask, param));

								for (coord_t i = 0; i < N - 1; i++)
									buffer(i) = buffer(i + 1);
								buffer(N - 1) = img(x + r + 1, y, z);
							}

							// Process end of line
							for (coord_t x = math::max((coord_t)0, img.width() - r - 1); x < img.width(); x++)
							{
								img(x, y, z) = math::pixelRound<pixel_t>(processNeighbourhood(buffer, mask, param));

								for (coord_t i = 0; i < N - 1; i++)
									buffer(i) = buffer(i + 1);
								buffer(N - 1) = edgeVal;
							}

						}

						showThreadProgress(counter, img.depth());
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
							pixel_t edgeVal = bc == Zero ? pixel_t() : img(x, 0, z);
							setValue(buffer, edgeVal);
							for (coord_t y = 0; y < math::min(r + 1, img.height()); y++)
							{
								buffer(r + y) = img(x, y, z);
							}

							edgeVal = bc == Zero ? pixel_t() : img(x, img.height() - 1, z);
							for (coord_t y = math::min(r + 1, img.height()); y < N - r; y++)
							{
								buffer(r + y) = edgeVal;
							}

							// Process
							for (coord_t y = 0; y < img.height() - r - 1; y++)
							{
								img(x, y, z) = math::pixelRound<pixel_t>(processNeighbourhood(buffer, mask, param));

								for (coord_t i = 0; i < N - 1; i++)
									buffer(i) = buffer(i + 1);
								buffer(N - 1) = img(x, y + r + 1, z);
							}

							// Process end of line
							for (coord_t y = img.height() - r - 1; y < img.height(); y++)
							{
								img(x, y, z) = math::pixelRound<pixel_t>(processNeighbourhood(buffer, mask, param));

								for (coord_t i = 0; i < N - 1; i++)
									buffer(i) = buffer(i + 1);
								buffer(N - 1) = edgeVal;
							}

						}

						showThreadProgress(counter, img.depth());
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
							pixel_t edgeVal = bc == Zero ? pixel_t() : img(x, y, 0);
							setValue(buffer, edgeVal);
							for (coord_t z = 0; z < math::min(r + 1, img.depth()); z++)
							{
								buffer(r + z) = img(x, y, z);
							}

							edgeVal = bc == Zero ? pixel_t() : img(x, y, img.depth() - 1);
							for (coord_t z = math::min(r + 1, img.depth()); z < N - r; z++)
							{
								buffer(r + z) = edgeVal;
							}

							// Process
							for (coord_t z = 0; z < img.depth() - r - 1; z++)
							{
								img(x, y, z) = math::pixelRound<pixel_t>(processNeighbourhood(buffer, mask, param));

								for (coord_t i = 0; i < N - 1; i++)
									buffer(i) = buffer(i + 1);
								buffer(N - 1) = img(x, y, z + r + 1);
							}

							// Process end of line
							for (coord_t z = img.depth() - r - 1; z < img.depth(); z++)
							{
								img(x, y, z) = math::pixelRound<pixel_t>(processNeighbourhood(buffer, mask, param));

								for (coord_t i = 0; i < N - 1; i++)
									buffer(i) = buffer(i + 1);
								buffer(N - 1) = edgeVal;
							}

						}

						showThreadProgress(counter, img.height());
					}
				}
			}
			else
			{
				throw ITLException("Invalid dimension.");
			}
		}

		template<typename pixel_t, double processNeighbourhood(const Image<pixel_t>& nb, const Image<pixel_t>& mask)> double paramRemover(const Image<pixel_t>& nb, const Image<pixel_t>& mask, int param)
		{
			return processNeighbourhood(nb, mask);
		}

		template<typename pixel_t, double processNeighbourhood(const Image<pixel_t>& nb, const Image<pixel_t>& mask)>
		void sepFilterOneDimension(Image<pixel_t>& img, coord_t r, size_t dim, BoundaryCondition bc)
		{
			internals::sepFilterOneDimension<pixel_t, int, internals::paramRemover<pixel_t, processNeighbourhood> >(img, r, dim, 0, bc);
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
	template<typename pixel_t, double processNeighbourhood(const Image<pixel_t>& nb, const Image<pixel_t>& mask)>
	void sepFilter(Image<pixel_t>& img, const math::Vec3c& nbRadius, BoundaryCondition bc)
	{
		for (size_t n = 0; n < img.dimensionality(); n++)
		{
			internals::sepFilterOneDimension<pixel_t, processNeighbourhood>(img, nbRadius[n], n, bc);
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
	template<typename pixel_t, typename param_t, double processNeighbourhood(const Image<pixel_t>& nb, const Image<pixel_t>& mask, param_t param)>
	void sepFilter(Image<pixel_t>& img, const math::Vec3c& nbRadius, const math::Vec3<param_t>& params, BoundaryCondition bc)
	{
		for (size_t n = 0; n < img.dimensionality(); n++)
		{
			internals::sepFilterOneDimension<pixel_t, param_t, processNeighbourhood>(img, nbRadius[n], n, params[n], bc);
		}
	}

	namespace internals
	{
		template<typename pixel_t> double meanOp(const Image<pixel_t>& nb, const Image<pixel_t>& mask)
		{
			double sum = 0;
			double count = 0;
			for (coord_t n = 0; n < nb.pixelCount(); n++)
			{
				double val = (double)nb(n);
				double w = (double)mask(n);

				sum += val * w;
				count += w;
			}
			return sum / count;
		}

		template<typename pixel_t> double varianceOp(const Image<pixel_t>& nb, const Image<pixel_t>& mask)
		{
			double count = 0;
			double total = 0;
			double total2 = 0;

			for (coord_t n = 0; n < nb.pixelCount(); n++)
			{
				double val = (double)nb(n);
				double w = (double)mask(n);
				count += w; // Note: w is either 0 or 1
				val *= w;
				total += val;
				total2 += val * val;
			}

			double var = (double)((total2 - (total * total / count)) / (count - 1));

			return var;
		}

		template<typename pixel_t> double medianOp(const Image<pixel_t>& nb, const Image<pixel_t>& mask)
		{
			vector<pixel_t> values;
			values.reserve(nb.pixelCount());

			for (coord_t n = 0; n < nb.pixelCount(); n++)
			{
				if (mask(n) != 0)
					values.push_back(nb(n));
			}

			size_t size = values.size();

			if (size == 0)
			{
				return pixel_t();
			}
			else if (size == 1)
			{
				return values[0];
			}
			else
			{
				std::sort(values.begin(), values.end());
				if (size % 2 == 0)
				{
					return (values[size / 2 - 1] + values[size / 2]) / 2.0;
				}
				else
				{
					return values[size / 2];
				}
			}
		}

		template<typename pixel_t> double maskedMedianOp(const Image<pixel_t>& nb, const Image<pixel_t>& mask, pixel_t badValue)
		{
			vector<pixel_t> values;
			values.reserve(nb.pixelCount());

			for (coord_t n = 0; n < nb.pixelCount(); n++)
			{
				if (mask(n) != 0 && nb(n) != badValue)
					values.push_back(nb(n));
			}

			size_t size = values.size();

			if (size == 0)
			{
				return pixel_t();
			}
			else if (size == 1)
			{
				return values[0];
			}
			else
			{
				std::sort(values.begin(), values.end());
				if (size % 2 == 0)
				{
					return (values[size / 2 - 1] + values[size / 2]) / 2.0;
				}
				else
				{
					return values[size / 2];
				}
			}
		}

		template<typename pixel_t> double minOp(const Image<pixel_t>& nb, const Image<pixel_t>& mask)
		{
			double res = numeric_limits<double>::max();
			for (coord_t n = 0; n < nb.pixelCount(); n++)
			{
				double val = (double)nb(n);
				if (mask(n) != 0 && val < res)
					res = val;

			}
			return res;
		}

		template<typename pixel_t> double maxOp(const Image<pixel_t>& nb, const Image<pixel_t>& mask)
		{
			double res = numeric_limits<double>::lowest();
			for (coord_t n = 0; n < nb.pixelCount(); n++)
			{
				double val = (double)nb(n);
				if (mask(n) != 0 && val > res)
					res = val;

			}
			return res;
		}

		/**
		Variance weighted mean filter.
		*/
		template<typename pixel_t> double vaweOp(const Image<pixel_t>& nb, const Image<pixel_t>& mask, double noiseStdDev)
		{
			double sum = 0.0;
			double sum_sqr = 0.0;
			double count = 0;
			for (coord_t n = 0; n < nb.pixelCount(); n++)
			{
				if (mask(n) != 0)
				{
					double pix = (double)nb(n);
					sum += pix;
					sum_sqr += pix * pix;
					count++;
				}
			}

			double mean = sum / count;
			double lvar = (sum_sqr - count * mean * mean) / (count - 1.0);
			double gvar = noiseStdDev * noiseStdDev;

			double multip;
			if (gvar <= lvar)
				multip = gvar / lvar;
			else
				multip = 1.0;

			double pixel = (double)nb(nb.pixelCount() / 2);
			return pixel - (multip * (pixel - mean));
		}

		/*
		static double normTable[] = { 0.3989422804014327,
					0.3969525474770118,
					0.3910426939754559,
					0.3813878154605241,
					0.3682701403033233,
					0.3520653267642995,
					0.3332246028917997,
					0.3122539333667613,
					0.2896915527614827,
					0.2660852498987548,
					0.2419707245191434,
					0.2178521770325505,
					0.1941860549832129,
					0.1713685920478074,
					0.1497274656357448,
					0.1295175956658917,
					0.1109208346794555,
					0.09404907737688691,
					0.07895015830089415,
					0.06561581477467658,
					0.05399096651318806,
					0.04398359598042719,
					0.03547459284623142,
					0.02832703774160116,
					0.02239453029484288,
					0.01752830049356854,
					0.01358296923368561,
					0.01042093481442259,
					0.007915451582979956,
					0.005952532419775849,
					0.004431848411938008,
					0.003266819056199918,
					0.00238408820146484,
					0.001722568939053678,
					0.001232219168473018,
					0.0008726826950457602,
					0.0006119019301137719,
					0.0004247802705507514,
					0.00029194692579146,
					0.0001986554713927724,
					0.0001338302257648854,
					8.926165717713262e-005,
					5.894306775653985e-005,
					3.853519674208713e-005,
					2.494247129005354e-005,
					1.598374110690547e-005,
					1.014085206548672e-005,
					6.36982517886709e-006,
					3.961299091032061e-006,
					2.438960745893352e-006,
					1.486719514734298e-006,
					8.972435162383307e-007,
					5.361035344697615e-007,
					3.171349216715964e-007,
					1.85736184455529e-007,
					1.076976004254328e-007,
					6.182620500165824e-008,
					3.513955094820434e-008,
					1.97731964062446e-008,
					1.101576362468231e-008,
					6.075882849823286e-009,
					3.317884243547281e-009,
					1.793783907964079e-009,
					9.601433370312266e-010,
					5.088140281645039e-010,
					2.669556614762852e-010,
					1.386679994165307e-010,
					7.131328123996076e-011,
					3.630961501791775e-011,
					1.830332217015571e-011,
					9.134720408364594e-012,
					4.513543677205485e-012,
					2.207989963137139e-012,
					1.069383787154157e-012,
					5.127753636796663e-013,
					2.43432053302901e-013,
					1.144156490180133e-013,
					5.324148372252943e-014,
					2.452855285696415e-014,
					1.118795621435182e-014,
					5.052271083536893e-015,
					2.258809403154303e-015,
					9.998378748497037e-016,
					4.381639435509327e-016,
					1.901081537907964e-016,
					8.16623563166955e-017,
					3.472962748566208e-017,
					1.462296357500637e-017,
					6.095758129562418e-018,
					2.515805776951405e-018,
					1.027977357166892e-018,
					4.15859897911516e-019,
					1.665588032379905e-019,
					6.604579860739308e-020,
					2.592864701100371e-020,
					1.007793539430001e-020,
					3.878111931746907e-021,
					1.477495492704244e-021,
					5.573000022720691e-022,
					2.081176820202825e-022,
					7.69459862670642e-023
				};

		inline double normPdfApprox(double x, double mu, double sigma)
		{
			x = fabs(x / sigma);
			
			if (x >= 10)
				return 0;

			int i = (int)x;
			double f = x - i;

			i *= 10;
			double  res = normTable[i] * (1 - f) + normTable[i + 1] * f;
			return res / sigma + mu;
		}
		*/

		/**
		Bilateral filter
		Mask is not used.
		*/
		template<typename pixel_t> double bilateralOp(const Image<pixel_t>& nb, const Image<pixel_t>& mask, math::Vec2d spatialandrangesigma)
		{
			double sigmas = spatialandrangesigma[0];
			double sigmat = spatialandrangesigma[1];


			double sum = 0.0;
			double wsum = 0.0;

			math::Vec3c center((nb.dimensions() - math::Vec3c(1, 1, 1)) / 2);
			double centerVal = (double)nb(center);

			for (coord_t z = 0; z < nb.depth(); z++)
			{
				for (coord_t y = 0; y < nb.height(); y++)
				{
					for (coord_t x = 0; x < nb.width(); x++)
					{
						double pix = (double)nb(x, y, z);
						
						double r2 = (double)(math::Vec3c(x, y, z) - center).normSquared();
						//double r2 = (x - cx) * (x - cx) + (y - cy) * (y - cy) + (z - cz) * (z - cz);
						double c = pix - centerVal;

						double w = (1 / (sigmas * sqrt(2 * math::PI)) * ::exp(-(r2) / (2 * sigmas * sigmas))) *	    // Spatial
							(1 / (sigmat * sqrt(2 * math::PI)) * ::exp(-(c * c) / (2 * sigmat * sigmat)));   // Range
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
		template<typename pixel_t> double convolution1DOp(const Image<pixel_t>& nb, const Image<pixel_t>& mask, const Image<float32_t>* kernel)
		{
			double sum = 0;
			coord_t N = kernel->pixelCount() - 1;
			for(coord_t n = 0; n <= N; n++)
			{
				sum += (double)nb(n) * (double)(*kernel)(N - n);
			}
			return sum;
		}

		/*
		Convolution.
		Mask is not used, so this works only with rectangular neighbourhood.
		Kernel size must match neighbourhood size.
		*/
		template<typename pixel_t> double convolution3DOp(const Image<pixel_t>& nb, const Image<pixel_t>& mask, const Image<float32_t>* kernel)
		{
			double sum = 0;
			Vec3c N = kernel->dimensions() - Vec3c(1, 1, 1);
			for (coord_t z = 0; z < nb.depth(); z++)
			{
				for (coord_t y = 0; y < nb.height(); y++)
				{
					for (coord_t x = 0; x < nb.width(); x++)
					{
						double pix = (double)nb(x, y, z);
						double weight = (double)(*kernel)(N.x - x, N.y - y, N.z - z);
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

				nbRadius = (coord_t)ceil(3.0 * sigma);

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
		inline void gaussianKernel3D(const math::Vec3d& sigma, math::Vec3c& nbRadius, Image<float32_t>& kernel, coord_t derivativeDimension1, coord_t derivativeDimension2)
		{
			if (derivativeDimension1 > 2 || derivativeDimension2 > 2 || (derivativeDimension1 < 0 && derivativeDimension2 >= 0))
				throw ITLException("Invalid derivative dimension.");

			nbRadius = componentwiseCeil(3.0 * sigma);

			math::Vec3c kernelSize = 2 * nbRadius + math::Vec3c(1, 1, 1);
			kernel.ensureSize(kernelSize);

			float32_t kernelSum = 0.0;
			for (coord_t z = 0; z < kernel.depth(); z++)
			{
				for (coord_t y = 0; y < kernel.height(); y++)
				{
					for (coord_t x = 0; x < kernel.width(); x++)
					{
						math::Vec3c X(x, y, z);
						X -= nbRadius;

						float32_t sum = 0.0;

						// Add contribution from each dimension.
						for (size_t i = 0; i < X.size(); i++)
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


		/**
		Non-separable Gaussian convolution.
		*/
		template<typename pixel_t> void gauss(const Image<pixel_t>& in, Image<pixel_t>& out, const math::Vec3d& sigma, coord_t derivativeDimension1, coord_t derivativeDimension2, BoundaryCondition bc)
		{
			math::Vec3c nbRadius;
			Image<float32_t> kernel;

			gaussianKernel3D(sigma, nbRadius, kernel, derivativeDimension1, derivativeDimension2);

			// Filter
			filter<pixel_t, pixel_t, const Image<float32_t>*, internals::convolution3DOp<pixel_t> >(in, out, nbRadius, &kernel, Rectangular, bc);
		}

		/**
		Separable Gaussian filtering in-place.
		*/
		template<typename pixel_t> void sepgauss(Image<pixel_t>& img, const math::Vec3d& sigma, coord_t derivativeDimension1, coord_t derivativeDimension2, BoundaryCondition bc)
		{
			math::Vec3<const Image<float32_t>* > kernels;
			math::Vec3c nbRadius;
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

			sepFilter<pixel_t, const Image<float32_t>*, internals::convolution1DOp<pixel_t> >(img, nbRadius, kernels, bc);
		}

		/**
		Separable Gaussian filtering.
		Use only if data type has good enough accuracy.
		*/
		template<typename input_t, typename output_t> void sepgauss(const Image<input_t>& in, Image<output_t>& out, const math::Vec3d& sigma, coord_t derivativeDimension1, coord_t derivativeDimension2, BoundaryCondition bc)
		{
			setValue(out, in);
			sepgauss(out, sigma, derivativeDimension1, derivativeDimension2, bc);
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
	gaussFilter(const Image<pixel_t>& in, Image<pixel_t>& out, const math::Vec3d& sigma, bool allowOpt = true, BoundaryCondition bc = Zero)
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
	template<typename pixel_t> void gaussFilter(const Image<pixel_t>& in, Image<pixel_t>& out, double sigma, bool allowOpt = true, BoundaryCondition bc = Zero)
	{
		gaussFilter(in, out, math::Vec3d(sigma, sigma, sigma), allowOpt, bc);
	}

	template<typename pixel_t>
	typename std::enable_if<std::is_integral<pixel_t>::value && (sizeof(pixel_t) >= 2)>::type
	gaussFilter(const Image<pixel_t>& in, Image<pixel_t>& out, const math::Vec3d& sigma, bool allowOpt, BoundaryCondition bc)
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

	//template<> inline void gaussFilter<uint16_t>(const Image<uint16_t>& in, Image<uint16_t>& out, const math::Vec3d& sigma, bool allowOpt, BoundaryCondition bc)
	//{
	//	// Perform separable filtering if allowOpt is true
	//	// Perform non-separable filtering if allowOpt is false
	//	if (allowOpt)
	//	{
	//		internals::sepgauss(in, out, sigma, -1, -1, bc);
	//	}
	//	else
	//	{
	//		internals::gauss(in, out, sigma, -1, -1, bc);
	//	}
	//}

	//template<> inline void gaussFilter<uint32_t>(const Image<uint32_t>& in, Image<uint32_t>& out, const math::Vec3d& sigma, bool allowOpt, BoundaryCondition bc)
	//{
	//	// Perform separable filtering if allowOpt is true
	//	// Perform non-separable filtering if allowOpt is false
	//	if (allowOpt)
	//	{
	//		internals::sepgauss(in, out, sigma, -1, -1, bc);
	//	}
	//	else
	//	{
	//		internals::gauss(in, out, sigma, -1, -1, bc);
	//	}
	//}

	//template<> inline void gaussFilter<uint64_t>(const Image<uint64_t>& in, Image<uint64_t>& out, const math::Vec3d& sigma, bool allowOpt, BoundaryCondition bc)
	//{
	//	// Perform separable filtering if allowOpt is true
	//	// Perform non-separable filtering if allowOpt is false
	//	if (allowOpt)
	//	{
	//		internals::sepgauss(in, out, sigma, -1, -1, bc);
	//	}
	//	else
	//	{
	//		internals::gauss(in, out, sigma, -1, -1, bc);
	//	}
	//}

	//template<> inline void gaussFilter<float32_t>(const Image<float32_t>& in, Image<float32_t>& out, const math::Vec3d& sigma, bool allowOpt, BoundaryCondition bc)
	inline void gaussFilter(const Image<float32_t>& in, Image<float32_t>& out, const math::Vec3d& sigma, bool allowOpt, BoundaryCondition bc)
	{
		// Perform FFT filtering if allowOpt is true
		// Perform separable filtering if allowOpt is false
		if (allowOpt && bc == Zero)
		{
			setValue(out, in);
			gaussFilter(out, sigma);
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
	inline void gaussFilter(Image<uint16_t>& img, const math::Vec3d& sigma, BoundaryCondition bc)
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
	inline void gaussFilter(Image<uint16_t>& img, double sigma, BoundaryCondition bc)
	{
		gaussFilter(img, math::Vec3d(sigma, sigma, sigma), bc);
	}

	/**
	Gaussian filtering in-place using separable algorithm.
	@param out Image to filter.
	@param sigma Standard deviation of the Gaussian kernel.
	@param derivativeDimension Set to negative value to calculate Gaussian filtering, and positive value less than image dimensionality to calculate image derivative in that direction using Gaussian convolution.
	@param bc Boundary condition.
	*/
	inline void gaussFilter(Image<uint32_t>& img, const math::Vec3d& sigma, BoundaryCondition bc)
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
	inline void gaussFilter(Image<uint32_t>& img, double sigma, BoundaryCondition bc)
	{
		gaussFilter(img, math::Vec3d(sigma, sigma, sigma), bc);
	}

	/**
	Gaussian filtering in-place using separable algorithm.
	@param out Image to filter.
	@param sigma Standard deviation of the Gaussian kernel.
	@param derivativeDimension Set to negative value to calculate Gaussian filtering, and positive value less than image dimensionality to calculate image derivative in that direction using Gaussian convolution.
	@param bc Boundary condition.
	*/
	inline void gaussFilter(Image<uint64_t>& img, const math::Vec3d& sigma, BoundaryCondition bc)
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
		gaussFilter(img, math::Vec3d(sigma, sigma, sigma), bc);
	}

	/**
	Gaussian filtering in-place using separable algorithm.
	@param out Image to filter.
	@param sigma Standard deviation of the Gaussian kernel.
	@param derivativeDimension Set to negative value to calculate Gaussian filtering, and positive value less than image dimensionality to calculate image derivative in that direction using Gaussian convolution.
	@param bc Boundary condition.
	*/
	inline void gaussFilter(Image<float32_t>& img, const math::Vec3d& sigma, BoundaryCondition bc)
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
	inline void gaussFilter(Image<float32_t>& img, double sigma, BoundaryCondition bc)
	{
		gaussFilter(img, math::Vec3d(sigma, sigma, sigma), bc);
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
		gaussDerivative(const Image<input_t>& in, Image<output_t>& out, const math::Vec3d& sigma, coord_t derivativeDimension1, coord_t derivativeDimension2, BoundaryCondition bc = Zero)
	{
		internals::sepgauss(in, out, sigma, derivativeDimension1, derivativeDimension2, bc);
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
		gaussDerivative(const Image<input_t>& in, Image<output_t>& out, double sigma, coord_t derivativeDimension1, coord_t derivativeDimension2, BoundaryCondition bc = Zero)
	{
		gaussDerivative(in, out, math::Vec3d(sigma, sigma, sigma), derivativeDimension1, derivativeDimension2, bc);
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
		gaussDerivative(const Image<pixel_t>& img, const math::Vec3d& sigma, coord_t derivativeDimension1, coord_t derivativeDimension2, BoundaryCondition bc = Zero)
	{
		internals::sepgauss(img, sigma, derivativeDimension1, derivativeDimension2, bc);
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
		gaussDerivative(const Image<pixel_t>& img, double sigma, coord_t derivativeDimension1, coord_t derivativeDimension2, BoundaryCondition bc = Zero)
	{
		gaussDerivative(img, math::Vec3d(sigma, sigma, sigma), derivativeDimension1, derivativeDimension2, bc);
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
	template<typename pixel_t> void highpassFilter(const Image<pixel_t>& in, Image<pixel_t>& out, const math::Vec3d& sigma, pixel_t shift = 0, bool allowOpt = true, BoundaryCondition bc = Zero)
	{
		gaussFilter(in, out, sigma, allowOpt, bc);
		if (shift == 0)
		{
			// Calculate out = in - out
			invSubtract(out, in);
		}
		else
		{
			invSubtractAdd(out, in, shift);
		}
	}


	/*
	Mean, variance, etc. filters
	*/
	// First define macro that creates two shorthand methods, first where neighbourhood radius is vector, and second where neighbourhood
	// radius is the same for all coordinate directions.
	// The comments are used to generate doxygen documentation.
	// This version includes separable filtering optimization for all data types.
#define DEFINE_FILTER_SEP_ALL(name, help) \
/** \
help \
\
In-place filtering supports only rectangular neighbourhood. \
@param in Image to filter. \
@param nbRadius Radius of filtering neighbourhood. \
@param bc Boundary condition (Zero or Nearest). \
*/ \
template<typename pixel_t> void name##Filter(Image<pixel_t>& img, const math::Vec3c& nbRadius, BoundaryCondition bc = Zero) \
{ \
	sepFilter<pixel_t, internals::name##Op<pixel_t> >(img, nbRadius, bc); \
} \
/** \
help \
\
Separable filtering is used for all pixel data types for rectangular neighbourhoods. \
@param in Input image. \
@param out Output image. \
@param nbRadius Radius of filtering neighbourhood. \
@param nbType Neighbourhood type (Rectangular or Ellipsoidal). \
@param bc Boundary condition (Zero or Nearest). \
*/ \
template<typename pixel_t, typename out_t> void name##Filter(const Image<pixel_t>& in, Image<out_t>& out, const math::Vec3c& nbRadius, NeighbourhoodType nbType = Ellipsoidal, BoundaryCondition bc = Zero) \
{ \
	if(nbType == Rectangular) \
	{ \
		out.ensureSize(in); \
		setValue<out_t, pixel_t>(out, in); \
		name##Filter<out_t>(out, nbRadius, bc); \
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
@param in Input image. \
@param out Output image. \
@param nbRadius Radius of filtering neighbourhood. \
@param nbType Neighbourhood type (Rectangular or Ellipsoidal). \
@param bc Boundary condition (Zero or Nearest). \
*/ \
template<typename pixel_t, typename out_t> void name##Filter(const Image<pixel_t>& in, Image<out_t>& out, coord_t nbRadius, NeighbourhoodType nbType = Ellipsoidal, BoundaryCondition bc = Zero) \
{ \
	name##Filter<pixel_t, out_t>(in, out, math::Vec3c(nbRadius, nbRadius, nbRadius), nbType, bc); \
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
@param nbType Neighbourhood type (Rectangular or Ellipsoidal). \
@param bc Boundary condition (Zero or Nearest). \
*/ \
template<typename pixel_t, typename out_t> void name##Filter(const Image<pixel_t>& in, Image<out_t>& out, const math::Vec3c& nbRadius, NeighbourhoodType nbType = Ellipsoidal, BoundaryCondition bc = Zero) \
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
@param nbType Neighbourhood type (Rectangular or Ellipsoidal). \
@param bc Boundary condition (Zero or Nearest). \
*/ \
template<typename pixel_t, typename out_t> void name##Filter(const Image<pixel_t>& in, Image<out_t>& out, coord_t nbRadius, NeighbourhoodType nbType = Ellipsoidal, BoundaryCondition bc = Zero) \
{ \
	name##Filter<pixel_t, out_t>(in, out, math::Vec3c(nbRadius, nbRadius, nbRadius), nbType, bc); \
} \
/** \
help \
\
In-place filtering supports only rectangular neighbourhood. \
@param img Image to process. \
@param nbRadius Radius of filtering neighbourhood. \
@param nbType Neighbourhood type (Rectangular or Ellipsoidal). \
@param bc Boundary condition (Zero or Nearest). \
*/ \
inline void name##Filter(Image<float32_t>& img, const math::Vec3c& nbRadius, BoundaryCondition bc = Zero) \
{ \
	sepFilter<float32_t, internals::name##Op<float32_t> >(img, nbRadius, bc); \
} \
/** \
help \
\
In-place filtering supports only rectangular neighbourhood. \
@param img Image to process. \
@param nbRadius Radius of filtering neighbourhood. \
@param nbType Neighbourhood type (Rectangular or Ellipsoidal). \
@param bc Boundary condition (Zero or Nearest). \
*/ \
inline void name##Filter(Image<float32_t>& img, coord_t nbRadius, BoundaryCondition bc = Zero) \
{ \
	sepFilter<float32_t, internals::name##Op<float32_t> >(img, math::Vec3c(nbRadius, nbRadius, nbRadius), bc); \
} \
/** \
help \
\
Separable filtering is used for floating point pixel data types for rectangular neighbourhoods. \
@param in Input image. \
@param out Output image. \
@param nbRadius Radius of filtering neighbourhood. \
@param nbType Neighbourhood type (Rectangular or Ellipsoidal). \
@param bc Boundary condition (Zero or Nearest). \
*/ \
template<> inline void name##Filter(const Image<float32_t>& in, Image<float32_t>& out, const math::Vec3c& nbRadius, NeighbourhoodType nbType, BoundaryCondition bc) \
{ \
	if(nbType == Rectangular) \
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
@param nbType Neighbourhood type (Rectangular or Ellipsoidal). \
@param bc Boundary condition (Zero or Nearest). \
*/ \
template<typename pixel_t, typename out_t> void name##Filter(const Image<pixel_t>& in, Image<out_t>& out, const math::Vec3c& nbRadius, NeighbourhoodType nbType = Ellipsoidal, BoundaryCondition bc = Zero) \
{ \
	filter<pixel_t, out_t, internals::name##Op<pixel_t> >(in, out, nbRadius, nbType, bc); \
} \
 \
/** \
help \
@param in Input image. \
@param out Output image. \
@param nbRadius Radius of filtering neighbourhood. \
@param nbType Neighbourhood type (Rectangular or Ellipsoidal). \
@param bc Boundary condition (Zero or Nearest). \
*/ \
template<typename pixel_t, typename out_t> void name##Filter(const Image<pixel_t>& in, Image<out_t>& out, coord_t nbRadius, NeighbourhoodType nbType = Ellipsoidal, BoundaryCondition bc = Zero) \
{ \
	name##Filter<pixel_t, out_t>(in, out, math::Vec3c(nbRadius, nbRadius, nbRadius), nbType, bc); \
}


// Version with one parameter and no separable optimization
#define DEFINE_FILTER_1PARAM(name, paramtype, help, paramhelp) \
/** \
help \
@param in Input image. \
@param out Output image. \
@param nbRadius Radius of filtering neighbourhood. \
@param parameter paramhelp \
@param nbType Neighbourhood type (Rectangular or Ellipsoidal). \
@param bc Boundary condition (Zero or Nearest). \
*/ \
template<typename pixel_t, typename out_t> void name##Filter(const Image<pixel_t>& in, Image<out_t>& out, const math::Vec3c& nbRadius, paramtype parameter, NeighbourhoodType nbType = Ellipsoidal, BoundaryCondition bc = Zero) \
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
@param nbType Neighbourhood type (Rectangular or Ellipsoidal). \
@param bc Boundary condition (Zero or Nearest). \
*/ \
template<typename pixel_t, typename out_t> void name##Filter(const Image<pixel_t>& in, Image<out_t>& out, coord_t nbRadius, paramtype parameter, NeighbourhoodType nbType = Ellipsoidal, BoundaryCondition bc = Zero) \
{ \
	name##Filter<pixel_t, out_t>(in, out, math::Vec3c(nbRadius, nbRadius, nbRadius), parameter, nbType, bc); \
}



	// Now define the filtering operations
	DEFINE_FILTER_SEP_ALL(min, Calculates minimum filtering.)
	DEFINE_FILTER_SEP_ALL(max, Calculates maximum filtering.)
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
	@param nbType Neighbourhood type (Rectangular or Ellipsoidal).
	@param bc Boundary condition (Zero or Nearest).
	*/
	template<typename pixel_t, typename out_t> void varianceFilter(const Image<pixel_t>& in, Image<out_t>& out, const math::Vec3c& nbRadius, NeighbourhoodType nbType = Ellipsoidal, BoundaryCondition bc = Zero)
	{
		filter<pixel_t, out_t, internals::varianceOp<pixel_t> >(in, out, nbRadius, nbType, bc);
	}
	
	/**
	Calculates variance filtering.

	Separable optimization is used for rectangular neighbourhoods.
	@param in Input image.
	@param out Output image.
	@param nbRadius Radius of filtering neighbourhood.
	@param nbType Neighbourhood type (Rectangular or Ellipsoidal).
	@param bc Boundary condition (Zero or Nearest).
	*/
	template<typename pixel_t, typename out_t> void varianceFilter(const Image<pixel_t>& in, Image<out_t>& out, coord_t nbRadius, NeighbourhoodType nbType = Ellipsoidal, BoundaryCondition bc = Zero)
	{
		varianceFilter<pixel_t, out_t>(in, out, math::Vec3c(nbRadius, nbRadius, nbRadius), nbType, bc);
	}

	namespace internals
	{
		inline float32_t rectVarianceFinalization(float32_t out, float32_t tmp)
		{
			return (float32_t)tmp - (float32_t)out * (float32_t)out;
		}
	}

	template<> inline void varianceFilter(const Image<float32_t>& in, Image<float32_t>& out, const math::Vec3c& nbRadius, NeighbourhoodType nbType, BoundaryCondition bc)
	{
		if(nbType == Rectangular)
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
			pointProcessImageImage<float32_t, float32_t, float32_t, internals::rectVarianceFinalization>(out, tmp);
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
	@param nbType Neighbourhood type (Rectangular or Ellipsoidal).
	@param bc Boundary condition (Zero or Nearest).
	*/
	template<typename pixel_t, typename out_t> void stddevFilter(const Image<pixel_t>& in, Image<out_t>& out, const math::Vec3c& nbRadius, NeighbourhoodType nbType = Ellipsoidal, BoundaryCondition bc = Zero)
	{
		filter<pixel_t, out_t, internals::varianceOp<pixel_t> >(in, out, nbRadius, nbType, bc);
		squareRoot(out);
	}

	/**
	Calculates standard deviation filtering.

	Separable optimization is used for rectangular neighbourhoods.
	@param in Input image.
	@param out Output image.
	@param nbRadius Radius of filtering neighbourhood.
	@param nbType Neighbourhood type (Rectangular or Ellipsoidal).
	@param bc Boundary condition (Zero or Nearest).
	*/
	template<typename pixel_t, typename out_t> void stddevFilter(const Image<pixel_t>& in, Image<out_t>& out, coord_t nbRadius, NeighbourhoodType nbType = Ellipsoidal, BoundaryCondition bc = Zero)
	{
		stddevFilter<pixel_t, out_t>(in, out, math::Vec3c(nbRadius, nbRadius, nbRadius), nbType, bc);
	}


	/**
	Opening operation.
	@param img Image that is to be processed.
	@param tmp Temporary image having the same size as the input image. If the size is wrong, it is changed by this method.
	@param nbType Neighbourhood type.
	@param nbRadius Neighbourhood radius.
	*/
	template<typename pixel_t> void openingFilter(Image<pixel_t>& img, Image<pixel_t>& tmp, const math::Vec3c& nbRadius, NeighbourhoodType nbType = Ellipsoidal, BoundaryCondition bc = Zero)
	{
		tmp.ensureSize(img);
		minFilter<pixel_t, pixel_t>(img, tmp, nbRadius, nbType, bc);
		maxFilter<pixel_t, pixel_t>(tmp, img, nbRadius, nbType, bc);
	}

	/**
	Opening operation.
	@param img Image that is to be processed.
	@param tmp Temporary image having the same size as the input image. If the size is wrong, it is changed by this method.
	@param nbType Neighbourhood type.
	@param nbRadius Neighbourhood radius.
	*/
	template<typename pixel_t> void openingFilter(Image<pixel_t>& img, Image<pixel_t>& tmp, coord_t nbRadius, NeighbourhoodType nbType = Ellipsoidal, BoundaryCondition bc = Zero)
	{
		openingFilter(img, tmp, math::Vec3c(nbRadius, nbRadius, nbRadius), nbType, bc);
	}

	/**
	Closing operation.
	@param img Image that is to be processed.
	@param tmp Temporary image having the same size as the input image. If the size is wrong, it is changed by this method.
	@param nbType Neighbourhood type.
	@param nbRadius Neighbourhood radius.
	*/
	template<typename pixel_t> void closingFilter(Image<pixel_t>& img, Image<pixel_t>& tmp, const math::Vec3c& nbRadius, NeighbourhoodType nbType = Ellipsoidal, BoundaryCondition bc = Zero)
	{
		tmp.ensureSize(img);
		maxFilter<pixel_t, pixel_t>(img, tmp, nbRadius, nbType, bc);
		minFilter<pixel_t, pixel_t>(tmp, img, nbRadius, nbType, bc);
	}

	/**
	Closing operation.
	@param img Image that is to be processed.
	@param tmp Temporary image having the same size as the input image. If the size is wrong, it is changed by this method.
	@param nbType Neighbourhood type.
	@param nbRadius Neighbourhood radius.
	*/
	template<typename pixel_t> void closingFilter(Image<pixel_t>& img, Image<pixel_t>& tmp, coord_t nbRadius, NeighbourhoodType nbType = Ellipsoidal, BoundaryCondition bc = Zero)
	{
		closingFilter(img, tmp, math::Vec3c(nbRadius, nbRadius, nbRadius), nbType, bc);
	}

	/**
	Calculates bilateral filtering.
	@param in Input image.
	@param out Output image.
	@param spatialSigma Standard deviation of gaussian kernel used for spatial smoothing.
	@param rangeSigma Standard deviation of gaussian kernel used to avoid smoothing edges of features. Order of magnitude must be similar to difference between gray levels of background and objects.
	@param bc Boundary condition (Zero or Nearest).
	*/
	template<typename pixel_t, typename out_t> void bilateralFilter(const Image<pixel_t>& in, Image<out_t>& out, double spatialSigma, double rangeSigma, BoundaryCondition bc = Zero)
	{
		coord_t nbRadius = (coord_t)ceil(3 * spatialSigma);
		filter<pixel_t, out_t, Vec2d, internals::bilateralOp<pixel_t> >(in, out, math::Vec3c(nbRadius, nbRadius, nbRadius), Vec2d(spatialSigma, rangeSigma), Rectangular, bc);
	}

	
	namespace tests
	{
		void filters();
		void separableOptimization();
		void gaussFilters();
		void bilateral();
	}
}
