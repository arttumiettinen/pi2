#pragma once

#include <omp.h>
#include <vector>
#include "image.h"
#include "math/vec3.h"
#include "interpolation.h"
#include "utilities.h"
#include "math/numberutils.h"
#include "math/vectoroperations.h"

using math::Vec3d;
using math::Vec3c;
using math::Vec3f;
using math::clamp;
using math::pixelRound;
using math::NumberUtils;

namespace itl2
{
	/**
	Shifts input image and stores the result in output image.
	@param in Image to be shifted.
	@param out Output image.
	@param shift Shift vector.
	@param interpolate Interpolation type.
	*/
	template<typename pixel_t, typename out_t> void translate(const Image<pixel_t>& in, Image<out_t>& out, const Vec3d& shift, const Interpolator<out_t, pixel_t>& interpolate = LinearInterpolator<out_t, pixel_t>(BoundaryCondition::Zero))
	{
		out.mustNotBe(in);
		out.ensureSize(in);

		#pragma omp parallel for if(out.pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		for (coord_t z = 0; z < out.depth(); z++)
		{
			for (coord_t y = 0; y < out.height(); y++)
			{
				for (coord_t x = 0; x < out.width(); x++)
				{
					typename NumberUtils<out_t>::RealFloatType xs = (typename NumberUtils<out_t>::RealFloatType)(x - shift.x);
					typename NumberUtils<out_t>::RealFloatType ys = (typename NumberUtils<out_t>::RealFloatType)(y - shift.y);
					typename NumberUtils<out_t>::RealFloatType zs = (typename NumberUtils<out_t>::RealFloatType)(z - shift.z);
					out(x, y, z) = interpolate(in, xs, ys, zs);
				}
			}
		}

	}

	/**
	Crops input image to size of output image. Left-top of output image is placed at the given position.
	*/
	template<typename pixel_t, typename out_t> void crop(const Image<pixel_t>& in, Image<out_t>& out, const Vec3c& outPos)
	{
		out.mustNotBe(in);

		#pragma omp parallel for if(out.pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		for (coord_t z = 0; z < out.depth(); z++)
		{
			for (coord_t y = 0; y < out.height(); y++)
			{
				for (coord_t x = 0; x < out.width(); x++)
				{
					Vec3c xi = Vec3c(x, y, z) + outPos;
					if(in.isInImage(xi))
						out(x, y, z) = in(xi);
				}
			}
		}
	}

	/**
	Converts input image to smaller scale by averaging binSize^dimensionality blocks.
	Does not use separable algorithm (that saves some memory).
	@param in Input image.
	@param out Output image. The image is automatically initialized to correct size.
	@param binSize Bin size. 2 makes the output image dimensions half of the input image dimensions, 3 makes them one third etc.
	@param indicateProgress Set to true to show a progress bar.
	*/
	template<typename pixel_t, typename out_t> void binning(const Image<pixel_t>& in, Image<out_t>& out, const Vec3c& binSize, bool indicateProgress = true)
	{
		if (binSize.min() <= 0)
			throw ITLException("Bins size must be positive.");

		out.mustNotBe(in);
		out.ensureSize(in.dimensions().componentwiseDivide(binSize));

		size_t counter = 0;
		#pragma omp parallel for if(out.pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		for (coord_t z = 0; z < out.depth(); z++)
		{
			coord_t inz = z * binSize.z;
			coord_t iny = 0;
			for (coord_t y = 0; y < out.height(); y++, iny += binSize.y)
			{
				coord_t inx = 0;
				for (coord_t x = 0; x < out.width(); x++, inx += binSize.x)
				{

					coord_t inzEnd = math::min(inz + binSize.x, in.depth());
					coord_t inyEnd = math::min(iny + binSize.y, in.height());
					coord_t inxEnd = math::min(inx + binSize.z, in.width());
					typename NumberUtils<out_t>::FloatType M = 0;
					for (coord_t zz = inz; zz < inzEnd; zz++)
					{
						for (coord_t yy = iny; yy < inyEnd; yy++)
						{
							for (coord_t xx = inx; xx < inxEnd; xx++)
							{
								M += in(xx, yy, zz);
							}
						}
					}

					M /= (typename NumberUtils<out_t>::RealFloatType)((inzEnd - inz) * (inyEnd - iny) * (inxEnd - inx));

					out(x, y, z) = pixelRound<out_t>(M);
				}
			}

			showThreadProgress(counter, out.depth(), indicateProgress);
		}
	}

	/**
	Converts input image to smaller scale by averaging binSize^dimensionality blocks.
	Does not use separable algorithm (that saves some memory).
	@param in Input image.
	@param out Output image. The image is automatically initialized to correct size.
	@param binSize Bin size. 2 makes the output image dimensions half of the input image dimensions, 3 makes them one third etc.
	@param indicateProgress Set to true to show a progress bar.
	*/
	template<typename pixel_t, typename out_t> void binning(const Image<pixel_t>& in, Image<out_t>& out, size_t binSize, bool indicateProgress = true)
	{
		binning(in, out, Vec3c(binSize, binSize, binSize), indicateProgress);
	}

	/**
	Converts input image to smaller scale by averaging binSize^dimensionality blocks.
	Does not average value specified as argument.
	@param in Input image.
	@param out Output image. The image is automatically initialized to correct size.
	@param amount Bin size. 2 makes the output image dimensions half of the input image dimensions, 3 makes them one third etc.
	@param badValue Value that should not be considered in the averaging calculations.
	@param undefinedValue Value that is placed to those pixels of the output image that do not correspond to any valid pixels in the input image.
	*/
	template<typename pixel_t, typename out_t> void maskedBinning(const Image<pixel_t>& in, Image<out_t>& out, size_t amount, pixel_t badValue, out_t undefinedValue, bool indicateProgress = true)
	{
		out.mustNotBe(in);
		coord_t binSize = (coord_t)amount;
		out.ensureSize(round(Vec3d(in.dimensions()) / (double)binSize));

		size_t counter = 0;
#pragma omp parallel for if(out.pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		for (coord_t z = 0; z < out.depth(); z++)
		{
			coord_t inz = z * binSize;
			coord_t iny = 0;
			for (coord_t y = 0; y < out.height(); y++, iny += binSize)
			{
				coord_t inx = 0;
				for (coord_t x = 0; x < out.width(); x++, inx += binSize)
				{

					coord_t inzEnd = math::min(inz + binSize, in.depth());
					coord_t inyEnd = math::min(iny + binSize, in.height());
					coord_t inxEnd = math::min(inx + binSize, in.width());
					typename NumberUtils<out_t>::FloatType M = 0;
					typename NumberUtils<out_t>::RealFloatType count = 0;
					for (coord_t zz = inz; zz < inzEnd; zz++)
					{
						for (coord_t yy = iny; yy < inyEnd; yy++)
						{
							for (coord_t xx = inx; xx < inxEnd; xx++)
							{
								pixel_t pix = in(xx, yy, zz);
								if (pix != badValue)
								{
									M += pix;
									count++;
								}
							}
						}
					}

					if (count > 0)
					{
						M /= count;
						out(x, y, z) = pixelRound<out_t>(M);
					}
					else
					{
						out(x, y, z) = undefinedValue;
					}

				}
			}

			showThreadProgress(counter, out.depth(), indicateProgress);
		}
	}

	/**
	Scales input image to the size of the output image and replaces output image by the scaled image.
	Does not try to suppress aliasing artifacts when downscaling.
	@param in Input image.
	@param out Output image.
	@param interpolate Interpolation type.
	@param indicateProgress Set to true to show a progress bar.
	@param factor Scaling factor. If zero or negative, determined from dimensions of input and output image.
	@param delta Shift that is added to coordinates of each input point. Used in distributed processing.
	*/
	template<typename in_t, typename out_t> void scale(const Image<in_t>& in, Image<out_t>& out, const Interpolator<out_t, in_t>& interpolate = LinearInterpolator<out_t, in_t>(BoundaryCondition::Zero), bool indicateProgress = true, Vec3d factor = Vec3d(0, 0, 0), const Vec3d& delta = Vec3d(0, 0, 0))
	{

		if (factor.max() <= 0)
		{
			factor = Vec3d(out.dimensions()).componentwiseDivide(Vec3d(in.dimensions()));
		}

		size_t counter = 0;
		#pragma omp parallel for if(out.pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		for (coord_t z = 0; z < out.depth(); z++)
		{
			typename NumberUtils<out_t>::RealFloatType sz = (typename NumberUtils<out_t>::RealFloatType)(z / factor.z) + (typename NumberUtils<out_t>::RealFloatType)delta.z;
			for (coord_t y = 0; y < out.height(); y++)
			{
				typename NumberUtils<out_t>::RealFloatType sy = (typename NumberUtils<out_t>::RealFloatType)(y / factor.y) + (typename NumberUtils<out_t>::RealFloatType)delta.y;
				for (coord_t x = 0; x < out.width(); x++)
				{
					typename NumberUtils<out_t>::RealFloatType sx = (typename NumberUtils<out_t>::RealFloatType)(x / factor.x) + (typename NumberUtils<out_t>::RealFloatType)delta.x;
					out(x, y, z) = interpolate(in, sx, sy, sz);
				}
			}

			showThreadProgress(counter, out.depth(), indicateProgress);
		}
	}

	namespace internals
	{
		/**
		Inverse distance interpolation.
		@param p Smoothing exponent.
		*/
		inline Vec3f inverseDistanceInterpolate(const vector<Vec3f>& refPoints, const vector<Vec3f>& values, const Vec3f& x, float p = 2.5)
		{
			if (refPoints.size() != values.size())
				throw ITLException("refPoints and values lists must have the same size.");

			if (refPoints.size() <= 0)
				return x;

			Vec3f sum = Vec3f(0, 0, 0);
			float wsum = 0;
			for (size_t n = 0; n < refPoints.size(); n++)
			{
				float dist = (x - refPoints[n]).norm();

				if (dist < 1e-7)
				{
					return values[n];
				}

				float w = 1 / ::pow(dist, p);

				sum += w * values[n];
				wsum += w;
			}

			return sum / wsum;
		}
	}

	/**
	Transforms img using a free-form point-to-point transformation.
	@param img Input image.
	@param out Output image. Size of this image must be set to the size of desired output.
	@param outPos Position of the output image relative to the origin of the input image.
	@param refPoints, defPoints List of point pairs. Element of refPoints defines a position in the input image and the corresponding element in defPoints gives the corresponding deformed position.
	*/
	template<typename pixel_t, typename out_t> void genericTransform(const Image<pixel_t>& img,
		Image<out_t>& out,
		const Vec3c& outPos,
		const std::vector<Vec3f>& refPoints,
		const std::vector<Vec3f>& defPoints,
		float exponent = 2.5f,
		const Interpolator<out_t, pixel_t>& interpolate = LinearInterpolator<out_t, pixel_t>(BoundaryCondition::Zero))
	{
		out.mustNotBe(img);

		std::vector<Vec3f> shifts = defPoints - refPoints;

		size_t counter = 0;
		#pragma omp parallel for if (!omp_in_parallel())
		for (coord_t z = 0; z < out.depth(); z++)
		{
			for (coord_t y = 0; y < out.height(); y++)
			{
				for (coord_t x = 0; x < out.width(); x++)
				{
					Vec3c oix(x, y, z);

					Vec3f ix = Vec3f(oix + outPos);
					Vec3f transformed = ix + internals::inverseDistanceInterpolate(refPoints, shifts, ix, exponent);

					out(oix) = interpolate(img, transformed.x, transformed.y, transformed.z);
				}
			}

			showThreadProgress(counter, out.pixelCount());
		}

	}

	namespace tests
	{
		void translate();
		void binning();
		void genericTransform();
	}

}
