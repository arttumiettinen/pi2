#pragma once

#include "image.h"
#include "math/mathutils.h"
#include "math/vec2.h"
#include "math/vec4.h"
#include "math/numberutils.h"

namespace itl2
{

	/**
	Process img in place.
	*/
	template<typename pixel_t, typename intermediate_t, intermediate_t process(pixel_t)> void pointProcess(Image<pixel_t>& img)
	{
		#pragma omp parallel for if(img.pixelCount() > PARALLELIZATION_THRESHOLD)
		for (coord_t n = 0; n < img.pixelCount(); n++)
		{
			img(n) = math::pixelRound<pixel_t, intermediate_t>(process(img(n)));

			// Showing progress info here would induce more processing than is done in the whole loop.
		}
	}


	/**
	Process corresponding pixels from l and r, place result to l.
	*/
	template<typename pixel1_t, typename pixel2_t, typename intermediate_t, intermediate_t process(pixel1_t, pixel2_t)> void pointProcessImageImage(Image<pixel1_t>& l, const Image<pixel2_t>& r)
	{
		l.checkSize(r);
		
		#pragma omp parallel for if(l.pixelCount() > PARALLELIZATION_THRESHOLD)
		for (coord_t n = 0; n < l.pixelCount(); n++)
		{
			l(n) = math::pixelRound<pixel1_t, intermediate_t>(process(l(n), r(n)));

			// Showing progress info here would induce more processing than is done in the whole loop.
		}
	}

	/**
	Process corresponding pixels from l and r, place result to l.
	*/
	template<typename pixel1_t, typename pixel2_t, typename param_t, typename intermediate_t, intermediate_t process(pixel1_t, pixel2_t, param_t)> void pointProcessImageImageParam(Image<pixel1_t>& l, const Image<pixel2_t>& r, param_t c)
	{
		l.checkSize(r);

		#pragma omp parallel for if(l.pixelCount() > PARALLELIZATION_THRESHOLD)
		for (coord_t n = 0; n < l.pixelCount(); n++)
		{
			l(n) = math::pixelRound<pixel1_t, intermediate_t>(process(l(n), r(n), c));

			// Showing progress info here would induce more processing than is done in the whole loop.
		}
	}

	/**
	Process pixel from img with parameter, place result to img.
	*/
	template<typename pixel_t, typename param_t, typename intermediate_t, intermediate_t process(pixel_t, param_t)> void pointProcessImageParam(Image<pixel_t>& img, param_t param)
	{
		#pragma omp parallel for if(img.pixelCount() > PARALLELIZATION_THRESHOLD)
		for (coord_t n = 0; n < img.pixelCount(); n++)
		{
			img(n) = math::pixelRound<pixel_t, intermediate_t>(process(img(n), param));

			// Showing progress info here would induce more processing than is done in the whole loop.
		}
	}

	/**
	Process pixel from img with parameter if pixel from img does not equal badValue, place result to img.
	*/
	template<typename pixel_t, typename param_t, typename intermediate_t, intermediate_t process(pixel_t, param_t)> void maskedPointProcessImageParam(Image<pixel_t>& img, param_t param, pixel_t badValue)
	{
		#pragma omp parallel for if(img.pixelCount() > PARALLELIZATION_THRESHOLD)
		for (coord_t n = 0; n < img.pixelCount(); n++)
		{
			pixel_t pix = img(n);
			if(pix != badValue)
				img(n) = math::pixelRound<pixel_t, intermediate_t>(process(pix, param));

			// Showing progress info here would induce more processing than is done in the whole loop.
		}
	}





	namespace internals
	{
		
		// Unary operators

		template<typename pixel_t, typename intermediate_t> intermediate_t negateOp(pixel_t val)
		{
			return -(intermediate_t)val;
		}

		template<typename pixel_t, typename intermediate_t> intermediate_t exponentiateOp(pixel_t val)
		{
			return exp((intermediate_t)val);
		}

		template<typename pixel_t, typename intermediate_t> intermediate_t squareOp(pixel_t val)
		{
			return (intermediate_t)val * (intermediate_t)val;
		}

		template<typename pixel_t, typename intermediate_t> intermediate_t squareRootOp(pixel_t val)
		{
			return sqrt((intermediate_t)val);
		}

		template<typename pixel_t, typename intermediate_t> intermediate_t absOp(pixel_t val)
		{
			return std::abs((intermediate_t)val);
		}

		template<typename pixel_t, typename intermediate_t> intermediate_t logOp(pixel_t val)
		{
			return std::log((intermediate_t)val);
		}

		template<typename pixel_t, typename intermediate_t> intermediate_t log10Op(pixel_t val)
		{
			return std::log10((intermediate_t)val);
		}

		template<typename pixel_t, typename intermediate_t> intermediate_t sinOp(pixel_t val)
		{
			return std::sin((intermediate_t)val);
		}

		template<typename pixel_t, typename intermediate_t> intermediate_t cosOp(pixel_t val)
		{
			return std::cos((intermediate_t)val);
		}

		template<typename pixel_t, typename intermediate_t> intermediate_t tanOp(pixel_t val)
		{
			return std::tan((intermediate_t)val);
		}

		template<typename pixel_t, typename intermediate_t> intermediate_t invOp(pixel_t val)
		{
			return (intermediate_t)1.0 / (intermediate_t)val;
		}

		template<typename pixel_t, typename intermediate_t> intermediate_t roundOp(pixel_t val)
		{
			return round((intermediate_t)val);
		}

		template<typename pixel_t, typename intermediate_t> intermediate_t ceilOp(pixel_t val)
		{
			return ceil((intermediate_t)val);
		}

		template<typename pixel_t, typename intermediate_t> intermediate_t floorOp(pixel_t val)
		{
			return floor((intermediate_t)val);
		}


		// Unary operators for complex data type

		template<typename pixel_t, typename intermediate_t> intermediate_t conjugateOp(pixel_t val)
		{
			return conj(val);
		}

		template<typename pixel_t, typename intermediate_t> intermediate_t normalizeOp(pixel_t val)
		{
			// NOTE: stl norm for complex type returns squared norm.
			intermediate_t L = norm(val);
			if (!math::NumberUtils<intermediate_t>::equals(L, 0))
				return val / sqrt(L);
			else
				return val;
		}

		template<typename pixel_t, typename intermediate_t> intermediate_t realOp(pixel_t val)
		{
			return real(val);
		}

		template<typename pixel_t, typename intermediate_t> intermediate_t imagOp(pixel_t val)
		{
			return imag(val);
		}

		template<typename pixel_t, typename intermediate_t> intermediate_t argOp(pixel_t val)
		{
			return arg(val);
		}

		template<typename pixel_t, typename intermediate_t> intermediate_t normSquaredOp(pixel_t val)
		{
			// NOTE: stl norm for complex type returns squared norm.
			return norm(val);
		}


		// Binary operators

		template<typename pixel_t, typename exp_t, typename intermediate_t> intermediate_t powOp(pixel_t val, exp_t exponent)
		{
			return (intermediate_t)std::pow(val, exponent);
		}

		template<typename t1, typename t2, typename intermediate_t> intermediate_t addOp(t1 a, t2 b)
		{
			return (intermediate_t)a + (intermediate_t)b;
		}

		template<typename t1, typename t2, typename intermediate_t> intermediate_t maxOp(t1 a, t2 b)
		{
			return std::max((intermediate_t)a, (intermediate_t)b);
		}

		template<typename t1, typename t2, typename intermediate_t> intermediate_t minOp(t1 a, t2 b)
		{
			return std::min((intermediate_t)a, (intermediate_t)b);
		}

		/**
		a-b
		*/
		template<typename t1, typename t2, typename intermediate_t> intermediate_t subtractOp(t1 a, t2 b)
		{
			return (intermediate_t)a - (intermediate_t)b;
		}

		/**
		b-a
		*/
		template<typename t1, typename t2, typename intermediate_t> intermediate_t invSubtractOp(t1 a, t2 b)
		{
			return (intermediate_t)b - (intermediate_t)a;
		}

		/**
		b-a+c
		*/
		template<typename t1, typename t2, typename t3, typename intermediate_t> intermediate_t invSubtractAddOp(t1 a, t2 b, t3 c)
		{
			return (intermediate_t)b - (intermediate_t)a + (intermediate_t)c;
		}

		template<typename t1, typename t2, typename intermediate_t> intermediate_t multiplyOp(t1 a, t2 b)
		{
			return (intermediate_t)a * (intermediate_t)b;
		}

		template<typename t1, typename t2, typename intermediate_t> intermediate_t divideOp(t1 a, t2 b)
		{
			return (intermediate_t)a / (intermediate_t)b;
		}

		template<typename pixel_t, typename out_t> out_t setValueOp(pixel_t val, out_t param)
		{
			return param;
		}

		template<typename pixel_t, typename pixel2_t> pixel2_t copyOp(pixel_t l, pixel2_t r)
		{
			return r;
		}


		/**
		Thresholds a with b.
		Returns 1 if a > b and 0 otherwise.
		*/
		template<typename pixel_t, typename th_t> pixel_t thresholdOp(pixel_t a, th_t b)
		{
			if (a > math::pixelRound<pixel_t>(b))
				return (pixel_t)1;
			else
				return (pixel_t)0;
		}

		/**
		Thresholds a with range (minmax.x, minmax.y].
		Returns 1 if minmax.x < a <= minmax.y, and 0 otherwise.
		*/
		template<typename pixel_t> pixel_t thresholdRangeOp(pixel_t a, const math::Vec2d minmax)
		{
			if (minmax.x < a && a <= minmax.y)
				return (pixel_t)1;
			else
				return (pixel_t)0;
		}

		/**
		Performs periodic thresholding.
		@param x The value to threshold.
		@param inputs 4-component vector containing (period start, period end, threshold range start, threshold range end).
		@return 1 if x is in (threshold range start, threshold range end] modulo period; 0 otherwise. If threshold range is larger than period,
		returns always 1. If threshold range start is larger than or equal to threshold range end, returns always 0.
		*/
		template<typename pixel_t> pixel_t thresholdPeriodicOp(pixel_t x, const math::Vec4d& inputs)
		{
			double periodStart = inputs[0];
			double periodEnd = inputs[1];
			double thresholdStart = inputs[2];
			double thresholdEnd = inputs[3];

			double periodLength = periodEnd - periodStart;
			double periodicX = math::realmod(x - periodStart, periodLength) + periodStart;

			if (thresholdStart >= thresholdEnd)
			{
				return (pixel_t)0;
			}
			else if (thresholdEnd - thresholdStart < periodLength)
			{
				double periodicStart = math::realmod(thresholdStart - periodStart, periodLength) + periodStart;
				double periodicEnd = math::realmod(thresholdEnd - periodStart, periodLength) + periodStart;

				if (math::NumberUtils<double>::lessThan(periodicStart, periodicEnd))
					return periodicStart < periodicX && periodicX <= periodicEnd ? (pixel_t)1 : (pixel_t)0;
				else
					return periodicX <= periodicEnd || periodicStart < periodicX ? (pixel_t)1 : (pixel_t)0;
			}
			else
			{
				return (pixel_t)1;
			}
		}

		
		/**
		Map values linearly.
		Maps position a in range [bounds.x, bounds.y] to position in range [bounds.z, bounds.w].
		i.e. if a = 0.5 and [bounds.x, bounds.y] = [0, 1] and [bounds.z, bounds.w] = [3, 4],
		return value will be 3 + (0.5 - 0) / (1 - 0) * (4 - 3) = 3.5.
		*/
		template<typename pixel_t, typename intermediate_t = typename math::NumberUtils<pixel_t>::FloatType> intermediate_t linearMapOp(pixel_t a, const math::Vec4d& bounds)
		{
			intermediate_t min = (intermediate_t)bounds.x;
			intermediate_t max = (intermediate_t)bounds.y;
			intermediate_t newmin = (intermediate_t)bounds.z;
			intermediate_t newmax = (intermediate_t)bounds.w;

			if (a < min)
				a = (pixel_t)min;
			else if (a > max)
				a = (pixel_t)max;

			return newmin + (a - min) / (max - min) * (newmax - newmin);
		}

	}


#define DEF_OPERATION_ALL(name, help) \
	/** \
	help \
	\
	@param in Image to process. \
	*/ \
	template<typename pixel_t> void name (Image<pixel_t>& img) \
	{ \
		pointProcess<pixel_t, typename math::NumberUtils<pixel_t>::FloatType, internals::name##Op<pixel_t, typename math::NumberUtils<pixel_t>::FloatType> >(img); \
	}

#define DEF_OPERATION_COMPLEX(name, help) \
	/** \
	help \
	\
	@param in Image to process. \
	*/ \
	inline void name (Image<complex32_t>& img) \
	{ \
		pointProcess<complex32_t, complex32_t, internals::name##Op<complex32_t, complex32_t> >(img); \
	}

	DEF_OPERATION_ALL(negate, "Negates pixel values.")
	DEF_OPERATION_ALL(exponentiate, "Exponentates pixel values.")
	DEF_OPERATION_ALL(square, "Calculates square of pixel values.")
	DEF_OPERATION_ALL(squareRoot, "Calculates square root of pixel values.")
	DEF_OPERATION_ALL(abs, "Calculates absolute value of pixel values.")
	DEF_OPERATION_ALL(log, "Calculates natural logarithm of pixel values.")
	DEF_OPERATION_ALL(log10, "Calculates base-10 logarithm of pixel values.")
	DEF_OPERATION_ALL(sin, "Calculates sine of pixel values.")
	DEF_OPERATION_ALL(cos, "Calculates cosine of pixel values.")
	DEF_OPERATION_ALL(tan, "Calculates tangent of pixel values.")
	DEF_OPERATION_ALL(inv, "Calculates inverse (1/x) of pixel values.")

	DEF_OPERATION_ALL(round, "Rounds pixel values.")
	DEF_OPERATION_ALL(ceil, "Calculates ceiling of pixel values.")
	DEF_OPERATION_ALL(floor, "Calculates floor of pixel values.")

	DEF_OPERATION_COMPLEX(conjugate, "Calculates complex conjugate of pixel values.")
	DEF_OPERATION_COMPLEX(normalize, "Makes norm of each pixel one.")
	DEF_OPERATION_COMPLEX(real, "Calculates real part of pixel values.")
	DEF_OPERATION_COMPLEX(imag, "Calculates imaginary part of pixel values.")
	DEF_OPERATION_COMPLEX(arg, "Calculates argument of pixel values.")
	DEF_OPERATION_COMPLEX(normSquared, "Calculates squared norm of pixel values.")

	// These make sure that we have normal non-image functions available, too.
	using std::floor;
	using std::ceil;
	using std::round;
	


	/**
	Sets pixel values of first image to those copied from the second image.
	@param l First image.
	@param r Second image.
	*/
	template<typename pixel_t, typename pixel2_t> void setValue(Image<pixel_t>& l, const Image<pixel2_t>& r)
	{
		l.ensureSize(r);
		pointProcessImageImage<pixel_t, pixel2_t, pixel2_t, internals::copyOp<pixel_t, pixel2_t> >(l, r);
	}

	/**
	Sets pixel values of image to constant.
	@param img Image.
	@param value Constant value.
	*/
	template<typename pixel_t, typename param_t> void setValue(Image<pixel_t>& img, param_t value)
	{
		pointProcessImageParam<pixel_t, param_t, param_t, internals::setValueOp<pixel_t, param_t> >(img, value);
	}


	/**
	Adds two images and places the result to the first image.
	@param l First image.
	@param r Second image.
	*/
	template<typename pixel1_t, typename pixel2_t> void add(Image<pixel1_t>& l, const Image<pixel2_t>& r)
	{
		pointProcessImageImage<pixel1_t, pixel2_t, typename math::NumberUtils<pixel1_t>::FloatType, internals::addOp<pixel1_t, pixel2_t, typename math::NumberUtils<pixel1_t>::FloatType> >(l, r);
	}

	/**
	Adds constant to pixel values.
	@param l Image.
	@param r Constant value.
	*/
	template<typename pixel_t, typename param_t> void add(Image<pixel_t>& l, param_t r)
	{
		pointProcessImageParam<pixel_t, param_t, typename math::NumberUtils<pixel_t>::FloatType, internals::addOp<pixel_t, param_t, typename math::NumberUtils<pixel_t>::FloatType> >(l, r);
	}

	/**
	Adds constant to all pixel values except those that correspond to a specific indicator value.
	@param l Image.
	@param r Constant value.
	@param badValue Value indicating pixels in the image that should not be processed.
	*/
	template<typename pixel_t, typename param_t> void maskedAdd(Image<pixel_t>& l, param_t r, pixel_t badValue)
	{
		maskedPointProcessImageParam<pixel_t, param_t, typename math::NumberUtils<pixel_t>::FloatType, internals::addOp<pixel_t, param_t, typename math::NumberUtils<pixel_t>::FloatType> >(l, r, badValue);
	}


	/**
	Subtracts second image from the first image and places the result to the first image (first image = first image - second image).
	@param l First image.
	@param r Second image.
	*/
	template<typename pixel_t, typename pixel2_t> void subtract(Image<pixel_t>& l, const Image<pixel2_t>& r)
	{
		pointProcessImageImage<pixel_t, pixel2_t, typename math::NumberUtils<pixel_t>::FloatType, internals::subtractOp<pixel_t, pixel2_t, typename math::NumberUtils<pixel_t>::FloatType> >(l, r);
	}

	/**
	Subtracts constant from pixel values.
	@param l Image.
	@param r Constant value.
	*/
	template<typename pixel_t, typename param_t> void subtract(Image<pixel_t>& l, param_t r)
	{
		pointProcessImageParam<pixel_t, param_t, typename math::NumberUtils<pixel_t>::FloatType, internals::subtractOp<pixel_t, param_t, typename math::NumberUtils<pixel_t>::FloatType> >(l, r);
	}

	/**
	Subtracts first image from the second image and places the result to the first image (first image = second image - first image).
	@param l First image.
	@param r Second image.
	*/
	template<typename pixel_t, typename pixel2_t> void invSubtract(Image<pixel_t>& l, const Image<pixel2_t>& r)
	{
		pointProcessImageImage<pixel_t, pixel2_t, typename math::NumberUtils<pixel_t>::FloatType, internals::invSubtractOp<pixel_t, pixel2_t, typename math::NumberUtils<pixel_t>::FloatType> >(l, r);
	}


	/**
	Subtracts first image from the second image, adds constant, and places the result to the first image (first image = second image - first image + constant).
	@param l First image.
	@param r Second image.
	@param constant Constant to add.
	*/
	template<typename pixel_t> void invSubtractAdd(Image<pixel_t>& l, const Image<pixel_t>& r, pixel_t constant)
	{
		pointProcessImageImageParam<pixel_t, pixel_t, pixel_t, typename math::NumberUtils<pixel_t>::FloatType, internals::invSubtractAddOp<pixel_t, pixel_t, pixel_t, typename math::NumberUtils<pixel_t>::FloatType> >(l, r, constant);
	}

	/**
	Multiplies two images and places the result to the first image.
	@param l First image.
	@param r Second image.
	*/
	template<typename pixel1_t, typename pixel2_t> void multiply(Image<pixel1_t>& l, const Image<pixel2_t>& r)
	{
		pointProcessImageImage<pixel1_t, pixel2_t, typename math::NumberUtils<pixel1_t>::FloatType, internals::multiplyOp<pixel1_t, pixel2_t, typename math::NumberUtils<pixel1_t>::FloatType> >(l, r);
	}

	/**
	Multiplies pixel values by a constant.
	@param l Image.
	@param r Constant value.
	*/
	template<typename pixel_t, typename param_t> void multiply(Image<pixel_t>& l, param_t r)
	{
		pointProcessImageParam<pixel_t, param_t, typename math::NumberUtils<pixel_t>::FloatType, internals::multiplyOp<pixel_t, param_t, typename math::NumberUtils<pixel_t>::FloatType> >(l, r);
	}


	/**
	Raises pixels of first image to power given by pixel value of the second image, and places the result to the first image.
	@param l First image.
	@param r Second image.
	*/
	template<typename pixel_t, typename pixel2_t> void pow(Image<pixel_t>& l, const Image<pixel2_t>& r)
	{
		pointProcessImageImage<pixel_t, pixel2_t, typename math::NumberUtils<pixel_t>::FloatType, internals::powOp<pixel_t, pixel2_t, typename math::NumberUtils<pixel_t>::FloatType> >(l, r);
	}

	/**
	Raises pixel values to power.
	@param l Image.
	@param r Exponent.
	*/
	template<typename pixel_t, typename param_t> void pow(Image<pixel_t>& l, param_t r)
	{
		pointProcessImageParam<pixel_t, param_t, typename math::NumberUtils<pixel_t>::FloatType, internals::powOp<pixel_t, param_t, typename math::NumberUtils<pixel_t>::FloatType> >(l, r);
	}


	/**
	Divides two images and places the result to the first image (first image = first image / second image).
	@param l First image.
	@param r Second image.
	*/
	template<typename pixel1_t, typename pixel2_t> void divide(Image<pixel1_t>& l, const Image<pixel2_t>& r)
	{
		pointProcessImageImage<pixel1_t, pixel2_t, typename math::NumberUtils<pixel1_t>::FloatType, internals::divideOp<pixel1_t, pixel2_t, typename math::NumberUtils<pixel1_t>::FloatType> >(l, r);
	}

	/**
	Divides pixel values by a constant.
	@param l Image.
	@param r Constant value.
	*/
	template<typename pixel_t,  typename param_t> void divide(Image<pixel_t>& l, param_t r)
	{
		pointProcessImageParam<pixel_t, param_t, typename math::NumberUtils<pixel_t>::FloatType, internals::divideOp<pixel_t, param_t, typename math::NumberUtils<pixel_t>::FloatType> >(l, r);
	}



	/**
	Calculates maximum of two images and places the result to the first image.
	@param l First image.
	@param r Second image.
	*/
	template<typename pixel1_t, typename pixel2_t> void max(Image<pixel1_t>& l, const Image<pixel2_t>& r)
	{
		pointProcessImageImage<pixel1_t, pixel2_t, pixel1_t, internals::maxOp<pixel1_t, pixel2_t, pixel1_t> >(l, r);
	}

	/**
	Calculates maximum of an image and a constant.
	@param l Image.
	@param r Constant value.
	*/
	template<typename pixel_t, typename param_t> void max(Image<pixel_t>& l, param_t r)
	{
		pointProcessImageParam<pixel_t, param_t, pixel_t, internals::maxOp<pixel_t, param_t, pixel_t> >(l, r);
	}

	/**
	Calculates minimum of two images and places the result to the first image.
	@param l First image.
	@param r Second image.
	*/
	template<typename pixel1_t, typename pixel2_t> void min(Image<pixel1_t>& l, const Image<pixel2_t>& r)
	{
		pointProcessImageImage<pixel1_t, pixel2_t, pixel1_t, internals::minOp<pixel1_t, pixel2_t, pixel1_t> >(l, r);
	}

	/**
	Calculates minimum of an image and a constant.
	@param l Image.
	@param r Constant value.
	*/
	template<typename pixel_t, typename param_t> void min(Image<pixel_t>& l, param_t r)
	{
		pointProcessImageParam<pixel_t, param_t, pixel_t, internals::minOp<pixel_t, param_t, pixel_t> >(l, r);
	}

	
	/**
	Thresholds img by another image.
	Sets pixel in img to 1 if pixel value > corresponding pixel value in the threshold image, and to 0 otherwise.
	*/
	template<typename pixel_t, typename pixel2_t> void threshold(Image<pixel_t>& img, const Image<pixel2_t>& threshold)
	{
		pointProcessImageImage<pixel_t, pixel2_t, pixel_t, internals::thresholdOp<pixel_t, pixel2_t> >(img, threshold);
	}

	/**
	Thresholds img by constant.
	Sets pixel to 1 if pixel value > threshold, and to 0 otherwise.
	*/
	template<typename pixel_t, typename param_t> void threshold(Image<pixel_t>& img, param_t threshold)
	{
		pointProcessImageParam<pixel_t, pixel_t, pixel_t, internals::thresholdOp<pixel_t, pixel_t> >(img, math::pixelRound<pixel_t>(threshold));
	}
	

	/**
	Thresholds image with range (minmax.x, minmax.y].
	Sets to 1 those pixels whose value a satisfies range.x < a <= range.y, and to 0 otherwise.
	*/
	template<typename pixel_t> void thresholdRange(Image<pixel_t>& img, const math::Vec2d& range)
	{
		pointProcessImageParam<pixel_t, math::Vec2d, pixel_t, internals::thresholdRangeOp<pixel_t> >(img, range);
	}

	/**
	Performs periodic thresholding.
	Sets pixel to 1 if its value is in range (threshold range start, threshold range end] modulo period; sets the pixel value to 
	0 otherwise. If threshold range is larger than period, sets all pixels to 1. If threshold range start is larger than or equal
	to threshold range end, sets all pixels to 0.
	@param img The image to threshold.
	@param inputs 4-component vector containing (period start, period end, threshold range start, threshold range end).
	*/
	template<typename pixel_t> void thresholdPeriodic(Image<pixel_t>& img, const math::Vec4d& inputs)
	{
		pointProcessImageParam<pixel_t, const math::Vec4d&, pixel_t, internals::thresholdPeriodicOp<pixel_t> >(img, inputs);
	}

	/**
	Maps pixel values linearly.
	Maps position a in range [bounds.x, bounds.y] to position in range [bounds.z, bounds.w].
	i.e. if a = 0.5 and [bounds.x, bounds.y] = [0, 1] and [bounds.z, bounds.w] = [3, 4],
	return value will be 3 + (0.5 - 0) / (1 - 0) * (4 - 3) = 3.5.
	@param bounds Bounds values: beginning of old range, end of old range, beginning of new range, end of new range.
	*/
	template<typename pixel_t> void linearMap(Image<pixel_t>& img, const math::Vec4d& bounds)
	{
		pointProcessImageParam<pixel_t, const math::Vec4d&, typename math::NumberUtils<pixel_t>::FloatType, internals::linearMapOp<pixel_t> >(img, bounds);
	}


	namespace tests
	{
		void pointProcess();
		void pointProcessComplex();
	}
}
