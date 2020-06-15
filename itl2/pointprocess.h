#pragma once

#include "image.h"
#include "math/mathutils.h"
#include "math/vec2.h"
#include "math/vec3.h"
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
			img(n) = pixelRound<pixel_t, intermediate_t>(process(img(n)));

			// Showing progress info here would induce more processing than is done in the whole loop.
		}
	}


	/**
	Process corresponding pixels from l and r, place result to l.
	*/
	template<typename pixel1_t, typename pixel2_t, typename intermediate_t, intermediate_t process(pixel1_t, pixel2_t)> void pointProcessImageImage(Image<pixel1_t>& l, const Image<pixel2_t>& r, bool allowBroadcast)
	{
		if (!allowBroadcast)
		{
			l.checkSize(r);

			#pragma omp parallel for if(l.pixelCount() > PARALLELIZATION_THRESHOLD)
			for (coord_t n = 0; n < l.pixelCount(); n++)
			{
				l(n) = pixelRound<pixel1_t, intermediate_t>(process(l(n), r(n)));

				// Showing progress info here would induce more processing than is done in the whole loop.
			}
		}
		else
		{
			Vec3c m(0, 0, 0);
			Vec3c M = r.dimensions() - Vec3c(1, 1, 1);

			#pragma omp parallel for if(l.pixelCount() > PARALLELIZATION_THRESHOLD)
			for (coord_t z = 0; z < l.depth(); z++)
			{
				for (coord_t y = 0; y < l.height(); y++)
				{
					for (coord_t x = 0; x < l.width(); x++)
					{
						Vec3c posl(x, y, z);
						Vec3c posr(x, y, z);
						clamp(posr, m, M);
						l(posl) = pixelRound<pixel1_t, intermediate_t>(process(l(posl), r(posr)));
					}
				}
			}
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
			l(n) = pixelRound<pixel1_t, intermediate_t>(process(l(n), r(n), c));

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
			img(n) = pixelRound<pixel_t, intermediate_t>(process(img(n), param));

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
				img(n) = pixelRound<pixel_t, intermediate_t>(process(pix, param));

			// Showing progress info here would induce more processing than is done in the whole loop.
		}
	}





	namespace internals
	{
		
		// Unary operators

		template<typename pixel_t> pixel_t swapByteOrderOp(pixel_t val)
		{
			return swapByteOrder(val);
		}

		template<typename pixel_t, typename intermediate_t> intermediate_t negateOp(pixel_t val)
		{
			return NumberUtils<intermediate_t>::saturatingSubtract(0, (intermediate_t)val);
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
			// std::abs is ambiguous (and unnecessary) for unsigned primitive types
			if constexpr (std::is_signed_v<pixel_t> || std::is_compound_v<pixel_t>)
				return std::abs((intermediate_t)val);
			else
				return (intermediate_t)val;
		}

		template<typename pixel_t, typename intermediate_t> intermediate_t logOp(pixel_t val)
		{
			return std::log((intermediate_t)val);
		}

		template<typename pixel_t, typename intermediate_t> intermediate_t negLogOp(pixel_t val)
		{
			return -std::log((intermediate_t)val);
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
			return (intermediate_t)std::round((intermediate_t)val);
		}

		template<typename pixel_t, typename intermediate_t> intermediate_t ceilOp(pixel_t val)
		{
			return (intermediate_t)std::ceil((intermediate_t)val);
		}

		template<typename pixel_t, typename intermediate_t> intermediate_t floorOp(pixel_t val)
		{
			return (intermediate_t)std::floor((intermediate_t)val);
		}


		// Unary operators for complex data type

		template<typename pixel_t, typename intermediate_t> intermediate_t conjugateOp(pixel_t val)
		{
			return conj(val);
		}

		template<typename pixel_t, typename intermediate_t> intermediate_t normalizeOp(pixel_t val)
		{
			intermediate_t L = std::abs((intermediate_t)val);
			if (!NumberUtils<intermediate_t>::equals(L, 0))
				return (intermediate_t)val / L;
			else
				return (intermediate_t)val;
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
			return NumberUtils<intermediate_t>::saturatingAdd((intermediate_t)a, (intermediate_t)b);
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
			return NumberUtils<intermediate_t>::saturatingSubtract((intermediate_t)a, (intermediate_t)b);
		}

		/**
		b-a
		*/
		template<typename t1, typename t2, typename intermediate_t> intermediate_t invsubtractOp(t1 a, t2 b)
		{
			return NumberUtils<intermediate_t>::saturatingSubtract((intermediate_t)b, (intermediate_t)a);
		}

		/**
		b-a+c
		*/
		template<typename t1, typename t2, typename t3, typename intermediate_t> intermediate_t invsubtractAddOp(t1 a, t2 b, t3 c)
		{
			return (intermediate_t)b - (intermediate_t)a + (intermediate_t)c;
		}

		template<typename t1, typename t2, typename intermediate_t> intermediate_t multiplyOp(t1 a, t2 b)
		{
			return NumberUtils<intermediate_t>::saturatingMultiply((intermediate_t)a, (intermediate_t)b);
		}

		template<typename t1, typename t2, typename intermediate_t> intermediate_t divideOp(t1 a, t2 b)
		{
			return NumberUtils<intermediate_t>::saturatingDivide((intermediate_t)a, (intermediate_t)b);
		}

		template<typename pixel_t> pixel_t setValueOp(pixel_t val, pixel_t param)
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
			//if (a > pixelRound<pixel_t>(b))
			if(intuitive::gt(a, b))
				return (pixel_t)1;
			else
				return (pixel_t)0;
		}

		/**
		Thresholds a with range (minmax.x, minmax.y].
		Returns 1 if minmax.x < a <= minmax.y, and 0 otherwise.
		*/
		template<typename pixel_t> pixel_t thresholdRangeOp(pixel_t a, const Vec2d minmax)
		{
			if (minmax.x < a && a <= minmax.y)
				return (pixel_t)1;
			else
				return (pixel_t)0;
		}

		/**
		Thresholds pixels to multiple classes.
		Returns 0 if pixel value is less than the first threshold.
		Returns 1 if pixel value is larger than or equal to the first threshold and less than the second threshold (if any)
		Returns 2 if pixel value is larger than or equal to the second threshold and less than the third threshold (if any)
		etc.
		*/
		template<typename pixel_t> pixel_t multiThresholdOp(pixel_t a, const std::vector<pixel_t>& thresholds)
		{
			size_t n = 0;
			for (; n < thresholds.size(); n++)
			{
				if (a < thresholds[n])
					break;
			}

			return pixelRound<pixel_t>(n);
		}

		/**
		Performs periodic thresholding.
		@param x The value to threshold.
		@param inputs 4-component vector containing (period start, period end, threshold range start, threshold range end).
		@return 1 if x is in (threshold range start, threshold range end] modulo period; 0 otherwise. If threshold range is larger than period,
		returns always 1. If threshold range start is larger than or equal to threshold range end, returns always 0.
		*/
		template<typename pixel_t> pixel_t thresholdPeriodicOp(pixel_t x, const Vec4d& inputs)
		{
			double periodStart = inputs[0];
			double periodEnd = inputs[1];
			double thresholdStart = inputs[2];
			double thresholdEnd = inputs[3];

			double periodLength = periodEnd - periodStart;
			double periodicX = realmod(x - periodStart, periodLength) + periodStart;

			if (thresholdStart >= thresholdEnd)
			{
				return (pixel_t)0;
			}
			else if (thresholdEnd - thresholdStart < periodLength)
			{
				double periodicStart = realmod(thresholdStart - periodStart, periodLength) + periodStart;
				double periodicEnd = realmod(thresholdEnd - periodStart, periodLength) + periodStart;

				if (NumberUtils<double>::lessThan(periodicStart, periodicEnd))
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
		template<typename pixel_t, typename intermediate_t = typename NumberUtils<pixel_t>::FloatType> intermediate_t linearMapOp(pixel_t a, const Vec4d& bounds)
		{
			intermediate_t min = (intermediate_t)bounds.x;
			intermediate_t max = (intermediate_t)bounds.y;
			intermediate_t newmin = (intermediate_t)bounds.z;
			intermediate_t newmax = (intermediate_t)bounds.w;

			//if (a < min)
			if(intuitive::lt(a, min))
				a = (pixel_t)min;
			//else if (a > max)
			else if(intuitive::gt(a, max))
				a = (pixel_t)max;

			return newmin + (a - min) / (max - min) * (newmax - newmin);
		}

		/**
		Replaces color data.x by data.y.
		*/
		template<typename pixel_t> pixel_t replaceOp(pixel_t a, const Vec2<pixel_t> data)
		{
			if (NumberUtils<pixel_t>::equals(a, data.x) || (NumberUtils<pixel_t>::isnan(a) && NumberUtils<pixel_t>::isnan(data.x)))
				return data.y;

			return a;
		}
	}


#define DEF_OPERATION_ALL_FLOAT_INTERMEDIATE_TYPE(name, help) \
	/** \
	help \
	\
	@param in Image to process. \
	*/ \
	template<typename pixel_t> void name (Image<pixel_t>& img) \
	{ \
		pointProcess<pixel_t, typename NumberUtils<pixel_t>::FloatType, internals::name##Op<pixel_t, typename NumberUtils<pixel_t>::FloatType> >(img); \
	}

#define DEF_OPERATION_ALL_NO_INTERMEDIATE_TYPE(name, help) \
	/** \
	help \
	\
	@param in Image to process. \
	*/ \
	template<typename pixel_t> void name (Image<pixel_t>& img) \
	{ \
		pointProcess<pixel_t, pixel_t, internals::name##Op<pixel_t> >(img); \
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



	DEF_OPERATION_ALL_NO_INTERMEDIATE_TYPE(swapByteOrder, "Swaps byte order of each pixel value. Use this command to convert images read in wrong endianness to the correct one, or to before saving and image if it should be saved in non-native byte order.")
	DEF_OPERATION_ALL_NO_INTERMEDIATE_TYPE(negate, "Negates pixel values. The processing is performed using saturation arithmetic.")
	DEF_OPERATION_ALL_NO_INTERMEDIATE_TYPE(abs, "Calculates absolute value of pixel values.")

	DEF_OPERATION_ALL_FLOAT_INTERMEDIATE_TYPE(exponentiate, "Exponentates pixel values.")
	DEF_OPERATION_ALL_FLOAT_INTERMEDIATE_TYPE(square, "Calculates square of pixel values.")
	DEF_OPERATION_ALL_FLOAT_INTERMEDIATE_TYPE(squareRoot, "Calculates square root of pixel values.")
	DEF_OPERATION_ALL_FLOAT_INTERMEDIATE_TYPE(log, "Calculates natural logarithm of pixel values.")
	DEF_OPERATION_ALL_FLOAT_INTERMEDIATE_TYPE(negLog, "Calculates negative of natural logarithm of pixel values.")
	DEF_OPERATION_ALL_FLOAT_INTERMEDIATE_TYPE(log10, "Calculates base-10 logarithm of pixel values.")
	DEF_OPERATION_ALL_FLOAT_INTERMEDIATE_TYPE(sin, "Calculates sine of pixel values.")
	DEF_OPERATION_ALL_FLOAT_INTERMEDIATE_TYPE(cos, "Calculates cosine of pixel values.")
	DEF_OPERATION_ALL_FLOAT_INTERMEDIATE_TYPE(tan, "Calculates tangent of pixel values.")
	DEF_OPERATION_ALL_FLOAT_INTERMEDIATE_TYPE(inv, "Calculates inverse (1/x) of pixel values.")

	DEF_OPERATION_ALL_FLOAT_INTERMEDIATE_TYPE(round, "Rounds pixel values.")
	DEF_OPERATION_ALL_FLOAT_INTERMEDIATE_TYPE(ceil, "Calculates ceiling of pixel values.")
	DEF_OPERATION_ALL_FLOAT_INTERMEDIATE_TYPE(floor, "Calculates floor of pixel values.")

	DEF_OPERATION_COMPLEX(conjugate, "Calculates complex conjugate of pixel values.")
	DEF_OPERATION_COMPLEX(normalize, "Makes norm of each pixel one.")
	DEF_OPERATION_COMPLEX(real, "Calculates real part of pixel values.")
	DEF_OPERATION_COMPLEX(imag, "Calculates imaginary part of pixel values.")
	DEF_OPERATION_COMPLEX(arg, "Calculates argument of pixel values.")
	DEF_OPERATION_COMPLEX(normSquared, "Calculates squared norm of pixel values.")

	// These make sure that we have normal non-image functions available, too.
	//using std::floor;
	//using std::ceil;
	//using std::round;
	using std::abs;
	using std::log;
	using std::sin;
	using std::cos;
	using std::tan;
	using std::real;
	using std::imag;
	using std::arg;
	


	/**
	Sets pixel values of first image to those copied from the second image.
	@param l First image.
	@param r Second image.
	*/
	template<typename pixel_t, typename pixel2_t> void setValue(Image<pixel_t>& l, const Image<pixel2_t>& r, bool allowBroadcast = false)
	{
		l.ensureSize(r);
		pointProcessImageImage<pixel_t, pixel2_t, pixel2_t, internals::copyOp<pixel_t, pixel2_t> >(l, r, allowBroadcast);
	}

	/**
	Sets pixel values in an image to a constant.
	@param img Image.
	@param value Constant value.
	*/
	template<typename pixel_t, typename param_t> void setValue(Image<pixel_t>& img, param_t value)
	{
		pointProcessImageParam<pixel_t, pixel_t, pixel_t, internals::setValueOp<pixel_t> >(img, pixelRound<pixel_t>(value));
	}



	// NOTE: This region contains structures that are used to find intermediate type for mathematical operations such as +, -, * and /.
	// The intermediate types are found roughly as follows:
	// Type1		Type2				Resulting intermediate type
	// class		class				class
	// class		not class			ERROR (void type)
	// double		any					double
	// float32_t	any not double		float32_t
	// larger int	smaller int			larger int
	// larger uint	smaller uint		larger uint
	// uint			int					smallest type that can represent values of both types


	/**
	Defines wider::type as the wider of the two argument types.
	*/
	template<class T, class U> struct wider {
		using type = typename std::conditional<sizeof(T) >= sizeof(U), T, U>::type;
	};


	/**
	Finds out a signed type that can contain values of both argument types.
	*/
	template<class T, class U> struct one_wider_signed {
		using type = typename wider<
			typename NumberUtils<T>::SignedType,
			typename NumberUtils<U>::SignedType
		>::type;
	};

	/**
	Finds out type suitable to be used as intermediate type in mathematical calculations, given two signed or unsigned integer argument types.
	*/
	template<class T, class U> struct math_intermediate_int {
		using type = typename std::conditional <
			std::is_signed_v<T>,
			// T is signed. Test if U is.
			typename std::conditional<
			std::is_signed_v<U>,
			// T is signed and U is signed. Return wider type.
			typename wider<T, U>::type,
			// T is signed and U is not signed. Return signed type wider than both arguments
			typename one_wider_signed<T, U>::type
			>::type,
			// T is not signed. Test if U is.
			typename std::conditional<
			std::is_signed_v<U>,
			// T is not signed, U is signed. Return signed type wider than both arguments.
			typename one_wider_signed<T, U>::type
			,
			// T is not signed, U is not signed. Return wider type.
			typename wider<T, U>::type
			>::type
		>::type;
	};

	/**
	Finds out type suitable to be used as intermediate type in mathematical calculations, given two argument types that are signed or unsigned float or int.
	*/
	template<class T, class U> struct math_intermediate_float_or_int {
		using type = typename std::conditional <
			std::is_floating_point_v<T>,

			// T is floating point. Test if U is
			typename std::conditional<
			std::is_floating_point_v<U>,

			// Both T and U are floating point. Return the larger type.
			typename wider<T, U>::type,
			// T is floating point, U is not. Return T
			T
			>::type,

			// T is not floating point. Test if U is
			typename std::conditional<
			std::is_floating_point_v<U>,
			// T is not floating point, U is. Return U
			U,
			// Both T and U are not floating point. Return the wider type.
			//typename wider<T, U>::type
			typename math_intermediate_int<T, U>::type
			>::type
		> ::type;
	};


	/**
	Finds out type suitable to be used as intermediate type in mathematical calculations, given two argument types.
	*/
	template<class T, class U> struct math_intermediate_type {
		using type = typename std::conditional <
			std::is_compound_v<T>,

			// T is class type. Test if U is
			typename std::conditional<
			std::is_compound_v<U>,

			// Both T and U are class types. Test if they are the same.
			typename std::conditional<
			std::is_same_v<T, U>,
			// T and U are the same, return T.
			T,
			// T and U are not the same. We don't know what is going on, so return void type to mark an error.
			std::void_t<T, U>
			>::type,

			// T is class, U is not. Return T
			T

			>::type,

			// T is not class. Test if U is
			typename std::conditional<
			std::is_compound_v<U>,

			// T is not class, U is. Return U.
			U,

			// Both T and U are not class types.
			typename math_intermediate_float_or_int<T, U>::type
			>::type
		> ::type;
	};

	



	/**
	Adds two images and places the result to the first image.
	@param l First image.
	@param r Second image.
	*/
	template<typename pixel1_t, typename pixel2_t> void add(Image<pixel1_t>& l, const Image<pixel2_t>& r, bool allowBroadcast = false)
	{
		using intermediate_t = typename math_intermediate_type<pixel1_t, pixel2_t>::type;
		pointProcessImageImage<pixel1_t, pixel2_t, intermediate_t, internals::addOp<pixel1_t, pixel2_t, intermediate_t> >(l, r, allowBroadcast);
	}

	/**
	Adds constant to pixel values.
	@param l Image.
	@param r Constant value.
	*/
	template<typename pixel_t, typename param_t> void add(Image<pixel_t>& l, param_t r)
	{
		using intermediate_t = typename math_intermediate_type<pixel_t, param_t>::type;
		pointProcessImageParam<pixel_t, param_t, intermediate_t, internals::addOp<pixel_t, param_t, intermediate_t> >(l, r);
	}

	/**
	Adds constant to all pixel values except those that correspond to a specific indicator value.
	@param l Image.
	@param r Constant value.
	@param badValue Value indicating pixels in the image that should not be processed.
	*/
	template<typename pixel_t, typename param_t> void maskedAdd(Image<pixel_t>& l, param_t r, pixel_t badValue)
	{
		using intermediate_t = typename math_intermediate_type<pixel_t, param_t>::type;
		maskedPointProcessImageParam<pixel_t, param_t, intermediate_t, internals::addOp<pixel_t, param_t, intermediate_t> >(l, r, badValue);
	}




	/**
	Subtracts second image from the first image and places the result to the first image (first image = first image - second image).
	@param l First image.
	@param r Second image.
	*/
	template<typename pixel_t, typename pixel2_t> void subtract(Image<pixel_t>& l, const Image<pixel2_t>& r, bool allowBroadcast = false)
	{
		using intermediate_t = typename math_intermediate_type<pixel_t, pixel2_t>::type;
		pointProcessImageImage<pixel_t, pixel2_t, intermediate_t, internals::subtractOp<pixel_t, pixel2_t, intermediate_t> >(l, r, allowBroadcast);
	}

	/**
	Subtracts constant from pixel values.
	@param l Image.
	@param r Constant value.
	*/
	template<typename pixel_t, typename param_t> void subtract(Image<pixel_t>& l, param_t r)
	{
		using intermediate_t = typename math_intermediate_type<pixel_t, param_t>::type;
		pointProcessImageParam<pixel_t, param_t, intermediate_t, internals::subtractOp<pixel_t, param_t, intermediate_t> >(l, r);
	}

	/**
	Subtracts first image from the second image and places the result to the first image (first image = second image - first image).
	@param l First image.
	@param r Second image.
	*/
	template<typename pixel_t, typename pixel2_t> void invsubtract(Image<pixel_t>& l, const Image<pixel2_t>& r, bool allowBroadcast = false)
	{
		using intermediate_t = typename math_intermediate_type<pixel_t, pixel2_t>::type;
		pointProcessImageImage<pixel_t, pixel2_t, intermediate_t, internals::invsubtractOp<pixel_t, pixel2_t, intermediate_t> >(l, r, allowBroadcast);
	}

	/**
	Subtracts pixel values from a constant.
	@param l Image.
	@param r Constant value.
	*/
	template<typename pixel_t, typename param_t> void invsubtract(Image<pixel_t>& l, param_t r)
	{
		using intermediate_t = typename math_intermediate_type<pixel_t, param_t>::type;
		pointProcessImageParam<pixel_t, param_t, intermediate_t, internals::invsubtractOp<pixel_t, param_t, intermediate_t> >(l, r);
	}

	/**
	Subtracts first image from the second image, adds constant, and places the result to the first image (first image = second image - first image + constant).
	@param l First image.
	@param r Second image.
	@param constant Constant to add.
	*/
	template<typename pixel_t> void invsubtractAdd(Image<pixel_t>& l, const Image<pixel_t>& r, pixel_t constant)
	{
		using intermediate_t = typename math_intermediate_type<pixel_t, pixel_t>::type;
		pointProcessImageImageParam<pixel_t, pixel_t, pixel_t, intermediate_t, internals::invsubtractAddOp<pixel_t, pixel_t, pixel_t, intermediate_t> >(l, r, constant);
	}

	/**
	Multiplies two images and places the result to the first image.
	@param l First image.
	@param r Second image.
	*/
	template<typename pixel1_t, typename pixel2_t> void multiply(Image<pixel1_t>& l, const Image<pixel2_t>& r, bool allowBroadcast = false)
	{
		using intermediate_t = typename math_intermediate_type<pixel1_t, pixel2_t>::type;
		pointProcessImageImage<pixel1_t, pixel2_t, intermediate_t, internals::multiplyOp<pixel1_t, pixel2_t, intermediate_t> >(l, r, allowBroadcast);
	}

	/**
	Multiplies pixel values by a constant.
	@param l Image.
	@param r Constant value.
	*/
	template<typename pixel_t, typename param_t> void multiply(Image<pixel_t>& l, param_t r)
	{
		using intermediate_t = typename math_intermediate_type<pixel_t, param_t>::type;
		pointProcessImageParam<pixel_t, param_t, intermediate_t, internals::multiplyOp<pixel_t, param_t, intermediate_t> >(l, r);
	}


	/**
	Raises pixels of first image to power given by pixel value of the second image, and places the result to the first image.
	@param l First image.
	@param r Second image.
	*/
	template<typename pixel_t, typename pixel2_t> void pow(Image<pixel_t>& l, const Image<pixel2_t>& r, bool allowBroadcast = false)
	{
		pointProcessImageImage<pixel_t, pixel2_t, typename NumberUtils<pixel_t>::FloatType, internals::powOp<pixel_t, pixel2_t, typename NumberUtils<pixel_t>::FloatType> >(l, r, allowBroadcast);
	}

	/**
	Raises pixel values to power.
	@param l Image.
	@param r Exponent.
	*/
	template<typename pixel_t, typename param_t> void pow(Image<pixel_t>& l, param_t r)
	{
		pointProcessImageParam<pixel_t, param_t, typename NumberUtils<pixel_t>::FloatType, internals::powOp<pixel_t, param_t, typename NumberUtils<pixel_t>::FloatType> >(l, r);
	}


	/**
	Divides two images and places the result to the first image (first image = first image / second image).
	@param l First image.
	@param r Second image.
	*/
	template<typename pixel1_t, typename pixel2_t> void divide(Image<pixel1_t>& l, const Image<pixel2_t>& r, bool allowBroadcast = false)
	{
		// NOTE: Here we use FloatType as intermediate as we want the division to be correctly rounded instead of e.g. integer division.
		pointProcessImageImage<pixel1_t, pixel2_t, typename NumberUtils<pixel1_t>::FloatType, internals::divideOp<pixel1_t, pixel2_t, typename NumberUtils<pixel1_t>::FloatType> >(l, r, allowBroadcast);
	}

	/**
	Divides pixel values by a constant.
	@param l Image.
	@param r Constant value.
	*/
	template<typename pixel_t,  typename param_t> void divide(Image<pixel_t>& l, param_t r)
	{
		// NOTE: Here we use FloatType as intermediate as we want the division to be correctly rounded instead of e.g. integer division.
		pointProcessImageParam<pixel_t, param_t, typename NumberUtils<pixel_t>::FloatType, internals::divideOp<pixel_t, param_t, typename NumberUtils<pixel_t>::FloatType> >(l, r);
	}



	/**
	Calculates maximum of two images and places the result to the first image.
	@param l First image.
	@param r Second image.
	*/
	template<typename pixel1_t, typename pixel2_t> void max(Image<pixel1_t>& l, const Image<pixel2_t>& r, bool allowBroadcast = false)
	{
		pointProcessImageImage<pixel1_t, pixel2_t, pixel1_t, internals::maxOp<pixel1_t, pixel2_t, pixel1_t> >(l, r, allowBroadcast);
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
	template<typename pixel1_t, typename pixel2_t> void min(Image<pixel1_t>& l, const Image<pixel2_t>& r, bool allowBroadcast = false)
	{
		pointProcessImageImage<pixel1_t, pixel2_t, pixel1_t, internals::minOp<pixel1_t, pixel2_t, pixel1_t> >(l, r, allowBroadcast);
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
	template<typename pixel_t, typename pixel2_t> void threshold(Image<pixel_t>& img, const Image<pixel2_t>& threshold, bool allowBroadcast = false)
	{
		pointProcessImageImage<pixel_t, pixel2_t, pixel_t, internals::thresholdOp<pixel_t, pixel2_t> >(img, threshold, allowBroadcast);
	}

	/**
	Thresholds img by constant.
	Sets pixel to 1 if pixel value > threshold, and to 0 otherwise.
	*/
	template<typename pixel_t, typename param_t> void threshold(Image<pixel_t>& img, param_t threshold)
	{
		pointProcessImageParam<pixel_t, pixel_t, pixel_t, internals::thresholdOp<pixel_t, pixel_t> >(img, pixelRound<pixel_t>(threshold));
	}
	

	/**
	Thresholds image with range (minmax.x, minmax.y].
	Sets to 1 those pixels whose value a satisfies range.x < a <= range.y, and to 0 otherwise.
	*/
	template<typename pixel_t> void thresholdRange(Image<pixel_t>& img, const Vec2d& range)
	{
		pointProcessImageParam<pixel_t, Vec2d, pixel_t, internals::thresholdRangeOp<pixel_t> >(img, range);
	}

	/**
	Thresholds img by set of thresholds.
	Sets pixel to 0 if its value is less than the first threshold.
	Sets pixel to 1 if pixel its is larger than or equal to the first threshold and less than the second threshold (if any).
	Sets pixel to 2 if pixel its is larger than or equal to the second threshold and less than the third threshold (if any).
	etc.
	*/
	template<typename pixel_t> void multiThreshold(Image<pixel_t>& img, const std::vector<pixel_t>& thresholds)
	{
		pointProcessImageParam<pixel_t, const std::vector<pixel_t>&, pixel_t, internals::multiThresholdOp<pixel_t> >(img, thresholds);
	}

	/**
	Performs periodic thresholding.
	Sets pixel to 1 if its value is in range (threshold range start, threshold range end] modulo period; sets the pixel value to 
	0 otherwise. If threshold range is larger than period, sets all pixels to 1. If threshold range start is larger than or equal
	to threshold range end, sets all pixels to 0.
	@param img The image to threshold.
	@param inputs 4-component vector containing (period start, period end, threshold range start, threshold range end).
	*/
	template<typename pixel_t> void thresholdPeriodic(Image<pixel_t>& img, const Vec4d& inputs)
	{
		pointProcessImageParam<pixel_t, const Vec4d&, pixel_t, internals::thresholdPeriodicOp<pixel_t> >(img, inputs);
	}

	/**
	Maps pixel values linearly.
	Maps position a in range [bounds.x, bounds.y] to position in range [bounds.z, bounds.w].
	i.e. if a = 0.5 and [bounds.x, bounds.y] = [0, 1] and [bounds.z, bounds.w] = [3, 4],
	return value will be 3 + (0.5 - 0) / (1 - 0) * (4 - 3) = 3.5.
	@param bounds Bounds values: beginning of old range, end of old range, beginning of new range, end of new range.
	*/
	template<typename pixel_t> void linearMap(Image<pixel_t>& img, const Vec4d& bounds)
	{
		// NOTE: Float intermediate type for good rounding performance.
		pointProcessImageParam<pixel_t, const Vec4d&, typename NumberUtils<pixel_t>::FloatType, internals::linearMapOp<pixel_t> >(img, bounds);
	}

	/**
	Sets to v.b those pixels whose value is v.a.
	*/
	template<typename pixel_t> void replace(Image<pixel_t>& img, const Vec2<pixel_t>& v)
	{
		pointProcessImageParam<pixel_t, Vec2<pixel_t>, pixel_t, internals::replaceOp<pixel_t> >(img, v);
	}


	namespace tests
	{
		void pointProcess();
		void pointProcessComplex();
		void broadcast();
		void byteOrder();
		void intermediateTypes();
	}
}
