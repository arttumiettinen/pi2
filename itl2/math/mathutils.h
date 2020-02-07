#pragma once

//#include <math.h>
#include <cmath>
#include <limits>
#include "datatypes.h"
#include "type.h"

#ifdef min
	#undef min
#endif
#ifdef max
	#undef max
#endif

namespace itl2
{
	const double PI = 3.14159265359;
	const float PIf = 3.14159265359f;
	
	/**
	Converts from degrees to radians.
	*/
	inline double degToRad(double deg)
	{
		return deg / 180.0 * PI;
	}

	/**
	Converts from radians to degrees.
	*/
	inline double radToDeg(double rad)
	{
		return rad / PI * 180.0;
	}

	
	/**
	Returns largest integer value whose square is less than given value.
	*/
	template<typename T> inline T largestIntWhoseSquareIsLessThan(T square)
	{
		// Initial guess using floating point math
		T result = (T)std::floor(std::sqrt(square));

		// Refine the result in the case there are floating point inaccuracies
		while (result * result < square)
			result++;
		result--;
		return result;
	}

	/**
	Function that should behave like MatLab's mod function.
	*/
	inline double realmod(double x, double y)
	{
		double result = std::fmod(x, y);
		return result >= 0 ? result : result + y;
	}

	/**
	Clamps the given value to range [lower, upper].
	*/
	template<typename T> void clamp(T& value, T lower, T upper)
	{
		if(value < lower)
			value = lower;
		else if(value > upper)
			value = upper;
	}

	///**
	//Calculates minimum of two values.
	//*/
	//template<typename T> T min(const T& a, const T& b)
	//{
	//	if(a < b)
	//		return a;
	//	return b;
	//}

	///**
	//Calculates maximum of two values.
	//*/
	//template<typename T> T max(const T& a, const T& b)
	//{
	//	if(a > b)
	//		return a;
	//	return b;
	//}

	/**
	Returns hypotenuse of a and b while avoiding underflow/overflow
	using (a * sqrt( 1 + (b/a) * (b/a))), rather than sqrt(a*a + b*b).
	*/
	template<typename T> T hypot(const T &a, const T &b)
	{
		if (a == 0)
		{
			return std::abs(b);
		}
		else
		{
			T c = b / a;
			return std::abs(a) * std::sqrt(1 + c * c);
		}
	}

// *** Conversion to/from spherical and polar coordinates.
	/**
	 * Convert from Cartesian to spherical coordinates.
	 * @param x, y, z The cartesian coordinates.
	 * @param r Radial coordinate. (from origo outwards. [0, inf])
	 * @param azimuthal Azimuthal angle. (from positive x-axis towards positive y-axis, [-pi, pi])
	 * @param polar Polar angle. (from positive z-axis towards xy-plane, [0, pi])
	 */
	inline void toSpherical(double x, double y, double z, double& r, double& azimuthal, double& polar)
	{
		r = std::sqrt(x * x + y * y + z * z);
		azimuthal = std::atan2(y, x);
		polar = std::acos(z / r);
	}

	/**
	 * Convert from spherical to Cartesian coordinates.
	 * @param r Radial coordinate.
	 * @param azimuthal Azimuthal angle. (from positive x-axis towards positive y-axis, [-pi, pi])
	 * @param polar Polar angle. (from positive z-axis towards xy-plane, [0, pi])
	 * @param x, y, z The cartesian coordinates.
	 */
	inline void toCartesian(double r, double azimuthal, double polar, double& x, double& y, double& z)
	{
		x = r * std::cos(azimuthal) * std::sin(polar);
		y = r * std::sin(azimuthal) * std::sin(polar);
		z = r * std::cos(polar);
	}

	/**
	 * Convert from Cartesian to polar coordinates.
	 * @param x, y The cartesian coordinates.
	 * @param r Radial coordinate. (from origo outwards. [0, inf])
	 * @param azimuthal Azimuthal angle. (from positive x-axis towards positive y-axis, [-pi, pi])
	 */
	inline void toPolar(double x, double y, double& r, double& azimuthal)
	{
		r = std::sqrt(x * x + y * y);
		azimuthal = std::atan2(y, x);
	}

	/**
	 * Convert from polar to Cartesian coordinates.
	 * @param r Radial coordinate.
	 * @param azimuthal Azimuthal angle. (from positive x-axis towards positive y-axis, [-pi, pi])
	 * @param x, y The cartesian coordinates.
	 */
	inline void toCartesian(double r, double azimuthal, double& x, double& y)
	{
		x = r * std::cos(azimuthal);
		y = r * std::sin(azimuthal);
	}

	
	// *** Rounding of pixel values

	/**
	* Use this to round pixel values from double to lower precision numbers.
	* Also clamps values correctly.
	*/
	template<typename Tout, typename Tin> Tout pixelRound(Tin value)
	{
		if constexpr (std::is_same<Tin, Tout>::value)
		{
			return value;
		}
		else if constexpr (std::is_integral<Tout>::value && std::is_integral<Tin>::value)
		{
			if (intuitive::ge(value, std::numeric_limits<Tout>::max()))
				return std::numeric_limits<Tout>::max();
			else if (intuitive::le(value, std::numeric_limits<Tout>::lowest()))
				return std::numeric_limits<Tout>::lowest();

			return (Tout)value;
		}
		else if constexpr(std::is_integral<Tout>::value)
		{
			if(value >= std::numeric_limits<Tout>::max())
				return std::numeric_limits<Tout>::max();
			else if(value <= std::numeric_limits<Tout>::lowest())
				return std::numeric_limits<Tout>::lowest();

			return (Tout)::round(value);
		}
		else if constexpr (std::is_same<Tout, complex32_t>::value)
		{
			// Convert real to complex
			complex32_t c(pixelRound<itl2::float32_t>(value), 0);
			return c;
		}
		else // Tout is floating point value
		{
			return (Tout)value;
		}
	}

	/**
	* Returns the parameter as-is.
	*/
	//template<> inline std::complex<itl2::float32_t> pixelRound(std::complex<itl2::float32_t> value)
	//{
	//	return value;
	//}

	/**
	This overload is needed so that some algorithms work with e.g. both Vec3<float32_t> and Vec3<coord_t> etc.
	*/
	inline itl2::coord_t round(itl2::coord_t value)
	{
		return value;
	}
	
	/**
	Round floating point value to coordinate type.
	*/
	inline itl2::coord_t round(double value)
	{
		return (itl2::coord_t)std::round(value);
	}

	/**
	Round floating point value to coordinate type.
	*/
	inline itl2::coord_t round(float value)
	{
		return (itl2::coord_t)std::round(value);
	}


	/**
	This overload is needed so that some algorithms work with e.g. both Vec3<float32_t> and Vec3<coord_t> etc.
	*/
	inline itl2::coord_t ceil(itl2::coord_t value)
	{
		return value;
	}

	/**
	Ceil floating point value to coordinate type.
	*/
	inline itl2::coord_t ceil(double value)
	{
		return (itl2::coord_t)std::ceil(value);
	}

	/**
	Ceil floating point value to coordinate type.
	*/
	inline itl2::coord_t ceil(float value)
	{
		return (itl2::coord_t)std::ceil(value);
	}


	/**
	This overload is needed so that some algorithms work with e.g. both Vec3<float32_t> and Vec3<coord_t> etc.
	*/
	inline itl2::coord_t floor(itl2::coord_t value)
	{
		return value;
	}

	/**
	Floor floating point value to coordinate type.
	*/
	inline itl2::coord_t floor(double value)
	{
		return (itl2::coord_t)std::floor(value);
	}

	/**
	Floor floating point value to coordinate type.
	*/
	inline itl2::coord_t floor(float value)
	{
		return (itl2::coord_t)std::floor(value);
	}



// *** Random number generation

	/**
	Returns pseudo-random number between 0 and 1.
	*/
	inline double frand()
	{
		return (double)rand() / (double)RAND_MAX;
	}

	/**
	Returns pseudo-random number between 0 and max.
	*/
	inline double frand(double max)
	{
		return frand() * max;
	}

	/**
	Returns pseudo-random number between min and max.
	*/
	inline double frand(double min, double max)
	{
		return min + frand() * (max - min);
	}

	/**
	Returns random coordinate in range [0, max-1] (== [0, max[).
	*/
	inline itl2::coord_t randc(itl2::coord_t max)
	{
		return round(frand() * (max - 1));
	}

	/**
	Returns random coordinate in range [min, max-1] (== [min, max[).
	*/
	inline itl2::coord_t randc(itl2::coord_t min, itl2::coord_t max)
	{
		return min + randc(max - min);
	}

}

