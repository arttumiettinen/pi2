#pragma once

#include <math.h>
#include <limits>
#include "datatypes.h"
#include "type.h"

#ifdef min
	#undef min
#endif
#ifdef max
	#undef max
#endif

using namespace itl2;

namespace math
{
	const double PI = 3.14159265359;
	const float PIf = 3.14159265359f;
		
#if defined(_MSC_VER)

	template <typename T> bool isnan(const T &x)
	{
		return _isnan(x) != 0;
	}

	template <typename T> bool isinf(const T &x)
	{
		return _finite(x) == 0;
	}
#else
	template <typename T> bool isnan(const T &x)
	{
		return std::isnan(x);
	}

	template <typename T> bool isinf(const T &x)
	{
		return std::isinf(x);
	}
#endif


		

	/**
	Function that should behave like MatLab's mod function.
	*/
	inline double realmod(double x, double y)
	{
		double result = fmod(x, y);
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

	/**
	Calculates minimum of two values.
	*/
	template<typename T> T min(const T& a, const T& b)
	{
		if(a < b)
			return a;
		return b;
	}

	/**
	Calculates maximum of two values.
	*/
	template<typename T> T max(const T& a, const T& b)
	{
		if(a > b)
			return a;
		return b;
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
		r = sqrt(x * x + y * y + z * z);
		azimuthal = atan2(y, x);
		polar = acos(z / r);
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
		x = r * cos(azimuthal) * sin(polar);
		y = r * sin(azimuthal) * sin(polar);
		z = r * cos(polar);
	}

	/**
	 * Convert from Cartesian to polar coordinates.
	 * @param x, y The cartesian coordinates.
	 * @param r Radial coordinate. (from origo outwards. [0, inf])
	 * @param azimuthal Azimuthal angle. (from positive x-axis towards positive y-axis, [-pi, pi])
	 */
	inline void toPolar(double x, double y, double& r, double& azimuthal)
	{
		r = sqrt(x * x + y * y);
		azimuthal = atan2(y, x);
	}

	/**
	 * Convert from polar to Cartesian coordinates.
	 * @param r Radial coordinate.
	 * @param azimuthal Azimuthal angle. (from positive x-axis towards positive y-axis, [-pi, pi])
	 * @param x, y The cartesian coordinates.
	 */
	inline void toCartesian(double r, double azimuthal, double& x, double& y)
	{
		x = r * cos(azimuthal);
		y = r * sin(azimuthal);
	}

	
	// *** Rounding of pixel values

	/**
	* Use this to round pixel values from double to lower precision numbers.
	* Also clamps values correctly.
	*/
	template<typename Tout, typename Tin> Tout pixelRound(Tin value)
	{
		if constexpr (std::is_integral<Tout>::value && std::is_integral<Tin>::value)
		{
			if (intuitive::ge(value, std::numeric_limits<Tout>::max()))
				return std::numeric_limits<Tout>::max();
			else if (intuitive::le(value, std::numeric_limits<Tout>::lowest()))
				return std::numeric_limits<Tout>::lowest();

			return (Tout)::round(value);
		}
		else if constexpr(std::is_integral<Tout>::value)
		{
			if(value >= std::numeric_limits<Tout>::max())
				return std::numeric_limits<Tout>::max();
			else if(value <= std::numeric_limits<Tout>::lowest())
				return std::numeric_limits<Tout>::lowest();

			return (Tout)::round(value);
		}
		else // Tout is floating point value
		{
			return (Tout)value;
		}
	}

	///**
	//* Returns the parameter as-is.
	//*/
	//template<> inline double pixelRound(double value)
	//{
	//	return value;
	//}

	///**
	//* Returns the parameter converted to float.
	//*/
	//template<> inline float pixelRound(double value)
	//{
	//	return (float)value;
	//}

	///**
	//* Returns the parameter as-is.
	//*/
	//template<> inline float pixelRound(float value)
	//{
	//	return value;
	//}

	///**
	//* Returns the parameter as-is.
	//*/
	//template<> inline double pixelRound(float value)
	//{
	//	return (double)value;
	//}

	/**
	* Returns the parameter as-is.
	*/
	template<> inline std::complex<float32_t> pixelRound(std::complex<float32_t> value)
	{
		return value;
	}

	
	/**
	Round double value to coordinate type.
	*/
	inline coord_t round(double value)
	{
		return (coord_t)::round(value);
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
	Returns random coordinate in range [0, max-1] (==[0, max[).
	*/
	inline coord_t randc(coord_t max)
	{
		return math::round(frand() * (max - 1));
	}

}

