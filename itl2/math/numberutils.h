#pragma once

/*
This file contains saturated arithmetic code adapted from
http://locklessinc.com/articles/sat_arithmetic/
The code can be found in stackoverflow etc. so I suppose it is free to use, although
there is no explicit license.
*/


#include <limits>

#include "datatypes.h"

#if defined(_WIN32)
// Enable intrinsics
#include <intrin.h>
#pragma intrinsic(_mul128, _umul128)
#endif

namespace itl2
{
	template<typename T> class NumberUtils
	{
	public:
		static inline T tolerance()
		{
			T::NumberUtils_undefined;
		}

		/**
		Type that can contain floating point values in the same range than T.
		*/
		typedef double FloatType;

		/**
		Type that can contain floating point values in the same range than T and is real.
		*/
		typedef double RealFloatType;

		/**
		Type that is signed and can contain values of T.
		*/
		typedef double SignedType;

		/**
		Returns typical maximum value of numbers of type T.
		(numeric_limits<T>::max() for integral types, 1.0 for floating point types)
		*/
		static inline T scale()
		{
			return 1;
		}

		/**
		Saturating addition.
		*/
		static inline T saturatingAdd(const T a, const T b)
		{
			return a + b;
		}

		/**
		Saturating subtraction.
		*/
		static inline T saturatingSubtract(const T a, const T b)
		{
			return a - b;
		}

		/**
		Saturating multiplication.
		*/
		static inline T saturatingMultiply(const T a, const T b)
		{
			return a * b;
		}

		/**
		Saturating division.
		*/
		static inline T saturatingDivide(const T a, const T b)
		{
			return a / b;
		}

		
		static inline bool isnan(const T a)
		{
			return false;
		}

		static inline bool equals(const T a, const T b, const T tol = tolerance())
		{
			return fabs(a - b) < tol;
		}

		static inline bool lessThan(const T a, const T b, const T tol = tolerance())
		{
			return a < b && !equals(a, b, tol);
		}

		static inline bool greaterThan(const T a, const T b, const T tol = tolerance())
		{
			return a > b && !equals(a, b, tol);
		}

		static inline bool lessThanOrEqual(const T a, const T b, const T tol = tolerance())
		{
			return a < b || equals(a, b, tol);
		}

		static inline bool greaterThanOrEqual(const T a, const T b, const T tol = tolerance())
		{
			return a > b || equals(a, b, tol);
		}
	};

	template<> class NumberUtils<double>
	{
	public:
		static inline double tolerance()
		{
			return 1e-15;
		}
		
		typedef double FloatType;
		typedef double RealFloatType;
		typedef double SignedType;

		static inline double scale()
		{
			return 1.0;
		}

		static inline double saturatingAdd(const double a, const double b)
		{
			return a + b;
		}

		static inline double saturatingSubtract(const double a, const double b)
		{
			return a - b;
		}

		static inline double saturatingMultiply(const double a, const double b)
		{
			return a * b;
		}

		static inline double saturatingDivide(const double a, const double b)
		{
			return a / b;
		}
		
		static inline bool isnan(const double a)
		{
			return std::isnan(a);
		}

		static inline bool equals(const double a, const double b, const double tol = tolerance())
		{
			return fabs(a - b) < tol;
		}

		static inline bool lessThan(const double a, const double b, const double tol = tolerance())
		{
			return a < b && !equals(a, b, tol);
		}

		static inline bool greaterThan(const double a, const double b, const double tol = tolerance())
		{
			return a > b && !equals(a, b, tol);
		}

		static inline bool lessThanOrEqual(const double a, const double b, const double tol = tolerance())
		{
			return a < b || equals(a, b, tol);
		}

		static inline bool greaterThanOrEqual(const double a, const double b, const double tol = tolerance())
		{
			return a > b || equals(a, b, tol);
		}
		
	};

	template<> class NumberUtils<float32_t>
	{
	public:
		static inline float32_t tolerance()
		{
			return 1e-10f;
		}

		typedef float32_t FloatType;
		typedef float32_t RealFloatType;
		typedef float32_t SignedType;
		
		static inline float32_t scale()
		{
			return 1.0f;
		}

		static inline bool isnan(const float32_t a)
		{
			return std::isnan(a);
		}

		static inline float32_t saturatingAdd(const float32_t a, const float32_t b)
		{
			return a + b;
		}

		static inline float32_t saturatingSubtract(const float32_t a, const float32_t b)
		{
			return a - b;
		}

		static inline float32_t saturatingMultiply(const float32_t a, const float32_t b)
		{
			return a * b;
		}

		static inline float32_t saturatingDivide(const float32_t a, const float32_t b)
		{
			return a / b;
		}
		
		static inline bool equals(const float32_t a, const float32_t b, const float32_t tol = tolerance())
		{
			return fabs(a - b) < tol;
		}

		static inline bool lessThan(const float32_t a, const float32_t b, const float32_t tol = tolerance())
		{
			return a < b && !equals(a, b, tol);
		}

		static inline bool greaterThan(const float32_t a, const float32_t b, const float32_t tol = tolerance())
		{
			return a > b && !equals(a, b, tol);
		}

		static inline bool lessThanOrEqual(const float32_t a, const float32_t b, const float32_t tol = tolerance())
		{
			return a < b || equals(a, b, tol);
		}

		static inline bool greaterThanOrEqual(const float32_t a, const float32_t b, const float32_t tol = tolerance())
		{
			return a > b || equals(a, b, tol);
		}
		
	};

	template<> class NumberUtils<uint8_t>
	{
	public:
		static inline uint8_t tolerance()
		{
			return 0;
		}

		typedef float FloatType;
		typedef float RealFloatType;
		typedef int16_t SignedType;
		
		static inline uint8_t scale()
		{
			return std::numeric_limits<uint8_t>::max();
		}

		static inline bool isnan(const uint8_t a)
		{
			return false;
		}

		static inline uint8_t saturatingAdd(const uint8_t x, const uint8_t y)
		{
			uint8_t res = x + y;
			res |= -(res < x);

			return res;
		}

		static inline uint8_t saturatingSubtract(const uint8_t x, const uint8_t y)
		{
			uint8_t res = x - y;
			res &= -(res <= x);

			return res;
		}

		static inline uint8_t saturatingMultiply(const uint8_t x, const uint8_t y)
		{
			uint16_t res = (uint16_t)x * (uint16_t)y;

			uint8_t hi = res >> 8 * sizeof(uint8_t);
			uint8_t lo = (uint8_t)res;

			return lo | -!!hi;
		}

		static inline uint8_t saturatingDivide(const uint8_t x, const uint8_t y)
		{
			return x / y;
		}

		static inline bool equals(const uint8_t a, const uint8_t b, const uint8_t dummy = 0)
		{
			return a == b;
		}

		static inline bool lessThan(const uint8_t a, const uint8_t b, const uint8_t dummy = 0)
		{
			return a < b;
		}

		static inline bool greaterThan(const uint8_t a, const uint8_t b, const uint8_t dummy = 0)
		{
			return a > b;
		}

		static inline bool lessThanOrEqual(const uint8_t a, const uint8_t b, const uint8_t dummy = 0)
		{
			return a <= b;
		}

		static inline bool greaterThanOrEqual(const uint8_t a, const uint8_t b, const uint8_t dummy = 0)
		{
			return a >= b;
		}

	};

	template<> class NumberUtils<uint16_t>
	{
	public:
		static inline uint16_t tolerance()
		{
			return 0;
		}

		typedef float FloatType;
		typedef float RealFloatType;
		typedef int32_t SignedType;

		static inline uint16_t scale()
		{
			return std::numeric_limits<uint16_t>::max();
		}

		static inline bool isnan(const uint16_t a)
		{
			return false;
		}

		static inline uint16_t saturatingAdd(const uint16_t x, const uint16_t y)
		{
			uint16_t res = x + y;
			res |= -(res < x);

			return res;
		}

		static inline uint16_t saturatingSubtract(const uint16_t x, const uint16_t y)
		{
			uint16_t res = x - y;
			res &= -(res <= x);

			return res;
		}

		static inline uint16_t saturatingMultiply(const uint16_t x, const uint16_t y)
		{
			uint32_t res = (uint32_t)x * (uint32_t)y;

			uint16_t hi = res >> 8 * sizeof(uint16_t);
			uint16_t lo = (uint16_t)res;

			return lo | -!!hi;
		}

		static inline uint16_t saturatingDivide(const uint16_t x, const uint16_t y)
		{
			return x / y;
		}

		static inline bool equals(const uint16_t a, const uint16_t b, const uint16_t dummy = 0)
		{
			return a == b;
		}

		static inline bool lessThan(const uint16_t a, const uint16_t b, const uint16_t dummy = 0)
		{
			return a < b;
		}

		static inline bool greaterThan(const uint16_t a, const uint16_t b, const uint16_t dummy = 0)
		{
			return a > b;
		}

		static inline bool lessThanOrEqual(const uint16_t a, const uint16_t b, const uint16_t dummy = 0)
		{
			return a <= b;
		}

		static inline bool greaterThanOrEqual(const uint16_t a, const uint16_t b, const uint16_t dummy = 0)
		{
			return a >= b;
		}

	};

	template<> class NumberUtils<uint32_t>
	{
	public:
		static inline uint32_t tolerance()
		{
			return 0;
		}

		typedef float FloatType;
		typedef float RealFloatType;
		typedef int64_t SignedType;

		static inline uint32_t scale()
		{
			return std::numeric_limits<uint32_t>::max();
		}

		static inline bool isnan(const uint32_t a)
		{
			return false;
		}

		static inline uint32_t saturatingAdd(const uint32_t x, const uint32_t y)
		{
			uint32_t res = x + y;
			res |= -(res < x);

			return res;
		}

		static inline uint32_t saturatingSubtract(const uint32_t x, const uint32_t y)
		{
			uint32_t res = x - y;
			res &= -(res <= x);

			return res;
		}

		static inline uint32_t saturatingMultiply(const uint32_t x, const uint32_t y)
		{
			uint64_t res = (uint64_t)x * (uint64_t)y;

			uint32_t hi = res >> 8 * sizeof(uint32_t);
			uint32_t lo = (uint32_t)res;

			return lo | -!!hi;
		}

		static inline uint32_t saturatingDivide(const uint32_t x, const uint32_t y)
		{
			return x / y;
		}

		static inline bool equals(const uint32_t a, const uint32_t b, const uint32_t dummy = 0)
		{
			return a == b;
		}

		static inline bool lessThan(const uint32_t a, const uint32_t b, const uint32_t dummy = 0)
		{
			return a < b;
		}

		static inline bool greaterThan(const uint32_t a, const uint32_t b, const uint32_t dummy = 0)
		{
			return a > b;
		}

		static inline bool lessThanOrEqual(const uint32_t a, const uint32_t b, const uint32_t dummy = 0)
		{
			return a <= b;
		}

		static inline bool greaterThanOrEqual(const uint32_t a, const uint32_t b, const uint32_t dummy = 0)
		{
			return a >= b;
		}

	};

	template<> class NumberUtils<uint64_t>
	{
	public:
		static inline uint64_t tolerance()
		{
			return 0;
		}

		typedef float FloatType;
		typedef float RealFloatType;

		// NOTE: This is not strictly to the spec as sizeof(int64_t) == sizeof(uint64_t).
		typedef int64_t SignedType;

		static inline uint64_t scale()
		{
			return std::numeric_limits<uint64_t>::max();
		}

		static inline bool isnan(const uint64_t a)
		{
			return false;
		}

		static inline uint64_t saturatingAdd(const uint64_t x, const uint64_t y)
		{
			uint64_t res = x + y;
			res |= -(res < x);

			return res;
		}

		static inline uint64_t saturatingSubtract(const uint64_t x, const uint64_t y)
		{
			uint64_t res = x - y;
			res &= -(res <= x);

			return res;
		}

		static inline uint64_t saturatingMultiply(const uint64_t x, const uint64_t y)
		{
#if defined(__linux__) || defined(__APPLE__)
			__uint128_t  res = (__uint128_t)x * (__uint128_t)y;

			uint64_t hi = res >> 8 * sizeof(uint64_t);
			uint64_t lo = (uint64_t)res;
#elif defined(_WIN32)
			uint64_t lo, hi;
			lo = _umul128(x, y, &hi);
#else

#error numberutils.h not configured for this platform.

#endif

			return lo | -!!hi;
		}

		static inline uint64_t saturatingDivide(const uint64_t x, const uint64_t y)
		{
			return x / y;
		}

		static inline bool equals(const uint64_t a, const uint64_t b, const uint64_t dummy = 0)
		{
			return a == b;
		}

		static inline bool lessThan(const uint64_t a, const uint64_t b, const uint64_t dummy = 0)
		{
			return a < b;
		}

		static inline bool greaterThan(const uint64_t a, const uint64_t b, const uint64_t dummy = 0)
		{
			return a > b;
		}

		static inline bool lessThanOrEqual(const uint64_t a, const uint64_t b, const uint64_t dummy = 0)
		{
			return a <= b;
		}

		static inline bool greaterThanOrEqual(const uint64_t a, const uint64_t b, const uint64_t dummy = 0)
		{
			return a >= b;
		}

	};


	template<> class NumberUtils<int8_t>
	{
	public:
		static inline int8_t tolerance()
		{
			return 0;
		}

		typedef float FloatType;
		typedef float RealFloatType;
		typedef int8_t SignedType;

		static inline int8_t scale()
		{
			return std::numeric_limits<int8_t>::max();
		}

		static inline bool isnan(const int8_t a)
		{
			return false;
		}

		static inline int8_t saturatingAdd(const int8_t x, const int8_t y)
		{
			uint8_t ux = x;
			uint8_t uy = y;
			uint8_t res = ux + uy;

			// Calculate overflowed result. (Don't change the sign bit of ux)
			ux = (ux >> (8 * sizeof(uint8_t) - 1)) + std::numeric_limits<int8_t>::max();

			// Force compiler to use cmovns instruction
			if ((int8_t)((ux ^ uy) | ~(uy ^ res)) >= 0)
			{
				res = ux;
			}

			return res;
		}

		static inline int8_t saturatingSubtract(const int8_t x, const int8_t y)
		{
			uint8_t ux = x;
			uint8_t uy = y;
			uint8_t res = ux - uy;

			// Calculate overflowed result. (Don't change the sign bit of ux)
			ux = (ux >> (8 * sizeof(uint8_t) - 1)) + std::numeric_limits<int8_t>::max();

			// Force compiler to use cmovns instruction
			if ((int8_t)((ux ^ uy) & (ux ^ res)) < 0)
			{
				res = ux;
			}

			return res;
		}

		static inline int8_t saturatingMultiply(const int8_t x, const int8_t y)
		{
			int16_t res = (int16_t)x * (int16_t)y;
			uint8_t res2 = ((uint8_t)(x ^ y) >> (8 * sizeof(uint8_t) - 1)) + std::numeric_limits<int8_t>::max();

			int8_t hi = (res >> 8 * sizeof(uint8_t));
			int8_t lo = (int8_t)res;

			if (hi != (lo >> (8 * sizeof(uint8_t) - 1)))
				res = res2;

			return (int8_t)res;
		}

		static inline int8_t saturatingDivide(const int8_t x, const int8_t y)
		{
			// Only one way to overflow, so test for and prevent it.
			int8_t xx = x + !((y + 1) | ((uint8_t)x + std::numeric_limits<int8_t>::lowest()));

			return xx / y;
		}

		static inline bool equals(const int8_t a, const int8_t b, const int8_t dummy = 0)
		{
			return a == b;
		}

		static inline bool lessThan(const int8_t a, const int8_t b, const int8_t dummy = 0)
		{
			return a < b;
		}

		static inline bool greaterThan(const int8_t a, const int8_t b, const int8_t dummy = 0)
		{
			return a > b;
		}

		static inline bool lessThanOrEqual(const int8_t a, const int8_t b, const int8_t dummy = 0)
		{
			return a <= b;
		}

		static inline bool greaterThanOrEqual(const int8_t a, const int8_t b, const int8_t dummy = 0)
		{
			return a >= b;
		}

	};



	template<> class NumberUtils<int16_t>
	{
	public:
		static inline int16_t tolerance()
		{
			return 0;
		}

		typedef float FloatType;
		typedef float RealFloatType;
		typedef int16_t SignedType;

		static inline int16_t scale()
		{
			return std::numeric_limits<int16_t>::max();
		}

		static inline bool isnan(const int16_t a)
		{
			return false;
		}

		static inline int16_t saturatingAdd(const int16_t x, const int16_t y)
		{
			uint16_t ux = x;
			uint16_t uy = y;
			uint16_t res = ux + uy;

			// Calculate overflowed result. (Don't change the sign bit of ux)
			ux = (ux >> (8 * sizeof(uint16_t) - 1)) + std::numeric_limits<int16_t>::max();

			// Force compiler to use cmovns instruction
			if ((int16_t)((ux ^ uy) | ~(uy ^ res)) >= 0)
			{
				res = ux;
			}

			return res;
		}

		static inline int16_t saturatingSubtract(const int16_t x, const int16_t y)
		{
			uint16_t ux = x;
			uint16_t uy = y;
			uint16_t res = ux - uy;

			// Calculate overflowed result. (Don't change the sign bit of ux)
			ux = (ux >> (8 * sizeof(uint16_t) - 1)) + std::numeric_limits<int16_t>::max();

			// Force compiler to use cmovns instruction
			if ((int16_t)((ux ^ uy) & (ux ^ res)) < 0)
			{
				res = ux;
			}

			return res;
		}

		static inline int16_t saturatingMultiply(const int16_t x, const int16_t y)
		{
			int32_t res = (int32_t)x * (int32_t)y;
			uint16_t res2 = ((uint16_t)(x ^ y) >> (8 * sizeof(uint16_t) - 1)) + std::numeric_limits<int16_t>::max();

			int16_t hi = (res >> 8 * sizeof(uint16_t));
			int16_t lo = (int16_t)res;

			if (hi != (lo >> (8 * sizeof(uint16_t) - 1)))
				res = res2;

			return (int16_t)res;
		}

		static inline int16_t saturatingDivide(const int16_t x, const int16_t y)
		{
			// Only one way to overflow, so test for and prevent it.
			int16_t xx = x + !((y + 1) | ((uint16_t)x + std::numeric_limits<int16_t>::lowest()));

			return xx / y;
		}

		static inline bool equals(const int16_t a, const int16_t b, const int16_t dummy = 0)
		{
			return a == b;
		}

		static inline bool lessThan(const int16_t a, const int16_t b, const int16_t dummy = 0)
		{
			return a < b;
		}

		static inline bool greaterThan(const int16_t a, const int16_t b, const int16_t dummy = 0)
		{
			return a > b;
		}

		static inline bool lessThanOrEqual(const int16_t a, const int16_t b, const int16_t dummy = 0)
		{
			return a <= b;
		}

		static inline bool greaterThanOrEqual(const int16_t a, const int16_t b, const int16_t dummy = 0)
		{
			return a >= b;
		}

	};


	template<> class NumberUtils<int32_t>
	{
	public:
		static inline int32_t tolerance()
		{
			return 0;
		}

		typedef float FloatType;
		typedef float RealFloatType;
		typedef int32_t SignedType;
		
		static inline int32_t scale()
		{
			return std::numeric_limits<int32_t>::max();
		}

		static inline bool isnan(const int32_t a)
		{
			return false;
		}

		static inline int32_t saturatingAdd(const int32_t x, const int32_t y)
		{
			uint32_t ux = x;
			uint32_t uy = y;
			uint32_t res = ux + uy;

			// Calculate overflowed result. (Don't change the sign bit of ux)
			ux = (ux >> (8 * sizeof(uint32_t) - 1)) + std::numeric_limits<int32_t>::max();

			// Force compiler to use cmovns instruction
			if ((int32_t)((ux ^ uy) | ~(uy ^ res)) >= 0)
			{
				res = ux;
			}

			return res;
		}

		static inline int32_t saturatingSubtract(const int32_t x, const int32_t y)
		{
			uint32_t ux = x;
			uint32_t uy = y;
			uint32_t res = ux - uy;

			// Calculate overflowed result. (Don't change the sign bit of ux)
			ux = (ux >> (8 * sizeof(uint32_t) - 1)) + std::numeric_limits<int32_t>::max();

			// Force compiler to use cmovns instruction
			if ((int32_t)((ux ^ uy) & (ux ^ res)) < 0)
			{
				res = ux;
			}

			return res;
		}

		static inline int32_t saturatingMultiply(const int32_t x, const int32_t y)
		{
			int64_t res = (int64_t)x * (int64_t)y;
			uint32_t res2 = ((uint32_t)(x ^ y) >> (8 * sizeof(uint32_t) - 1)) + std::numeric_limits<int32_t>::max();

			int32_t hi = (res >> 8 * sizeof(uint32_t));
			int32_t lo = (int32_t)res;

			if (hi != (lo >> (8 * sizeof(uint32_t) - 1)))
				res = res2;

			return (int32_t)res;
		}

		static inline int32_t saturatingDivide(const int32_t x, const int32_t y)
		{
			// Only one way to overflow, so test for and prevent it.
			int32_t xx = x + !((y + 1) | ((uint32_t)x + std::numeric_limits<int32_t>::lowest()));

			return xx / y;
		}

		static inline bool equals(const int32_t a, const int32_t b, const int32_t dummy = 0)
		{
			return a == b;
		}

		static inline bool lessThan(const int32_t a, const int32_t b, const int32_t dummy = 0)
		{
			return a < b;
		}

		static inline bool greaterThan(const int32_t a, const int32_t b, const int32_t dummy = 0)
		{
			return a > b;
		}

		static inline bool lessThanOrEqual(const int32_t a, const int32_t b, const int32_t dummy = 0)
		{
			return a <= b;
		}

		static inline bool greaterThanOrEqual(const int32_t a, const int32_t b, const int32_t dummy = 0)
		{
			return a >= b;
		}

	};

	template<> class NumberUtils<int64_t>
	{
	public:
		static inline int64_t tolerance()
		{
			return 0;
		}

		typedef float FloatType;
		typedef float RealFloatType;
		typedef int64_t SignedType;

		static inline int64_t scale()
		{
			return std::numeric_limits<int64_t>::max();
		}
		
		static inline bool isnan(const int64_t a)
		{
			return false;
		}

		static inline int64_t saturatingAdd(const int64_t x, const int64_t y)
		{
			uint64_t ux = x;
			uint64_t uy = y;
			uint64_t res = ux + uy;

			// Calculate overflowed result. (Don't change the sign bit of ux)
			ux = (ux >> (8 * sizeof(uint64_t) - 1)) + std::numeric_limits<int64_t>::max();

			// Force compiler to use cmovns instruction
			if ((int64_t)((ux ^ uy) | ~(uy ^ res)) >= 0)
			{
				res = ux;
			}

			return res;
		}

		static inline int64_t saturatingSubtract(const int64_t x, const int64_t y)
		{
			uint64_t ux = x;
			uint64_t uy = y;
			uint64_t res = ux - uy;

			// Calculate overflowed result. (Don't change the sign bit of ux)
			ux = (ux >> (8 * sizeof(uint64_t) - 1)) + std::numeric_limits<int64_t>::max();

			// Force compiler to use cmovns instruction
			if ((int64_t)((ux ^ uy) & (ux ^ res)) < 0)
			{
				res = ux;
			}

			return res;
		}

		static inline int64_t saturatingMultiply(const int64_t x, const int64_t y)
		{

#if defined(__linux__) || defined(__APPLE__)
			__int128_t  res = (__int128_t)x * (__int128_t)y;
			int64_t hi = (res >> 8 * sizeof(uint64_t));
			int64_t lo = (int64_t)res;
#elif defined(_WIN32)
			int64_t lo, hi;
			lo = _mul128(x, y, &hi);
#else

#error numberutils.h not configured for this platform.

#endif

			uint64_t res2 = ((uint64_t)(x ^ y) >> (8 * sizeof(uint64_t) - 1)) + std::numeric_limits<int64_t>::max();

			if (hi != (lo >> (8 * sizeof(uint64_t) - 1)))
				lo = res2;
				//return (int64_t)res2;
				//res = res2;

			//return (int64_t)res;
			return lo;
		}

		static inline int64_t saturatingDivide(const int64_t x, const int64_t y)
		{
			// Only one way to overflow, so test for and prevent it.
			int64_t xx = x + !((y + 1) | ((uint64_t)x + std::numeric_limits<int64_t>::lowest()));

			return xx / y;
		}

		static inline bool equals(const int64_t a, const int64_t b, const int64_t dummy = 0)
		{
			return a == b;
		}

		static inline bool lessThan(const int64_t a, const int64_t b, const int64_t dummy = 0)
		{
			return a < b;
		}

		static inline bool greaterThan(const int64_t a, const int64_t b, const int64_t dummy = 0)
		{
			return a > b;
		}

		static inline bool lessThanOrEqual(const int64_t a, const int64_t b, const int64_t dummy = 0)
		{
			return a <= b;
		}

		static inline bool greaterThanOrEqual(const int64_t a, const int64_t b, const int64_t dummy = 0)
		{
			return a >= b;
		}

	};

	template<> class NumberUtils<complex32_t>
	{
	public:
		static inline complex32_t tolerance()
		{
			return 0;
		}

		typedef complex32_t FloatType;
		typedef float RealFloatType;
		typedef complex32_t SignedType;

		static inline bool isnan(const complex32_t a)
		{
			return false;
		}

		/**
		Saturating addition.
		*/
		static inline complex32_t saturatingAdd(const complex32_t a, const complex32_t b)
		{
			return a + b;
		}

		/**
		Saturating subtraction.
		*/
		static inline complex32_t saturatingSubtract(const complex32_t a, const complex32_t b)
		{
			return a - b;
		}

		/**
		Saturating multiplication.
		*/
		static inline complex32_t saturatingMultiply(const complex32_t a, const complex32_t b)
		{
			return a * b;
		}

		/**
		Saturating division.
		*/
		static inline complex32_t saturatingDivide(const complex32_t a, const complex32_t b)
		{
			return a / b;
		}

		static inline bool equals(const complex32_t a, const complex32_t b)
		{
			return a == b;
		}

	};

	namespace tests
	{
		void equals();
		void saturatingArithmetic();
	}
}

