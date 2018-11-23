#pragma once

#include <limits>

#include "test.h"
#include "datatypes.h"

using itl2::complex32_t;
using std::numeric_limits;

namespace math
{
	template<typename T> class NumberUtils
	{
	public:
		static inline const T tolerance()
		{
			T::NumberUtils_undefined;
		}

		typedef double FloatType;
		typedef double RealFloatType;

		/**
		Returns typical maximum value of numbers of type T.
		(numeric_limits<T>::max() for integral types, 1.0 for floating point types)
		*/
		static inline T scale()
		{
			return 1;
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
		static inline const double tolerance()
		{
			return 1e-15;
		}
		
		typedef double FloatType;
		typedef double RealFloatType;

		static inline double scale()
		{
			return 1.0;
		}
		
		static inline bool isnan(const double a)
		{
			return ::isnan(a);
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

	template<> class NumberUtils<float>
	{
	public:
		static inline const float tolerance()
		{
			return 1e-10f;
		}

		typedef float FloatType;
		typedef float RealFloatType;
		
		static inline float scale()
		{
			return 1.0f;
		}

		static inline bool isnan(const float a)
		{
			return ::isnan(a);
		}
		
		static inline bool equals(const float a, const float b, const float tol = tolerance())
		{
			return fabs(a - b) < tol;
		}

		static inline bool lessThan(const float a, const float b, const float tol = tolerance())
		{
			return a < b && !equals(a, b, tol);
		}

		static inline bool greaterThan(const float a, const float b, const float tol = tolerance())
		{
			return a > b && !equals(a, b, tol);
		}

		static inline bool lessThanOrEqual(const float a, const float b, const float tol = tolerance())
		{
			return a < b || equals(a, b, tol);
		}

		static inline bool greaterThanOrEqual(const float a, const float b, const float tol = tolerance())
		{
			return a > b || equals(a, b, tol);
		}
		
	};

	template<> class NumberUtils<uint8_t>
	{
	public:
		static inline const uint8_t tolerance()
		{
			return 0;
		}

		typedef float FloatType;
		typedef float RealFloatType;
		
		static inline uint8_t scale()
		{
			return numeric_limits<uint8_t>::max();
		}

		static inline bool isnan(const uint8_t a)
		{
			return false;
		}

		static inline bool equals(const uint8_t a, const uint8_t b)
		{
			return a == b;
		}

		static inline bool lessThan(const uint8_t a, const uint8_t b)
		{
			return a < b;
		}

		static inline bool greaterThan(const uint8_t a, const uint8_t b)
		{
			return a > b;
		}

		static inline bool lessThanOrEqual(const uint8_t a, const uint8_t b)
		{
			return a <= b;
		}

		static inline bool greaterThanOrEqual(const uint8_t a, const uint8_t b)
		{
			return a >= b;
		}

	};

	template<> class NumberUtils<uint16_t>
	{
	public:
		static inline const uint16_t tolerance()
		{
			return 0;
		}

		typedef float FloatType;
		typedef float RealFloatType;

		static inline uint16_t scale()
		{
			return numeric_limits<uint16_t>::max();
		}

		static inline bool isnan(const uint16_t a)
		{
			return false;
		}

		static inline bool equals(const uint16_t a, const uint16_t b)
		{
			return a == b;
		}

		static inline bool lessThan(const uint16_t a, const uint16_t b)
		{
			return a < b;
		}

		static inline bool greaterThan(const uint16_t a, const uint16_t b)
		{
			return a > b;
		}

		static inline bool lessThanOrEqual(const uint16_t a, const uint16_t b)
		{
			return a <= b;
		}

		static inline bool greaterThanOrEqual(const uint16_t a, const uint16_t b)
		{
			return a >= b;
		}

	};

	template<> class NumberUtils<uint32_t>
	{
	public:
		static inline const uint32_t tolerance()
		{
			return 0;
		}

		typedef float FloatType;
		typedef float RealFloatType;

		static inline uint32_t scale()
		{
			return numeric_limits<uint32_t>::max();
		}

		static inline bool isnan(const uint32_t a)
		{
			return false;
		}

		static inline bool equals(const uint32_t a, const uint32_t b)
		{
			return a == b;
		}

		static inline bool lessThan(const uint32_t a, const uint32_t b)
		{
			return a < b;
		}

		static inline bool greaterThan(const uint32_t a, const uint32_t b)
		{
			return a > b;
		}

		static inline bool lessThanOrEqual(const uint32_t a, const uint32_t b)
		{
			return a <= b;
		}

		static inline bool greaterThanOrEqual(const uint32_t a, const uint32_t b)
		{
			return a >= b;
		}

	};

	template<> class NumberUtils<uint64_t>
	{
	public:
		static inline const uint64_t tolerance()
		{
			return 0;
		}

		typedef float FloatType;
		typedef float RealFloatType;

		static inline uint64_t scale()
		{
			return numeric_limits<uint64_t>::max();
		}

		static inline bool isnan(const uint64_t a)
		{
			return false;
		}

		static inline bool equals(const uint64_t a, const uint64_t b)
		{
			return a == b;
		}

		static inline bool lessThan(const uint64_t a, const uint64_t b)
		{
			return a < b;
		}

		static inline bool greaterThan(const uint64_t a, const uint64_t b)
		{
			return a > b;
		}

		static inline bool lessThanOrEqual(const uint64_t a, const uint64_t b)
		{
			return a <= b;
		}

		static inline bool greaterThanOrEqual(const uint64_t a, const uint64_t b)
		{
			return a >= b;
		}

	};

	template<> class NumberUtils<int>
	{
	public:
		static inline const int tolerance()
		{
			return 0;
		}

		typedef float FloatType;
		typedef float RealFloatType;
		
		static inline int scale()
		{
			return numeric_limits<int>::max();
		}

		static inline bool isnan(const int a)
		{
			return false;
		}

		static inline bool equals(const int a, const int b)
		{
			return a == b;
		}

		static inline bool lessThan(const int a, const int b)
		{
			return a < b;
		}

		static inline bool greaterThan(const int a, const int b)
		{
			return a > b;
		}

		static inline bool lessThanOrEqual(const int a, const int b)
		{
			return a <= b;
		}

		static inline bool greaterThanOrEqual(const int a, const int b)
		{
			return a >= b;
		}

	};

	template<> class NumberUtils<long long>
	{
	public:
		static inline const long long tolerance()
		{
			return 0;
		}

		typedef float FloatType;
		typedef float RealFloatType;

		static inline long long scale()
		{
			return numeric_limits<long long>::max();
		}
		
		static inline bool isnan(const long long a)
		{
			return false;
		}

		static inline bool equals(const long long a, const long long b)
		{
			return a == b;
		}

		static inline bool lessThan(const long long a, const long long b)
		{
			return a < b;
		}

		static inline bool greaterThan(const long long a, const long long b)
		{
			return a > b;
		}

		static inline bool lessThanOrEqual(const long long a, const long long b)
		{
			return a <= b;
		}

		static inline bool greaterThanOrEqual(const long long a, const long long b)
		{
			return a >= b;
		}

	};

	template<> class NumberUtils<complex32_t>
	{
	public:
		static inline const complex32_t tolerance()
		{
			return 0;
		}

		typedef complex32_t FloatType;
		typedef float RealFloatType;

		static inline bool isnan(const complex32_t a)
		{
			return false;
		}

		static inline bool equals(const complex32_t a, const complex32_t b)
		{
			return a == b;
		}

	};

	namespace tests
	{
		inline void numberUtils()
		{
			itl2::testAssert(NumberUtils<double>::equals(1.0, 1.0 + 0.5 * NumberUtils<double>::tolerance()), "Equals");
		}
	}
}

