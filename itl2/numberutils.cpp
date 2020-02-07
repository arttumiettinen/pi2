
#include "math/numberutils.h"
#include "test.h"

using namespace std;

namespace itl2
{
	namespace tests
	{
		void equals()
		{
			itl2::testAssert(NumberUtils<double>::equals(1.0, 1.0 + 0.5 * NumberUtils<double>::tolerance()), "Equals");
		}

		template<typename T> void testSaturatingUnsigned()
		{
			T M = numeric_limits<T>::max();
			T m = numeric_limits<T>::lowest();

			itl2::testAssert(NumberUtils<T>::saturatingAdd(M-10, M-15) == M, "saturating add");
			itl2::testAssert(NumberUtils<T>::saturatingAdd(M/2-5, M/2-5) == M/2+M/2-10, "add");
			itl2::testAssert(NumberUtils<T>::saturatingSubtract(M/2, M-10) == m, "saturating subtract");
			itl2::testAssert(NumberUtils<T>::saturatingSubtract(M / 2 + 10, M / 2) == 10, "subtract");
			itl2::testAssert(NumberUtils<T>::saturatingMultiply(M - 10, M - 15) == M, "saturating mul");
			itl2::testAssert(NumberUtils<T>::saturatingMultiply(M/3, 2) == M/3*2, "mul");

			itl2::testAssert(NumberUtils<T>::saturatingDivide(M, 10) == M/10, "div");
		}

		template<typename T> void testSaturatingSigned()
		{
			T M = numeric_limits<T>::max();
			T m = numeric_limits<T>::lowest();

			itl2::testAssert(NumberUtils<T>::saturatingAdd(M - 10, M - 15) == M, "saturating add");
			itl2::testAssert(NumberUtils<T>::saturatingAdd(-10, -M) == m, "saturating add negative");
			itl2::testAssert(NumberUtils<T>::saturatingAdd(M / 2 - 5, M / 2 - 5) == M / 2 + M / 2 - 10, "add");

			itl2::testAssert(NumberUtils<T>::saturatingSubtract(-1, M) == m, "saturating subtract");
			itl2::testAssert(NumberUtils<T>::saturatingSubtract(1, -M) == M, "saturating subtract negative");
			itl2::testAssert(NumberUtils<T>::saturatingSubtract(M / 2 + 10, M / 2) == 10, "subtract");

			itl2::testAssert(NumberUtils<T>::saturatingMultiply(M - 10, M - 15) == M, "saturating mul");
			itl2::testAssert(NumberUtils<T>::saturatingMultiply(M - 10, -M + 15) == m, "saturating mul");
			itl2::testAssert(NumberUtils<T>::saturatingMultiply(M / 3, 2) == M / 3 * 2, "mul");

			itl2::testAssert(NumberUtils<T>::saturatingDivide(m, -1) == M, "saturating div");
			itl2::testAssert(NumberUtils<T>::saturatingDivide(M, 10) == M / 10, "div");
		}

		void saturatingArithmetic()
		{
			testSaturatingUnsigned<uint8_t>();
			testSaturatingUnsigned<uint16_t>();
			testSaturatingUnsigned<uint32_t>();
			testSaturatingUnsigned<uint64_t>();

			testSaturatingSigned<int8_t>();
			testSaturatingSigned<int16_t>();
			testSaturatingSigned<int32_t>();
			testSaturatingSigned<int64_t>();

			// These test do not work as floating point math produces infinities instead of saturation (and that's how we want it to work!)
			//testSaturatingSigned<float32_t>();
			//testSaturatingSigned<double>();
		}
	}
}