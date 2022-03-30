
#include "aabox.h"
#include "test.h"
#include "utilities.h"
#include <string>

using namespace std;

namespace itl2
{

	namespace tests
	{
		void boxTestExclusive(const AABoxc& a, const AABoxc& b, bool overlapResult, int testIndex)
		{
			testAssert(a.overlapsExclusive(b) == overlapResult, string("exclusive overlap test left-right ") + toString(testIndex));
			testAssert(b.overlapsExclusive(a) == overlapResult, string("exclusive overlap test right-left ") + toString(testIndex));
		}

		void boxTestInclusive(const AABoxc& a, const AABoxc& b, bool overlapResult, int testIndex)
		{
			testAssert(a.overlapsInclusive(b) == overlapResult, string("inclusive overlap test left-right ") + toString(testIndex));
			testAssert(b.overlapsInclusive(a) == overlapResult, string("inclusive overlap test right-left ") + toString(testIndex));
		}

		void aabox()
		{
			AABoxc l = AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(1, 1, 1));

			boxTestExclusive(AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(10, 10, 10)), AABoxc::fromMinMax(Vec3c(5, 0, 0), Vec3c(15, 10, 10)), true, 1);
			boxTestExclusive(AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(10, 10, 10)), AABoxc::fromMinMax(Vec3c(10, 0, 0), Vec3c(15, 10, 10)), false, 2);
			boxTestExclusive(AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(10, 10, 10)), AABoxc::fromMinMax(Vec3c(11, 0, 0), Vec3c(15, 10, 10)), false, 3);
			boxTestExclusive(AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(10, 10, 10)), AABoxc::fromMinMax(Vec3c(-10, 0, 0), Vec3c(15, 10, 10)), true, 4);
			boxTestExclusive(AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(10, 10, 10)), AABoxc::fromMinMax(Vec3c(-10, 0, 0), Vec3c(5, 10, 10)), true, 5);
			boxTestExclusive(AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(10, 10, 10)), AABoxc::fromMinMax(Vec3c(-10, 0, 0), Vec3c(0, 10, 10)), false, 6);
			boxTestExclusive(AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(10, 10, 10)), AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(10, 10, 10)), true, 7);
			boxTestExclusive(AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(10, 10, 10)), AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(5, 10, 10)), true, 8);
			boxTestExclusive(AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(10, 10, 10)), AABoxc::fromMinMax(Vec3c(-10, 0, 0), Vec3c(-5, 10, 10)), false, 9);

			boxTestInclusive(AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(10, 10, 10)), AABoxc::fromMinMax(Vec3c(5, 0, 0), Vec3c(15, 10, 10)), true, 1);
			boxTestInclusive(AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(10, 10, 10)), AABoxc::fromMinMax(Vec3c(10, 0, 0), Vec3c(15, 10, 10)), true, 2);
			boxTestInclusive(AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(10, 10, 10)), AABoxc::fromMinMax(Vec3c(11, 0, 0), Vec3c(15, 10, 10)), false, 3);
			boxTestInclusive(AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(10, 10, 10)), AABoxc::fromMinMax(Vec3c(-10, 0, 0), Vec3c(15, 10, 10)), true, 4);
			boxTestInclusive(AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(10, 10, 10)), AABoxc::fromMinMax(Vec3c(-10, 0, 0), Vec3c(5, 10, 10)), true, 5);
			boxTestInclusive(AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(10, 10, 10)), AABoxc::fromMinMax(Vec3c(-10, 0, 0), Vec3c(0, 10, 10)), true, 6);
			boxTestInclusive(AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(10, 10, 10)), AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(10, 10, 10)), true, 7);
			boxTestInclusive(AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(10, 10, 10)), AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(5, 10, 10)), true, 8);
			boxTestInclusive(AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(10, 10, 10)), AABoxc::fromMinMax(Vec3c(-10, 0, 0), Vec3c(-5, 10, 10)), false, 9);
		}
	}
}