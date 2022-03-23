
#include "aabox.h"
#include "test.h"
#include "utilities.h"
#include <string>

using namespace std;

namespace itl2
{

	namespace tests
	{
		void boxTest(const AABoxc& a, const AABoxc& b, bool overlapResult, int testIndex)
		{
			testAssert(a.overlaps(b) == overlapResult, string("overlap test left-right ") + toString(testIndex));
			testAssert(b.overlaps(a) == overlapResult, string("overlap test right-left ") + toString(testIndex));
		}

		void aabox()
		{
			AABoxc l = AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(1, 1, 1));

			boxTest(AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(10, 10, 10)), AABoxc::fromMinMax(Vec3c(5, 0, 0), Vec3c(15, 10, 10)), true, 1);
			boxTest(AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(10, 10, 10)), AABoxc::fromMinMax(Vec3c(10, 0, 0), Vec3c(15, 10, 10)), false, 2);
			boxTest(AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(10, 10, 10)), AABoxc::fromMinMax(Vec3c(11, 0, 0), Vec3c(15, 10, 10)), false, 3);
			boxTest(AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(10, 10, 10)), AABoxc::fromMinMax(Vec3c(-10, 0, 0), Vec3c(15, 10, 10)), true, 4);
			boxTest(AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(10, 10, 10)), AABoxc::fromMinMax(Vec3c(-10, 0, 0), Vec3c(5, 10, 10)), true, 5);
			boxTest(AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(10, 10, 10)), AABoxc::fromMinMax(Vec3c(-10, 0, 0), Vec3c(0, 10, 10)), false, 6);
			boxTest(AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(10, 10, 10)), AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(10, 10, 10)), true, 7);
			boxTest(AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(10, 10, 10)), AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(5, 10, 10)), true, 8);
			boxTest(AABoxc::fromMinMax(Vec3c(0, 0, 0), Vec3c(10, 10, 10)), AABoxc::fromMinMax(Vec3c(-10, 0, 0), Vec3c(-5, 10, 10)), false, 9);
		}
	}
}