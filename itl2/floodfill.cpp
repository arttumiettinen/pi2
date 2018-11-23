
#include "floodfill.h"
#include "io/raw.h"

using namespace math;

namespace itl2
{
	namespace tests
	{

		void floodfill()
		{
			// NOTE: No asserts!

			Image<uint8_t> head;
			raw::readd(head, "t1-head_bin_256x256x129.raw");

			itl2::floodfill(head, Vec3c(110, 110, 25), (uint8_t)128, (uint8_t)128);

			raw::writed(head, "./floodfill/filled");
		}

	}
}
