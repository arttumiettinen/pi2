
#include "pathopening.h"
#include "generation.h"

namespace itl2
{

	namespace tests
	{
		void pathopening()
		{
			// Create original
			Image<uint8_t> orig(100, 100, 100);
			draw(orig, Capsule<double>(Vec3d(10, 10, 10), Vec3d(70, 70, 70), 5), (uint8_t)255);
			draw(orig, Capsule<double>(Vec3d(20, 80, 15), Vec3d(50, 50, 70), 8), (uint8_t)255);

			// Calculate
			Image<float32_t> length;
			itl2::pathLength2Binary3dNormalOrChamferMemorySave(orig, length);

			// Save results
			raw::writed(orig, "./pathopening/orig");
			raw::writed(length, "./pathopening/length");
		}
	}
}