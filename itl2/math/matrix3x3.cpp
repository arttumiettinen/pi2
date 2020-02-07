
#include "matrix3x3.h"
#include "test.h"

namespace itl2
{
	namespace tests
	{
		void matrix3x3()
		{
			Matrix3x3f mat(1, 2, 3, 2, 4, 5, 3, 5, 9);

			Vec3f v1, v2, v3;
			float32_t l1, l2, l3;
			mat.eigsym(v1, v2, v3, l1, l2, l3);

			itl2::testAssert(NumberUtils<float32_t>::equals(l1, 13.1582f, 0.0001f), "lambda1");
			itl2::testAssert(NumberUtils<float32_t>::equals(l2, 0.924002f, 0.0001f), "lambda2");
			itl2::testAssert(NumberUtils<float32_t>::equals(l3, -0.0822487f, 0.0001f), "lambda3");

			itl2::testAssert(NumberUtils<float32_t>::equals(v1.x / v1.z, 0.349095f, 0.0001f), "v1");
			itl2::testAssert(NumberUtils<float32_t>::equals(v1.y / v1.z, 0.622192f, 0.0001f), "v1");

			itl2::testAssert(NumberUtils<float32_t>::equals(v2.x / v2.z, -0.204981f, 0.0001f), "v2");
			itl2::testAssert(NumberUtils<float32_t>::equals(v2.y / v2.z, -1.49221f, 0.0001f), "v2");

			itl2::testAssert(NumberUtils<float32_t>::equals(v3.x / v3.z, -5.37488f, 0.0001f), "v3");
			itl2::testAssert(NumberUtils<float32_t>::equals(v3.y / v3.z, 1.40848f, 0.0001f), "v3");
		}
	}
}