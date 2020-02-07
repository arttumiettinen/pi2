
#include "generation.h"
#include "io/itltiff.h"

namespace itl2
{
	namespace tests
	{
		void ellipsoid()
		{

			Vec3d u1(1, 1, 1);
			Vec3d u2(-1, 1, -1);
			Ellipsoid ell(Vec3d(50, 60, 70), Vec3d(30, 20, 10), u1, u2);

			AABox<double> bounds = ell.boundingBox();

			Image<uint8_t> result(100, 120, 140);

			draw(result, bounds, (uint8_t)128);
			draw(result, ell, (uint8_t)255);

			tiff::writed(result, "./generation/ellipsoid");
			
		}
	}
}