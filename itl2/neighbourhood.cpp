
#include "neighbourhood.h"

#include "io/raw.h"

using namespace math;

namespace itl2
{
	namespace tests
	{
		void neighbourhoodTools()
		{
			// TODO: Is there any way to add asserts here?

			Image<uint8_t> mask;
			createNeighbourhoodMask(NeighbourhoodType::Rectangular, Vec3c(1, 2, 1), mask);
			raw::writed(mask, "./neighbourhood/rect_mask");

			createNeighbourhoodMask(NeighbourhoodType::Ellipsoidal, Vec3c(1, 2, 1), mask);
			raw::writed(mask, "./neighbourhood/ell_mask");

			createNeighbourhoodMask(NeighbourhoodType::Ellipsoidal, Vec3c(10, 20, 5), mask);
			raw::writed(mask, "./neighbourhood/ell_mask");

			createNeighbourhoodMask(NeighbourhoodType::Rectangular, Vec3c(10, 20, 0), mask);
			raw::writed(mask, "./neighbourhood/rect_mask_2D");

			createNeighbourhoodMask(NeighbourhoodType::Ellipsoidal, Vec3c(10, 20, 0), mask);
			raw::writed(mask, "./neighbourhood/ell_mask_2D");
		}
	}
}