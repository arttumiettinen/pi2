

#include "io/raw.h"
#include "generation.h"

#include "surfacecurvature.h"

using namespace std;

namespace itl2
{
	namespace tests
	{


		void surfaceCurvature()
		{
			Image<uint8_t> img;
			Image<float32_t> kappa1, kappa2;
			//raw::read(img, "./structure/simple_structures_128x128x128.raw");
			//raw::read(img, "./structure/simple_structures_20x20x20.raw");

			//img.ensureSize(Vec3c(20, 20, 20));
			//draw(img, Capsule(Vec3d(10, 10, 0), Vec3d(10, 10, 20), 10.0), (uint8_t)128);
			//draw(img, Capsule(Vec3d(10, 10, 0), Vec3d(10, 10, 20), 5.0), (uint8_t)0);
			//linearMap(img, Vec4d(0, 128, 128, 0));

			// Crossing bars
			img.ensureSize(Vec3c(40, 40, 40));
			draw(img, Capsule(Vec3d(20, 10, 0), Vec3d(20, 10, 40), 13.0), (uint8_t)128);
			draw(img, Capsule(Vec3d(0, 30, 20), Vec3d(40, 30, 20), 13.0), (uint8_t)128);
			// Add this to make the crossing bars hollow:
			draw(img, Capsule(Vec3d(20, 10, 0), Vec3d(20, 10, 40), 10.0), (uint8_t)0);
			draw(img, Capsule(Vec3d(0, 30, 20), Vec3d(40, 30, 20), 10.0), (uint8_t)0);

			// Very small sphere
			//img.ensureSize(Vec3c(40, 40, 40));
			//draw(img, Sphere(Vec3d(20, 20, 20), 3.0), (uint8_t)128);

			// Hollow cylinder
			//img.ensureSize(Vec3c(40, 40, 40));
			//draw(img, Capsule(Vec3d(20, 20, 0), Vec3d(20, 20, 40), 10.0), (uint8_t)128);
			//draw(img, Capsule(Vec3d(20, 20, 0), Vec3d(20, 20, 40), 5.0), (uint8_t)0);

			// Hollow sphere
			//img.ensureSize(Vec3c(40, 40, 40));
			//draw(img, Sphere(Vec3d(20, 20, 20), 10.0), (uint8_t)128);
			//draw(img, Sphere(Vec3d(20, 20, 20), 5.0), (uint8_t)0);

			raw::writed(img, "./surface_curvature/geometry");

			itl2::surfaceCurvature(img, 10, &kappa1, &kappa2, nullptr, nullptr, BoundaryCondition::Zero, numeric_limits<float32_t>::signaling_NaN());

			raw::writed(kappa1, "./surface_curvature/kappa1");
			raw::writed(kappa2, "./surface_curvature/kappa2");
		}
	}
}