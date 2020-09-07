
#include "csa.h"
#include "io/raw.h"
#include "structure.h"
#include "conversions.h"

namespace itl2
{
	namespace tests
	{
		void csa()
		{
			// Create original
			Image<uint8_t> orig(100, 100, 100);
			draw(orig, Capsule<double>(Vec3d(10, 10, 10), Vec3d(70, 70, 70), 5), (uint8_t)255);
			draw(orig, Capsule<double>(Vec3d(20, 80, 15), Vec3d(50, 50, 70), 8), (uint8_t)255);

			// Calculate orientation
			Image<float32_t> energy, phi, theta;
			//convert(orig, energy);
			//cylinderOrientation(energy, phi, theta, 3, 3);
			raw::read(energy, "./csa/energy");
			raw::read(phi, "./csa/phi");
			raw::read(theta, "./csa/theta");

			// Calculate length
			// TODO

			// Calculate CSA
			Results results;
			Image<uint8_t> slices;
			Image<uint8_t> vis;
			convert(orig, vis);
			itl2::csa<uint8_t>(orig, energy, phi, theta, nullptr, allCrossSectionAnalyzers<uint8_t>(), results, 20, 300, 123, &slices, nullptr, &vis);

			// Save results
			raw::writed(orig, "./csa/orig");
			raw::writed(energy, "./csa/energy");
			raw::writed(phi, "./csa/phi");
			raw::writed(theta, "./csa/theta");
			raw::writed(slices, "./csa/slices");
			raw::writed(vis, "./csa/vis");
		}
	}
}