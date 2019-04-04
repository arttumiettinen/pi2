
#include "transform.h"
#include "io/raw.h"
#include "math/vec3.h"
#include "pointprocess.h"

using namespace math;

namespace itl2
{
	namespace tests
	{
		void translate()
		{
			// NOTE: No asserts!

			Image<uint16_t> head16;
			raw::read(head16, "./t1-head_noisy_256x256x129.raw");

			add(head16, 1);

			Image<uint16_t> shifted(head16.dimensions());

			translate(head16, shifted, Vec3d(10.5, 0.7, -3.2), NearestNeighbourInterpolator<uint16_t, uint16_t>(BoundaryCondition::Zero));
			raw::writed(shifted, "./transform/shift_nearest");

			translate(head16, shifted, Vec3d(10.5, 0.7, -3.2), LinearInterpolator<uint16_t, uint16_t>(BoundaryCondition::Zero));
			raw::writed(shifted, "./transform/shift_linear");

			translate(head16, shifted, Vec3d(10.5, 0.7, -3.2), LinearInvalidValueInterpolator<uint16_t, uint16_t>(BoundaryCondition::Zero, 0, 0));
			raw::writed(shifted, "./transform/shift_linear_invalidvalue");

			translate(head16, shifted, Vec3d(10.5, 0.7, -3.2), CubicInterpolator<uint16_t, uint16_t>(BoundaryCondition::Zero));
			raw::writed(shifted, "./transform/shift_cubic");

			translate(head16, shifted, Vec3d(10.5, 0.7, -3.2), CubicInvalidValueInterpolator<uint16_t, uint16_t>(BoundaryCondition::Zero, 0, 0));
			raw::writed(shifted, "./transform/shift_cubic_invalidvalue");
		}

		void binning()
		{
			// NOTE: No asserts!

			Image<uint16_t> head16;
			raw::read(head16, "./t1-head_noisy_256x256x129.raw");

			Image<uint16_t> headb;
			binning(head16, headb, 2);
			raw::writed(headb, "./transform/binning_2");

			binning(head16, headb, 3);
			raw::writed(headb, "./transform/binning_3");

			binning(head16, headb, 4);
			raw::writed(headb, "./transform/binning_4");
		}

		void genericTransform()
		{
			// NOTE: No asserts!

			Image<uint16_t> head16, headb;
			raw::read(head16, "./t1-head_noisy_256x256x129.raw");

			vector<Vec3f> refPoints, defPoints;
			defPoints.push_back(Vec3f(0, 0, 64));
			refPoints.push_back(Vec3f(0, 50, 64));

			
			defPoints.push_back(Vec3f(256, 256, 64));
			refPoints.push_back(Vec3f(256, 256-50, 64));
			
			headb.ensureSize(head16.dimensions());
			itl2::genericTransform(head16, headb, Vec3c(0, 0, 0), refPoints, defPoints);

			raw::writed(headb, "./transform/generic_transform");
		}
	}
}
