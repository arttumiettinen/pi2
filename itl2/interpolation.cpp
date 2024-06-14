
#include "interpolation.h"
#include "pointprocess.h"
#include "transform.h"
#include "io/raw.h"

namespace itl2
{

	namespace tests
	{
		void invalidValueInterpolation()
		{
			//Image<uint16_t> input(20, 20, 1);
			//setValue(input, 100);
			//input(10, 10, 0) = 0;
			//input(11, 10, 0) = 0;
			
			Image<uint16_t> input;
			raw::read(input, "./interpolation/piece");
			
			raw::writed(input, "./interpolation/original");

			CubicInvalidValueInterpolator<float32_t, uint16_t, float32_t> interpCub = CubicInvalidValueInterpolator<float32_t, uint16_t, float32_t>(BoundaryCondition::Zero, 0, 0);
			LinearInvalidValueInterpolator<float32_t, uint16_t, float32_t> interpLin = LinearInvalidValueInterpolator<float32_t, uint16_t, float32_t>(BoundaryCondition::Zero, 0, 0);

			Vec3f point(10.1f, 10.0f, 10);
			float32_t rlin = interpLin(input, point);
			float32_t rcub = interpCub(input, point);

			std::cout << "Linear: " << point << " => " << rlin << std::endl;
			std::cout << "Cubic: " << point << " => " << rcub << std::endl;

			double angle = 1 / 180.0 * 3.14;

			Image<float32_t> outLin(input.dimensions());
			itl2::rotate(input, outLin, angle, Vec3d(0, 0, 1), interpLin);
			raw::writed(outLin, "./interpolation/linear");

			Image<float32_t> outCub(input.dimensions());
			itl2::rotate(input, outCub, angle, Vec3d(0, 0, 1), interpCub);
			raw::writed(outLin, "./interpolation/cubic");
		}
	}

}