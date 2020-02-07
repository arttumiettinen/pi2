
#include "io/vol.h"
#include "pointprocess.h"
#include "testutils.h"
#include "math/vec3.h"

namespace itl2
{

	namespace vol
	{
		namespace tests
		{
			void volio()
			{
				Vec3c dims;
				ImageDataType dt;
				string end;
				size_t hs;
				string reason;
				vol::getInfo("./input_data/simple_structures.vol", dims, dt, end, hs, reason);

				Image<uint8_t> img;
				vol::read(img, "./input_data/simple_structures.vol");

				Image<float32_t> gt;
				raw::read(gt, "./input_data/simple_structures_128x128x128.raw");
				multiply(gt, 255);

				itl2::checkDifference(img, gt, "Same structure from .vol and .raw file.");
			}
		}
	}

}
