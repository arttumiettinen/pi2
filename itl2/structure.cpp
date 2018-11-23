
#include "structure.h"
#include "io/raw.h"

namespace itl2
{
	namespace tests
	{
		void curvature()
		{
			// NOTE: No asserts!

			Image<float32_t> img, curv;
			raw::readd(img, "./simple_structures_256x256x256.raw");

			meanCurvature(img, 10, curv);

			raw::writed(curv, "./structure/mean_curvature");
		}

		void structureTensor()
		{
			// NOTE: No asserts!

			Image<float32_t> img, cyl, plan, energy, l1, l2, l3;
			raw::readd(img, "./simple_structures_128x128x128.raw");

			itl2::structureTensor<float32_t>(img, 3, 3, &l1, &l2, &l3, 0, 0, 0, 0, 0, 0, &cyl, &plan, &energy);

			raw::writed(l1, "./structure/l1");
			raw::writed(l2, "./structure/l2");
			raw::writed(l3, "./structure/l3");
			raw::writed(cyl, "./structure/cylindricality");
			raw::writed(plan, "./structure/planarity");
			raw::writed(energy, "./structure/energy");
		}

		void lineFilter()
		{

			// NOTE: No asserts!
			Image<float32_t> img, lambda123, Vo;
			raw::readd(img, "./simple_structures_128x128x128.raw");

			//multiply(img, 255);

			itl2::lineFilter<float32_t>(img, 3, &lambda123, 1.0, 1.0, 0.5, &Vo, 0.25, 0.5, 0.5);

			raw::writed(lambda123, "./structure/lambda123_32bit");
			raw::writed(Vo, "./structure/Vo_32bit");
			
			Image<uint8_t> img8, lambda1238, Vo8;
			setValue(img8, img);
			itl2::lineFilter<uint8_t>(img8, 3, &lambda1238, 1.0, 1.0, 0.5, &Vo8, 0.25, 0.5, 0.5);

			raw::writed(lambda1238, "./structure/lambda123_8bit");
			raw::writed(Vo8, "./structure/Vo_8bit");
		}
	}
}