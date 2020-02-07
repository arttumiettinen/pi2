
#include "structure.h"
#include "io/raw.h"
#include "conversions.h"

namespace itl2
{
	namespace tests
	{

		void structureTensor()
		{
			// NOTE: No asserts!

			Image<float32_t> img, cyl, plan, energy, l1, l2, l3;
			raw::read(img, "./input_data/simple_structures_128x128x128.raw");

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
			raw::read(img, "./input_data/simple_structures_128x128x128.raw");

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

		void canny()
		{
			Image<uint16_t> head16;
			raw::read(head16, "./input_data/t1-head_256x256x129.raw");

			//// 'Manual' Canny edge detection logic
			//Image<float32_t> head;
			//convert(head16, head);

			//// 1. Pre-smoothing
			////gaussFilter(head, 1.0, BoundaryCondition::Nearest);

			//// 2. Gradient
			//Image<float32_t> dx, dy, dz;
			//gradient(head, dx, dy, dz, 1.0);

			//// 3. Non-maximum suppression - find local maxima of gradient
			//Image<float32_t> out;
			//nonMaximumSuppression(dx, dy, dz, out);
			//raw::writed(out, "canny/non-maximum_suppression");

			//// 4. Dual threshold - edges with gradient value above upper threshold are "surely" edges,
			//// and edges with gradient value between lower and upper threshold are edges only if
			//// they touch "sure" edge.
			//vector<float32_t> th = { 10, 100 };
			//multiThreshold(out, th);
			//raw::writed(out, "canny/dual_threshold");

			//// 5. Edge tracking - Convert all those edges to "sure" that touch a "sure" edge.
			//grow(out, 2.0f, 1.0f);
			//raw::writed(out, "canny/edge_tracking");
			//
			//// 6. Remove all other non-sure edges.
			//threshold(out, 1);
			//raw::writed(out, "canny/final");

			canny(head16, 1.0, 10.0, 100.0);
			raw::writed(head16, "canny/final_auto");
		}
	}
}