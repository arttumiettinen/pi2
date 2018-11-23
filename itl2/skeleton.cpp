
#include "hybridskeleton.h"
#include "lineskeleton.h"

#include "io/raw.h"
#include "pointprocess.h"
#include "projections.h"

namespace itl2
{
	namespace tests
	{
		void hybridSkeleton()
		{
			Image<uint8_t> head;
			raw::readd(head, "t1-head_bin_256x256x129.raw");

			hybridSkeleton(head);

			raw::writed(head, "./skeleton/head_hybrid_skeleton");

			Image<uint8_t> gt;
			raw::readd(gt, "./fiji_hybrid_skeleton_256x256x129.raw");

			subtract(head, gt);
			double diff = max(head);

			testAssert(diff == 0, "skeleton compared to Fiji skeleton");
		}

		void lineSkeleton()
		{
			// NOTE: No asserts!

			Image<uint8_t> head;
			raw::readd(head, "t1-head_bin_256x256x129.raw");
			
			lineSkeleton(head);

			raw::writed(head, "./skeleton/head_line_skeleton");
		}

		//void cavities()
		//{
		//	Image<uint8_t> head;
		//	raw::readd(head, "./skeleton/cavities/geometry_with_cavity_100x100x100.raw");

		//	hybridSkeleton(head);

		//	raw::writed(head, "./skeleton/cavities/cavity_hybrid_skeleton");
		//}
	}
}