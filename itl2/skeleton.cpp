
#include "surfaceskeleton.h"
#include "surfaceskeleton2.h"
#include "lineskeleton.h"

#include "io/raw.h"
#include "pointprocess.h"
#include "projections.h"

namespace itl2
{
	namespace experimental
	{
		namespace tests
		{
			void surfaceSkeleton2()
			{
				// NOTE: No asserts!

				Image<uint8_t> head;
				raw::read(head, "./input_data/t1-head_bin_256x256x129.raw");

				experimental::surfaceSkeleton2(head);

				raw::writed(head, "./skeleton/head_surface_skeleton2");

				//Image<uint8_t> img;
				//raw::read(img, "./skeleton/in/planes");

				//surfaceSkeleton2(img);

				//raw::writed(img, "./skeleton/skele");
			}
		}
	}

	namespace tests
	{

		void surfaceSkeleton()
		{
			// NOTE: No asserts!

			Image<uint8_t> img;
			raw::read(img, "./input_data/t1-head_bin_256x256x129.raw");
			//raw::read(img, "./skeleton/in/planes");

			surfaceSkeleton(img);

			raw::writed(img, "./skeleton/surface_skeleton");

			for (size_t n = 0; n < 10; n++)
			{
				thin(img, false);
				raw::writed(img, "./skeleton/surface_skeleton_to_line_skeleton_iteration_" + toString(n));
			}
		}

		void lineSkeleton()
		{
			// NOTE: No asserts!

			Image<uint8_t> head;
			raw::read(head, "./input_data/t1-head_bin_256x256x129.raw");
			
			lineSkeleton(head);

			raw::writed(head, "./skeleton/head_line_skeleton");
		}

		//void cavities()
		//{
		//	Image<uint8_t> head;
		//	raw::read(head, "./skeleton/cavities/geometry_with_cavity_100x100x100.raw");

		//	hybridSkeleton(head);

		//	raw::writed(head, "./skeleton/cavities/cavity_hybrid_skeleton");
		//}
	}
}