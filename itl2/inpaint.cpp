
#include "inpaint.h"
#include "io/raw.h"
#include "projections.h"

namespace itl2
{
	namespace tests
	{
		void inpaintNearest()
		{
			// NOTE: No asserts!

			Image<uint16_t> head;
			raw::read(head, "./t1-head_256x256x129.raw");

			inpaintNearest(head);

			raw::writed(head, "./inpaint/head_inpaint_nearest");
		}

		void inpaintGarcia()
		{
			// NOTE: No asserts!

			Image<uint8_t> img;
			raw::read(img, "./orig_20x20x1.raw");

			for (coord_t y = 5; y < 10; y++)
			{
				for (coord_t x = 5; x < 10; x++)
				{
					img(x, y) = 0;
				}
			}

			raw::writed(img, "./inpaint/data_missing");

			itl2::inpaintGarcia<uint8_t>(img, 0, true, 0.5);

			raw::writed(img, "./inpaint/inpaint_garcia");

		}

		void garciaToleranceTest(const string& imgname, float32_t tolerance)
		{
			Image<float32_t> img;
			Image<float32_t> inpaintedZeroTolerance;
			raw::read(img, "./" + imgname + "_53x53x53.raw");

			for (coord_t z = 0; z < img.depth(); z++)
			{
				for (coord_t y = 0; y < img.height(); y++)
				{
					for (coord_t x = 0; x < img.width(); x++)
					{
						if (x == 15 || x == 16 || x == 17 || y == 15 || y == 16 || y == 17 || z == 15 || z == 16 || z == 17)
							img(x, y, z) = 0;
					}
				}
			}

			
			raw::writed(img, "./inpaint/" + imgname + "_data_missing");

			setValue(inpaintedZeroTolerance, img);
			itl2::inpaintGarcia<float32_t>(inpaintedZeroTolerance, 0, true, 0);
			itl2::inpaintGarcia<float32_t>(img, 0, true, tolerance);

			raw::writed(img, "./inpaint/" + imgname + "_inpaint_garcia_tol" + toString(tolerance));

			subtract(img, inpaintedZeroTolerance);
			abs(img);
			double diff = max(img);

			cout << "Difference between zero tolerance and " << tolerance << " = " << diff << endl;
		}

		void inpaintGarcia2()
		{
			// NOTE: No asserts!

			for (float32_t tolerance = 0; tolerance < 1; tolerance += 0.1f)
			{
				garciaToleranceTest("dx", tolerance);
				garciaToleranceTest("dy", tolerance);
				garciaToleranceTest("dz", tolerance);
			}
			
		}
	}
}
