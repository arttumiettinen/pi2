
#include "misc.h"
#include "io/raw.h"
#include "iteration.h"
#include "testutils.h"

namespace itl2
{
	namespace tests
	{
		void edges()
		{
			Image<uint8_t> img(10);
			setEdges(img, 255);
			raw::writed(img, "./misc/edges1d");

			testAssert(img(0) == 255 && img(img.width() - 1) == 255, "1D edges");
			for(coord_t x = 1; x < img.width()-1; x++)
				testAssert(img(x) == 0, "1D center");


			Image<uint8_t> img_test(img.dimensions());
			forEdges(img_test.bounds(), img.dimensionality(), [&](coord_t x, coord_t y, coord_t z)
				{
					img_test(x, y, z) = 255;
				});
			raw::writed(img_test, "./misc/edges1d_test");
			checkDifference(img, img_test, "set edges of 1D image");




			Image<uint8_t> img2(10, 10);
			setEdges(img2, 255);

			raw::writed(img2, "./misc/edges2");

			Image<uint8_t> img2_test(img2.dimensions());
			forEdges(img2_test.bounds(), img2.dimensionality(), [&](coord_t x, coord_t y, coord_t z)
				{
					img2_test(x, y, z) = 255;
				});
			raw::writed(img2_test, "./misc/edges2_test");
			checkDifference(img2, img2_test, "set edges of 2D image");




			Image<uint8_t> img3(10, 10, 10);
			setEdges(img3, 255);

			raw::writed(img3, "./misc/edges3");

			Image<uint8_t> img3_test(img3.dimensions());
			forEdges(img3_test.bounds(), img3.dimensionality(), [&](coord_t x, coord_t y, coord_t z)
				{
					img3_test(x, y, z) = 255;
				});
			raw::writed(img3_test, "./misc/edges3_test");
			checkDifference(img3, img3_test, "set edges of 3D image");




			Image<uint8_t> img4(10, 10, 10);
			setEdges(img4, 128, 3);

			raw::writed(img4, "./misc/thick_edges");

			Image<uint8_t> img4_test(img4.dimensions());
			for (coord_t u = 0; u < 3; u++)
			{
				AABox bounds = img4_test.bounds();
				bounds.inflate(-u);
				forEdges(bounds, img4.dimensionality(), [&](coord_t x, coord_t y, coord_t z)
					{
						img4_test(x, y, z) = 128;
					});
			}
			raw::writed(img4_test, "./misc/thick_edges_test");
			checkDifference(img4, img4_test, "set edges of 3D image (thick)");

		}

		void normalizeZ()
		{
			// NOTE: No asserts!

			Image<uint8_t> img(10, 11, 12);

			for (coord_t z = 0; z < img.depth(); z++)
			{
				for (coord_t y = 0; y < img.height(); y++)
				{
					for (coord_t x = 0; x < img.width(); x++)
					{
						img(x, y, z) += (uint8_t)z;
					}
				}
			}

			raw::writed(img, "misc/normalizez_before");
			normalizeZ(img);
			raw::writed(img, "misc/normalizez_after");
		}
	}

}