
#include "misc.h"
#include "io/raw.h"

namespace itl2
{
	namespace tests
	{
		void edges()
		{
			Image<uint8_t> img(10);
			setEdges(img, 255);

			testAssert(img(0) == 255 && img(img.width() - 1) == 255, "1D edges");
			for(coord_t x = 1; x < img.width()-1; x++)
				testAssert(img(x) == 0, "1D center");

			Image<uint8_t> img2(10, 10);
			setEdges(img2, 255);

			raw::writed(img2, "./misc/edges");

			// NOTE: No assert :(

			Image<uint8_t> img3(10, 10, 10);
			setEdges(img3, 255);

			raw::writed(img3, "./misc/edges");

			// NOTE: No assert :(


			Image<uint8_t> img4(10, 10, 10);
			setEdges(img4, 128, 3);

			raw::writed(img4, "./misc/thick_edges");

			// NOTE: No assert :(

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