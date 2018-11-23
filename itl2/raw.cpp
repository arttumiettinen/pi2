
#include "io/raw.h"
#include "pointprocess.h"
#include "projections.h"
#include "conversions.h"

namespace itl2
{
	namespace raw
	{
		namespace tests
		{
			void raw()
			{
				Vec3c dim;
				ImageDataType dt;

				testAssert(internals::parseDimensions("filexx_ddd_10x20.raw", dim, dt) == true, "parseDimension return value");
				testAssert(dim == Vec3c(10, 20, 1), "dimensions 2D");
				testAssert(dt == Unknown, "data type");

				testAssert(internals::parseDimensions("filexx_ddd_100x2000x3456.raw", dim, dt) == true, "parseDimension return value");
				testAssert(dim == Vec3c(100, 2000, 3456), "dimensions");
				testAssert(dt == Unknown, "data type");
			}


			template<typename pixel_t> void printDifference(const Image<pixel_t>& a, const Image<pixel_t>& b)
			{
				a.checkSize(b);

				Image<float32_t> af(a.dimensions());
				convert(a, af);

				Image<float32_t> bf(b.dimensions());
				convert(b, bf);

				subtract(af, bf);
				abs(af);
				double diff = sum(af);

				cout << "|a-b| = " << diff << endl;

				testAssert(NumberUtils<double>::equals(diff, 0.0), "difference");
			}

			void readWriteBlock()
			{
				Image<uint16_t> head(256, 256, 129);
				raw::read(head, "./t1-head_256x256x129.raw");

				//raw::writed(head, "./raw/head");

				Vec3c outputDimensions = round(1.5 * Vec3d(head.dimensions()));
				string outFile = raw::internals::concatDimensions("./raw/head_3D_montage", outputDimensions);

				Vec3c blockStart(50, 50, 0);
				Vec3c blockSize = head.dimensions() - blockStart - Vec3c(0, 10, 0);

				for (coord_t z = 0; z < 4; z++)
				{
					for (coord_t y = 0; y < 4; y++)
					{
						for (coord_t x = 0; x < 3; x++)
						{
							Vec3c pos(x * blockSize.x, y * blockSize.y, z * blockSize.z);
							raw::writeBlock(head, outFile, pos, outputDimensions, blockStart, blockSize, true);
						}
					}
				}

			}
		}
	}
}