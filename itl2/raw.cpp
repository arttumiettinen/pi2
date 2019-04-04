
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

				testAssert(raw::getInfo("nonexisting_filexx_ddd_10x20.raw", dim, dt) == false, "parseDimension return value");

				testAssert(raw::getInfo("nonexisting_filexx_ddd_100x2000x3456.raw", dim, dt) == false, "parseDimension return value");
			}

			void writeBlock()
			{
				Image<uint16_t> head(256, 256, 129);
				raw::read(head, "./t1-head_256x256x129.raw");

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

			void writeBlockFast()
			{
				Image<uint16_t> head(256, 256, 129);
				raw::read(head, "./t1-head_256x256x129.raw");

				Vec3c outputDimensions = head.dimensions();
				string outFile = raw::internals::concatDimensions("./raw/head_fast_blocks", outputDimensions);

				for (coord_t z = 0; z < 7; z++)
				{
					Vec3c blockSize = head.dimensions();
					blockSize.z = 20;

					Vec3c pos(0, 0, z * blockSize.z);
					
					if (pos.z + blockSize.z > head.dimensions().z)
						blockSize.z = head.dimensions().z - pos.z;
					raw::writeBlock(head, outFile, pos, outputDimensions, pos, blockSize, true);
				}

				Image<uint16_t> headBlocks(256, 256, 129);
				raw::read(headBlocks, outFile);

				testAssert(equals(head, headBlocks), "by block-written image differs from original");
			}

			void expandFilename()
			{
				string f1 = "t1-head";
				string f2 = "t1-head_bin";
				string f3 = "t1-head_bin_dmap";
				string f4 = "t1-head_";

				raw::internals::expandRawFilename(f1);
				raw::internals::expandRawFilename(f2);
				raw::internals::expandRawFilename(f3);
				raw::internals::expandRawFilename(f4);

				testAssert(f1 == ".\\t1-head_256x256x129.raw", "filename 1");
				testAssert(f2 == ".\\t1-head_bin_256x256x129.raw", "filename 2");
				testAssert(f3 == ".\\t1-head_bin_dmap_256x256x129.raw", "filename 3");
				testAssert(f4 == ".\\t1-head_256x256x129.raw", "filename 4");
			}
		}
	}
}