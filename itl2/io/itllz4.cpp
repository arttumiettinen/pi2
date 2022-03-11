
#include "io/itllz4.h"
#include "projections.h"
#include "generation.h"
#include "transform.h"

namespace itl2
{
	namespace lz4
	{
		namespace internals
		{
			bool getInfo(std::ifstream& in, Vec3c& dimensions, ImageDataType& dataType, string& reason)
			{
				coord_t w = internals::readSafe<coord_t>(in);
				coord_t h = internals::readSafe<coord_t>(in);
				coord_t d = internals::readSafe<coord_t>(in);
				dimensions = Vec3c(w, h, d);

				if (w <= 0 || h <= 0 || d <= 0)
				{
					reason = "Invalid image dimensions.";
					return false;
				}

				dataType = (ImageDataType)internals::readSafe<int32_t>(in);
				return true;
			}
		}

		bool getInfo(const std::string& filename, Vec3c& dimensions, ImageDataType& dataType, string& reason)
		{
			std::ifstream in(filename.c_str(), std::ios_base::in | std::ios_base::binary);
			if (!in)
			{
				reason = std::string("Unable to open ") + filename + std::string(", ") + getStreamErrorMessage();
				return false;
			}

			return internals::getInfo(in, dimensions, dataType, reason);
		}

		bool isFile(const std::string& filename)
		{
			Vec3c dimensions;
			ImageDataType dataType;
			string reason;
			return getInfo(filename, dimensions, dataType, reason);
		}


		namespace tests
		{
			void singleLZ4Test(const Vec3c& dimensions)
			{
				Image<uint16_t> img(dimensions);
				ramp3(img);

				lz4::writed(img, "./lz4/image");

				Image<uint16_t> fromDisk;
				lz4::readd(fromDisk, "./lz4/image");

				testAssert(equals(img, fromDisk), "LZ4 compressed raw image");
			}

			void lz4io()
			{
				singleLZ4Test(Vec3c(1, 1, 1));
				singleLZ4Test(Vec3c(10, 10, 10));
				singleLZ4Test(Vec3c(100, 10, 10));
				singleLZ4Test(Vec3c(100, 100, 100));
				singleLZ4Test(Vec3c(100, 100, 100));

				singleLZ4Test(Vec3c(10, 10, 1));
				singleLZ4Test(Vec3c(10, 1, 10));
				singleLZ4Test(Vec3c(1, 10, 10));

				singleLZ4Test(Vec3c(1000, 1000, 1));
				singleLZ4Test(Vec3c(1000, 1, 1000));
				singleLZ4Test(Vec3c(1, 1000, 1000));
			}

			void lz4blockIo()
			{
				Vec3c dimensions(100, 200, 300);

				Image<uint16_t> img(dimensions);
				ramp3(img);

				Vec3c blockStart(10, 20, 30);
				Vec3c blockSize(50, 60, 70);

				// Write entire image and a block.
				lz4::write(img, "./lz4block/entire_image");
				lz4::writeBlock(img, "./lz4block/block", Vec3c(0, 0, 0), Vec3c(0, 0, 0), blockStart, blockSize);

				// Generate ground truth block by cropping the image.
				Image<uint16_t> gtBlock(blockSize);
				crop(img, gtBlock, blockStart);

				// Read entire file written using writeBlock and check against ground truth.
				Image<uint16_t> fileBlock;
				lz4::readd(fileBlock, "./lz4block/block");

				testAssert(equals(fileBlock, gtBlock), "LZ4raw writeBlock");


				// Read using readBlock and compare to the ground truth.
				Image<uint16_t> readBlockResult(blockSize);
				lz4::readBlock(readBlockResult, "./lz4block/entire_image", blockStart);
				testAssert(equals(readBlockResult, gtBlock), "LZ4raw readBlock");

				// Now write entire image in multiple blocks, read back, and compare to original
				fs::remove_all("./lz4_multiple_blocks_full_file");
				Vec3c blockCount = img.dimensions().componentwiseDivide(blockSize);
				for (coord_t z = 0; z < blockCount.z+1; z++)
				{
					for (coord_t y = 0; y < blockCount.y+1; y++)
					{
						for (coord_t x = 0; x < blockCount.x+1; x++)
						{
							Vec3c pos(x * blockSize.x, y * blockSize.y, z * blockSize.z);
							if (img.isInImage(pos))
							{
								Vec3c blockEnd = pos + blockSize;
								blockEnd = min(img.dimensions(), blockEnd);
								Vec3c realBlockSize = blockEnd - pos;
								lz4::writeBlock(img, "./lz4_multiple_blocks_full_file/image", pos, img.dimensions(), pos, realBlockSize);
							}
						}
					}
				}

				Image<uint16_t> multiBlockResult;
				lz4::read(multiBlockResult, "./lz4_multiple_blocks_full_file/image");

				raw::writed(multiBlockResult, "./lz4_multiple_blocks_result");

				testAssert(equals(img, multiBlockResult), "LZ4 file written in multiple blocks compared to the original.");
			}
		}
	}
}
