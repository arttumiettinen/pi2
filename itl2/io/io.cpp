
#include "io.h"
#include "projections.h"

namespace itl2
{
	namespace io
	{
	
	
    	bool getInfo(const std::string& filename, Vec3c& dimensions, ImageDataType& dataType, string& reason)
		{
			string volReason, tiffReason, nrrdReason, sequenceReason, rawReason, pcrReason, nn5Reason, lz4Reason;
			if (vol::getInfo(filename, dimensions, dataType, volReason))
			{
				return true;
			}
			else if (tiff::getInfo(filename, dimensions, dataType, tiffReason))
			{
				return true;
			}
			else if (nrrd::getInfo(filename, dimensions, dataType, nrrdReason))
			{
				return true;
			}
			else if (sequence::getInfo(filename, dimensions, dataType, sequenceReason))
			{
				return true;
			}
			else if (pcr::getInfo(filename, dimensions, dataType, pcrReason))
			{
				return true;
			}
			else if (raw::getInfo(filename, dimensions, dataType, rawReason))
			{
				return true;
			}
			else if (nn5::getInfo(filename, dimensions, dataType, nn5Reason))
			{
				return true;
			}
			else if (lz4::getInfo(filename, dimensions, dataType, lz4Reason))
			{
				return true;
			}
			else
			{
				reason = internals::combineReasons(rawReason, tiffReason, sequenceReason, volReason, nrrdReason, pcrReason, nn5Reason, lz4Reason);
				return false;
			}
		}
	
	
		namespace tests
		{
			void readWrite()
			{
				Image<uint8_t> img8;
				Image<uint16_t> img;
				
				io::read(img8, "../test_input_data/uint8.png");
				io::read(img8, "../test_input_data/t1-head_bin_");
				
				io::read(img, "../test_input_data/t1-head");
				io::read(img, "../test_input_data/t1-head_");
				io::read(img, "../test_input_data/t1-head_256x256x129");
				io::read(img, "../test_input_data/t1-head_256x256x129.raw");
				
				io::read(img, "../test_input_data/t1-head.tif");

				sequence::write(img, "./sequence/head/head_@(5)_test.tif");

				Image<uint16_t> seq;
				io::read(seq, "./sequence/head/head_@_test.tif");

				testAssert(equals(img, seq), "read and written sequence are not equal.");
			}

			void badnn5()
			{
				string path = "./bad_nn5/0/0/0";
				fs::create_directories(path);

				Vec3c dimensions;
				ImageDataType dt;
				string reason;
				bool result = io::getInfo(path, dimensions, dt, reason);
				std::cout << reason << std::endl;
				testAssert(result == false, "reading bad nn5");

			}
		}
	}
}
