
#include "io.h"
#include "projections.h"

namespace itl2
{
	namespace io
	{
	
	
    	bool getInfo(const std::string& filename, Vec3c& dimensions, ImageDataType& dataType, string& reason)
		{
			string volReason, tiffReason, nrrdReason, sequenceReason, rawReason, pcrReason, nn5Reason;
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
			else
			{
				reason = internals::combineReasons(rawReason, tiffReason, sequenceReason, volReason, nrrdReason, pcrReason, nn5Reason);
				return false;
			}
		}
	
	
		namespace tests
		{
			void readWrite()
			{
				Image<uint8_t> img8;
				Image<uint16_t> img;
				io::read(img8, "./input_data/uint8.png");
				io::read(img8, "./input_data/t1-head_bin_");
				io::read(img, "./input_data/t1-head.tif");

				sequence::write(img, "./sequence/head/head_@(5)_test.tif");

				Image<uint16_t> seq;
				io::read(seq, "./sequence/head/head_@_test.tif");

				testAssert(equals(img, seq), "read and written sequence are not equal.");
			}
		}
	}
}
