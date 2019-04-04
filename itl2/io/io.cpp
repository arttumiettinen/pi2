
#include "io.h"
#include "projections.h"

namespace itl2
{
	namespace io
	{
		namespace tests
		{
			void readWrite()
			{
				Image<uint8_t> img8;
				Image<uint16_t> img;
				io::read(img8, "uint8.png");
				io::read(img8, "t1-head_bin_");
				io::read(img, "t1-head.tif");

				sequence::write(img, "./sequence/head/head_@(5)_test.tif");

				Image<uint16_t> seq;
				io::read(seq, "./sequence/head/head_@_test.tif");

				testAssert(equals(img, seq), "read and written sequence are not equal.");
			}
		}
	}
}