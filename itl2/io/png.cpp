
#include "io/itlpng.h"
#include "io/raw.h"

namespace itl2
{
	namespace png
	{
		namespace tests
		{
			void png()
			{
				coord_t w, h;
				ImageDataType dt;

				png::getInfo("./uint8.png", w, h, dt);
				testAssert(w == 100, "png width");
				testAssert(h == 200, "png height");
				testAssert(dt == UInt8, "png data type (uint8)");

				Image<uint8_t> img1(w, h);
				png::read(img1, "./uint8.png", 0);
				raw::writed(img1, "./png/uint8");
				png::writed(img1, "./png/uint8_out", 0);


				png::getInfo("./uint16.png", w, h, dt);
				testAssert(w == 100, "png width");
				testAssert(h == 200, "png height");
				testAssert(dt == UInt16, "png data type (uint16)");

				Image<uint16_t> img2(w, h);
				png::read(img2, "./uint16.png", 0);
				raw::writed(img2, "./png/uint16");
				png::writed(img2, "./png/uint16_out", 0);

			}
		}
	}
}