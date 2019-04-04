
#include "io/itlpng.h"
#include "io/raw.h"

namespace itl2
{
	namespace png
	{
		namespace internals
		{
			string lastPngError = "";

			void pngErrorFunc(png_structp png_ptr, png_const_charp error_msg)
			{
				lastPngError = error_msg;
				longjmp(png_jmpbuf(png_ptr), 1);
			}

			void pngWarningFunc(png_structp png_ptr, png_const_charp error_msg)
			{
#if defined(_DEBUG)
				cout << error_msg << endl;
#endif
				lastPngError = error_msg;
			}

			string pngLastError()
			{
				return lastPngError;
			}
		}

		bool getInfo(const string& filename, coord_t& width, coord_t& height, ImageDataType& dataType)
		{
			width = 0;
			height = 0;
			dataType = ImageDataType::Unknown;

			FILE* f;
			if (fopen_s(&f, filename.c_str(), "rb") != 0)
				return false;

			// Check that the file is .png
			uint8_t buf[8];
			if (fread(buf, 1, 8, f) != 8)
			{
				fclose(f);
				return false;
			}

			if (png_sig_cmp(buf, 0, 8) != 0)
			{
				fclose(f);
				return false;
			}

			png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, internals::pngErrorFunc, internals::pngErrorFunc);
			if (!png)
			{
				fclose(f);
				return false;
			}

			png_infop pngInfo = png_create_info_struct(png);
			if (!pngInfo)
			{
				fclose(f);
				png_destroy_read_struct(&png, NULL, NULL);
				return false;
			}

			png_set_sig_bytes(png, 8);

			if (setjmp(png_jmpbuf(png)) == 0)
			{
				png_init_io(png, f);
				png_read_info(png, pngInfo);

				png_uint_32 pngWidth, pngHeight;
				int bitDepth;
				int colorType;
				png_get_IHDR(png, pngInfo, &pngWidth, &pngHeight, &bitDepth, &colorType, NULL, NULL, NULL);

				width = pngWidth;
				height = pngHeight;

				if (bitDepth <= 8)
					dataType = ImageDataType::UInt8;
				else
					dataType = ImageDataType::UInt16;

				if (colorType != PNG_COLOR_TYPE_GRAY)
					dataType = ImageDataType::Unknown;
			}
			else
			{
				// libPng error
				dataType = ImageDataType::Unknown;
			}

			png_destroy_read_struct(&png, &pngInfo, NULL);
			fclose(f);

			if (dataType != ImageDataType::Unknown)
			{
				//cout << "png::getInfo returns true for " << filename << endl;
				return true;
			}

			return false;
		}

		namespace tests
		{
			void png()
			{
				coord_t w, h;
				ImageDataType dt;

				png::getInfo("./uint8.png", w, h, dt);
				testAssert(w == 100, "png width");
				testAssert(h == 200, "png height");
				testAssert(dt == ImageDataType::UInt8, "png data type (uint8)");

				Image<uint8_t> img1(w, h);
				png::read(img1, "./uint8.png", 0);
				raw::writed(img1, "./png/uint8");
				png::writed(img1, "./png/uint8_out", 0);


				png::getInfo("./uint16.png", w, h, dt);
				testAssert(w == 100, "png width");
				testAssert(h == 200, "png height");
				testAssert(dt == ImageDataType::UInt16, "png data type (uint16)");

				Image<uint16_t> img2(w, h);
				png::read(img2, "./uint16.png", 0);
				raw::writed(img2, "./png/uint16");
				png::writed(img2, "./png/uint16_out", 0);

			}
		}
	}
}