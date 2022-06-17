
#include "io/itlpng.h"
#include "io/raw.h"

#include <iostream>
using namespace std;

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
    			string msg = error_msg;
    			if(startsWith(msg, "Application built with"))
    			    throw ITLException(error_msg);
	    
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

		bool getInfo(const string& filename, coord_t& width, coord_t& height, ImageDataType& dataType, string& reason)
		{
			width = 0;
			height = 0;
			dataType = ImageDataType::Unknown;

			// This is required in some Linux systems to differentiate files from directories.
			if (!fs::is_regular_file(filename))
			{
				reason = "Not a file.";
				return false;
			}

			FILE* f;
			if (fopen_s(&f, filename.c_str(), "rb") != 0)
			{
				reason = "Unable to open file.";
				return false;
			}
				
			// Check that the file is .png
			uint8_t buf[8];
			if (fread(buf, 1, 8, f) != 8)
			{
				fclose(f);
				reason = "Unable to read header from the file.";
				return false;
			}

			if (png_sig_cmp(buf, 0, 8) != 0)
			{
				fclose(f);
				reason = "The file does not contain .png header.";
				return false;
			}

			png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, nullptr, internals::pngErrorFunc, internals::pngWarningFunc);
			if (!png)
			{
				fclose(f);
				reason = "Unable to initialize .png reader.";
				return false;
			}

			png_infop pngInfo = png_create_info_struct(png);
			if (!pngInfo)
			{
				fclose(f);
				png_destroy_read_struct(&png, nullptr, nullptr);
				reason = "Unable to initialize .png reader.";
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
				
				png_get_IHDR(png, pngInfo, &pngWidth, &pngHeight, &bitDepth, &colorType, nullptr, nullptr, nullptr);

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

			png_destroy_read_struct(&png, &pngInfo, nullptr);
			
			fclose(f);

			if (dataType != ImageDataType::Unknown)
			{
				//cout << "png::getInfo returns true for " << filename << endl;
				return true;
			}

			reason = "The .png file contains data in unsupported pixel format.";
			return false;
		}

		namespace tests
		{
			void png()
			{
				coord_t w, h;
				ImageDataType dt;
				string reason;

				png::getInfo("../test_input_data/uint8.png", w, h, dt, reason);
				testAssert(w == 100, "png width");
				testAssert(h == 200, "png height");
				testAssert(dt == ImageDataType::UInt8, "png data type (uint8)");

				Image<uint8_t> img1(w, h);
				png::read(img1, "../test_input_data/uint8.png", 0);
				raw::writed(img1, "./png/uint8");
				png::writed(img1, "./png/uint8_out", 0);


				png::getInfo("../test_input_data/uint16.png", w, h, dt, reason);
				testAssert(w == 100, "png width");
				testAssert(h == 200, "png height");
				testAssert(dt == ImageDataType::UInt16, "png data type (uint16)");

				Image<uint16_t> img2(w, h);
				png::read(img2, "../test_input_data/uint16.png", 0);
				raw::writed(img2, "./png/uint16");
				png::writed(img2, "./png/uint16_out", 0);

			}
		}
	}
}
