
#include "itljpeg.h"

#include "io/raw.h"
#include "io/itltiff.h"
#include "projections.h"

namespace itl2
{
	namespace jpeg
	{
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

			FILE* infile;
			if (fopen_s(&infile, filename.c_str(), "rb") != 0)
			{
				reason = "Unable to open file.";
				return false;
			}

			// Check that the file is .jpg
			uint16_t buf;
			if (fread(&buf, 1, 2, infile) != 2)
			{
				fclose(infile);
				reason = "Unable to read header from the file.";
				return false;
			}

			if (buf != 0xd8ff)
			{
				fclose(infile);
				reason = "The file does not contain a JPEG header.";
				return false;
			}

			// Put file read pointer back to the beginning of the file.
			fseek(infile, 0, SEEK_SET);

			// Read header
			jpeg_decompress_struct cinfo;
			internals::my_error_mgr jerr;
			cinfo.err = jpeg_std_error(&jerr.pub);
			jerr.pub.error_exit = internals::my_error_exit;
			// Establish the setjmp return context for my_error_exit to use.
			if (setjmp(jerr.setjmp_buffer))
			{
				// If we get here, the JPEG code has signaled an error.
				jpeg_destroy_decompress(&cinfo);
				fclose(infile);
				reason = "Unable to decompress jpeg data.";
				return false;
			}

			jpeg_create_decompress(&cinfo);
			jpeg_stdio_src(&cinfo, infile);
			jpeg_read_header(&cinfo, TRUE);
			jpeg_start_decompress(&cinfo);
			if (cinfo.output_components != 1)
			{
				reason = "The jpeg file contains more than 1 color component.";
				jpeg_destroy_decompress(&cinfo);
				fclose(infile);
				return false;
			}
			jpeg_destroy_decompress(&cinfo);
			fclose(infile);

			width = cinfo.output_width;
			height = cinfo.output_height;
			dataType = ImageDataType::UInt8;
			return true;
		}

		namespace tests
		{
			void jpeg()
			{
				coord_t w, h;
				ImageDataType dt;
				string reason;
				testAssert(jpeg::getInfo("../test_input_data/uint8.jpg", w, h, dt, reason), "jpeg getinfo return value");
				testAssert(w == 100, "jpeg width");
				testAssert(h == 200, "jpeg height");
				testAssert(dt == ImageDataType::UInt8, "jpeg datatype");

				testAssert(jpeg::getInfo("../test_input_data/non existent file.jpg", w, h, dt, reason) == false, "jpeg getinfo return value");


				Image<uint8_t> img;
				jpeg::read(img, "../test_input_data/uint8.jpg");
				raw::writed(img, "./jpeg/from_jpeg_uint8");

				Image<uint8_t> img2;
				tiff::read(img2, "../test_input_data/uint8.tif");

				testAssert(equalsTol(img, img2, (uint8_t)25), "Tiff and Jpeg images differ.");
			}
		}
	}
}