#pragma once

#include <string>
#include "image.h"

#include "jpeglib.h"

namespace itl2
{

	namespace jpeg
	{
		namespace internals
		{
			struct my_error_mgr {
				struct jpeg_error_mgr pub;
				jmp_buf setjmp_buffer;
			};
			
			inline void my_error_exit(j_common_ptr cinfo)
			{
				// cinfo->err really points to a my_error_mgr struct, so coerce pointer
				my_error_mgr* myerr = (my_error_mgr*)cinfo->err;

				// Return control to the setjmp point
				longjmp(myerr->setjmp_buffer, 1);
			}

			/*
			Reads a .jpg file, does not throw exceptions but returns success/failure and error message.
			*/
			template<typename pixel_t> bool readNoThrow(Image<pixel_t>& img, const string& filename, coord_t z, string& errorMessage)
			{
				errorMessage = "";

				if (imageDataType<pixel_t>() != ImageDataType::UInt8)
					throw ITLException("Jpeg format supports uint8 images only.");

				if (z < 0 || z >= img.depth())
				{
					errorMessage = "Invalid z-coordinate.";
					return false;
				}
				
				// Source file
				FILE* infile = NULL;
				if (fopen_s(&infile, filename.c_str(), "rb") != 0)
				{
					errorMessage = string("Unable to open file.");
					return false;
				}

				// We set up the normal JPEG error routines, then override error_exit.
				jpeg_decompress_struct cinfo;

				// We use our private extension JPEG error handler.
				my_error_mgr jerr;

				cinfo.err = jpeg_std_error(&jerr.pub);
				jerr.pub.error_exit = my_error_exit;
				// Establish the setjmp return context for my_error_exit to use.
				if (setjmp(jerr.setjmp_buffer))
				{
					// If we get here, the JPEG code has signaled an error.
					jpeg_destroy_decompress(&cinfo);
					fclose(infile);
					errorMessage = string("Unable to decompress jpeg data.");
					return false;
				}

				jpeg_create_decompress(&cinfo);

				// specify data source (eg, a file)
				jpeg_stdio_src(&cinfo, infile);

				// read file parameters with jpeg_read_header()
				jpeg_read_header(&cinfo, TRUE);
				// We can ignore the return value from jpeg_read_header since
				//   (a) suspension is not possible with the stdio data source, and
				//   (b) we passed TRUE to reject a tables-only JPEG file as an error.

				// Start decompressor
				jpeg_start_decompress(&cinfo);

				if (cinfo.output_components == 1)
				{

					if (img.width() != cinfo.output_width || img.height() != cinfo.output_height)
						img.ensureSize(cinfo.output_width, cinfo.output_height, img.depth());

					// Check that image size is correct.
					//if (img.width() == cinfo.output_width && img.height() == cinfo.output_height)
					//{
						// Read scanlines
						// Here we use the library's state variable cinfo.output_scanline as the
						// loop counter, so that we don't have to keep track ourselves.
						while (cinfo.output_scanline < cinfo.output_height)
						{
							// jpeg_read_scanlines expects an array of pointers to scanlines.
							JSAMPROW rowPtr = (JSAMPLE*)&img(0, cinfo.output_scanline, 0);
							jpeg_read_scanlines(&cinfo, &rowPtr, 1);
						}

						// Finish decompression
						jpeg_finish_decompress(&cinfo);
					//}
					//else
					//{
					//	errorMessage = "Image width and height do not match to file width and height.";
					//}
				}
				else
				{
					errorMessage = string("Only images with 1 color component are supported, but the image contains ") + toString(cinfo.output_components) + " color components.";
				}

				// Release JPEG decompression object
				jpeg_destroy_decompress(&cinfo);
				fclose(infile);

				// At this point you may want to check to see whether any corrupt-data
				// warnings occurred (test whether jerr.pub.num_warnings is nonzero).

				if (errorMessage != "")
					return false;

				return true;
			}
		}


		/*
		Get information from a .jpg image file.
		Supports 8- and 16-bit grayscale images only.
		@param width, height Dimensions of the image
		@param dataType Pixel data type of the image.
		@return True if the file seems to be an existing, valid .jpeg file with supported pixel data type.
		*/
		bool getInfo(const string& filename, coord_t& width, coord_t& height, ImageDataType& dataType, string& reason);

		inline bool getInfo(const std::string& filename, Vec3c& dimensions, ImageDataType& dataType, string& reason)
		{
			coord_t w, h;
			bool result = getInfo(filename, w, h, dataType, reason);
			if (result)
				dimensions = Vec3c(w, h, 1);
			return true;
		}

		/*
		Read a .jpg file.
		Supports only 8- and 16-bit grayscale images.
		The image will be resized in x and y so that it has the same size than the .png image.
		@param z Z-coordinate where the read data will be placed.
		*/
		template<typename pixel_t> void read(Image<pixel_t>& img, const string& filename, coord_t z = 0)
		{
			string errorMessage;
			if (!internals::readNoThrow(img, filename, z, errorMessage))
				throw ITLException(errorMessage);
		}

		namespace tests
		{
			void jpeg();
		}
	}
}