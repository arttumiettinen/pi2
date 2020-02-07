#pragma once

#include <stdio.h>
#include "png.h"
#include "image.h"
#include "io/imagedatatype.h"

namespace itl2
{

	namespace png
	{

#ifndef _WIN32	
		inline int fopen_s(FILE **f, const char *name, const char *mode)
		{
			int ret = 0;
			*f = fopen(name, mode);
			if (!*f)
				ret = errno;
			return ret;
		}
#endif

		namespace internals
		{
			void pngErrorFunc(png_structp png_ptr, png_const_charp error_msg);

			void pngWarningFunc(png_structp png_ptr, png_const_charp error_msg);

			string pngLastError();

			/*
			Reads .png file, does not throw exceptions but returns success/failure and error message.
			*/
			template<typename pixel_t> bool readNoThrow(Image<pixel_t>& img, const string& filename, coord_t z, string& errorMessage)
			{
				
				if (z < 0 || z >= img.depth())
				{
					errorMessage = "Invalid z-coordinate.";
					return false;
				}

				FILE* f;
				if (fopen_s(&f, filename.c_str(), "rb") != 0)
				{
					errorMessage = string("Unable to open file ") + filename + ".";
					return false;
				}

				png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, nullptr, internals::pngErrorFunc, internals::pngErrorFunc);
				if (!png)
				{
					errorMessage = "Unable to create png support data structure.";
					fclose(f);
					return false;
				}

				png_infop pngInfo = png_create_info_struct(png);
				if (!pngInfo)
				{
					errorMessage = "Unable to create png info data structure.";
					fclose(f);
					png_destroy_read_struct(&png, nullptr, nullptr);
					return false;
				}

				volatile png_bytepp rowPointers = 0;
				volatile ImageDataType dataType = ImageDataType::Unknown;

				if (setjmp(png_jmpbuf(png)) == 0)
				{
					png_init_io(png, f);

					png_read_info(png, pngInfo);

					png_uint_32 pngWidth, pngHeight;
					int bitDepth;
					int colorType;
					png_get_IHDR(png, pngInfo, &pngWidth, &pngHeight, &bitDepth, &colorType, nullptr, nullptr, nullptr);

					// Check that image size is correct.
					if (img.width() == pngWidth && img.height() == pngHeight)
					{

						// Make sure that image data type and png file data type match.
						if (colorType == PNG_COLOR_TYPE_GRAY)
						{
							if ((bitDepth <= 8 && imageDataType<pixel_t>() == ImageDataType::UInt8)||
								(bitDepth == 16 && imageDataType<pixel_t>() == ImageDataType::UInt16))
							{
								// Expand < 8 bit depths to 8
								if (bitDepth < 8)
									png_set_expand(png);

								// Convert to little endian
								if (bitDepth == 16)
									png_set_swap(png);

								png_read_update_info(png, pngInfo);

								// Make sure rows have the correct size.
								png_size_t rowBytes = png_get_rowbytes(png, pngInfo);
								if (rowBytes == pngWidth * sizeof(pixel_t))
								{
									// Construct row pointers to libpng.
									rowPointers = new png_bytep[img.height()];
									for (coord_t y = 0; y < img.height(); y++)
										rowPointers[y] = (png_bytep)&img(0, y, z);

									// Final image reading is done here.
									png_read_image(png, rowPointers);

									// Set dataType to signal succesful read.
									if (bitDepth == 16)
										dataType = ImageDataType::UInt16;
									else
										dataType = ImageDataType::UInt8;
								}
								else
								{
									errorMessage = "Invalid rowBytes value.";
								}
							}
							else
							{
								errorMessage = "Pixel data type of the .png file does not match to the pixel data type requested.";
							}
						}
						else
						{
							errorMessage = "Png file does not contain a grayscale image.";
						}
					}
					else
					{
						errorMessage = "Image width and height do not match to file width and height.";
					}
				}
				else
				{
					errorMessage = internals::pngLastError();
				}

				if(rowPointers)
					delete[] rowPointers;
				png_destroy_read_struct(&png, &pngInfo, nullptr);
				fclose(f);

				if (dataType != ImageDataType::Unknown)
					return true;

				return false;
			}

			/*
			Writes a .png file, does not throw exceptions but returns success/failure and error message.
			*/
			template<typename pixel_t> bool writeNoThrow(const Image<pixel_t>& img, const string& filename, coord_t z, string& errorMessage)
			{
				if (sizeof(pixel_t) > 16)
				{
					errorMessage = "Png supports 8-bit and 16-bit images only.";
					return false;
				}

				if (z < 0 || z >= img.depth())
				{
					errorMessage = "Invalid z-coordinate.";
					return false;
				}

				createFoldersFor(filename);

				FILE *f;
				if (fopen_s(&f, filename.c_str(), "wb") != 0)
				{
					errorMessage = string("Unable to open file ") + filename + " for writing.";
					return false;
				}

				png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, internals::pngErrorFunc, internals::pngErrorFunc);
				if (!png)
				{
					fclose(f);
					errorMessage = "Unable to create png structure.";
					return false;
				}

				png_infop pngInfo = png_create_info_struct(png);
				if (!png)
				{
					fclose(f);
					png_destroy_write_struct(&png, nullptr);
					errorMessage = "Unable to create png info structure.";
					return false;
				}

				volatile bool result = false;
				volatile png_bytepp rowPointers = 0;
				if (setjmp(png_jmpbuf(png)) == 0)
				{
					png_init_io(png, f);

					png_set_IHDR(png, pngInfo, (png_uint_32)img.width(), (png_uint_32)img.height(), sizeof(pixel_t) * 8, PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
					
					// Construct row pointers to libpng.
					rowPointers = new png_bytep[img.height()];
					for (coord_t y = 0; y < img.height(); y++)
						rowPointers[y] = (png_bytep)&img(0, y, z);

					png_set_rows(png, pngInfo, rowPointers);

					png_write_png(png, pngInfo, PNG_TRANSFORM_SWAP_ENDIAN, nullptr);

					result = true;
				}
				else
				{
					errorMessage = internals::pngLastError();
				}

				if (rowPointers)
					delete[] rowPointers;
				png_destroy_write_struct(&png, &pngInfo);
				fclose(f);

				return result;
			}
		}

		/*
		Get information from .png image file.
		Supports only 8- and 16-bit grayscale images.
		@param width, height Dimensions of the image
		@param dataType Pixel data type of the image.
		@return True if the file seems to be an existing, valid .png file with supported pixel data type.
		*/
		bool getInfo(const string& filename, coord_t& width, coord_t& height, ImageDataType& dataType, string& reason);

		/*
		Read a .png file.
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

		/*
		Write a .png file.
		Supports only 8- and 16-bit grayscale images.
		@param z Z-coordinate of the slice that will be written.
		*/
		template<typename pixel_t> void write(const Image<pixel_t>& img, const string& filename, coord_t z = 0)
		{
			string errorMessage;
			if (!internals::writeNoThrow(img, filename, z, errorMessage))
				throw ITLException(errorMessage);
		}

		/*
		Write a .png file, adds .png to the file name.
		Supports only 8- and 16-bit grayscale images.
		@param z Z-coordinate of the slice that will be written.
		*/
		template<typename pixel_t> void writed(const Image<pixel_t>& img, const string& filename, coord_t z = 0)
		{
			if(endsWithIgnoreCase(filename, ".png"))
				write(img, filename, z);
			else
				write(img, filename + ".png", z);
		}


		namespace tests
		{
			void png();
		}
	}

}
