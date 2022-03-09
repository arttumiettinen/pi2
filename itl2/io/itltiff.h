#pragma once

#include "tiffio.h"

#include "math/vec3.h"
#include "io/imagedatatype.h"
#include "image.h"
#include "math/mathutils.h"
#include "transform.h"

#include <string>
#include <memory>

namespace itl2
{
	namespace tiff
	{
		namespace internals
		{
			bool getInfo(TIFF* tif, Vec3c& dimensions, ImageDataType& dataType, size_t& pixelSizeBytes, uint64_t& rawDataOffset, std::string& reason);

			/**
			Initialize .tiff reading library.
			This function can be called multiple times.
			*/
			void initTIFF();

			/**
			Gets the last error that occured in .tiff processing.
			*/
			std::string tiffLastError();

			template<typename pixel_t> void readDirectories(Image<pixel_t>& img, TIFF* tif, size_t z, const Vec3c& dimensions, const std::string& filename, uint64_t rawDataOffset)
			{
				if (z >= (size_t)img.depth())
					throw ITLException("Invalid target z coordinate.");


				if (rawDataOffset != 0)
				{
					// This is an ImageJ fake tiff. Read it as raw data file.
					raw::readBlockNoParse(img, filename, dimensions, Vec3c(0, 0, 0), false, (size_t)rawDataOffset);
					return;
				}

				// Read all directories
				do
				{
					if (z >= (size_t)img.depth())
						throw ITLException("TIFF file contains more frames than expected. This is a bug in the software. Please report it to the authors.");

					bool isTiled = TIFFIsTiled(tif) != 0;
					if (isTiled)
					{
						// Read tiles
						// TODO: This is not very well tested.

						uint32_t tileWidth = 0;
						uint32_t tileHeight = 0;
						uint32_t tileDepth = 0;

						TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tileWidth);
						TIFFGetField(tif, TIFFTAG_TILELENGTH, &tileHeight);
						TIFFGetField(tif, TIFFTAG_TILEDEPTH, &tileDepth);

						if (tileDepth < 1)
							tileDepth = 1;

						tmsize_t tiffTileSize = TIFFTileSize(tif);
						auto buf = std::unique_ptr<pixel_t, decltype(_TIFFfree)*>( (pixel_t*)_TIFFmalloc(tiffTileSize), _TIFFfree );

						for (coord_t z0 = z; z0 < dimensions.z; z0 += tileDepth)
						{
							for (coord_t y0 = 0; y0 < dimensions.y; y0 += tileHeight)
							{
								for (coord_t x0 = 0; x0 < dimensions.x; x0 += tileWidth)
								{
									TIFFReadTile(tif, buf.get(), (uint32_t)x0, (uint32_t)y0, (uint32_t)z0, 0);

									// Plot tile to image
									size_t n = 0;
									for (coord_t zz = z0; zz < std::min(z0 + tileDepth, img.depth()); zz++)
									{
										for (coord_t yy = y0; yy < std::min(y0 + tileHeight, img.height()); yy++)
										{
											for (coord_t xx = x0; xx < std::min(x0 + tileWidth, img.width()); xx++)
											{
												img(xx, yy, zz) = buf.get()[n];
												n++;
											}
										}
									}
								}
							}
						}

					}
					else
					{
						// Read strips

						tmsize_t stripSize = TIFFStripSize(tif);
						tstrip_t stripCount = TIFFNumberOfStrips(tif);

						uint32_t rowsPerStrip = 0;
						TIFFGetFieldDefaulted(tif, TIFFTAG_ROWSPERSTRIP, &rowsPerStrip);

						size_t imgPos = img.getLinearIndex(0, 0, z);
						for (tstrip_t strip = 0; strip < stripCount; strip++)
						{
							if (TIFFReadEncodedStrip(tif, strip, &(img.getData()[imgPos]), (tsize_t)-1) < 0)
								throw ITLException(string("TIFF read error: ") + internals::tiffLastError());

							imgPos += rowsPerStrip * img.width();
						}
					}

					z++;
				} while (TIFFReadDirectory(tif) == 1);

			}

			/*
			Read a .tif file.
			@param z Z-coordinate where the read data will be placed.
			*/
			template<typename pixel_t> void read(Image<pixel_t>& img, const std::string& filename, size_t z, bool is2D, bool allowResize)
			{		
				internals::initTIFF();

				auto tifObj = std::unique_ptr<TIFF, decltype(TIFFClose)*>(TIFFOpen(filename.c_str(), "r"), TIFFClose);
				TIFF* tif = tifObj.get();

				if (tif)
				{	
					Vec3c dimensions;
					ImageDataType dataType;
					size_t pixelSizeBytes;
					string reason;
					uint64_t rawDataOffset;
					if (internals::getInfo(tif, dimensions, dataType, pixelSizeBytes, rawDataOffset, reason))
					{
						if (dataType != imageDataType<pixel_t>() && pixelSizeBytes != sizeof(pixel_t))
							throw ITLException(string("Pixel data type in .tiff file is ") + toString(dataType) + " (" + toString(pixelSizeBytes) + " bytes per pixel), but image data type is " + toString(imageDataType<pixel_t>()) + " (" + toString(sizeof(pixel_t)) + " bytes per pixel).");

						if (is2D)
						{
							if (dimensions.z > 1)
								throw ITLException(string("Trying to read a 3D tiff as a 2D tiff: ") + filename);
							dimensions.z = 1;
							
							if(allowResize)
								img.ensureSize(dimensions.x, dimensions.y, img.depth());

							if(img.dimensions() == Vec3c(dimensions.x, dimensions.y, img.depth()))
							{
								internals::readDirectories(img, tif, z, dimensions, filename, rawDataOffset);
							}
							else
							{
								//if(img.dimensions() != Vec3c(dimensions.x, dimensions.y, img.depth()))
								//    throw ITLException(string("Dimensions of the data in the 2D .tif file ") + filename + string(" do not match to the dimensions of the image."));

								Image<pixel_t> tempImage(dimensions.x, dimensions.y, 1);
								internals::readDirectories(tempImage, tif, 0, dimensions, filename, rawDataOffset);
								
								// Calculate shift such that the on-disk image becomes centered in the space available.
								Vec3c target = (img.dimensions() - Vec3c(dimensions.x, dimensions.y, img.depth())) / 2;
								target.z = z;
								
								// Copy pixel values to final image.
								copyValues(img, tempImage, target);
							}
						}
						else
						{
							if(allowResize)
								img.ensureSize(dimensions);
								
							if(img.dimensions() == Vec3c(dimensions.x, dimensions.y, img.depth()))
							{
								internals::readDirectories(img, tif, 0, dimensions, filename, rawDataOffset);
							}
							else
							{
								//if(img.dimensions() != Vec3c(dimensions.x, dimensions.y, img.depth()))
								//    throw ITLException(string("Dimensions of the data in the .tif file ") + filename + string(" do not match to the dimensions of the image."));

								Image<pixel_t> tempImage(dimensions);
								internals::readDirectories(tempImage, tif, 0, dimensions, filename, rawDataOffset);
								
								// Calculate shift such that the on-disk image becomes centered in the space available.
								Vec3c target = (img.dimensions() - dimensions) / 2;
								
								// Copy pixel values to final image.
								copyValues(img, tempImage, target);
							}
						}
						
						return;
					}
					else
					{
						throw ITLException(reason);
					}
				}

				throw ITLException(string("Error while reading .tif image: ") + internals::tiffLastError());
			}
		}

		/*
		Get information of .tif image file.
		@param dimensions Dimensions of the image
		@param dataType Pixel data type of the image.
		@return True if the file seems to be an existing, valid .tiff file with supported pixel data type.
		*/
		bool getInfo(const std::string& filename, Vec3c& dimensions, ImageDataType& dataType, string& reason);
		
		template<typename pixel_t> void read2D(Image<pixel_t>& img, const std::string& filename, coord_t z, bool allowResize)
		{
			internals::read(img, filename, z, true, allowResize);
		}

		/*
		Read a .tif file.
		*/
		template<typename pixel_t> void read(Image<pixel_t>& img, const std::string& filename)
		{
			internals::read(img, filename, 0, false, true);
		}



		/**
		Reads part of a .tif file to given image.
		NOTE: Does not support out of bounds start position.
		@param img Image where the data is placed. The size of the image defines the size of the block that is read.
		@param filename The name of the file to read.
		@param start Start location of the read. The size of the image defines the size of the block that is read.
		*/
		template<typename pixel_t> void readBlock(Image<pixel_t>& img, const std::string& filename, const Vec3c& start, bool showProgressInfo = false)
		{
			internals::initTIFF();

			if (img.width() >= std::numeric_limits<uint32_t>::max() ||
				img.height() >= std::numeric_limits<uint32_t>::max() ||
				img.depth() >= std::numeric_limits<uint16_t>::max())
				throw ITLException("The image is too large to be read as a block of .tif file.");


			if ((start.x + img.width()) >= std::numeric_limits<uint32_t>::max() ||
				(start.y + img.height()) >= std::numeric_limits<uint32_t>::max() ||
				(start.z + img.depth()) >= std::numeric_limits<uint16_t>::max())
				throw ITLException("The block to be read is too far from origin to be located from a .tif file.");

			auto tifObj = std::unique_ptr<TIFF, decltype(TIFFClose)*>(TIFFOpen(filename.c_str(), "r"), TIFFClose);
			TIFF* tif = tifObj.get();

			if (tif)
			{
				Vec3c dimensions;
				ImageDataType dataType;
				size_t pixelSizeBytes;
				string reason;
				uint64_t rawDataOffset;
				if (internals::getInfo(tif, dimensions, dataType, pixelSizeBytes, rawDataOffset, reason))
				{
					if (dataType != imageDataType<pixel_t>() && pixelSizeBytes != sizeof(pixel_t))
						throw ITLException(string("Pixel data type in .tiff file is ") + toString(dataType) + " (" + toString(pixelSizeBytes) + " bytes per pixel), but image data type is " + toString(imageDataType<pixel_t>()) + " (" + toString(sizeof(pixel_t)) + " bytes per pixel).");

					if (rawDataOffset != 0)
					{
						// This is an ImageJ fake tiff. Read it as a .raw file.
						raw::readBlockNoParse(img, filename, dimensions, start, showProgressInfo, rawDataOffset);
						return;
					}

					// Read all directories
					coord_t tifz = start.z;
					do
					{
						if(TIFFSetDirectory(tif, (uint16_t)tifz) != 1)
							throw ITLException(string("Unable to read expected number of slices from .tif file ") + filename);

						coord_t z = tifz - start.z;

						if (z >= img.depth())
							break; // All read!

						bool isTiled = TIFFIsTiled(tif) != 0;
						if (isTiled)
						{
							// Read tiles
							
							throw ITLException("Block-by-block reading of tiled .tif files is not implemented due to lack of tiled example .tif files. Please send your tiled .tif file to the developers of the software to get this feature implemented. Please use the info() command or the website to find the latest contact details.");
						}
						else
						{
							// Read strips

							tmsize_t stripSize = TIFFStripSize(tif);
							tstrip_t stripCount = TIFFNumberOfStrips(tif);
							auto buf = std::unique_ptr<pixel_t, decltype(_TIFFfree)*>((pixel_t*)_TIFFmalloc(stripSize), _TIFFfree);
							uint32_t rowsPerStrip = 0;
							TIFFGetFieldDefaulted(tif, TIFFTAG_ROWSPERSTRIP, &rowsPerStrip);

							coord_t y = 0; // y-coordinate in full image coordinates
							for (tstrip_t strip = 0; strip < stripCount; strip++)
							{
								if (TIFFReadEncodedStrip(tif, strip, buf.get(), (tsize_t)-1) < 0)
									throw ITLException(string("TIFF read error: ") + internals::tiffLastError());


								// stripy = y-coordinate in strip coordinates
								// y in image coordinates = y - start.y
								for (coord_t stripy = 0; y < start.y + img.height() && stripy < rowsPerStrip; y++, stripy++)
								{
									if (y >= start.y)
									{
										// Put the data to the image
										for (coord_t x = start.x; x < start.x + img.width(); x++)
											img(x - start.x, y - start.y, z) = buf.get()[x + stripy * dimensions.x];
									}
								}
							}
						}

						tifz++;
					} while (TIFFReadDirectory(tif) == 1);

					return;
				}
				else
				{
					throw ITLException(reason);
				}
			}

			throw ITLException(string("Error while reading .tif image: ") + internals::tiffLastError());
		}



		/**
		Writes a .tif file.
		*/
		template<typename pixel_t> void write(const Image<pixel_t>& img, const std::string& filename)
		{
			createFoldersFor(filename);

			internals::initTIFF();

			if (img.width() >= std::numeric_limits<uint32_t>::max() ||
				img.height() >= std::numeric_limits<uint32_t>::max() ||
				img.depth() >= std::numeric_limits<uint16_t>::max())
				throw ITLException("The image is too large to be written to a .tif file.");

			string mode = "w"; // Normal tif file
			if (img.pixelCount() * img.pixelSize() > (size_t)4 * (size_t)1000 * (size_t)1000 * (size_t)1000)
				mode = "w8"; // BigTIFF file
			auto tifObj = std::unique_ptr<TIFF, decltype(TIFFClose)*>(TIFFOpen(filename.c_str(), mode.c_str()), TIFFClose);
			TIFF* tif = tifObj.get();

			if (tif)
			{
				for (coord_t z = 0; z < img.depth(); z++)
				{
					// Write header information
					uint32_t w = (uint32_t)img.width();
					uint32_t h = (uint32_t)img.height();
					TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, w);
					TIFFSetField(tif, TIFFTAG_IMAGELENGTH, h);

					TIFFSetField(tif, TIFFTAG_IMAGEDEPTH, 1);
					TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);


					uint16_t bitsPerSample = sizeof(pixel_t) * 8;
					TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, bitsPerSample);

					uint16_t sampleFormat = 0;
					uint16_t tiffDataType = 0;
					switch (imageDataType<pixel_t>())
					{
					case ImageDataType::UInt8:
						tiffDataType = TIFF_BYTE;
						sampleFormat = SAMPLEFORMAT_UINT;
						break;
					case ImageDataType::UInt16:
						tiffDataType = TIFF_SHORT;
						sampleFormat = SAMPLEFORMAT_UINT;
						break;
					case ImageDataType::UInt32:
						tiffDataType = TIFF_LONG;
						sampleFormat = SAMPLEFORMAT_UINT;
						break;
					case ImageDataType::UInt64:
						tiffDataType = TIFF_LONG8;
						sampleFormat = SAMPLEFORMAT_UINT;
						break;
					case ImageDataType::Int8:
						tiffDataType = TIFF_BYTE;
						sampleFormat = SAMPLEFORMAT_INT;
						break;
					case ImageDataType::Int16:
						tiffDataType = TIFF_SHORT;
						sampleFormat = SAMPLEFORMAT_INT;
						break;
					case ImageDataType::Int32:
						tiffDataType = TIFF_LONG;
						sampleFormat = SAMPLEFORMAT_INT;
						break;
					case ImageDataType::Int64:
						tiffDataType = TIFF_LONG8;
						sampleFormat = SAMPLEFORMAT_INT;
						break;
					case ImageDataType::Float32:
						tiffDataType = TIFF_FLOAT;
						sampleFormat = SAMPLEFORMAT_IEEEFP;
						break;
					case ImageDataType::Complex32:
						tiffDataType = TIFF_NOTYPE;
						sampleFormat = SAMPLEFORMAT_COMPLEXIEEEFP;
						break;
					default:
						tiffDataType = TIFF_NOTYPE;
						sampleFormat = SAMPLEFORMAT_VOID;
						break;
					}
					TIFFSetField(tif, TIFFTAG_DATATYPE, tiffDataType);
					TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, sampleFormat);

					TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, h);
					TIFFSetField(tif, TIFFTAG_PLANARCONFIG, 1);

					TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);

					string softName = "pi2";
					TIFFSetField(tif, TIFFTAG_SOFTWARE, softName.c_str());

					// Write image data
					void* ptr = (void*)&img(0, 0, z);
					if (TIFFWriteEncodedStrip(tif, 0, ptr, img.width() * img.height() * sizeof(pixel_t)) < 0)
					{
						throw ITLException(string("Unable to write data to .tif file ") + filename + ": " + internals::tiffLastError());
					}


					// Go to next directory if we are going to write more slices
					if (z < img.depth() - 1)
					{
						if(TIFFWriteDirectory(tif) == 0)
							throw ITLException(string("Error while writing to .tif file ") + filename + ": " + internals::tiffLastError());
					}
				}
			}
			else
			{
				throw ITLException(string("Error while opening .tif image file ") + filename + " for writing: " + internals::tiffLastError());
			}
		}

		/*
		Writes a 2D .tif file.
		@param z Z-coordinate of the slice that will be written.
		*/
		template<typename pixel_t> void write2D(const Image<pixel_t>& img, const std::string& filename, coord_t z)
		{
			Image<pixel_t> view(img, z, z);
			write(view, filename);
		}

		/*
		Write a .tif file, adds .tif to the file name if it does not end with .tif or .tiff.
		*/
		template<typename pixel_t> void writed(const Image<pixel_t>& img, const std::string& filename)
		{
			if (endsWithIgnoreCase(filename, ".tif") || endsWithIgnoreCase(filename, ".tiff"))
				write(img, filename);
			else
				write(img, filename + ".tif");
		}


		namespace tests
		{
			void readWrite();
			void imageJLargeTiff();
		}
	}
}
