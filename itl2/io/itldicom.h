#pragma once

#include "image.h"
#include "io/raw.h"

// Header path: $(SolutionDir)gdcm\include\gdcm-3.0;
// Libs: gdcmMSFF.lib;gdcmDSED.lib
//#pragma warning(disable: 996)
//#include "gdcmImageReader.h"
//#pragma warning(default: 996)

namespace itl2
{
	namespace dicom
	{
		namespace internals
		{
			bool getInfo(const std::string& filename, Vec3c& dimensions, ImageDataType& dataType, size_t& pixelSizeBytes, uint64_t& rawDataOffset, bool& littleEndian, std::string& reason);

			/**
			Read a .dcm file.
			@param z Z-coordinate where the read data will be placed.
			*/
			template<typename pixel_t> void read(Image<pixel_t>& img, const std::string& filename, size_t z, bool is2D, bool allowResize)
			{
				Vec3c dimensions;
				ImageDataType dataType;
				size_t pixelSizeBytes;
				uint64_t offset;
				string reason;
				bool littleEndian;
				if (!getInfo(filename, dimensions, dataType, pixelSizeBytes, offset, littleEndian, reason))
					throw ITLException(reason);

				//cout << "Reason = " << reason << endl;
				//cout << "Dimensions = " << dimensions << endl;
				//cout << "Data type = " << toString(dt) << endl;
				//cout << "Offset = " << offset << endl;
				//cout << "Little endian = " << littleEndian << endl;
				//cout << "Pixel size (bytes) = " << pixSize << endl;

				if (dataType != imageDataType<pixel_t>() && pixelSizeBytes != sizeof(pixel_t))
					throw ITLException(string("Pixel data type in DICOM file is ") + toString(dataType) + " (" + toString(pixelSizeBytes) + " bytes per pixel), but image data type is " + toString(imageDataType<pixel_t>()) + " (" + toString(sizeof(pixel_t)) + " bytes per pixel).");

				if (is2D)
				{
					if (dimensions.z > 1)
						throw ITLException(string("Trying to read a 3D DICOM as a 2D DICOM: ") + filename);
					dimensions.z = 1;

					if (allowResize)
						img.ensureSize(dimensions.x, dimensions.y, img.depth());

					if (z == 0)
					{
						itl2::raw::readBlockNoParse(img, filename, dimensions, Vec3c(0, 0, 0), (size_t)offset);
					}
					else
					{
						Image<pixel_t> view(img, z, z);
						itl2::raw::readBlockNoParse(view, filename, dimensions, Vec3c(0, 0, 0), (size_t)offset);
					}
				}
				else
				{
					if (allowResize)
						img.ensureSize(dimensions);
					itl2::raw::readBlockNoParse(img, filename, dimensions, Vec3c(0, 0, 0), (size_t)offset);
				}
			}
		}

		/**
		Get information of a .dcm DICOM image file.
		@param dimensions Dimensions of the image
		@param dataType Pixel data type of the image.
		@return True if the file seems to be an existing, valid DICOM file with supported pixel data type.
		*/
		bool getInfo(const std::string& filename, Vec3c& dimensions, ImageDataType& dataType, string& reason);

		/**
		Read a 2D .dcm file.
		@param img Target image.
		@param filename Name of the file to read.
		@param z Z-position of the 2D slice in the image.
		@param allowResize Resize img to match file dimensions.
		*/
		template<typename pixel_t> void read2D(Image<pixel_t>& img, const std::string& filename, coord_t z, bool allowResize)
		{
			internals::read(img, filename, z, true, allowResize);
		}

		/**
		Read a .dcm file into the given image.
		*/
		template<typename pixel_t> void read(Image<pixel_t>& img, const std::string& filename)
		{
			internals::read(img, filename, 0, false, true);
		}

		namespace tests
		{
			void read();
		}
	}
}