#pragma once

#include <string>

#include "io/sequence.h"
#include "io/raw.h"
#include "io/vol.h"

namespace itl2
{
	namespace io
	{
		/**
		Finds data format of given file and reads it to the given image.
		The data type of the target image must be correct but its size is set automatically.
		*/
		template<typename pixel_t> void read(Image<pixel_t>& img, const std::string& filename)
		{
			math::Vec3c dimensions;
			ImageDataType dt;
			
			if (vol::getInfo(filename, dimensions, dt))
			{
				vol::read(img, filename);
			}
			else if (tiff::getInfo(filename, dimensions, dt))
			{
				tiff::read(img, filename);
			}
			else if (sequence::getInfo(filename, dimensions, dt))
			{
				sequence::read(img, filename);
			}
			else if (raw::getInfo(filename, dimensions, dt))
			{
				raw::read(img, filename);
			}
			else
			{
				throw ITLException(string("Unsupported file type, file not found, or cannot be read: ") + filename);
			}
		}

		/**
		Parses the given file and finds the dimensions and pixel data type of image stored in the file.
		@return True if the parsing succeeded, false if the file cannot be loaded.
		*/
		inline bool getInfo(const std::string& filename, math::Vec3c& dimensions, ImageDataType& dataType)
		{
			if (vol::getInfo(filename, dimensions, dataType))
			{
				return true;
			}
			else if (tiff::getInfo(filename, dimensions, dataType))
			{
				return true;
			}
			else if (sequence::getInfo(filename, dimensions, dataType))
			{
				return true;
			}
			else if (raw::getInfo(filename, dimensions, dataType))
			{
				return true;
			}
			else
			{
				return false;
			}
		}

		/**
		Finds data format of given file and reads part of it to the given image.
		The data type of the target image must be correct but its size is set automatically.
		The size of the block to read is determined by the size of the image.
		@param blockStart The start position of the block.
		*/
		template<typename pixel_t> void readBlock(Image<pixel_t>& img, const std::string& filename, const math::Vec3c& blockStart, bool showProgressInfo = false)
		{
			math::Vec3c dimensions;
			ImageDataType dt;
			if (vol::getInfo(filename, dimensions, dt))
			{
				vol::readBlock(img, filename, blockStart, showProgressInfo);
			}
			else if (tiff::getInfo(filename, dimensions, dt))
			{
				tiff::readBlock(img, filename, blockStart, showProgressInfo);
			}
			else if (sequence::getInfo(filename, dimensions, dt))
			{
				sequence::readBlock(img, filename, blockStart, showProgressInfo);
			}
			else if (raw::getInfo(filename, dimensions, dt))
			{
				raw::readBlock(img, filename, blockStart, showProgressInfo);
			}
			else
			{
				throw ITLException(string("Unsupported file type, file not found, or cannot be read: ") + filename);
			}
		}

		namespace tests
		{
			void readWrite();
		}
	}
}
