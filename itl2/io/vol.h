#pragma once

#include "image.h"
#include "stringutils.h"
#include "math/vec3.h"
#include "io/raw.h"
#include "byteorder.h"

#include <string>
#include <array>

namespace itl2
{
	namespace vol
	{
		bool getInfo(const std::string& filename, Vec3c& dimensions, ImageDataType& dataType, bool& isBigEndian, size_t& headerSize, string& reason);

		/**
		Gets information about .vol file in the given path.
		*/
		inline bool getInfo(const std::string& filename, Vec3c& dimensions, ImageDataType& dataType, std::string& reason)
		{
			bool isBigEndian;
			size_t headerSize;
			return getInfo(filename, dimensions, dataType, isBigEndian, headerSize, reason);
		}

		template<typename pixel_t> void getInfoAndCheck(const std::string& filename, Vec3c& dimensions, ImageDataType& dataType, size_t& headerSize, bool& isBigEndian)
		{
			std::string reason;
			if (!getInfo(filename, dimensions, dataType, isBigEndian, headerSize, reason))
				throw ITLException(std::string("Not a .vol file: ") + filename + ". " + reason);

			// Check data
			if (dimensions.x <= 0 || dimensions.y <= 0 || dimensions.z <= 0)
				throw ITLException(std::string(".vol file contains invalid size specification: ") + toString(dimensions));

			if (dataType != imageDataType<pixel_t>())
				throw ITLException(std::string("Expected data type is ") + toString(imageDataType<pixel_t>()) + " but the .vol file contains data of type " + toString(dataType) + ".");
		}

		/**
		Reads .vol file from the given path.
		*/
		template<typename pixel_t> void read(Image<pixel_t>& img, const std::string& filename)
		{
			Vec3c dimensions;
			ImageDataType dataType;
			size_t headerSize;
			bool isBigEndianFile;
			getInfoAndCheck<pixel_t>(filename, dimensions, dataType, headerSize, isBigEndianFile);

			// Read data
			img.ensureSize(dimensions);
			raw::readNoParse(img, filename, headerSize);

			if (isBigEndianFile != isBigEndian())
				swapByteOrder(img);
		}

		/**
		Reads a block of .vol file from the given path.
		*/
		template<typename pixel_t> void readBlock(Image<pixel_t>& img, const std::string& filename, const Vec3c& start, bool showProgressInfo = false)
		{
			Vec3c dimensions;
			ImageDataType dataType;
			size_t headerSize;
			bool isBigEndianFile;
			getInfoAndCheck<pixel_t>(filename, dimensions, dataType, headerSize, isBigEndianFile);

			raw::readBlockNoParse(img, filename, dimensions, start, showProgressInfo, headerSize);

			if (isBigEndianFile != isBigEndian())
				swapByteOrder(img);
		}

		namespace tests
		{
			void volio();
		}
	}

}
