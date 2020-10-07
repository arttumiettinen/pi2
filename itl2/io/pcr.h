#pragma once

#include <string>
#include <iostream>
#include "filesystem.h"

#include "math/vec3.h"
#include "io/imagedatatype.h"
#include "image.h"
#include "io/raw.h"
#include "pointprocess.h"

namespace itl2
{
	namespace pcr
	{

		/**
		Gets information about file in the given path.
		*/
		bool getInfo(const std::string& filename, Vec3c& dimensions, ImageDataType& dataType, string& failReason, string& datafile, bool& isBigEndian);

		/**
		Gets information about file in the given path.
		*/
		inline bool getInfo(const std::string& filename, Vec3c& dimensions, ImageDataType& dataType, string& reason)
		{
			std::string datafile;
			bool isBigEndian;
			return getInfo(filename, dimensions, dataType, reason, datafile, isBigEndian);
		}

		/**
		Reads file from the given path.
		*/
		template<typename pixel_t> void read(Image<pixel_t>& img, const std::string& filename)
		{
			Vec3c dimensions;
			ImageDataType dataType;
			std::string reason;
			std::string datafile;
			bool isBigEndian;
			if (!getInfo(filename, dimensions, dataType, reason, datafile, isBigEndian))
				throw ITLException(reason);

			if (dataType != imageDataType<pixel_t>())
				throw ITLException(string("Pixel data type in the PCR file is ") + toString(dataType) + ", but image data type is " + toString(imageDataType<pixel_t>()) + ".");

			img.ensureSize(dimensions);
			raw::readNoParse(img, datafile);

			if (isBigEndian)
				swapByteOrder(img);
		}

		/**
		Reads part of an PCR file to the given image.
		NOTE: Does not support out-of-bounds start position.
		@param img Image where the data is placed. The size of the image defines the size of the block that is read.
		@param filename Path to the file to read.
		@param start Start location of the read. The size of the image defines the size of the block that is read.
		*/
		template<typename pixel_t> void readBlock(Image<pixel_t>& img, const std::string& filename, const Vec3c& start, bool showProgressInfo = false)
		{
			Vec3c dimensions;
			ImageDataType dataType;
			std::string reason;
			std::string datafile;
			bool isBigEndian;
			if (!getInfo(filename, dimensions, dataType, reason, datafile, isBigEndian))
				throw ITLException(reason);

			if (dataType != imageDataType<pixel_t>())
				throw ITLException(string("Pixel data type in the PCR file is ") + toString(dataType) + ", but image data type is " + toString(imageDataType<pixel_t>()) + ".");

			raw::readBlockNoParse(img, datafile, dimensions, start, showProgressInfo);

			if (isBigEndian)
				swapByteOrder(img);
		}

		namespace tests
		{
			void read();
		}
	}
}
