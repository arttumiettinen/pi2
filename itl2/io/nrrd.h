#pragma once

#include <string>
#include <iostream>
#include <filesystem>

#include "math/vec3.h"
#include "io/imagedatatype.h"
#include "image.h"
#include "io/raw.h"
#include "pointprocess.h"

namespace fs = std::filesystem;

namespace itl2
{
	namespace nrrd
	{
		namespace internals
		{
			/**
			Converts image data type to corresponding NRRD data type string.
			*/
			std::string toNRRDType(ImageDataType dt, bool& writePixelSize);
		}

		/**
		Gets information about NRRD file in the given path.
		*/
		bool getInfo(const std::string& filename, Vec3c& dimensions, ImageDataType& dataType, string& failReason, size_t& headerSize, string& datafile, bool& isBigEndian);

		/**
		Gets information about NRRD file in the given path.
		*/
		inline bool getInfo(const std::string& filename, Vec3c& dimensions, ImageDataType& dataType, string& reason)
		{
			std::size_t headerSize;
			std::string datafile;
			bool isBigEndian;
			return getInfo(filename, dimensions, dataType, reason, headerSize, datafile, isBigEndian);
		}

		/**
		Reads NRRD file from the given path.
		*/
		template<typename pixel_t> void read(Image<pixel_t>& img, const std::string& filename)
		{
			// Read header
			Vec3c dimensions;
			ImageDataType dataType;
			std::string reason;
			size_t headerSize;
			std::string datafile;
			bool isBigEndian;
			if (!getInfo(filename, dimensions, dataType, reason, headerSize, datafile, isBigEndian))
				throw ITLException(reason);

			if (dataType != imageDataType<pixel_t>())
			{
				if (!(dataType == ImageDataType::Unknown && imageDataType<pixel_t>() == ImageDataType::Complex32))
					throw ITLException(string("Pixel data type in NRRD file is ") + toString(dataType) + ", but image data type is " + toString(imageDataType<pixel_t>()) + ".");
			}

			// Read data
			img.ensureSize(dimensions);

			if (datafile == "")
			{
				raw::readNoParse(img, filename, headerSize);
			}
			else
			{
				if (datafile == "LIST")
					throw ITLException("LIST type data file set is not supported at the moment.");

				datafile = fs::path(filename).replace_filename(fs::path(datafile)).string();

				if (!fileExists(datafile))
					throw ITLException(std::string("Data file not found or specifies unsupported image sequence: ") + datafile);

				raw::readNoParse(img, datafile);
			}

			if (isBigEndian)
				swapByteOrder(img);
		}

		/*
		Writes a NRRD file.
		*/
		template<typename pixel_t> void write(const Image<pixel_t>& img, const std::string& filename)
		{
			createFoldersFor(filename);

			std::ofstream out(filename.c_str(), std::ios_base::out | std::ios_base::trunc);

			if (!out)
				throw ITLException(std::string("Unable to write to ") + filename + std::string(", ") + getStreamErrorMessage());

			size_t dim = img.dimensionality();

			// Write header. We use the simplest NRRD version as that contains everything we need.
			bool writePixelSize;
			out << "NRRD0001" << std::endl;
			out << "type: " << internals::toNRRDType(imageDataType<pixel_t>(), writePixelSize) << std::endl;
			if (writePixelSize)
				out << "blocksize: " << sizeof(pixel_t) << std::endl;
			out << "dimension: " << dim << std::endl;
			out << "sizes: " << img.dimension(0);
			if (dim > 1)
				out << " " << img.dimension(1);
			if (dim > 2)
				out << " " << img.dimension(2);
			out << std::endl;
			out << "encoding: raw" << std::endl;
			out << "endian: little" << std::endl;
			out << "" << std::endl;

			out.close();

			// Write data
			raw::write(img, filename, false);
		}

		/**
		Write a NRRD file, adds .nrrd to the file name if it does not end with .nnrd.
		*/
		template<typename pixel_t> void writed(const Image<pixel_t>& img, const std::string& filename)
		{
			if (endsWithIgnoreCase(filename, ".nrrd"))
				write(img, filename);
			else
				write(img, filename + ".nrrd");
		}

		namespace tests
		{
			void readWrite();
		}
	}
}
