#pragma once

#include "image.h"
#include "stringutils.h"
#include "math/vec3.h"
#include "io/raw.h"

#include <string>
#include <array>

namespace itl2
{
	namespace vol
	{
		namespace internals
		{
			/**
			Tests if the string consists of characters that are allowed in .vol header.
			*/
			inline bool isOkVolHeaderCharacters(const string& s)
			{
				const string OK_CHARS = " abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890:-";
				for (size_t i = 0; i < s.length(); i++)
				{
					char c = s[i];
					if (OK_CHARS.find(c) == string::npos)
						return false;
				}

				return true;
			}
		}

		/**
		Gets information about .vol file in the given path.
		*/
		inline bool getInfo(const std::string& filename, Vec3c& dimensions, ImageDataType& dataType, string& endianness, size_t& headerSize, string& reason)
		{
			std::ifstream in(filename.c_str(), std::ios_base::in | std::ios_base::binary);

			if (!in)
			{
				reason = "Unable to open file.";
				return false;
				//throw ITLException(std::string("Unable to open ") + filename + string(", ") + getStreamErrorMessage());
			}

			coord_t voxelSize = 0;

			string line;
			while (true)
			{
				// Don't use bare getline here. If the file is NOT a .vol file, there might be tens of gigabytes of data read before the first
				// newline is encountered!
				//std::getline(in, line);

				std::array<char, 101> tmp{};
				in.getline(&tmp[0], tmp.size()-1);
				line = string(&tmp[0]);

				if (line.length() >= 100)
				{
					reason = "Excessively long line in file header.";
					return false; // .vol file does not PROBABLY contain this long lines!
				}

				trim(line);

				if (line.length() <= 0 || line == ".")
					break;

				if (!internals::isOkVolHeaderCharacters(line))
				{
					reason = "Invalid characters in file header.";
					return false;
				}

				toLower(line);

				std::vector<string> items = split(line, false, ':', true);
				if (items.size() == 2)
				{
					std::string key = items[0];
					std::string value = items[1];

					if (key == "x")
						dimensions.x = fromString<coord_t>(value);
					else if (key == "y")
						dimensions.y = fromString<coord_t>(value);
					else if (key == "z")
						dimensions.z = fromString<coord_t>(value);
					else if (key == "voxel-size" || key == "lvoxel-size")
						voxelSize = fromString<coord_t>(value);
					//else if (key == "voxel-endian" || key == "lvoxel-endian")
					//	endianness = value;
					else if(key == "int-endian")
						endianness = value;
					//else if (key == "alpha-color")
					//	cout << "Warning: .vol file specifies alpha color but it is not supported." << endl;
				}
				else
				{
					reason = string("Unexpected key in file header: ") + line;
					return false;
					//cout << "Warning: Unexpected key in .vol file header: " << line << endl;
				}
			}

			headerSize = in.tellg();
			
			if (voxelSize == 1)
				dataType = ImageDataType::UInt8;
			else if (voxelSize == 2)
				dataType = ImageDataType::UInt16;
			else if (voxelSize == 4)
				dataType = ImageDataType::UInt32;
			else
				dataType = ImageDataType::Unknown;

			return dataType != ImageDataType::Unknown;
		}

		/**
		Gets information about .vol file in the given path.
		*/
		inline bool getInfo(const std::string& filename, Vec3c& dimensions, ImageDataType& dataType, std::string& reason)
		{
			std::string endianness;
			size_t headerSize;
			return getInfo(filename, dimensions, dataType, endianness, headerSize, reason);
		}

		template<typename pixel_t> void getInfoAndCheck(const std::string& filename, Vec3c& dimensions, ImageDataType& dataType, size_t& headerSize)
		{
			std::string endianness;
			std::string reason;
			if (!getInfo(filename, dimensions, dataType, endianness, headerSize, reason))
				throw ITLException(std::string("Not a .vol file: ") + filename + ". " + reason);

			// Check data
			if (dimensions.x <= 0 || dimensions.y <= 0 || dimensions.z <= 0)
				throw ITLException(std::string(".vol file contains invalid size specification: ") + toString(dimensions));

			if (endianness != "0123")
				throw ITLException(std::string(".vol file contains unsupported endianness value: ") + toString(endianness));

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
			getInfoAndCheck<pixel_t>(filename, dimensions, dataType, headerSize);

			// Read data
			img.ensureSize(dimensions);
			raw::readNoParse(img, filename, headerSize);
		}

		/**
		Reads a block of .vol file from the given path.
		*/
		template<typename pixel_t> void readBlock(Image<pixel_t>& img, const std::string& filename, const Vec3c& start, bool showProgressInfo = false)
		{
			Vec3c dimensions;
			ImageDataType dataType;
			size_t headerSize;
			getInfoAndCheck<pixel_t>(filename, dimensions, dataType, headerSize);

			raw::readBlockNoParse(img, filename, dimensions, start, showProgressInfo, headerSize);
		}

		namespace tests
		{
			void volio();
		}
	}

}
