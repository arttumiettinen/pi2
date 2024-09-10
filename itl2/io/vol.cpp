
#include "io/vol.h"
#include "pointprocess.h"
#include "testutils.h"
#include "math/vec3.h"
#include "inireader.h"

namespace itl2
{

	namespace vol
	{
		namespace internals
		{
			/**
			Tests if the string consists of characters that are allowed in .vol header.
			*/
			bool isOkVolHeaderCharacters(const string& s)
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

			bool getInfoVol(const std::string& filename, Vec3c& dimensions, ImageDataType& dataType, bool& isBigEndian, size_t& headerSize, string& reason)
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
					in.getline(&tmp[0], tmp.size() - 1);
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
						else if (key == "int-endian")
						{
							//endianness = value;
							if (value == "0123")
								isBigEndian = false;
							else if (value == "1032")
								isBigEndian = true;
							else
							{
								reason = string(".vol file contains unsupported endianness value: ") + toString(value);
								return false;
							}
						}
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

				if (dataType == ImageDataType::Unknown)
				{
					reason = "Unknown image data type.";
					return false;
				}

				return true;
			}

			bool getInfoPyHSTVol(const string& imagefile, const string& volInfo, Vec3c& dimensions, ImageDataType& dataType, bool& isBigEndian, size_t& headerSize, string& reason)
			{
				headerSize = 0;

				if (fileSize(volInfo) > 1 * 1024 * 1024)
				{
					reason = "File is too large to be a PyHST .vol.info header.";
					return false;
				}

				INIReader reader(volInfo);

				if (reader.parseError() < 0)
				{
					reason = "Unable to open .vol.info file.";
					return false;
				}

				if (reader.parseError() > 0)
				{
					reason = string("Unable to parse contents of file ") + volInfo + string(". Problematic line number is ") + toString(reader.parseError());
					return false;
				}

				coord_t width = reader.get("", "NUM_X", (coord_t)0);
				coord_t height = reader.get("", "NUM_Y", (coord_t)0);
				coord_t depth = reader.get("", "NUM_Z", (coord_t)0);
				string byteorder = reader.get<string>("", "BYTEORDER", "LOWBYTEFIRST");

				if (byteorder == "LOWBYTEFIRST")
					isBigEndian = false;
				else
					isBigEndian = true;

				dimensions = Vec3c(width, height, depth);

				if (dimensions.min() <= 0)
				{
					reason = string("Invalid image dimensions: ") + toString(dimensions);
					return false;
				}

				size_t psize;
				dataType = raw::internals::estimateDataType(imagefile, dimensions, psize);

				if (dataType == ImageDataType::Unknown)
				{
					reason = "Unknown image data type.";
					return false;
				}

				return true;
			}
		}

		/**
		Gets information about .vol file in the given path.
		*/
		bool getInfo(const std::string& filename, Vec3c& dimensions, ImageDataType& dataType, bool& isBigEndian, size_t& headerSize, string& reason)
		{
			// This is required in some Linux systems to differentiate files from directories.
			if (!fs::is_regular_file(filename))
			{
				reason = "Not a file.";
				return false;
			}

			// Check if this is Coeurjolly's .vol file or PyHST .vol file.
			// The latter is accompanied with .vol.info filem and the former not.
			string volInfo = filename + ".info";
			if (fs::is_regular_file(volInfo))
			{
				// This is PyHST .vol file.
				return internals::getInfoPyHSTVol(filename, volInfo, dimensions, dataType, isBigEndian, headerSize, reason);
			}
			else
			{
				// This is Coeurjolly's .vol file.
				return internals::getInfoVol(filename, dimensions, dataType, isBigEndian, headerSize, reason);
			}

		}

		namespace tests
		{
			void volio()
			{
				Vec3c dims;
				ImageDataType dt;
				bool end;
				size_t hs;
				string reason;
				vol::getInfo("../test_input_data/simple_structures.vol", dims, dt, end, hs, reason);

				Image<uint8_t> img;
				vol::read(img, "../test_input_data/simple_structures.vol");

				Image<float32_t> gt;
				raw::read(gt, "../test_input_data/simple_structures_128x128x128.raw");
				multiply(gt, 255);

				itl2::checkDifference(img, gt, "Same structure from .vol and .raw file.");

				// Try PyHST .vol files, too.
				Image<uint16_t> img16;
				vol::read(img16, "../test_input_data/t1-head.vol");
				Image<uint16_t> gt16;
				raw::read(gt16, "../test_input_data/t1-head");
				itl2::checkDifference(img16, gt16, "Same structure from PyHST .vol and .raw file.");
			}
		}
	}

}
