
#include "nrrd.h"

#include "projections.h"

#include <array>

using namespace std;

namespace itl2
{
	namespace nrrd
	{
		namespace internals
		{
			/**
			Converts NRRD image data type string to ImageDataType value.
			*/
			ImageDataType fromNRRDType(const string& type)
			{
				if (type == "signed char" || type == "int8" || type == "int8_t")
					return ImageDataType::Int8;

				if(type == "uchar" || type == "unsigned char" || type == "uint8" || type == "uint8_t")
					return ImageDataType::UInt8;

				if (type == "short" || type == "short int" || type == "signed short" || type == "signed short int" || type == "int16" || type == "int16_t")
					return ImageDataType::Int16;

				if (type == "ushort" || type == "unsigned short" || type == "unsigned short int" || type == "uint16" || type == "uint16_t")
					return ImageDataType::UInt16;

				if (type == "int" || type == "signed int" || type == "int32" || type == "int32_t")
					return ImageDataType::Int32;

				if (type == "uint" || type == "unsigned int" || type == "uint32" || type == "uint32_t")
					return ImageDataType::UInt32;

				if (type == "longlong" || type == "long long" || type == "long long int" || type == "signed long long" || type == "signed long long int" || type == "int64" || type == "int64_t")
					return ImageDataType::Int64;

				if (type == "ulonglong" || type == "unsigned long long" || type == "unsigned long long int" || type == "uint64" || type == "uint64_t")
					return ImageDataType::UInt64;

				if (type == "float")
					return ImageDataType::Float32;

				// These are also possible
				//"double"
				//"block"

				return ImageDataType::Unknown;
			}

			std::string toNRRDType(ImageDataType dt, bool& writePixelSize)
			{
				writePixelSize = false;
				switch (dt)
				{
				case ImageDataType::UInt8: return "uint8";
				case ImageDataType::UInt16: return "uint16";
				case ImageDataType::UInt32: return "uint32";
				case ImageDataType::UInt64: return "uint64";
				case ImageDataType::Int8: return "int8";
				case ImageDataType::Int16: return "int16";
				case ImageDataType::Int32: return "int32";
				case ImageDataType::Int64: return "int64";
				case ImageDataType::Float32: return "float";
				default:
					// complex32_t, other data types
					writePixelSize = true;
					return "block";
				}
			}
		}

		bool getInfo(const std::string& filename, Vec3c& dimensions, ImageDataType& dataType, string& failReason, size_t& headerSize, string& datafile, bool& isBigEndian)
		{
			headerSize = 0;
			failReason = "";
			datafile = "";
			isBigEndian = false;

			ifstream in(filename.c_str(), ios_base::in | ios_base::binary);

			if (!in)
			{
				failReason = string("Unable to open file ") + filename;
				return false;
			}

			string line;

			// Check header
			std::array<char, 101> tmp{};
			in.getline(&tmp[0], tmp.size() - 1);
			line = string(&tmp[0]);

			if (!startsWithIgnoreCase(line, "NRRD"))
			{
				failReason = "The file does not contain NRRD header.";
				return false;
			}

			bool dimensionsFound = false;
			bool dataTypeFound = false;
			while (true)
			{
				std::getline(in, line);

				trim(line);

				if (line.length() <= 0)
					break; // End of header

				if (startsWithIgnoreCase(line, "#"))
					continue; // Comment

				if (contains(line, ":="))
					continue; // Key-value pair, we don't handle these

				
				char delim;
				string name = getToken(line, ":", delim);
				
				if (delim == 0)
				{
					failReason = string("NRRD header contains invalid line: ") + line;
					return false;
				}

				trim(line);
				toLower(name);

				if (name == "sizes")
				{
					vector<string> parts = split(line, false, ' ', true);
					if (parts.size() > 3)
					{
						failReason = string("Only 1-, 2- and 3-dimensional images are supported, but the file contains ") + toString(parts.size()) + "-dimensional image.";
						return false;
					}

					dimensions = Vec3c(1, 1, 1);

					for(size_t n = 0; n < parts.size(); n++)
					{
						dimensions[n] = fromString<coord_t>(parts[n]);
					}

					if (dimensions.min() < 0)
					{
						failReason = string("The file specifies invalid image dimensions: ") + toString(dimensions);
						return false;
					}

					dimensions = max(dimensions, Vec3c(1, 1, 1));

					dimensionsFound = true;
				}
				else if (name == "type")
				{
					dataType = internals::fromNRRDType(line);

					dataTypeFound = true;
				}
				else if (name == "encoding")
				{
					if (line != "raw")
					{
						failReason = string("Only raw encoding is currently supported, but the file contains data in ") + line + " encoding.";
						return false;
					}
				}
				else if (name == "data file" || name == "datafile")
				{
					datafile = line;
				}
				else if (name == "endian")
				{
					if (line == "big")
					{
						isBigEndian = true;
					}
					else if (line != "little")
					{
						failReason = string("The file specifies unsupported endianness value: ") + line;
						return false;
					}
				}
			}

			headerSize = in.tellg();

			return dimensionsFound && dataTypeFound;
		}

		namespace tests
		{
			void readWrite()
			{
				Image<uint16_t> img;
				raw::read(img, "./input_data/t1-head");

				nrrd::writed(img, "./nrrd/t1-head");

				Image<uint16_t> img2;
				nrrd::read(img2, "./nrrd/t1-head.nrrd");

				testAssert(equals(img, img2), "NRRD read and write");
			}
		}
	}
}