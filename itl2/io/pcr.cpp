
#include "pcr.h"

#include "inireader.h"
#include "projections.h"

using namespace std;

namespace itl2
{
	namespace pcr
	{
		namespace internals
		{
			/**
			Converts NRRD image data type string to ImageDataType value.
			*/
			ImageDataType fromPCRType(int type, const string& datafile, const Vec3c& dimensions, bool& isBigEndian)
			{
				isBigEndian = false;

				if (type == 5)
					return ImageDataType::UInt16;

				// TODO: I don't know any other data types at the moment, so we just try to parse them from the data file size.
				size_t pixelSizeBytes = 0;
				return raw::internals::estimateDataType(datafile, dimensions, pixelSizeBytes);
			}
		}

		bool getInfo(const std::string& filename, Vec3c& dimensions, ImageDataType& dataType, string& failReason, string& datafile, bool& isBigEndian)
		{
			failReason = "";
			datafile = "";
			isBigEndian = false;
			dataType = ImageDataType::Unknown;

			if (fs::path(filename).extension().string() != ".pcr")
			{
				failReason = "Not a .pcr file.";
				return false;
			}

			if (fileSize(filename) > 1 * 1024 * 1024)
			{
				failReason = "File is too large to be a PCR header.";
				return false;
			}

			INIReader reader(filename);

			if (reader.parseError() < 0)
			{
				failReason = string("Unable to open file ") + filename;
				return false;
			}

			if (reader.parseError() > 0)
			{
				failReason = string("Unable to parse contents of file ") + filename + string(". Problematic line number is ") + toString(reader.parseError());
				return false;
			}

			coord_t width = reader.get("VolumeData", "Volume_SizeX", (coord_t)0);
			coord_t height = reader.get("VolumeData", "Volume_SizeY", (coord_t)0);
			coord_t depth = reader.get("VolumeData", "Volume_SizeZ", (coord_t)0);
			int format = reader.get("VolumeData", "Format", -1);
			
			// One data file candidate is the file read from the header file.
			fs::path datafile1 = reader.get("VolumeData", "VOL_File", string());
			// Second candidate is file name read from the header file + header file path
			fs::path datafile2 = fs::path(filename).parent_path() / datafile1.filename();
			// Third candidate is header file name but with .vol extension
			fs::path datafile3 = fs::path(filename).replace_extension(".vol");

			if (fs::exists(datafile2))
			{
				datafile = datafile2.string();
			}
			else if (fs::exists(datafile1))
			{
				datafile = datafile1.string();
			}
			else if (fs::exists(datafile3))
			{
				datafile = datafile3.string();
			}
			else
			{
				failReason = string("Image data file not found. It was searched from ") + datafile2.string() + string(", ") + datafile1.string() + string(" and ") + datafile3.string() + string(".");
				return false;
			}
			
			dimensions = Vec3c(width, height, depth);

			if (dimensions.min() <= 0)
			{
				failReason = string("Invalid image dimensions: ") + toString(dimensions);
				return false;
			}

			dataType = internals::fromPCRType(format, datafile, dimensions, isBigEndian);

			if (dataType == ImageDataType::Unknown)
			{
				failReason = "Unknown image data type.";
				return false;
			}

			return true;
		}

		namespace tests
		{
			void read()
			{
				Image<uint16_t> img;
				raw::read(img, "../test_input_data/t1-head");

				Image<uint16_t> img2;
				pcr::read(img2, "../test_input_data/t1-head.pcr");

				testAssert(equals(img, img2), "PCR read");
			}
		}
	}
}