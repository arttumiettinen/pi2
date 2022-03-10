
#include <iostream>

#include "nn5.h"
#include "json.h"
#include "byteorder.h"
#include "generation.h"

using namespace std;

namespace itl2
{
	namespace nn5
	{
		bool getInfo(const std::string& path, Vec3c& dimensions, bool& isNativeByteOrder, ImageDataType& dataType, Vec3c& chunkSize, NN5Compression& compression, std::string& reason)
		{
			dimensions = Vec3c();
			isNativeByteOrder = true;
			dataType = ImageDataType::Unknown;
			chunkSize = Vec3c();
			compression = NN5Compression::Raw;

			// Check that metadata file exists.
			string metadataFilename = internals::nn5MetadataFilename(path);
			if (!fs::exists(metadataFilename))
			{
				reason = "Metadata file does not exist.";
				return false;
			}

			// Read metadata and populate dimensions and pixel data type.
			ifstream in(metadataFilename);
			nlohmann::json j;

			try
			{
				in >> j;
			}
			catch (nlohmann::json::exception ex)
			{
				reason = string("Unable to parse nn5 metadata: ") + ex.what();
				return false;
			}

			if (!j.contains("Dimensions"))
			{
				reason = "Dimensions are missing in nn5 metadata.";
				return false;
			}

			auto dims = j["Dimensions"];
			if (dims.size() > 3)
			{
				reason = "This nn5 implementation supports only 1-, 2-, or 3-dimensional datasets.";
				return false;
			}

			dimensions = Vec3c(1, 1, 1);
			dimensions[0] = dims[0].get<size_t>();
			if (dims.size() >= 2)
				dimensions[1] = dims[1].get<size_t>();
			if (dims.size() >= 3)
				dimensions[2] = dims[2].get<size_t>();


			if (!j.contains("Data type"))
			{
				reason = "Data type is missing in nn5 metadata.";
				return false;
			}

			try
			{
				dataType = fromString<ImageDataType>(j["Data type"].get<string>());
			}
			catch (ITLException& e)
			{
				reason = e.message();
				return false;
			}

			isNativeByteOrder = true;
			if (!j.contains("Byte order"))
			{
				// Assume native byte order
			}
			else
			{
				string boStr = j["Byte order"].get<string>();

				Endianness bo;
				try
				{
					bo = fromString<Endianness>(boStr);
				}
				catch (ITLException& e)
				{
					reason = e.message();
					return false;
				}

				if (bo != nativeByteOrder())
					isNativeByteOrder = false;
			}



			if (!j.contains("Chunk dimensions"))
			{
				// Assume one chunk
				chunkSize = dimensions;
			}
			else
			{
				auto chunkDims = j["Chunk dimensions"];
				if (chunkDims.size() != dims.size())
				{
					reason = "Chunk dimensions and dataset dimensions contain different number of elements.";
					return false;
				}

				chunkSize = Vec3c(1, 1, 1);
				chunkSize[0] = chunkDims[0].get<size_t>();
				if (chunkDims.size() >= 2)
					chunkSize[1] = chunkDims[1].get<size_t>();
				if (chunkDims.size() >= 3)
					chunkSize[2] = chunkDims[2].get<size_t>();
			}



			if (!j.contains("Compression"))
			{
				// Assume raw
				compression = NN5Compression::Raw;
			}
			else
			{
				string compStr = j["Compression"].get<string>();

				try
				{
					compression = fromString<NN5Compression>(compStr);
				}
				catch (ITLException& e)
				{
					reason = e.message();
					return false;
				}
			}

			return true;
		}

		namespace internals
		{
			void writeMetadata(const std::string& path, const Vec3c& dimensions, ImageDataType dataType, const Vec3c& chunkSize, NN5Compression compression)
			{
				nlohmann::json j;
				j["Dimensions"][0] = dimensions[0];
				j["Dimensions"][1] = dimensions[1];
				j["Dimensions"][2] = dimensions[2];
				j["Data type"] = toString(dataType);
				j["Byte order"] = toString(nativeByteOrder());
				j["Chunk dimensions"][0] = chunkSize[0];
				j["Chunk dimensions"][1] = chunkSize[1];
				j["Chunk dimensions"][2] = chunkSize[2];
				j["Compression"] = toString(compression);

				string metadataFilename = internals::nn5MetadataFilename(path);
				ofstream of(metadataFilename, ios_base::trunc | ios_base::out);
				of << std::setw(4) << j << endl;
			}

			vector<string> getFileList(const string& dir)
			{
				vector<string> filenames;

				if (fs::is_directory(dir)) // Note: This is required in Linux, or otherwise we get an exception for non-existing directories.
				{
					for (auto& p : fs::directory_iterator(dir))
					{
						if (p.is_regular_file())
						{
							filenames.push_back(p.path().filename().string());
						}
					}
				}

				return filenames;
			}
		}


		namespace tests
		{

			void nn5Metadata()
			{
				string data =
					R"END({
"Dimensions": [10, 20, 30],
"Data type": "uint8",
"Byte order": "Little endian",
"Chunk dimensions": [8, 9, 11],
"Compression": "Raw",
"extra data": "doh doh boh boh"
})END";
				string path = "./nn5_metadata";
				writeText(path + "/metadata.json", data);

				Vec3c dimensions;
				Vec3c chunkSize;
				ImageDataType datatype;
				string reason;
				bool isNativeByteOrder;
				nn5::NN5Compression compression;
				bool result = nn5::getInfo(path, dimensions, isNativeByteOrder, datatype, chunkSize, compression, reason);

				testAssert(result == true, "nn5 getinfo result");
				testAssert(dimensions == Vec3c(10, 20, 30), "nn5 dimensions");
				testAssert(datatype == ImageDataType::UInt8, "nn5 datatype");
				testAssert(isNativeByteOrder == true, "nn5 native byte order");
				testAssert(chunkSize == Vec3c(8, 9, 11), "nn5 chunk size");
				testAssert(reason == "", "nn5 reason");

				nn5::internals::writeMetadata(path, dimensions, datatype, chunkSize, compression);
				result = nn5::getInfo(path, dimensions, isNativeByteOrder, datatype, chunkSize, compression, reason);
				testAssert(result == true, "nn5 getinfo result from written file");
				testAssert(dimensions == Vec3c(10, 20, 30), "nn5 dimensions from written file");
				testAssert(datatype == ImageDataType::UInt8, "nn5 datatype from written file");
				testAssert(isNativeByteOrder == true, "nn5 native byte order from written file");
				testAssert(chunkSize == Vec3c(8, 9, 11), "nn5 chunk size from written file");
				testAssert(reason == "", "nn5 reason from written file");

			}

			void onenn5ioTest(const Image<uint16_t>& img, const Vec3c& chunkSize, NN5Compression compression)
			{
				cout << "Testing chunk size " << chunkSize << ", compression = " << toString(compression) << endl;

				// Write
				nn5::write(img, "./nn5_testimage", chunkSize, compression);

				// Read
				Image<uint16_t> read;
				nn5::read(read, "./nn5_testimage");

				raw::writed(read, "./nn5_results/read_from_disk");

				testAssert(equals(img, read), string("nn5 read image compared to written image, chunk size = ") + toString(chunkSize));
			}

			void nn5io()
			{
				coord_t w = 100;
				coord_t h = 200;
				coord_t d = 300;

				Image<uint16_t> img(w, h, d);
				ramp3(img);
				raw::writed(img, "./nn5_results/true");

				srand(1212);

				for (coord_t n = 0; n < 100; n++)
				{
					Vec3c chunkSize(randc(10, w + 10), randc(10, h + 10), randc(10, d + 10));

					onenn5ioTest(img, chunkSize, NN5Compression::Raw);
					onenn5ioTest(img, chunkSize, NN5Compression::LZ4);
				}
			}

			void nn5BlockIoOneTest(NN5Compression compression, const Vec3c& chunkSize)
			{
				cout << "Chunk size = " << chunkSize << ", compression = " << toString(compression) << endl;

				Vec3c dimensions(100, 200, 300);

				Image<uint16_t> img(dimensions);
				ramp3(img);

				Vec3c blockStart(10, 20, 30);
				Vec3c blockSize(50, 60, 70);

				// Write entire image and read
				nn5::write(img, "./nn5_block_io/entire_image", chunkSize, compression);
				Image<uint16_t> entireFromDisk;
				nn5::read(entireFromDisk, "./nn5_block_io/entire_image");
				testAssert(equals(entireFromDisk, img), "NN5 entire image read/write cycle");


				// Write block and read
				fs::remove_all("./nn5_block_io/block");
				nn5::writeBlock(img, "./nn5_block_io/block", chunkSize, compression, Vec3c(0, 0, 0), blockSize, blockStart, blockSize);

				// Generate ground truth block by cropping the image.
				Image<uint16_t> gtBlock(blockSize);
				crop(img, gtBlock, blockStart);

				// Read entire file written using writeBlock and check against ground truth.
				Image<uint16_t> fileBlock;
				nn5::read(fileBlock, "./nn5_block_io/block");
				testAssert(equals(fileBlock, gtBlock), "NN5 writeBlock");


				// Read using readBlock and compare to the ground truth.
				Image<uint16_t> readBlockResult(blockSize);
				nn5::readBlock(readBlockResult, "./nn5_block_io/entire_image", blockStart);
				testAssert(equals(readBlockResult, gtBlock), "NN5 readBlock");
			}

			void nn5BlockIo()
			{
				nn5BlockIoOneTest(NN5Compression::Raw, Vec3c(50, 102, 99));
				nn5BlockIoOneTest(NN5Compression::LZ4, Vec3c(50, 102, 99));

				nn5BlockIoOneTest(NN5Compression::Raw, Vec3c(40, 30, 20));
				nn5BlockIoOneTest(NN5Compression::LZ4, Vec3c(40, 30, 20));
			}
		}
	}
}