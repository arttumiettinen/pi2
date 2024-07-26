
#include <iostream>
#include <tiff.h>

#include "filesystem.h"
#include "zarr.h"
#include "json.h"

using namespace std;

namespace itl2
{
	namespace zarr
	{
		//TODO: remove isNativeByteOrder
		bool getInfo(const std::string& path,
			Vec3c& shape,
			bool& isNativeByteOrder,
			ImageDataType& dataType,
			Vec3c& chunkSize,
			std::list<codecs::ZarrCodec>& codecs,
			int& fillValue,
			std::string& reason)
		{
			shape = Vec3c();
			isNativeByteOrder = true;
			dataType = ImageDataType::Unknown;
			chunkSize = Vec3c();
			fillValue = 43;

			// Check that metadata file exists.
			string metadataFilename = zarr::internals::zarrMetadataFilename(path);
			if (!fs::exists(metadataFilename))
			{
				reason = "Metadata file does not exist.";
				return false;
			}

			// Read metadata and populate shape and pixel data type.
			ifstream in(metadataFilename);
			nlohmann::json j;

			try
			{
				in >> j;
			}

			catch (nlohmann::json::exception ex)
			{
				reason = string("Unable to parse zarr metadata: ") + ex.what();
				return false;
			}
			if (!j.contains("zarr_format"))
			{
				throw ITLException("zarr_format is missing in zarr metadata.");
			}
			else
			{
				int zarr_format = j["zarr_format"].get<int>();
				if (zarr_format != 3)
				{
					reason = "This zarr implementation supports only zarr_format 3.";
					return false;
				}
			}
			if (!j.contains("node_type"))
			{
				throw ITLException("node_type is missing in zarr metadata.");
			}
			else
			{
				string node_type = j["node_type"].get<string>();
				if (node_type != "array")
				{
					throw ITLException("This zarr implementation supports only node_type array.");
				}
			}
			if (!j.contains("shape"))
			{
				throw ITLException("shape is missing in zarr metadata.");
			}

			auto dims = j["shape"];
			if (dims.size() > 3)
			{
				reason = "This zarr implementation supports only 1-, 2-, or 3-dimensional datasets.";
				return false;
			}

			//transpose xyz -> yxz
			shape = Vec3c(1, 1, 1);
			shape[0] = dims[0].get<size_t>();
			if (dims.size() >= 2)
				shape[1] = dims[1].get<size_t>();
			if (dims.size() >= 3)
				shape[2] = dims[2].get<size_t>();

			cout << "getInfo shape: " << shape << endl;

			if (!j.contains("data_type"))
			{
				throw ITLException("data_type is missing in zarr metadata.");
			}
			try
			{
				string dataTypeString = j["data_type"].get<string>();
				dataType = fromString<ImageDataType>(dataTypeString);
				if (dataType == ImageDataType::Unknown)
				{
					throw ITLException("Could not identify datatype: " + dataTypeString);
				}
			}
			catch (ITLException& e)
			{
				reason = e.message();
				return false;
			}

			if (!j.contains("chunk_grid")
				|| !j["chunk_grid"].contains("configuration")
				|| !j["chunk_grid"]["configuration"].contains("chunk_shape")
				|| !j["chunk_grid"].contains("name")
				|| j["chunk_grid"]["name"].get<string>() != "regular"
				)
			{
				throw ITLException("chunk_grid is missing or malformed in zarr metadata.");
			}
			else
			{
				auto chunkDims = j["chunk_grid"]["configuration"]["chunk_shape"];
				if (chunkDims.size() != dims.size())
				{
					throw ITLException("Chunk shape and dataset shape contain different number of elements.");
				}

				chunkSize = Vec3c(1, 1, 1);
				chunkSize[0] = chunkDims[0].get<size_t>();
				if (chunkDims.size() >= 2)
					chunkSize[1] = chunkDims[1].get<size_t>();
				if (chunkDims.size() >= 3)
					chunkSize[2] = chunkDims[2].get<size_t>();

				cout << "getInfo chunkSize: " << chunkSize << endl;
			}

			if (!j.contains("chunk_key_encoding"))
			{
				throw ITLException("chunk_key_encoding is missing in zarr metadata.");
			}
			else
			{
				if (!j["chunk_key_encoding"].contains("name"))
				{
					throw ITLException("chunk_key_encoding name is missing in zarr metadata.");
					return false;
				}
				if (j["chunk_key_encoding"]["name"].get<string>() != "default")
				{
					reason = "This zarr implementation supports only default chunk_key_encoding.";
					return false;
				}
				if (j["chunk_key_encoding"].contains("configuration")
					&& !j["chunk_key_encoding"]["configuration"].contains("separator")
					&& j["chunk_key_encoding"]["configuration"]["separator"].get<string>() != "/")
				{
					reason = "This zarr implementation supports only default chunk_key_encoding with separator \"/\".";
					return false;
				}
			}

			if (!j.contains("fill_value"))
			{
				throw ITLException("fill_value is missing in zarr metadata.");
			}
			else
			{
				try
				{
					fillValue = j["fill_value"].get<int>();
				}
				catch (nlohmann::json::exception ex)
				{
					reason =
						string("Unable to parse fill_value in zarr metadata (this implementation only supports integers): ") +
							ex.what();
					return false;
				}
			}

			if (!j.contains("codecs"))
			{
				throw ITLException("codecs is missing in zarr metadata.");
			}
			else
			{
				if (!codecs::fromJSON(codecs, j["codecs"], reason))
				{
					return false;
				}
			}
			return true;
		}

		namespace internals
		{
			void writeMetadata(const std::string& path, const Vec3c& shape, ImageDataType dataType, const Vec3c& chunkSize, int fillValue, const std::list<codecs::ZarrCodec>& codecs)
			{
				nlohmann::json j =
					{
						{ "zarr_format", 3 },//zarr updated
						{ "node_type", "array" },//zarr updated
						{ "shape", { shape[0], shape[1], shape[2] }},//zarr updated
						{ "data_type", toString(dataType) },//zarr updated
						{ "chunk_grid", {{ "name", "regular" }, { "configuration", {{ "chunk_shape", { chunkSize[0], chunkSize[1], chunkSize[2] }}}}}},
						{ "chunk_key_encoding", {{ "name", "default" }, { "configuration", {{ "separator", "/" }}}}},
						{ "fill_value", fillValue },
						{ "codecs", {}}
					};
				for (auto& codec : codecs)
				{
					j["codecs"].push_back(codec.toJSON());
				}
				// TODO: allow other chunk_key_encoding
				// TODO: optional parameters
				string metadataFilename = zarr::internals::zarrMetadataFilename(path);
				ofstream of(metadataFilename, ios_base::trunc | ios_base::out);
				of << std::setw(4) << j << endl;
			}

			void handleExisting(const Vec3c& imageDimensions,
				ImageDataType imageDataType,
				const std::string& path,
				const Vec3c& chunkSize,
				int fillValue,
				const std::list<codecs::ZarrCodec>& codecs,
				bool deleteOldData)
			{
				// Delete old dataset if it exists.
				if (fs::exists(path))
				{
					bool oldIsNativeByteOrder;
					Vec3c oldDimensions;
					ImageDataType oldDataType;
					Vec3c oldChunkSize;
					std::list<codecs::ZarrCodec> oldCodecs;
					string dummyReason;
					int oldFillValue;
					if (!zarr::getInfo(path, oldDimensions, oldIsNativeByteOrder, oldDataType, oldChunkSize, oldCodecs, oldFillValue, dummyReason))
					{
						// The path does not contain a Zarr dataset.
						// If it is no known image, do not delete it.
						if (!io::getInfo(path, oldDimensions, oldDataType, dummyReason))
							throw ITLException(string("Unable to write a Zarr as the output folder already exists but cannot be verified to be an image: ") + path
								+ " Consider removing the existing dataset manually.");

						// Here the path does not contain Zarr but contains an image of known type.
						// Delete the old image.
						fs::remove_all(path);
					}
					else if (oldDimensions == imageDimensions &&
						oldIsNativeByteOrder &&
						oldDataType == imageDataType &&
						oldChunkSize == chunkSize &&
						oldCodecs == codecs &&
						oldFillValue == fillValue)
					{
						// The path contains a compatible Zarr dataset.
						// Delete it if we are not continuing a concurrent write.
						if (!fs::exists(internals::concurrentTagFile(path)))
							if (deleteOldData)
								fs::remove_all(path);
					}
					else
					{
						// The path contains an incompatible Zarr dataset. Delete it.
						if (fs::exists(internals::concurrentTagFile(path)) && !deleteOldData)
							throw ITLException(string("The output folder contains an incompatible zarr dataset that is currently being processed concurrently."));
						fs::remove_all(path);
					}
				}
			}

			std::vector<char> readBytesOfFile(std::string& filename)
			{
				std::ifstream ifs(filename, std::ios_base::binary | std::ios::ate);
				if (!ifs)
				{
					throw ITLException(std::string("Unable to open ") + filename + std::string(", ") + getStreamErrorMessage());
				}
				std::ifstream::pos_type pos = ifs.tellg();
				std::vector<char> result(pos);

				ifs.seekg(0, std::ios::beg);
				ifs.read(&result[0], pos);

				return result;
			}

			void writeBytesToFile(std::vector<char>& buffer, const std::string& filename, size_t startInFilePos)
			{
				createFoldersFor(filename);
				size_t fileSize = buffer.size();
				setFileSize(filename, fileSize);
				std::ofstream out(filename.c_str(), std::ios_base::in | std::ios_base::out | std::ios_base::binary);
				if (!out)
					throw ITLException(std::string("Unable to open ") + filename + std::string(", ") + getStreamErrorMessage());
				out.write(&buffer[startInFilePos], fileSize - startInFilePos);
				if (!out)
					throw ITLException(std::string("Unable to write to ") + filename + std::string(", ") + getStreamErrorMessage());
			}
		}

		namespace tests
		{
			void read()
			{
				Image<int32> fromDisk;
				zarr::read(fromDisk, "./testoutput/zarrita.zarr");
			}
			void write()
			{
				string path = "./testoutput/test_write.zarr";

				Image<uint16_t> img(Vec3c(2, 3, 4));
				ramp3(img);
				add(img, 10);
				zarr::write(img,
					path,
					zarr::DEFAULT_CHUNK_SIZE);
				Image<uint16_t> fromDisk;
				zarr::read(fromDisk, path);

				testAssert(equals(img, fromDisk), string("zarr test read and write"));
			}
			void writeBlock()
			{
				string path = "./testoutput/writeBlock.zarr";
				Vec3c size = Vec3c(10, 10, 10);
				Vec3c startBlock(2, 2, 2);
				Vec3c endBlock(3, 4, 5);

				Image<uint16_t> img(size, 0);
				zarr::write(img,path,zarr::DEFAULT_CHUNK_SIZE);
				add(img, 10);
				zarr::writeBlock(img, path,Vec3c(0,0,0),  size, startBlock, endBlock);

				Image<uint16_t> fromDisk;
				zarr::read(fromDisk, path);

				Image<uint16_t> expected(size, 0);
				draw(expected, AABoxsc::fromMinMax(Vec3<int>(startBlock), Vec3<int>(endBlock)), (uint16_t)10);

				testAssert(equals(img, fromDisk), string("zarr test read and write"));
			}
			void transpose()
			{
				string path = "./testoutput/test_transpose.zarr";

				Image<uint16_t> img(Vec3c(2, 3, 4));
				ramp3(img);
				add(img, 10);
				string transposeCodecConfig = R"({"order": [0, 2, 1]})";
				zarr::write(img,
					path,
					zarr::DEFAULT_CHUNK_SIZE,
					{ codecs::ZarrCodec(codecs::Name::Transpose, nlohmann::json::parse(transposeCodecConfig)), codecs::ZarrCodec(codecs::Name::Bytes), });

				Image<uint16_t> fromDisk;
				zarr::read(fromDisk, path);

				testAssert(equals(img, fromDisk), string("zarr test write transpose"));
			}

			void blosc()
			{
				string path = "./testoutput/test_blosc.zarr";

				Image<uint16_t> img(Vec3c(2, 5, 10));
				ramp(img, 0);
				add(img, 10);
				string bloscCodecConfig = R"({"cname": "lz4", "clevel": 1, "shuffle": "shuffle", "typesize": 4, "blocksize": 0})";
				zarr::write(img,
					path,
					zarr::DEFAULT_CHUNK_SIZE,
					{ codecs::ZarrCodec(codecs::Name::Bytes), codecs::ZarrCodec(codecs::Name::Blosc, nlohmann::json::parse(bloscCodecConfig)) });

				Image<uint16_t> fromDisk;
				zarr::read(fromDisk, path);

				testAssert(equals(img, fromDisk), string("zarr test write transpose"));
			}

		}
	}
}