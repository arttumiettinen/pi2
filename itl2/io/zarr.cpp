
#include <iostream>

#include "zarr.h"
#include "json.h"
#include "nn5.h"

using namespace std;

namespace itl2
{
	namespace zarr
	{
		//TODO: remove isNativeByteOrder
		bool getInfo(const std::string& path, Vec3c& shape, bool& isNativeByteOrder, ImageDataType& dataType, Vec3c& chunkSize, std::list<ZarrCodec>& codecs, int& fillValue, std::string& reason)
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
				reason = "zarr_format is missing in zarr metadata.";
				return false;
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
				reason = "node_type is missing in zarr metadata.";
				return false;
			}
			else
			{
				string node_type = j["node_type"].get<string>();
				if (node_type != "array")
				{
					reason = "This zarr implementation supports only node_type array.";
					return false;
				}
			}
			if (!j.contains("shape"))
			{
				reason = "shape is missing in zarr metadata.";
				return false;
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
				reason = "data_type is missing in zarr metadata.";
				return false;
			}
			try
			{
				dataType = fromString<ImageDataType>(j["data_type"].get<string>());
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
				reason = "chunk_grid is missing or malformed in zarr metadata.";
				return false;
			}
			else
			{
				auto chunkDims = j["chunk_grid"]["configuration"]["chunk_shape"];
				if (chunkDims.size() != dims.size())
				{
					reason = "Chunk shape and dataset shape contain different number of elements.";
					return false;
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
				reason = "chunk_key_encoding is missing in zarr metadata.";
				return false;
			}
			else
			{
				if (!j["chunk_key_encoding"].contains("name"))
				{
					reason = "chunk_key_encoding name is missing in zarr metadata.";
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
				reason = "fill_value is missing in zarr metadata.";
				return false;
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
				reason = "codecs is missing in zarr metadata.";
				return false;
			}
			else
			{
				int numberArrayBytesCodecs = 0;
				for (auto& codec : j["codecs"])
				{
					if (!codec.contains("name"))
					{
						reason = "codec name is missing in zarr metadata.";
						return false;
					}
					try
					{
						nlohmann::json codecConfig = {};
						if (codec.contains("configuration"))
						{
							codecConfig = codec["configuration"];
						}
						zarr::ZarrCodecName zarrCodecName = fromString<zarr::ZarrCodecName>(codec["name"].get<string>());
						zarr::ZarrCodec zarrCodec = zarr::ZarrCodec(zarrCodecName, codecConfig);
						codecs.push_back(zarrCodec);
						switch (zarrCodec.type)
						{
						case ZarrCodecType::ArrayArrayCodec:
							if (numberArrayBytesCodecs > 0)
							{
								reason = "ArrayArrayCodec cannot be used after ArrayBytesCodec.";
								return false;
							}
							break;
						case ZarrCodecType::ArrayBytesCodec:
							numberArrayBytesCodecs++;
							break;
						case ZarrCodecType::BytesBytesCodec:
							if (numberArrayBytesCodecs < 1)
							{
								reason = "ArrayBytesCodec must be used before BytesBytesCodec.";
								return false;
							}
							break;
						default:
							reason = "Unknown codec type.";
							return false;
						}
					}
					catch (ITLException& e)
					{
						reason = e.message();
						return false;
					}
				}
				if (numberArrayBytesCodecs != 1)
				{
					reason = "Exactly one ArrayBytesCodec was expected in the codecs list, got " + to_string(numberArrayBytesCodecs) + ".";
					return false;
				}
			}
			return true;
		}

		namespace internals
		{
			void writeMetadata(const std::string& path, const Vec3c& shape, ImageDataType dataType, const Vec3c& chunkSize, int fillValue, const std::list<ZarrCodec>& codecs)
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
				const std::list<ZarrCodec>& codecs,
				bool deleteOldData)
			{
				// Delete old dataset if it exists.
				if (fs::exists(path))
				{
					bool oldIsNativeByteOrder;
					Vec3c oldDimensions;
					ImageDataType oldDataType;
					Vec3c oldChunkSize;
					std::list<ZarrCodec> oldCodecs;
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

			/**
			Reads block position [X, Y, Z] from a filename in format 'chunk_X-Y-Z_something' or 'chunk_X-Y-Z.something'.
			*/
			Vec3c parsePosition(const string& filename)
			{
				string name = fs::path(filename).filename().string();
				vector<string> parts = split(name, false, '_');
				if (parts.size() < 2)
					throw ITLException(string("Invalid chunk writes name: ") + filename);
				if (parts[0] != "chunk")
					throw ITLException(string("Chunk filename does not begin with 'chunk_': ") + filename);
				string part = parts[1];
				parts = split(part, true, '-');
				if (parts.size() < 3)
					throw ITLException(string("Block coordinates do not contain three elements: ") + filename);
				coord_t x = fromString<coord_t>(parts[0]);
				coord_t y = fromString<coord_t>(parts[1]);
				coord_t z = fromString<coord_t>(parts[2]);
				return Vec3c(x, y, z);
			}
		}
	}
}