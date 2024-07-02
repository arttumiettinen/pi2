
#include <iostream>
#include <list>

#include "zarr.h"
#include "json.h"
#include "byteorder.h"
#include "generation.h"
#include "nn5.h"

using namespace std;

namespace itl2
{
	namespace zarr
	{
        //TODO: remove isNativeByteOrder
		bool getInfo(const std::string& path, Vec3c& dimensions, bool& isNativeByteOrder, ImageDataType& dataType, Vec3c& chunkSize, std::list<ZarrCodec>& codecs, int& fillValue, std::string& reason)
		{
			dimensions = Vec3c();
			isNativeByteOrder = true;
			dataType = ImageDataType::Unknown;
			chunkSize = Vec3c();
            fillValue = 0;

			// Check that metadata file exists.
			string metadataFilename = zarr::internals::zarrMetadataFilename(path);
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
				reason = string("Unable to parse zarr metadata: ") + ex.what();
				return false;
			}
            if(!j.contains("zarr_format"))
            {
                reason = "zarr_format is missing in zarr metadata.";
                return false;
            }
            else
            {
                int zarr_format = j["zarr_format"].get<int>();
                if(zarr_format != 3)
                {
                    reason = "This zarr implementation supports only zarr_format 3.";
                    return false;
                }
            }
            if(!j.contains("node_type"))
            {
                reason = "node_type is missing in zarr metadata.";
                return false;
            }
            else
            {
                string node_type = j["node_type"].get<string>();
                if(node_type != "array")
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

			dimensions = Vec3c(1, 1, 1);
			dimensions[0] = dims[0].get<size_t>();
			if (dims.size() >= 2)
				dimensions[1] = dims[1].get<size_t>();
			if (dims.size() >= 3)
				dimensions[2] = dims[2].get<size_t>();

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

            if(!j.contains("chunk_key_encoding"))
            {
                reason = "chunk_key_encoding is missing in zarr metadata.";
                return false;
            }
            else
            {
                if(!j["chunk_key_encoding"].contains("name"))
                {
                    reason = "chunk_key_encoding name is missing in zarr metadata.";
                    return false;
                }
                if(j["chunk_key_encoding"]["name"].get<string>() != "default")
                {
                    reason = "This zarr implementation supports only default chunk_key_encoding.";
                    return false;
                }
                if(j["chunk_key_encoding"].contains("configuration")
                    && !j["chunk_key_encoding"]["configuration"].contains("separator")
                    && j["chunk_key_encoding"]["configuration"]["separator"].get<string>() != "/")
                {
                    reason = "This zarr implementation supports only default chunk_key_encoding with separator \"/\".";
                    return false;
                }
            }

            if(!j.contains("fill_value"))
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
                catch(nlohmann::json::exception ex)
                {
                    reason = string("Unable to parse fill_value in zarr metadata (this implementation only supports integers): ") + ex.what();
                    return false;
                }

			if (!j.contains("codecs"))
			{
				reason = "codecs is missing in zarr metadata.";
                return false;
            }
			{
                int numberArrayBytesCodecs = 0;
				for (auto& codec : j["codecs"])
                {
                    if(!codec.contains("name"))
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
                        zarr::ZarrCodec zarrCodec = fromString<zarr::ZarrCodec>(codec["name"].get<string>());
                        zarrCodec.readConfig(codecConfig);
                        codecs.push_back(zarrCodec);
                        switch(zarrCodec.type)
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

			}

			return true;
		}

		namespace internals
		{
		    //zarr updated
            //TODO: write codecs
			void writeMetadata(const std::string& path, const Vec3c& dimensions, ImageDataType dataType, const Vec3c& chunkSize, int fillValue, std::list<ZarrCodec>& codecs)
			{
				nlohmann::json j =
                        {
                            {"zarr_format", 3},//zarr updated
                            {"node_type", "array"},//zarr updated
                            {"shape", {dimensions[0], dimensions[1], dimensions[2]}},//zarr updated
                            {"data_type", toString(dataType)},//zarr updated
                            {"chunk_grid", {{"name", "regular"}, {"configuration", {{"chunk_shape", {chunkSize[0], chunkSize[1], chunkSize[2]}}}}}},
                            {"chunk_key_encoding", {{"name", "default"}, {"configuration", {{"separator", "/"}}}}},
                            {"fill_value", fillValue},
                            {"codecs", {}}
                        };
                for (auto& codec : codecs)
                {
                    j["codecs"].push_back({codec.toJSON()});
                }
                // TODO: allow other chunk_key_encoding
                // TODO: optional parameters

				string metadataFilename = zarr::internals::zarrMetadataFilename(path);
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

			void check(const Vec3c& chunkSize)
			{
				if (chunkSize.min() <= 0)
					throw ITLException(string("NN5 chunk size must be positive, but it is ") + toString(chunkSize));
			}

			void beginWrite(const Vec3c& imageDimensions, ImageDataType imageDataType, const std::string& path, const Vec3c& chunkSize, int fillValue, std::list<ZarrCodec>& codecs, bool deleteOldData)
			{
				check(chunkSize);

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
							throw ITLException(string("Unable to write a Zarr as the output folder already exists but cannot be verified to be an image: ") + path + " Consider removing the existing dataset manually.");

						// Here the path does not contain Zarr but contains an image of known type.
						// Delete the old image.
						fs::remove_all(path);
					}
					else if (oldDimensions == imageDimensions &&
						oldIsNativeByteOrder == true &&
						oldDataType == imageDataType &&
						oldChunkSize == chunkSize &&
                        oldCodecs == codecs &&
                        oldFillValue == fillValue)
					{
						// The path contains a compatible NN5 dataset.
						// Delete it if we are not continuing a concurrent write.
						if (!fs::exists(internals::concurrentTagFile(path)))
							if (deleteOldData)
								fs::remove_all(path);
					}
					else
					{
						// The path contains an incompatible NN5 dataset. Delete it.
						if (fs::exists(internals::concurrentTagFile(path)) && !deleteOldData)
							throw ITLException(string("The output folder contains an incompatible NN5 dataset that is currently being processed concurrently."));
						fs::remove_all(path);
					}
				}

				fs::create_directories(path);

				// Write metadata
				internals::writeMetadata(path, imageDimensions, imageDataType, chunkSize, fillValue, codecs);
			}
		}

		//size_t countWritersAt(const AABoxc& box, const std::vector<NN5Process>& processes)
		//{
		//	size_t count = 0;
		//	for (const NN5Process& process : processes)
		//	{
		//		if (process.writeBlock.overlaps(box))
		//			count++;
		//	}
		//	return count;
		//}

		/**
		Finds out if a chunk (given its bounding box) is 'safe' or not.
		Safe chunks can be written to without any synchronization or post-processing of the results.
		*/
		bool isChunkSafe(const AABoxc& box, const std::vector<ZarrProcess>& processes)
		{
			vector<size_t> readerIndices;
			vector<size_t> writerIndices;
			for (size_t n = 0; n < processes.size(); n++)
			{
				const auto& process = processes[n];
				if (process.readBlock.overlapsExclusive(box))
					readerIndices.push_back(n);
				if (process.writeBlock.overlapsExclusive(box))
					writerIndices.push_back(n);
			}

			if (writerIndices.size() <= 0)
			{
				// No writers, the chunk is never written to, so it is safe.
				return true;
			}
			if (writerIndices.size() > 1)
			{
				// Multiple writers, the chunk is unsafe as the writers can write simultaneously.
				return false;
			}
			else
			{
				// One writer.
				if (readerIndices.size() <= 0)
				{
					// No readers, one writer, the chunk is safe.
					return true;
				}
				else if (readerIndices.size() > 1)
				{
					// Multiple readers, one writer, the chunk is not safe as it can be read from and written to simultaneously.
					return false;
				}
				else
				{
					// One reader, one writer.
					
					if (readerIndices[0] == writerIndices[0])
					{
						// Reader and writer are the same process.
						// The chunk is safe as the reader/writer process should control its possibly overlapping reads and writes internally.
						return true;
					}
					else
					{
						// Reader and writer are different processes.
						// The chunk is not safe as the reader and the writer might access the chunk simultaneously.
						return false;
					}
				}
			}
		}

		size_t startConcurrentWrite(const Vec3c& imageDimensions, ImageDataType imageDataType, const std::string& path, const Vec3c& chunkSize, int fillValue, std::list<ZarrCodec>& codecs, const std::vector<ZarrProcess>& processes)
		{
			// Find chunks that are
			// * written to by separate processes, or
			// * read from and written to by at least two separate processes,
			// and tag those unsafe by creating writes folder into the chunk folder.

			zarr::internals::beginWrite(imageDimensions, imageDataType, path, chunkSize,  fillValue, codecs, false);

			// Tag the image as concurrently processed
			ofstream out(internals::concurrentTagFile(path), ios_base::out | ios_base::trunc | ios_base::binary);

			size_t unsafeChunkCount = 0;
			internals::forAllChunks(imageDimensions, chunkSize, false, [&](const Vec3c& chunkIndex, const Vec3c& chunkStart)
				{
					string chunkFolder = internals::chunkFolder(path, getDimensionality(imageDimensions), chunkIndex);
					fs::create_directories(chunkFolder);

					string writesFolder = internals::writesFolder(chunkFolder);

					AABoxc chunkBox = AABoxc::fromPosSize(chunkStart, chunkSize);

					if (!isChunkSafe(chunkBox, processes))
					{
						// Mark the chunk as unsafe by creating writes folder.
						fs::create_directories(writesFolder);
						unsafeChunkCount++;
					}
				});
			return unsafeChunkCount;
		}

		bool needsEndConcurrentWrite(const std::string& path, size_t dimensionality, const Vec3c& chunkIndex)
		{
			string chunkFolder = internals::chunkFolder(path, dimensionality, chunkIndex);
			string writesFolder = internals::writesFolder(chunkFolder);
			return fs::exists(writesFolder);
		}

		bool needsEndConcurrentWrite(const std::string& path, const Vec3c& chunkIndex)
		{
			bool isNativeByteOrder;
			Vec3c imageDimensions;
			ImageDataType dataType;
			Vec3c chunkSize;
			int fillValue;
            std::list<ZarrCodec> codecs;
			string reason;
			if (!zarr::getInfo(path, imageDimensions, isNativeByteOrder, dataType, chunkSize, fillValue, codecs, reason))
				throw ITLException(string("Unable to read nn5 dataset: ") + reason);

			return needsEndConcurrentWrite(path, getDimensionality(imageDimensions), chunkIndex);
		}

		vector<Vec3c> getChunksThatNeedEndConcurrentWrite(const std::string& path)
		{
			bool isNativeByteOrder;
			Vec3c imageDimensions;
			ImageDataType dataType;
			Vec3c chunkSize;
            int fillValue;
            std::list<ZarrCodec> codecs;
			string reason;
			if (!zarr::getInfo(path, imageDimensions, isNativeByteOrder, dataType, chunkSize, compression, reason))
				throw ITLException(string("Unable to read nn5 dataset: ") + reason);

			size_t dimensionality = getDimensionality(imageDimensions);
			vector<Vec3c> result;
			internals::forAllChunks(imageDimensions, chunkSize, false, [&](const Vec3c& chunkIndex, const Vec3c& chunkStart)
				{
					if (needsEndConcurrentWrite(path, dimensionality, chunkIndex))
						result.push_back(chunkIndex);
				});
			return result;
		}

		namespace internals
		{
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

			template<typename pixel_t> void readAndAdd(Image<pixel_t>& img, const string& filename, int fillValue, std::list<ZarrCodec>& codecs)
			{
				Vec3c blockPos = parsePosition(filename);

				Image<pixel_t> block;
				readChunkFile(block, filename, fillValue, codecs);

				copyValues(img, block, blockPos);
			}

			template<typename pixel_t> struct CombineChunkWrites
			{
			public:
				static void run(const string& path, const Vec3c& datasetSize, const Vec3c& chunkSize, int fillValue, std::list<ZarrCodec>& codecs, const Vec3c& chunkIndex)
				{
					string chunkFolder = internals::chunkFolder(path, getDimensionality(datasetSize), chunkIndex);
					string writesFolder = internals::writesFolder(chunkFolder);

					if (fs::exists(writesFolder))
					{
						// Find all files in the writes folder
						vector<string> writesFiles = buildFileList(writesFolder + "/");

						if (writesFiles.size() > 0)
						{
							Vec3c realChunkSize = internals::clampedChunkSize(chunkIndex, chunkSize, datasetSize);
							Image<pixel_t> img(realChunkSize);

							// Read old data if it exists.
							// TODO: No need to read if the written blocks overwrite the chunk completely.
							std::vector<string> originalFiles = getFileList(chunkFolder);
							if (originalFiles.size() <= 0)
							{
								// No file => all pixels in the block are zeroes.
								setValue(img, (pixel_t)0);
							}
							else if (originalFiles.size() == 1)
							{
								string filename = chunkFolder + "/" + originalFiles[0];
								readChunkFile(img, filename, fillValue, codecs);
							}
							else
							{
								throw ITLException(string("Multiple image files found in block directory ") + chunkFolder + " while combining chunk writes.");
							}

							// Modify data with the new writes.
							for (const string& file : writesFiles)
							{
								readAndAdd(img, file, fillValue, codecs);
							}

							// Write back to disk.
							switch (compression)
							{
								case NN5Compression::Raw:
								{
									string filename = concatDimensions(chunkFolder + "/chunk", img.dimensions());
									raw::write(img, filename);
									break;
								}
								case NN5Compression::LZ4:
								{
									string filename = chunkFolder + "/chunk.lz4raw";
									lz4::write(img, filename);
									break;
								}
								default:
								{
									throw ITLException(string("Unsupported nn5 compression algorithm: ") + toString(compression));
								}
							}
						}

						// Remove the writes folder in order to mark this chunk processed.
						fs::remove_all(writesFolder);
					}
				}
			};

			void endConcurrentWrite(const std::string& path, const Vec3c& imageDimensions, ImageDataType dataType, const Vec3c& chunkSize, int fillValue, std::list<ZarrCodec>& codecs, const Vec3c& chunkIndex)
			{
				pick<zarr::internals::CombineChunkWrites>(dataType, path, imageDimensions, chunkSize, fillValue, codecs, chunkIndex);
			}
		}

		

		void endConcurrentWrite(const std::string& path, const Vec3c& chunkIndex)
		{
			bool isNativeByteOrder;
			Vec3c imageDimensions;
			ImageDataType dataType;
			Vec3c chunkSize;
            int fillValue;
            std::list<ZarrCodec> codecs;
			string reason;
			if (!zarr::getInfo(path, imageDimensions, isNativeByteOrder, dataType, chunkSize, fillValue, codecs, reason))
				throw ITLException(string("Unable to read nn5 dataset: ") + reason);

			zarr::internals::endConcurrentWrite(path, imageDimensions, dataType, chunkSize, fillValue, codecs, chunkIndex);
		}

		void endConcurrentWrite(const std::string& path, bool showProgressInfo)
		{
			bool isNativeByteOrder;
			Vec3c fileDimensions;
			ImageDataType dataType;
			Vec3c chunkSize;
            int fillValue;
            std::list<ZarrCodec> codecs;
			string reason;
			if (!zarr::getInfo(path, fileDimensions, isNativeByteOrder, dataType, chunkSize, fillValue, codecs, reason))
				throw ITLException(string("Unable to read nn5 dataset: ") + reason);
			size_t dimensionality = getDimensionality(fileDimensions);

			internals::forAllChunks(fileDimensions, chunkSize, showProgressInfo, [&](const Vec3c& chunkIndex, const Vec3c& chunkStart)
				{
					if(needsEndConcurrentWrite(path, dimensionality, chunkIndex))
						zarr::internals::endConcurrentWrite(path, fileDimensions, dataType, chunkSize, fillValue, codecs, chunkIndex);
				});
			
			// Remove concurrent tag file after all blocks are processed such that if exception is thrown during processing,
			// the endConcurrentWrite can continue simply by re-running it.
			fs::remove_all(internals::concurrentTagFile(path));
		}
	}
}