#pragma once
#include <list>
#include <variant>
#include <iostream>
#include <string>

#include "filesystem.h"
#include "image.h"
#include "utilities.h"
#include "io/raw.h"
#include "io/itllz4.h"
#include "math/aabox.h"
#include "io/zarrcodecs.h"
#include "generation.h"
#include "math/vec3.h"
#include "io/distributedimageprocess.h"
#include "json.h"
#include "datatypes.h"


namespace itl2
{
	using std::cout, std::endl, std::ifstream;
	namespace io
	{
		// Forward declaration of io::getInfo to avoid header loop.
		bool getInfo(const std::string& filename, Vec3c& dimensions, ImageDataType& dataType, std::string& reason);
	}
	namespace zarr
	{
		typedef std::variant<uint8_t, uint16_t, uint32_t, uint64_t, int8_t, int16_t, int32_t, int64_t, float32_t, complex32_t> ImageValue;

#define JSON_UNSIGNED_t uint16_t
#define JSON_SIGNED_t int64_t
#define JSON_FLOAT_t double
		typedef std::variant<JSON_UNSIGNED_t, JSON_SIGNED_t, JSON_FLOAT_t> JsonNumberType;

		struct ZarrMetadata
		{
			ImageDataType dataType;
			Vec3c chunkSize;
			codecs::Pipeline codecs;
			fillValue_t fillValue;
			std::string separator;
		};

		inline bool operator==(const ZarrMetadata& lhs, const ZarrMetadata& rhs)
		{
			return lhs.chunkSize == rhs.chunkSize &&
				lhs.dataType == rhs.dataType &&
				lhs.codecs == rhs.codecs &&
				lhs.fillValue == rhs.fillValue &&
				lhs.separator == rhs.separator;
		}

		namespace internals
		{

			inline string zarrMetadataFilename(const string& path)
			{
				return path + "/zarr.json";
			}

			inline string concurrentTagFile(const string& path)
			{
				return path + "/concurrent";
			}



			/**
			Writes Zarr metadata file.
			*/
			inline void writeMetadata(const std::string& path, const Vec3c& shape, const ZarrMetadata& metadata)
			{
				nlohmann::json j =
					{
						{ "zarr_format", 3 },
						{ "node_type", "array" },
						{ "shape", { shape[0], shape[1], shape[2] }},
						{ "data_type", toString(metadata.dataType) },
						{ "chunk_grid", {{ "name", "regular" }, { "configuration", {{ "chunk_shape", { metadata.chunkSize[0], metadata.chunkSize[1], metadata.chunkSize[2] }}}}}},
						{ "chunk_key_encoding", {{ "name", "default" }, { "configuration", {{ "separator", metadata.separator }}}}},
						{ "codecs", {}}
					};

				j["fill_value"] = metadata.fillValue;

				for (auto& codec : metadata.codecs)
				{
					j["codecs"].push_back(codec.toJSON());
				}
				// TODO: allow other chunk_key_encoding
				// TODO: optional parameters
				string metadataFilename = zarr::internals::zarrMetadataFilename(path);
				std::ofstream of(metadataFilename, std::ios_base::trunc | std::ios_base::out);
				of << std::setw(4) << j << endl;
			}

			inline std::string chunkFile(const std::string& path, size_t dimensionality, const Vec3c& chunkIndex, const string& separator)
			{
				string filename = path + string("/c");
				for (size_t n = 0; n < dimensionality; n++)
					filename += separator + toString(chunkIndex[n]);
				return filename;
			}

			inline std::string writesFolder(const std::string& chunkFile)
			{
				// the file is a flag
				// if it exists, the chunk is unsafe to write to because it is currently being written to.
				return chunkFile + "_writes";
			}

			inline Vec3c clampedChunkSize(const Vec3c& chunkIndex, const Vec3c& chunkSize, const Vec3c& datasetSize)
			{
				Vec3c datasetChunkStart = chunkIndex.componentwiseMultiply(chunkSize);
				Vec3c datasetChunkEnd = datasetChunkStart + chunkSize;
				for (size_t n = 0; n < datasetChunkEnd.size(); n++)
				{
					if (datasetChunkEnd[n] > datasetSize[n])
						datasetChunkEnd[n] = datasetSize[n];
				}
				Vec3c realChunkSize = datasetChunkEnd - datasetChunkStart;

				return realChunkSize;
			}

			std::vector<char> readBytesOfFile(std::string& filename);
			void writeBytesToFile(std::vector<char>& buffer, const std::string& filename, size_t startInFilePos = 0);


			/**
			Reads zarr chunk files.
			*/
			template<typename pixel_t>
			void readChunksInRange(Image<pixel_t>& outputBlock, const std::string& path,
				const ZarrMetadata& metadata,
				const Vec3c& datasetShape,
				const Vec3c& start, const Vec3c& end,
				bool showProgressInfo)
			{
				Vec3c transposedChunkShape = metadata.chunkSize;
				for (auto codec : metadata.codecs)
				{
					if (codec.name == codecs::Name::Transpose)
					{
						Vec3c order;
						codec.getTransposeConfiguration(order);
						transposedChunkShape = transposedChunkShape.transposed(order);
					}
				}
				const AABoxc selectedBlock = AABoxc::fromMinMax(start, end);
				forAllChunks(datasetShape, metadata.chunkSize, showProgressInfo, [&](const Vec3c& chunkIndex, const Vec3c& chunkStart)
				{
				  AABox<coord_t> currentChunk = AABox<coord_t>::fromPosSize(chunkStart, metadata.chunkSize);
				  if (selectedBlock.overlapsExclusive(currentChunk))
				  {
					  AABox<coord_t> readRegion = selectedBlock.intersection(currentChunk);
					  string filename = chunkFile(path, getDimensionality(datasetShape), chunkIndex, metadata.separator);
					  if (fs::is_directory(filename))
					  {
						  throw ITLException(filename + string(" is a directory, but it should be a file."));
					  }
					  if (fs::is_regular_file(filename))
					  {
						  std::vector<char> buffer = readBytesOfFile(filename);
						  Image<pixel_t> imgChunk(transposedChunkShape);
						  decodePipeline(metadata.codecs, imgChunk, buffer, metadata.fillValue);

						  //write all pixels of chunk back to img
						  forAllInBox(readRegion, [&](coord_t x, coord_t y, coord_t z)
						  {
							Vec3c pos = Vec3c(x, y, z);
							outputBlock(pos - selectedBlock.position()) = imgChunk(pos - currentChunk.position());
						  });
					  }
					  else
					  {
						  draw<pixel_t>(outputBlock, readRegion, static_cast<pixel_t>(metadata.fillValue));
					  }
				  }
				});
			}

			/*
			 * getInfo will return false if the provided data is either no zarr data or has a configuration not supported by this implementation
			 * getInfo will throw ITLException if zarr data is provided but does not match zarr specifications
			 */
			inline bool getInfo(const std::string& path,
				Vec3c& shape,
				ZarrMetadata& metadata,
				std::string& reason)
			{
				shape = Vec3c();
				metadata.dataType = ImageDataType::Unknown;

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

				if (!j.contains("data_type"))
				{
					throw ITLException("data_type is missing in zarr metadata.");
				}
				try
				{
					string dataTypeString = j["data_type"].get<string>();
					metadata.dataType = fromString<ImageDataType>(dataTypeString);
					if (metadata.dataType == ImageDataType::Unknown)
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

					metadata.chunkSize = Vec3c(1, 1, 1);
					metadata.chunkSize[0] = chunkDims[0].get<size_t>();
					if (chunkDims.size() >= 2)
						metadata.chunkSize[1] = chunkDims[1].get<size_t>();
					if (chunkDims.size() >= 3)
						metadata.chunkSize[2] = chunkDims[2].get<size_t>();
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
					}
					if (j["chunk_key_encoding"]["name"].get<string>() != "default")
					{
						reason = "This zarr implementation supports only default chunk_key_encoding.";
						return false;
					}
					if (j["chunk_key_encoding"].contains("configuration")
						&& !j["chunk_key_encoding"]["configuration"].contains("separator"))
					{
						throw ITLException("chunk_key_encoding configuration separator is missing in zarr metadata");
					}
					else
					{
						metadata.separator = j["chunk_key_encoding"]["configuration"]["separator"].get<string>();
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
						metadata.fillValue = static_cast<fillValue_t>(j["fill_value"]);
					}
					catch (nlohmann::json::exception ex)
					{
						reason =
							string("Unable to parse fill_value in zarr metadata (this implementation only supports floats and integers): ") +
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
					if (!codecs::fromJSON(metadata.codecs, j["codecs"], reason))
					{
						return false;
					}
				}
				return true;
			}

			/**
			Checks that provided zarr information is correct and writes metadata.
			@param deleteOldData Contents of the dataset are usually deleted when a write process is started. Set this to false to delete only if the path contains an incompatible dataset, but keep the dataset if it seems to be the same than the current image.
			*/
			inline void handleExisting(const Vec3c& imageDimensions,
				const std::string& path,
				const ZarrMetadata& metadata,
				bool deleteOldData)
			{
				// Delete old dataset if it exists.
				if (fs::exists(path))
				{
					Vec3c oldDimensions;
					ZarrMetadata oldMetadata;
					string dummyReason;
					string zarrReason;

					bool isZarrImage = false;
					try
					{
						isZarrImage = getInfo(path, oldDimensions, oldMetadata, zarrReason);
					}
					catch (ITLException& e)
					{
						zarrReason = e.message();
					}
					if (!isZarrImage)
					{
						// The path does not contain a Zarr dataset.
						// If it is no known image, do not delete it.
						bool isKnownImage = false;
						try
						{
							isKnownImage = io::getInfo(path, oldDimensions, oldMetadata.dataType, dummyReason);
						}
						catch (ITLException& e)
						{}
						if (!isKnownImage)
							throw ITLException(string("Unable to write a Zarr as the output folder already exists but cannot be verified to be an image: ") + path
								+ " Consider removing the existing dataset manually. Reason why could not read as zarr: " + zarrReason);

						// Here the path does not contain Zarr but contains an image of known type.
						// Delete the old image.
						fs::remove_all(path);
					}
					else if (oldDimensions == imageDimensions &&
						oldMetadata == metadata)
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
		}

		inline bool getInfo(const std::string& path, Vec3c& shape, ImageDataType& dataType, std::string& reason)
		{
			ZarrMetadata dummyMetadata;
			return zarr::internals::getInfo(path, shape, dummyMetadata, reason);
		}

		inline const Vec3c DEFAULT_CHUNK_SIZE = Vec3c(256, 256, 256);
		inline const string DEFAULT_SEPARATOR = "/";
		inline const fillValue_t DEFAULT_FILLVALUE = 0;
		inline const ImageDataType DEFAULT_DATATYPE = ImageDataType::Int32;
		inline const codecs::Pipeline DEFAULT_CODECS = { codecs::ZarrCodec(codecs::Name::Bytes) };
		inline const nlohmann::json DEFAULT_CODECS_JSON = { codecs::ZarrCodec(codecs::Name::Bytes).toJSON() };
		inline const ZarrMetadata DEFAULT_METADATA = {DEFAULT_DATATYPE, DEFAULT_CHUNK_SIZE, DEFAULT_CODECS, DEFAULT_FILLVALUE, DEFAULT_SEPARATOR};



		/**
		Reads a part of a .zarr dataset to the given image.
		NOTE: Does not support out of bounds start position.
		@param targetImg Image where the data is placed. The size of the image defines the size of the block that is read.
		@param filename The name of the dataset to read.
		@param fileStart Start location of the read in the file. The size of the image defines the size of the block that is read.
		*/
		template<typename pixel_t>
		void readBlock(Image<pixel_t>& img, const std::string& path, const Vec3c& fileStart, bool showProgressInfo = false)
		{
			Vec3c datasetShape;
			ZarrMetadata metadata;
			string reason;
			if (!itl2::zarr::internals::getInfo(path, datasetShape, metadata, reason))
				throw ITLException(string("Unable to read zarr dataset: ") + reason);

			if (metadata.dataType != img.dataType())
				throw ITLException(string("Expected data type is ") + toString(img.dataType()) + " but the zarr dataset contains data of type " + toString(metadata.dataType) + ".");

			if (fileStart.x < 0 || fileStart.y < 0 || fileStart.z < 0 || fileStart.x >= datasetShape.x || fileStart.y >= datasetShape.y || fileStart.z >= datasetShape.z)
				throw ITLException("Out of bounds start position in readBlock.");

			Vec3c cStart = fileStart;
			clamp(cStart, Vec3c(0, 0, 0), datasetShape);
			Vec3c cEnd = fileStart + img.dimensions();
			clamp(cEnd, Vec3c(0, 0, 0), datasetShape);

			internals::readChunksInRange(img, path, metadata, datasetShape, cStart, cEnd, showProgressInfo);
		}
		/**
		Reads a zarr dataset file to the given image.
		@param targetImg Image where the data is placed. The size of the image will be set based on the dataset contents.
		@param path Path to the root of the nn5 dataset.
		*/
		template<typename pixel_t>
		void read(Image<pixel_t>& img, const std::string& path, bool showProgressInfo = false)
		{
			Vec3c dimensions;
			ImageDataType dataType;
			ZarrMetadata metadata;
			string reason;
			if (!itl2::zarr::internals::getInfo(path, dimensions, metadata, reason))
				throw ITLException(string("Unable to read zarr dataset: ") + reason);
			img.ensureSize(dimensions);

			readBlock(img, path, Vec3c(0, 0, 0), showProgressInfo);
		}

		namespace internals
		{
			/**
				Writes single Zarr chunk file.
				@param chunkIndex Index of the chunk to write. This is used to determine the correct output folder.
				@param chunkSize Chunk size of the dataset.
				@param startInChunkCoords Start position of the block to be written in the coordinates of the chunk.
				@param chunkStart Start position of the block to be written in the coordinates of image targetImg.
				@param writeSize Size of block to be written.
				*/
			template<typename pixel_t>
			void writeSingleChunk(const Image<pixel_t>& img, const std::string& path, const Vec3c& chunkStart, const Vec3c& chunkIndex,
				AABoxc updateRegion, AABoxc currentChunk, const ZarrMetadata& metadata)
			{
				Image<pixel_t> imgChunk(metadata.chunkSize, metadata.fillValue);
				readBlock(imgChunk, path, chunkStart);
				string filename = chunkFile(path, getDimensionality(img.dimensions()), chunkIndex, metadata.separator);
				//write all pixels of chunk in img to imgChunk
				forAllInBox(updateRegion, [&](coord_t x, coord_t y, coord_t z)
				{
				  Vec3c pos = Vec3c(x, y, z);
				  imgChunk(pos - currentChunk.position()) = img(pos);
				});


				if (allEquals(imgChunk, static_cast<pixel_t>(metadata.fillValue)))
				{
					//chunk is empty so we can delete the chunk file
					//TODO: what if writesFolder exists?
					if (fs::exists(filename))
					{
						fs::remove(filename);
					}
				}
				else
				{

					// Check if we are in an unsafe chunk where writing to the chunk file is prohibited.
					// Chunk is unsafe if its folder contains writes folder.
					//todo: when will writesFolder be created?
					string writesFolder = internals::writesFolder(filename);
					if (fs::exists(writesFolder))
					{
						std::cout << "writeSingleChunk unsafe, writesFolder exist" << std::endl;
						// Unsafe chunk: write to separate writes folder.
						filename = writesFolder + toString(updateRegion);
					}
					std::vector<char> buffer;
					assert(imgChunk.dimensions() == currentChunk.size());
					encodePipeline(metadata.codecs, imgChunk, buffer, metadata.fillValue);
					writeBytesToFile(buffer, filename);
				}
			}

			template<typename pixel_t>
			void writeChunksInRange(const Image<pixel_t>& img, const std::string& path, const ZarrMetadata& metadata,
				const Vec3c& blockPosition,
				const Vec3c& blockDimensions,
				bool showProgressInfo)
			{
				const AABoxc selectedBlock = AABoxc::fromPosSize(blockPosition, blockDimensions);
				forAllChunks(img.dimensions(), metadata.chunkSize, showProgressInfo, [&](const Vec3c& chunkIndex, const Vec3c& chunkStart)
				{
				  AABoxc currentChunk = AABoxc::fromPosSize(chunkStart, metadata.chunkSize);
				  if (selectedBlock.overlapsExclusive(currentChunk))
				  {
					  AABoxc updateRegion = selectedBlock.intersection(currentChunk).intersection(img.bounds());
					  writeSingleChunk(img, path, chunkStart, chunkIndex, updateRegion, currentChunk, metadata);
				  }
				});
			}
		}
		/**
		Writes a block of an image to the specified location in an .zarr dataset.
		The output dataset is not truncated if it exists.
		If the output file does not exist, it is created and values outside the block
		are filled with fillValue
		@param targetImg Image to write.
		@param filename Name of file to write.
		@param filePosition Position in the file to write to.
		@param fileDimension Total dimensions of the entire output file.
		@param imagePosition Position in the image where the block to be written starts.
		@param blockDimensions Dimensions of the block of the source image to write.
		*/
		template<typename pixel_t>
		void writeBlock(const Image<pixel_t>& img, const std::string& path, const Vec3c& blockPosition,
			const Vec3c& blockDimensions,
			const Vec3c& chunkSize = DEFAULT_CHUNK_SIZE,
			codecs::Pipeline codecs = DEFAULT_CODECS,
			fillValue_t fillValue = DEFAULT_FILLVALUE,
			const string& separator = DEFAULT_SEPARATOR,
			bool showProgressInfo = false)
		{
			Vec3c clampedChunkSize = chunkSize;
			if (chunkSize.min()<=0)
				throw ITLException("Illegal chunk size: " + toString(chunkSize));
			fs::create_directories(path);
			ZarrMetadata metadata = {img.dataType(), clampedChunkSize, codecs, fillValue, separator};
			internals::writeMetadata(path, img.dimensions(), metadata);
			internals::writeChunksInRange(img, path, metadata, blockPosition, blockDimensions, showProgressInfo);
		}

		/**
		Write an image to a zarr dataset.
		@param targetImg Image to write.
		@param path Name of the top directory of the nn5 dataset.
		*/
		template<typename pixel_t>
		void write(const Image<pixel_t>& img,
			const std::string& path,
			const Vec3c& chunkSize = DEFAULT_CHUNK_SIZE,
			codecs::Pipeline codecs = DEFAULT_CODECS,
			const std::string& separator = DEFAULT_SEPARATOR,
			fillValue_t fillValue = DEFAULT_FILLVALUE,
			bool showProgressInfo = false)
		{
			bool deleteOldData = false;
			Vec3c dimensions = img.dimensions();
			if (chunkSize.min()<=0)
				throw ITLException("Illegal chunk size: " + toString(chunkSize));
			ZarrMetadata metadata = { img.dataType(), chunkSize, codecs, fillValue, separator };
			internals::handleExisting(dimensions, path, metadata, deleteOldData);
			writeBlock(img, path, Vec3c(0, 0, 0), dimensions, chunkSize, codecs, fillValue, separator, showProgressInfo);
		}

		/**
		Enables concurrent access from multiple processes for an existing or a new NN5 dataset.
		This function should be called before the processes are started.
		@param targetImg Image that is to be saved into the NN5 dataset by the processes.
		@param path Path to the NN5 dataset.
		@param chunkSize Chunk size for the NN5 dataset.
		@param compression Compression method to be used.
		@param processes A list of DistributedImageProcess objects that define the block that where each process will have read and write access. The blocks may overlap.
		*/
		template<typename pixel_t> void startConcurrentWrite(const Image<pixel_t>& img, const std::string& path, const Vec3c& chunkSize, const std::vector<io::DistributedImageProcess>& processes)
		{
			startConcurrentWrite(img.dimensions(), img.dataType(), path, chunkSize, processes);
		}

		/**
		Enables concurrent access from multiple processes for an existing or a new NN5 dataset.
		This function should be called before the processes are started.
		@param imageDimensions Dimensions of the image to be saved into the NN5 dataset.
		@param imageDataType Data type of the image.
		@param path Path to the NN5 dataset.
		@param chunkSize Chunk size for the NN5 dataset.
		@param compression Compression method to be used.
		@param processes A list of DistributedImageProcess objects that define the block that where each process will have read and write access. The blocks may overlap.
		@return Number of chunks that require special processing in endConcurrentWrite.
		*/
		size_t startConcurrentWrite(const Vec3c& imageDimensions, ImageDataType imageDataType, const std::string& path, const Vec3c& chunkSize, const std::vector<io::DistributedImageProcess>& processes);


		/**
		Finalizes concurrent access from multiple processes.
		This function should be called after all the processes have finished accessing the dataset.
		This function calls endConcurrentWrite(path, chunkIndex) for all chunks in the dataset for which needsEndConcurrentWrite return true,
		and removes concurrent tag file from the dataset root folder.
		@param path Path to the NN5 dataset.
		*/
		void endConcurrentWrite(const std::string& path, bool showProgressInfo = false);

		/**
		Used to test if a block in the given NN5 dataset requires calling endConcurrentWrite after concurrent access by multiple processes.
		@param path Path to the NN5 dataset.
		@param chunkIndex Index of the chunk to finalize.
		*/
		bool needsEndConcurrentWrite(const std::string& path, const Vec3c& chunkIndex);

		/**
		Get a list of chunks for which needsEndConcurrentWrite must be called.
		*/
		std::vector<Vec3c> getChunksThatNeedEndConcurrentWrite(const std::string& path);

		/**
		Finalizes concurrent access from multiple processes for a single chunk in an NN5 dataset.
		This function can be used instead of the other endConcurrentWrite overload, but it must be called for
		all chunks in the dataset.
		The function can be called concurrently for different chunk indices.
		Processing might involve doing nothing, or reading and re-writing the chunk.
		@param path Path to the NN5 dataset.
		@param chunkIndex Index of the chunk to finalize.
		*/
		void endConcurrentWrite(const std::string& path, const Vec3c& chunkIndex);

		namespace tests
		{
			void read();
			void write();
			void transpose();
			void blosc();
			void writeBlock();
			void readBlock();
			void zarrMetadataEquals();
			void separator();
			void sharding();
			void emptyChunks();

		}
	}
}