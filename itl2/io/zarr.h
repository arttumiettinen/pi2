#pragma once
#include <list>

#include "filesystem.h"
#include "image.h"
#include "utilities.h"
#include "io/raw.h"
#include "io/itllz4.h"
#include "math/aabox.h"
#include "io/zarrcodecs.h"
#include "generation.h"
#include "math/vec3.h"

namespace itl2
{
	using std::cout, std::endl;
	namespace io
	{
		// Forward declaration of io::getInfo to avoid header loop.
		bool getInfo(const std::string& filename, Vec3c& dimensions, ImageDataType& dataType, std::string& reason);
	}
	namespace zarr
	{
		typedef struct
		{
			Vec3c chunkSize;
			std::list<codecs::ZarrCodec> codecs;
			int fillValue;
			std::string separator;
		} ZarrMetadata;

		//todo: check if this is correctly checking for equality (not only identity)
		bool operator==(const ZarrMetadata& lhs, const ZarrMetadata& rhs)
		{
			return lhs.chunkSize == rhs.chunkSize &&
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
			void writeMetadata(const std::string& path, const Vec3c& shape, ImageDataType dataType, const ZarrMetadata& metadata);

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
			Writes single Zarr chunk file.
			@param chunkIndex Index of the chunk to write. This is used to determine the correct output folder.
			@param chunkSize Chunk size of the dataset.
			@param startInChunkCoords Start position of the block to be written in the coordinates of the chunk.
			@param chunkStart Start position of the block to be written in the coordinates of image targetImg.
			@param writeSize Size of block to be written.
			*/
			template<typename pixel_t>
			void writeSingleChunk(const Image<pixel_t>& img, const std::string& path, const Vec3c& chunkIndex, const Vec3c& datasetSize,
				Vec3c startInChunkCoords, const Vec3c& startInImageCoords, const Vec3c& writeSize, const ZarrMetadata& metadata)
			{
				//todo: read this chunk and only update the requested region by const startInImageCoords,startInChunkCoords, writeSize
				Image<pixel_t> imgChunk(metadata.chunkSize);
				const Vec3c chunkStartPos = chunkIndex.componentwiseMultiply(metadata.chunkSize);
				//write all pixels of chunk in img to imgChunk
				forAllPixels(imgChunk, [&](coord_t x, coord_t y, coord_t z)
				{
				  Vec3c pos = Vec3c(x, y, z) + chunkStartPos;
				  imgChunk(x, y, z) = img(pos);
				});

				string filename = chunkFile(path, getDimensionality(datasetSize), chunkIndex, metadata.separator);

				// Clamp write size to the size of the image.
				Vec3c imageChunkEnd = startInImageCoords + writeSize;
				for (size_t n = 0; n < imageChunkEnd.size(); n++)
				{
					if (imageChunkEnd[n] > img.dimension(n))
						imageChunkEnd[n] = img.dimension(n);
				}
				Vec3c realWriteSize = imageChunkEnd - startInImageCoords;
				Vec3c realChunkSize = clampedChunkSize(chunkIndex, metadata.chunkSize, datasetSize);
				realWriteSize = min(realWriteSize, realChunkSize);

				// Check if we are in an unsafe chunk where writing to the chunk file is prohibited.
				// Chunk is unsafe if its folder contains writes folder.
				//todo: when will writesFolder be created?
				string writesFolder = internals::writesFolder(filename);
				if (fs::exists(writesFolder))
				{
					std::cout << "writeSingleChunk unsafe, writesFolder exist" << std::endl;
					// Unsafe chunk: write to separate writes folder.
					filename = writesFolder + toString(startInChunkCoords.x) + string("-") +
						toString(startInChunkCoords.y) + string("-") + toString(startInChunkCoords.z);
					//print all writeblock parameters in one line
					startInChunkCoords = Vec3c(0, 0, 0);
					realChunkSize = realWriteSize;
				}
				std::list<codecs::ZarrCodec>::const_iterator codec = metadata.codecs.begin();
				for (; codec->type == codecs::Type::ArrayArrayCodec; ++codec)
				{
					codecs::encodeArrayArrayCodec(*codec, imgChunk, metadata.fillValue);
				}
				std::vector<char> buffer;
				codecs::encodeArrayBytesCodec(*codec, imgChunk, buffer);
				++codec;
				for (; codec != metadata.codecs.end(); ++codec)
				{
					codecs::encodeBytesBytesCodec(*codec, buffer);
				}
				writeBytesToFile(buffer, filename);
			}

			inline size_t countChunks(const Vec3c& imageDimensions, const Vec3c& chunkSize)
			{
				Vec3c chunkStart(0, 0, 0);
				size_t chunkCount = 0;
				while (chunkStart.z < imageDimensions.z)
				{
					while (chunkStart.y < imageDimensions.y)
					{
						while (chunkStart.x < imageDimensions.x)
						{
							chunkCount++;
							chunkStart.x += chunkSize.x;
						}
						chunkStart.x = 0;
						chunkStart.y += chunkSize.y;
					}
					chunkStart.x = 0;
					chunkStart.y = 0;
					chunkStart.z += chunkSize.z;
				}
				return chunkCount;
			}
			/**
			Call lambda(chunkIndex, chunkStart) for all chunks in an image of given dimensions and chunk size.
			*/
			template<typename F>
			void forAllChunks(const Vec3c& imageDimensions, const Vec3c& chunkSize, bool showProgressInfo, F&& lambda)
			{

				size_t maxSteps = 0;
				if (showProgressInfo)
					maxSteps = countChunks(imageDimensions, chunkSize);
				ProgressIndicator progress(maxSteps, showProgressInfo);

				Vec3c chunkStart(0, 0, 0);
				Vec3c chunkIndex(0, 0, 0);
				while (chunkStart.z < imageDimensions.z)
				{
					while (chunkStart.y < imageDimensions.y)
					{
						while (chunkStart.x < imageDimensions.x)
						{
							lambda(chunkIndex, chunkStart);
							progress.step();

							chunkIndex.x++;
							chunkStart.x += chunkSize.x;
						}

						chunkIndex.x = 0;
						chunkIndex.y++;
						chunkStart.x = 0;
						chunkStart.y += chunkSize.y;
					}
					chunkIndex.x = 0;
					chunkIndex.y = 0;
					chunkIndex.z++;
					chunkStart.x = 0;
					chunkStart.y = 0;
					chunkStart.z += chunkSize.z;
				}
			}

			template<typename pixel_t>
			void writeChunksInRange(const Image<pixel_t>& img, const std::string& path, const ZarrMetadata& metadata,
				const Vec3c& filePosition, const Vec3c& fileDimensions,
				const Vec3c& imagePosition,
				const Vec3c& blockDimensions,
				bool showProgressInfo)
			{
				const AABoxc fileTargetBlock = AABoxc::fromPosSize(filePosition, blockDimensions);

				forAllChunks(fileDimensions, metadata.chunkSize, showProgressInfo, [&](const Vec3c& chunkIndex, const Vec3c& chunkStart)
				{
				  AABoxc currentChunk = AABoxc::fromPosSize(chunkStart, metadata.chunkSize);
				  if (fileTargetBlock.overlapsExclusive(currentChunk))
				  {
					  AABoxc chunkUpdateRegion = fileTargetBlock.intersection(currentChunk);

					  Vec3c imageUpdateStartPosition = chunkUpdateRegion.position() - filePosition + imagePosition;
					  AABoxc imageUpdateRegion = AABoxc::fromPosSize(imageUpdateStartPosition, chunkUpdateRegion.size());
					  chunkUpdateRegion = chunkUpdateRegion.translate(-chunkStart);
					  if (chunkUpdateRegion.size().min() > 0)
						  writeSingleChunk(img, path, chunkIndex, fileDimensions, chunkUpdateRegion.position(), imageUpdateRegion.position(), chunkUpdateRegion.size(), metadata);
				  }
				});
			}

			/**
			Reads zarr chunk files.
			*/
			template<typename pixel_t>
			void readChunksInRange(Image<pixel_t>& img, const std::string& path,
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
				const AABoxc imageBox = AABoxc::fromMinMax(start, end);
				forAllChunks(datasetShape, metadata.chunkSize, showProgressInfo, [&](const Vec3c& chunkIndex, const Vec3c& chunkStart)
				{
				  AABox<coord_t> currentChunk = AABox<coord_t>::fromPosSize(chunkStart, metadata.chunkSize);
				  const Vec3c chunkStartInTarget = chunkStart - start;
				  const Vec3c readSize = currentChunk.intersection(imageBox).size();

				  if (currentChunk.overlapsExclusive(imageBox))
				  {
					  string filename = chunkFile(path, getDimensionality(datasetShape), chunkIndex, metadata.separator);
					  if (fs::is_directory(filename))
					  {
						  throw ITLException(filename + string(" is a directory, but it should be a file."));
					  }
					  if (fs::is_regular_file(filename))
					  {
						  std::vector<char> buffer = readBytesOfFile(filename);
						  //todo: does this work with const codecs
						  std::list<codecs::ZarrCodec>::const_reverse_iterator codec = metadata.codecs.rbegin();
						  for (; codec->type == codecs::Type::BytesBytesCodec; ++codec)
						  {
							  codecs::decodeBytesBytesCodec(*codec, buffer);
						  }
						  Image<pixel_t> imgChunk(transposedChunkShape);
						  codecs::decodeArrayBytesCodec(*codec, imgChunk, buffer);
						  ++codec;
						  for (; codec != metadata.codecs.rend(); ++codec)
						  {
							  codecs::decodeArrayArrayCodec(*codec, imgChunk, metadata.fillValue);
						  }

						  //write all pixels of chunk back to img
						  forAllPixels(imgChunk, [&](coord_t x, coord_t y, coord_t z)
						  {
							Vec3c pos = Vec3c(x, y, z) + chunkStartInTarget;
							img(pos) = imgChunk(x, y, z);
						  });
					  }
					  else
					  {
						  cout << "no file: " << filename << endl;
						  draw<pixel_t>(img, AABoxc::fromPosSize(chunkStartInTarget, readSize), (pixel_t)metadata.fillValue);
					  }
				  }
				});
			}

			/**
			Checks that provided zarr information is correct and writes metadata.
			@param deleteOldData Contents of the dataset are usually deleted when a write process is started. Set this to false to delete only if the path contains an incompatible dataset, but keep the dataset if it seems to be the same than the current image.
			*/
			void handleExisting(const Vec3c& imageDimensions,
				ImageDataType imageDataType,
				const std::string& path,
				const ZarrMetadata& metadata,
				bool deleteOldData);

		}

		bool getInfo(const std::string& path,
			Vec3c& shape,
			ImageDataType& dataType,
			ZarrMetadata& metadata,
			std::string& reason);

		inline bool getInfo(const std::string& path, Vec3c& shape, ImageDataType& dataType, std::string& reason)
		{
			ZarrMetadata dummyMetadata;
			return getInfo(path, shape, dataType, dummyMetadata, reason);
		}

		inline const Vec3c DEFAULT_CHUNK_SIZE = Vec3c(1536, 1536, 1536);
		inline const string DEFAULT_SEPARATOR = "/";
		inline const int DEFAULT_FILLVALUE = 42;
		inline const std::list<codecs::ZarrCodec> DEFAULT_CODECS = { codecs::ZarrCodec(codecs::Name::Bytes) };
		inline const nlohmann::json DEFAULT_CODECS_JSON = { codecs::ZarrCodec(codecs::Name::Bytes).toJSON() };

		/**
		Writes a block of an image to the specified location in an .zarr dataset.
		The output dataset is not truncated if it exists.
		If the output file does not exist, it is created.
		@param targetImg Image to write.
		@param filename Name of file to write.
		@param filePosition Position in the file to write to.
		@param fileDimension Total dimensions of the entire output file.
		@param imagePosition Position in the image where the block to be written starts.
		@param blockDimensions Dimensions of the block of the source image to write.
		*/
		template<typename pixel_t>
		void writeBlock(const Image<pixel_t>& img, const std::string& path, const Vec3c& filePosition, const Vec3c& fileDimensions,
			const Vec3c& imagePosition,
			const Vec3c& blockDimensions,
			const Vec3c& chunkSize = DEFAULT_CHUNK_SIZE,
			std::list<codecs::ZarrCodec> codecs = DEFAULT_CODECS,
			int fillValue = DEFAULT_FILLVALUE,
			const string& separator = DEFAULT_SEPARATOR,
			bool showProgressInfo = false)
		{
			Vec3c clampedChunkSize = chunkSize;
			clamp(clampedChunkSize, Vec3c(1, 1, 1), img.dimensions());
			fs::create_directories(path);
			ZarrMetadata metadata = { clampedChunkSize, codecs, fillValue, separator };
			internals::writeMetadata(path, fileDimensions, img.dataType(), metadata);
			internals::writeChunksInRange(img, path, metadata, filePosition, fileDimensions, imagePosition, blockDimensions, showProgressInfo);
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
			const std::string& separator = DEFAULT_SEPARATOR,
			std::list<codecs::ZarrCodec> codecs = DEFAULT_CODECS,
			int fillValue = DEFAULT_FILLVALUE,
			bool showProgressInfo = false)
		{
			bool deleteOldData = false;
			Vec3c dimensions = img.dimensions();
			Vec3c clampedChunkSize = chunkSize;
			clamp(clampedChunkSize, Vec3c(1, 1, 1), dimensions);
			ZarrMetadata metadata = { clampedChunkSize, codecs, fillValue, separator };
			internals::handleExisting(dimensions, img.dataType(), path, metadata, deleteOldData);
			writeBlock(img, path, Vec3c(0, 0, 0), dimensions, Vec3c(0, 0, 0), dimensions, clampedChunkSize, codecs, fillValue, separator, showProgressInfo);
		}

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
			bool isNativeByteOrder;
			Vec3c datasetShape;
			ImageDataType dataType;
			ZarrMetadata metadata;
			string reason;
			std::cout << "readBlock" << std::endl;
			if (!itl2::zarr::getInfo(path, datasetShape, dataType, metadata, reason))
				throw ITLException(string("Unable to read zarr dataset: ") + reason);

			if (dataType != img.dataType())
				throw ITLException(string("Expected data type is ") + toString(img.dataType()) + " but the zarr dataset contains data of type " + toString(dataType) + ".");

			if (fileStart.x < 0 || fileStart.y < 0 || fileStart.z < 0 || fileStart.x >= datasetShape.x || fileStart.y >= datasetShape.y || fileStart.z >= datasetShape.z)
				throw ITLException("Out of bounds start position in readBlock.");

			Vec3c cStart = fileStart;
			clamp(cStart, Vec3c(0, 0, 0), datasetShape);
			Vec3c cEnd = fileStart + img.dimensions();
			clamp(cEnd, Vec3c(0, 0, 0), datasetShape);

			internals::readChunksInRange(img, path, metadata, datasetShape, cStart, cEnd, showProgressInfo);

			if (!isNativeByteOrder)
				swapByteOrder(img);
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
			if (!itl2::zarr::getInfo(path, dimensions, dataType, metadata, reason))
				throw ITLException(string("Unable to read zarr dataset: ") + reason);
			img.ensureSize(dimensions);

			readBlock(img, path, Vec3c(0, 0, 0), showProgressInfo);
		}
		namespace tests
		{
			void read();
			void write();
			void transpose();
			void blosc();
			void writeBlock();
			void readBlock();
		}
	}
}