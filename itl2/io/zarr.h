#pragma once

#include "filesystem.h"
#include "image.h"
#include "utilities.h"
#include "io/raw.h"
#include "io/itllz4.h"
#include "math/aabox.h"
#include "io/nn5compression.h"
#include "generation.h"

namespace itl2
{

	namespace zarr
	{
	    using nn5::NN5Compression;
		namespace internals
		{
			/**
			Constructs NN5 dataset metadata filename.
			*/
			inline string nn5MetadataFilename(const string& path)
			{
				return path + "/metadata.json";
			}

			inline string concurrentTagFile(const string& path)
			{
				return path + "/concurrent";
			}

			/**
			Writes NN5 metadata file.
			*/
			void writeMetadata(const std::string& path, const Vec3c& dimensions, ImageDataType dataType, const Vec3c& chunkSize, NN5Compression compression);

			/**
			Retrieve a list of all files in the given directory.
			*/
			std::vector<string> getFileList(const std::string& dir);

			inline std::string chunkFolder(const std::string& path, size_t dimensionality, const Vec3c& chunkIndex)
			{
				// 0-dimensional images get a special treatment.
				if (dimensionality <= 0)
					return path + string("/0");

				string filename = path;
				for (size_t n = 0; n < dimensionality; n++)
					filename += string("/") + toString(chunkIndex[n]);
				return filename;
			}

			inline std::string writesFolder(const std::string& chunkFolder)
			{
				return chunkFolder + "/writes";
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

			/**
			Writes single NN5 chunk file.
			@param chunkIndex Index of the chunk to write. This is used to determine the correct output folder.
			@param chunkSize Chunk size of the dataset.
			@param startInChunkCoords Start position of the block to be written in the coordinates of the chunk.
			@param chunkStart Start position of the block to be written in the coordinates of image targetImg.
			@param writeSize Size of block to be written.
			*/
			template<typename pixel_t> void writeSingleChunk(const Image<pixel_t>& img, const std::string& path, const Vec3c& chunkIndex, const Vec3c& chunkSize, const Vec3c& datasetSize,
				const Vec3c& startInChunkCoords, const Vec3c& startInImageCoords, const Vec3c& writeSize, NN5Compression compression)
			{
				// Build path to chunk folder.
				string filename = chunkFolder(path, getDimensionality(datasetSize), chunkIndex);

				// Check if we are in an unsafe chunk where writing to the chunk file is prohibited.
				// Chunk is unsafe if its folder contains writes folder.
				string writesFolder = internals::writesFolder(filename);
				bool unsafe;
				if (fs::exists(writesFolder))
				{
					// Unsafe chunk, write to separate writes folder.
					unsafe = true;
					filename = writesFolder + string("/chunk_") + toString(startInChunkCoords.x) + string("-") + toString(startInChunkCoords.y) + string("-") + toString(startInChunkCoords.z);
				}
				else
				{
					// Safe chunk, write directly to the chunk file.
					unsafe = false;
					filename += "/chunk";
				}

				// Clamp write size to the size of the image.
				Vec3c imageChunkEnd = startInImageCoords + writeSize;
				for (size_t n = 0; n < imageChunkEnd.size(); n++)
				{
					if (imageChunkEnd[n] > img.dimension(n))
						imageChunkEnd[n] = img.dimension(n);
				}
				Vec3c realWriteSize = imageChunkEnd - startInImageCoords;

				// Determine if this is the last chunk, and reduce chunk size accordingly so that image size stays correct.
				//Vec3c datasetChunkStart = chunkIndex.componentwiseMultiply(chunkSize);
				//Vec3c datasetChunkEnd = datasetChunkStart + chunkSize;
				//for (size_t n = 0; n < datasetChunkEnd.size(); n++)
				//{
				//	if (datasetChunkEnd[n] > datasetSize[n])
				//		datasetChunkEnd[n] = datasetSize[n];
				//}
				//Vec3c realChunkSize = datasetChunkEnd - datasetChunkStart;
				Vec3c realChunkSize = clampedChunkSize(chunkIndex, chunkSize, datasetSize);

				realWriteSize = min(realWriteSize, realChunkSize);

				if (!unsafe)
				{
					switch (compression)
					{
					case NN5Compression::Raw:
					{
						filename = concatDimensions(filename, realWriteSize);
						raw::writeBlock(img, filename, startInChunkCoords, realChunkSize, startInImageCoords, realWriteSize, false);
						break;
					}
					case NN5Compression::LZ4:
					{
						filename += ".lz4raw";
						lz4::writeBlock(img, filename, startInChunkCoords, realChunkSize, startInImageCoords, realWriteSize);
						break;
					}
					default:
					{
						throw ITLException(string("Unsupported nn5 compression algorithm: ") + toString(compression));
					}
					}
				}
				else
				{
					switch (compression)
					{
					case NN5Compression::Raw:
					{
						filename = concatDimensions(filename, realWriteSize);
						raw::writeBlock(img, filename, Vec3c(0, 0, 0), realWriteSize, startInImageCoords, realWriteSize, false);
						break;
					}
					case NN5Compression::LZ4:
					{
						filename += ".lz4raw";
						lz4::writeBlock(img, filename, Vec3c(0, 0, 0), realWriteSize, startInImageCoords, realWriteSize);
						break;
					}
					default:
					{
						throw ITLException(string("Unsupported nn5 compression algorithm: ") + toString(compression));
					}
					}
				}
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
				if(showProgressInfo)
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

			/**
			Writes NN5 chunk files.
			*/
			template<typename pixel_t> void writeChunks(const Image<pixel_t>& img, const std::string& path, const Vec3c& chunkSize, NN5Compression compression, const Vec3c& datasetSize, bool showProgressInfo)
			{
				forAllChunks(img.dimensions(), chunkSize, showProgressInfo, [&](const Vec3c& chunkIndex, const Vec3c& chunkStart)
				{
					writeSingleChunk(img, path, chunkIndex, chunkSize, datasetSize, Vec3c(0, 0, 0), chunkStart, chunkSize, compression);
				});
			}

			template<typename pixel_t> void writeChunksInRange(Image<pixel_t>& img, const std::string& path, const Vec3c& chunkSize, NN5Compression compression,
				const Vec3c& filePosition, const Vec3c& fileDimensions,
				const Vec3c& imagePosition,
				const Vec3c& blockDimensions,
				bool showProgressInfo)
			{
				// Writes block of image defined by (imagePosition, blockDimensions) to file (defined by path),
				// to location defined by filePosition.

				AABoxc imageBlock = AABoxc::fromPosSize(imagePosition, blockDimensions);

				AABoxc fileTargetBlock = AABoxc::fromPosSize(filePosition, blockDimensions);

				forAllChunks(fileDimensions, chunkSize, showProgressInfo, [&](const Vec3c& chunkIndex, const Vec3c& chunkStart)
					{
						// This is done for all chunks in the output file.
						// We will need to update the chunk if the file block to be written (fileTargetBlock)
						// overlaps the current chunk.
						AABoxc currentChunk = AABoxc::fromPosSize(chunkStart, chunkSize);
						if (fileTargetBlock.overlapsExclusive(currentChunk))
						{
							// The region of the current chunk that will be updated is the intersection of the current chunk and the
							// file block to be written.
							// Here, the chunkUpdateRegion is in file coordinates.
							AABoxc chunkUpdateRegion = fileTargetBlock.intersection(currentChunk);
							
							// Convert chunkUpdateRegion to image coordinates
							Vec3c imageUpdateStartPosition = chunkUpdateRegion.position() - filePosition + imagePosition;
							AABoxc imageUpdateRegion = AABoxc::fromPosSize(imageUpdateStartPosition, chunkUpdateRegion.size());

							// Convert chunkUpdateRegion to chunk coordinates by subtracting chunk start position.
							chunkUpdateRegion = chunkUpdateRegion.translate(-chunkStart);

							// Bug check:
							if (chunkUpdateRegion.minc.x < 0 ||
								chunkUpdateRegion.minc.y < 0 ||
								chunkUpdateRegion.minc.z < 0 ||
								chunkUpdateRegion.maxc.x > chunkSize.x ||
								chunkUpdateRegion.maxc.y > chunkSize.y ||
								chunkUpdateRegion.maxc.z > chunkSize.z)
								throw ITLException(string("Sanity check failed: Chunk update region ") + toString(chunkUpdateRegion) + string(" is inconsistent with chunk size ") + toString(chunkSize) + string("."));

							if (imageUpdateRegion.minc.x < imagePosition.x ||
								imageUpdateRegion.minc.y < imagePosition.y ||
								imageUpdateRegion.minc.z < imagePosition.z ||
								imageUpdateRegion.maxc.x > imagePosition.x + blockDimensions.x ||
								imageUpdateRegion.maxc.y > imagePosition.y + blockDimensions.y ||
								imageUpdateRegion.maxc.z > imagePosition.z + blockDimensions.z
								)
								throw ITLException(string("Sanity check failed: Image update region ") + toString(imageUpdateRegion) + string(" is inconsistent with image position ") + toString(imagePosition) + string(" and update block dimensions") + toString(blockDimensions) + string("."));

							// Write chunk data only if the update region is non-empty.
							if(chunkUpdateRegion.size().min() > 0)
								writeSingleChunk(img, path, chunkIndex, chunkSize, fileDimensions, chunkUpdateRegion.position(), imageUpdateRegion.position(), chunkUpdateRegion.size(), compression);
						}

						//// We need to write the chunk only if the current output chunk overlaps with the region to be written = imageBlock
						//AABox<coord_t> currentChunk = AABox<coord_t>::fromPosSize(chunkStart, chunkSize);
						//if (fileTargetBlock.overlaps(currentChunk))
						//{
						//	if (fileTargetBlock.contains(chunkStart) && fileTargetBlock.contains(chunkStart + chunkSize - Vec3c(1, 1, 1)))
						//	{
						//		// Chunk start is inside file target block, so write start in chunk coords is 0.
						//		// Write start in image coordinates is different from imagePosition.
						//		// Chunk end is in file target block, so we write the entire chunk.
						//		Vec3c startInChunkCoords = Vec3c(0, 0, 0);
						//		Vec3c startInImageCoords = chunkStart - filePosition + imagePosition;
						//		writeSingleChunk(targetImg, path, chunkIndex, chunkSize, fileDimensions, Vec3c(0, 0, 0), startInImageCoords, chunkSize, compression);
						//	}
						//	else if (fileTargetBlock.contains(chunkStart))
						//	{
						//		// Chunk start is inside file target block, so write start in chunk coords is 0.
						//		// Write start in image coordinates is different from imagePosition.
						//		Vec3c startInChunkCoords = Vec3c(0, 0, 0);
						//		Vec3c startInImageCoords = chunkStart - filePosition + imagePosition;
						//		Vec3c endInImageCoords = imagePosition + blockDimensions;
						//		Vec3c writeSize = endInImageCoords - startInImageCoords;
						//		writeSingleChunk(targetImg, path, chunkIndex, chunkSize, fileDimensions, startInChunkCoords, startInImageCoords, writeSize, compression);
						//	}
						//	else
						//	{
						//		// Chunk start is outside of file target block, so write start in chunk coords is not 0.
						//		// Write start in image coordinates is imagePosition.
						//		Vec3c startInChunkCoords = filePosition - chunkStart;
						//		Vec3c endInChunkCoords = chunkSize;
						//		Vec3c startInImageCoords = imagePosition;
						//		Vec3c writeSize = endInChunkCoords - startInChunkCoords;
						//		writeSingleChunk(targetImg, path, chunkIndex, chunkSize, fileDimensions, startInChunkCoords, startInImageCoords, writeSize, compression);
						//	}
						//}
					});
			}


			template<typename pixel_t> void readChunkFile(Image<pixel_t>& img, const string& filename, NN5Compression compression)
			{
				switch (compression)
				{
					case NN5Compression::Raw:
					{
						raw::read(img, filename);
						break;
					}
					case NN5Compression::LZ4:
					{
						lz4::read(img, filename);
						break;
					}
					default:
					{
						throw ITLException(string("Unsupported nn5 decompression algorithm: ") + toString(compression));
					}
				}
			}

			
			template<typename pixel_t> void readFileIntoImageBlock(Image<pixel_t>& img, const string& filename, const Vec3c& imagePosition, NN5Compression compression, Image<pixel_t>& temp)
			{
				// TODO: This is not very efficient due to the copying of the block, and memory allocation, improve?
				//			Note that this could be easily improved using an image view to the desired block in the targetImg.

				readChunkFile(temp, filename, compression);

				copyValues(img, temp, imagePosition);
			}

			/**
			Reads single NN5 chunk file.
			*/
			template<typename pixel_t> void readSingleChunk(Image<pixel_t>& target, const std::string& path, const Vec3c& datasetDimensions, const Vec3c& chunkIndex, const Vec3c& chunkStartInTarget, const Vec3c& readSize, NN5Compression compression, Image<pixel_t>& temp)
			{
				string dir = chunkFolder(path, getDimensionality(datasetDimensions), chunkIndex);

				//Vec3c chunkEnd = chunkStartInTarget + readSize;
				//for (size_t n = 0; n < chunkEnd.size(); n++)
				//{
				//	if (chunkEnd[n] > target.dimension(n))
				//		chunkEnd[n] = target.dimension(n);
				//}
				//Vec3c realReadSize = chunkEnd - chunkStartInTarget;

				// Search for files in the directory
				std::vector<string> files = getFileList(dir);

				if (files.size() <= 0)
				{
					// No file => all pixels in the block are zeroes.
					draw<pixel_t>(target, AABoxc::fromPosSize(chunkStartInTarget, readSize), (pixel_t)0);
				}
				else if (files.size() == 1)
				{
					string filename = dir + "/" + files[0];
					readFileIntoImageBlock(target, filename, chunkStartInTarget, compression, temp);
				}
				else
				{
					throw ITLException(string("Multiple image files found in block directory ") + dir);
				}
				
			}

			/**
			Reads NN5 chunk files.
			*/
			template<typename pixel_t> void readChunks(Image<pixel_t>& img, const std::string& path, const Vec3c& chunkSize, NN5Compression compression, bool showProgressInfo)
			{
				Image<pixel_t> temp;

				forAllChunks(img.dimensions(), chunkSize, showProgressInfo, [&](const Vec3c& chunkIndex, const Vec3c& chunkStart)
					{
						readSingleChunk(img, path, img.dimensions(), chunkIndex, chunkStart, chunkSize, compression, temp);
					});
			}


			/**
			Reads NN5 chunk files.
			*/
			template<typename pixel_t> void readChunksInRange(Image<pixel_t>& img, const std::string& path, const Vec3c& chunkSize, NN5Compression compression,
				const Vec3c& datasetDimensions,
				const Vec3c& start, const Vec3c& end,
				bool showProgressInfo)
			{
				Image<pixel_t> temp;

				AABoxc imageBox = AABoxc::fromMinMax(start, end);

				// This is a check-all-chunks algoritm. Alternatively, we could calculate the required chunk range.
				forAllChunks(datasetDimensions, chunkSize, showProgressInfo, [&](const Vec3c& chunkIndex, const Vec3c& chunkStart)
					{
						AABox<coord_t> currentChunk = AABox<coord_t>::fromPosSize(chunkStart, chunkSize);
						if (currentChunk.overlapsExclusive(imageBox))
							readSingleChunk(img, path, datasetDimensions, chunkIndex, chunkStart - start, currentChunk.intersection(imageBox).size(), compression, temp);
					});
			}

			/**
			Checks that chunk size is valid and if not, throws an exception.
			*/
			void check(const Vec3c& chunkSize);

			/**
			Checks that provided NN5 information is correct and writes metadata.
			@param deleteOldData Contents of the dataset are usually deleted when a write process is started. Set this to false to delete only if the path contains an incompatible dataset, but keep the dataset if it seems to be the same than the current image.
			*/
			void beginWrite(const Vec3c& imageDimensions, ImageDataType imageDataType, const std::string& path, const Vec3c& chunkSize, NN5Compression compression, bool deleteOldData);
		

			template<typename pixel_t> void write(const Image<pixel_t>& img, const std::string& path, const Vec3c& chunkSize, NN5Compression compression, bool deleteOldData, bool showProgressInfo)
			{
				internals::beginWrite(img.dimensions(), img.dataType(), path, chunkSize, compression, deleteOldData);

				// Write data
				internals::writeChunks(img, path, chunkSize, compression, img.dimensions(), showProgressInfo);
			}
		}

		bool getInfo(const std::string& path, Vec3c& dimensions, bool& isNativeByteOrder, ImageDataType& dataType, Vec3c& chunkSize, NN5Compression& compression, std::string& reason);

		inline bool getInfo(const std::string& path, Vec3c& dimensions, ImageDataType& dataType, std::string& reason)
		{
			bool dummyIsNative;
			Vec3c dummyChunkSize;
			NN5Compression dummyCompression;
			return itl2::zarr::getInfo(path, dimensions, dummyIsNative, dataType, dummyChunkSize, dummyCompression, reason);
		}

		/**
		Write an image to an nn5 dataset.
		@param targetImg Image to write.
		@param path Name of the top directory of the nn5 dataset.
		*/
		template<typename pixel_t> void write(const Image<pixel_t>& img, const std::string& path, const Vec3c& chunkSize, NN5Compression compression, bool showProgressInfo = false)
		{
			internals::write(img, path, chunkSize, compression, false, showProgressInfo);
		}


		/**
		Write an image to an nn5 dataset.
		@param targetImg Image to write.
		@param path Name of the top directory of the nn5 dataset.
		*/
		template<typename pixel_t> void write(const Image<pixel_t>& img, const std::string& path, const Vec3c& chunkSize, bool showProgressInfo = false)
		{
			zarr::write(img, path, chunkSize, NN5Compression::LZ4, showProgressInfo);
		}

		/**
		Default chunk size for NN5 dataset.
		*/
		inline const Vec3c DEFAULT_CHUNK_SIZE = Vec3c(1536, 1536, 1536);

		/**
		Write an image to an nn5 dataset.
		@param targetImg Image to write.
		@param path Name of the top directory of the nn5 dataset.
		*/
		template<typename pixel_t> void write(const Image<pixel_t>& img, const std::string& path, bool showProgressInfo = false)
		{
			write(img, path, DEFAULT_CHUNK_SIZE, showProgressInfo);
		}

		/**
		Writes a block of an image to the specified location in an .nn5 dataset.
		The output dataset is not truncated if it exists.
		If the output file does not exist, it is created.
		@param targetImg Image to write.
		@param filename Name of file to write.
		@param filePosition Position in the file to write to.
		@param fileDimension Total dimensions of the entire output file.
		@param imagePosition Position in the image where the block to be written starts.
		@param blockDimensions Dimensions of the block of the source image to write.
		*/
		template<typename pixel_t> void writeBlock(Image<pixel_t>& img, const std::string& path, const Vec3c& chunkSize, NN5Compression compression,
			const Vec3c& filePosition, const Vec3c& fileDimensions,
			const Vec3c& imagePosition,
			const Vec3c& blockDimensions,
			bool showProgressInfo = false)
		{
			internals::check(chunkSize);

			fs::create_directories(path);

			// Write metadata
			internals::writeMetadata(path, fileDimensions, img.dataType(), chunkSize, compression);

			// Write data
			internals::writeChunksInRange(img, path, chunkSize, compression, filePosition, fileDimensions, imagePosition, blockDimensions, showProgressInfo);
		}



		/**
		Reads an nn5 dataset file to the given image.
		@param targetImg Image where the data is placed. The size of the image will be set based on the dataset contents.
		@param path Path to the root of the nn5 dataset.
		*/
		template<typename pixel_t> void read(Image<pixel_t>& img, const std::string& path, bool showProgressInfo = false)
		{
			bool isNativeByteOrder;
			Vec3c dimensions;
			ImageDataType dataType;
			Vec3c chunkSize;
			NN5Compression compression;
			string reason;
			if (!itl2::zarr::getInfo(path, dimensions, isNativeByteOrder, dataType, chunkSize, compression, reason))
				throw ITLException(string("Unable to read nn5 dataset: ") + reason);

			if (dataType != img.dataType())
				throw ITLException(string("Expected data type is ") + toString(img.dataType()) + " but the nn5 dataset contains data of type " + toString(dataType) + ".");

			img.ensureSize(dimensions);

			internals::readChunks(img, path, chunkSize, compression, showProgressInfo);

			if (!isNativeByteOrder)
				swapByteOrder(img);
		}

		/**
		Reads a part of a .nn5 dataset to the given image.
		NOTE: Does not support out of bounds start position.
		@param targetImg Image where the data is placed. The size of the image defines the size of the block that is read.
		@param filename The name of the dataset to read.
		@param fileStart Start location of the read in the file. The size of the image defines the size of the block that is read.
		*/
		template<typename pixel_t> void readBlock(Image<pixel_t>& img, const std::string& path, const Vec3c& fileStart, bool showProgressInfo = false)
		{
			bool isNativeByteOrder;
			Vec3c fileDimensions;
			ImageDataType dataType;
			Vec3c chunkSize;
			NN5Compression compression;
			string reason;
			if (!itl2::zarr::getInfo(path, fileDimensions, isNativeByteOrder, dataType, chunkSize, compression, reason))
				throw ITLException(string("Unable to read nn5 dataset: ") + reason);

			if (dataType != img.dataType())
				throw ITLException(string("Expected data type is ") + toString(img.dataType()) + " but the nn5 dataset contains data of type " + toString(dataType) + ".");


			if (fileStart.x < 0 || fileStart.y < 0 || fileStart.z < 0 || fileStart.x >= fileDimensions.x || fileStart.y >= fileDimensions.y || fileStart.z >= fileDimensions.z)
				throw ITLException("Out of bounds start position in readBlock.");

			Vec3c cStart = fileStart;
			clamp(cStart, Vec3c(0, 0, 0), fileDimensions);
			Vec3c cEnd = fileStart + img.dimensions();
			clamp(cEnd, Vec3c(0, 0, 0), fileDimensions);

			if (cStart == Vec3c(0, 0, 0) && cEnd == fileDimensions && fileDimensions == img.dimensions())
			{
				// Reading whole dataset, use the whole dataset reading function.
				read(img, path);
				return;
			}

			internals::readChunksInRange(img, path, chunkSize, compression, fileDimensions, cStart, cEnd, showProgressInfo);

			if (!isNativeByteOrder)
				swapByteOrder(img);
		}

		

		struct NN5Process
		{
			AABoxc readBlock;
			AABoxc writeBlock;
		};

		/**
		Enables concurrent access from multiple processes for an existing or a new NN5 dataset.
		This function should be called before the processes are started.
		@param imageDimensions Dimensions of the image to be saved into the NN5 dataset.
		@param imageDataType Data type of the image.
		@param path Path to the NN5 dataset.
		@param chunkSize Chunk size for the NN5 dataset.
		@param compression Compression method to be used.
		@param processes A list of NN5Process objects that define the block that where each process will have read and write access. The blocks may overlap.
		@return Number of chunks that require special processing in endConcurrentWrite.
		*/
		size_t startConcurrentWrite(const Vec3c& imageDimensions, ImageDataType imageDataType, const std::string& path, const Vec3c& chunkSize, NN5Compression compression, const std::vector<NN5Process>& processes);

		/**
		Enables concurrent access from multiple processes for an existing or a new NN5 dataset.
		This function should be called before the processes are started.
		@param targetImg Image that is to be saved into the NN5 dataset by the processes.
		@param path Path to the NN5 dataset.
		@param chunkSize Chunk size for the NN5 dataset.
		@param compression Compression method to be used.
		@param processes A list of NN5Process objects that define the block that where each process will have read and write access. The blocks may overlap.
		*/
		template<typename pixel_t> void startConcurrentWrite(const Image<pixel_t>& img, const std::string& path, const Vec3c& chunkSize, NN5Compression compression, const std::vector<NN5Process>& processes)
		{
			startConcurrentWrite(img.dimensions(), img.dataType(), path, chunkSize, compression, processes);
		}

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
			void nn5Metadata();
			void nn5io();
			void nn5BlockIo();
			void concurrency();
			void concurrencyLong();
		}
	}

}