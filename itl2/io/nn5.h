#pragma once

#include "filesystem.h"
#include "image.h"
#include "utilities.h"
#include "io/raw.h"
#include "io/itllz4.h"
#include "aabox.h"

namespace itl2
{
	namespace io
	{
		// Forward declaration of io::getInfo to avoid header loop.
		bool getInfo(const std::string& filename, Vec3c& dimensions, ImageDataType& dataType, std::string& reason);
	}

	namespace nn5
	{
		enum class NN5Compression
		{
			Raw,
			LZ4
		};

		namespace internals
		{
			/**
			Constructs NN5 dataset metadata filename.
			*/
			inline std::string nn5MetadataFilename(const string& path)
			{
				return path + "/metadata.json";
			}

			/**
			Writes NN5 metadata file.
			*/
			void writeMetadata(const std::string& path, const Vec3c& dimensions, ImageDataType dataType, const Vec3c& chunkSize, NN5Compression compression);

			/**
			Retrieve a list of all files in the given directory.
			*/
			std::vector<string> getFileList(const std::string& dir);

			/**
			Writes single NN5 chunk file.
			@param chunkIndex Index of the chunk to write. This is used to determine the correct output folder.
			@param chunkSize Chunk size of the dataset.
			@param startInChunkCoords Start position of the block to be written in the coordinates of the chunk.
			@param startInImageCoords Start position of the block to be written in the coordinates of image img.
			@param writeSize Size of block to be written.
			*/
			template<typename pixel_t> void writeSingleChunk(const Image<pixel_t>& img, const std::string& path, const Vec3c& chunkIndex, const Vec3c& chunkSize, const Vec3c& datasetSize,
				const Vec3c& startInChunkCoords, const Vec3c& startInImageCoords, const Vec3c& writeSize, NN5Compression compression)
			{
				string filename = path;
				for (size_t n = 0; n < img.dimensionality(); n++)
					filename += string("/") + toString(chunkIndex[n]);
				filename += "/chunk";

				// Clamp write size to the size of the image.
				Vec3c imageChunkEnd = startInImageCoords + writeSize;
				for (size_t n = 0; n < imageChunkEnd.size(); n++)
				{
					if (imageChunkEnd[n] > img.dimension(n))
						imageChunkEnd[n] = img.dimension(n);
				}
				Vec3c realWriteSize = imageChunkEnd - startInImageCoords;

				// Determine if this is the last chunk, and reduce chunk size accordingly so that image size stays correct.
				Vec3c datasetChunkStart = chunkIndex.componentwiseMultiply(chunkSize);
				Vec3c datasetChunkEnd = datasetChunkStart + chunkSize;
				for (size_t n = 0; n < datasetChunkEnd.size(); n++)
				{
					if (datasetChunkEnd[n] > datasetSize[n])
						datasetChunkEnd[n] = datasetSize[n];
				}
				Vec3c realChunkSize = datasetChunkEnd - datasetChunkStart;

				realWriteSize = min(realWriteSize, realChunkSize);

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
						filename += ".lz4";
						lz4::writeBlock(img, filename, startInChunkCoords, realChunkSize, startInImageCoords, realWriteSize);
						break;
					}
					default:
					{
						throw ITLException(string("Unsupported nn5 compression algorithm: ") + toString(compression));
					}
				}

			}

			/**
			Writes NN5 chunk files.
			*/
			template<typename pixel_t> void writeChunks(const Image<pixel_t>& img, const std::string& path, const Vec3c& chunkSize, NN5Compression compression, const Vec3c& datasetSize)
			{
				Vec3c chunkStart(0, 0, 0);
				Vec3c chunkIndex(0, 0, 0);
				while (chunkStart.z < img.depth())
				{
					while (chunkStart.y < img.height())
					{
						while (chunkStart.x < img.width())
						{
							writeSingleChunk(img, path, chunkIndex, chunkSize, datasetSize, Vec3c(0, 0, 0), chunkStart, chunkSize, compression);

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

			template<typename pixel_t> void writeChunksInRange(Image<pixel_t>& img, const std::string& path, const Vec3c& chunkSize, NN5Compression compression,
				const Vec3c& filePosition, const Vec3c& fileDimensions,
				const Vec3c& imagePosition,
				const Vec3c& blockDimensions)
			{
				// Writes block of image defined by (imagePosition, blockDimensions) to file (defined by path),
				// to location defined by filePosition.

				AABox<coord_t> imageBlock = AABox<coord_t>::fromPosSize(imagePosition, blockDimensions);

				AABox<coord_t> fileTargetBlock = AABox<coord_t>::fromPosSize(filePosition, blockDimensions);

				Vec3c chunkStart(0, 0, 0);
				Vec3c chunkIndex(0, 0, 0);
				while (chunkStart.z < fileDimensions.z)
				{
					while (chunkStart.y < fileDimensions.y)
					{
						while (chunkStart.x < fileDimensions.x)
						{
							// We need to write the chunk only if the current output chunk overlaps with the region to be written = imageBlock
							AABox<coord_t> currentChunk = AABox<coord_t>::fromPosSize(chunkStart, chunkSize);
							if (fileTargetBlock.overlaps(currentChunk))
							{
								if (fileTargetBlock.contains(chunkStart) && fileTargetBlock.contains(chunkStart + chunkSize))
								{
									// Chunk start is inside file target block, so write start in chunk coords is 0.
									// Write start in image coordinates is different from imagePosition.
									// Chunk end is in file target block, so we write the entire chunk.
									Vec3c startInChunkCoords = Vec3c(0, 0, 0);
									Vec3c startInImageCoords = chunkStart - filePosition + imagePosition;
									writeSingleChunk(img, path, chunkIndex, chunkSize, fileDimensions, Vec3c(0, 0, 0), startInImageCoords, chunkSize, compression);
								}
								else if (fileTargetBlock.contains(chunkStart))
								{
									// Chunk start is inside file target block, so write start in chunk coords is 0.
									// Write start in image coordinates is different from imagePosition.
									Vec3c startInChunkCoords = Vec3c(0, 0, 0);
									Vec3c startInImageCoords = chunkStart - filePosition + imagePosition;
									writeSingleChunk(img, path, chunkIndex, chunkSize, fileDimensions, startInChunkCoords, startInImageCoords, chunkSize, compression);
								}
								else
								{
									// Chunk start is outside of file target block, so write start in chunk coords is not 0.
									// Write start in image coordinates is imagePosition.
									Vec3c startInChunkCoords = filePosition - chunkStart;
									Vec3c startInImageCoords = imagePosition;
									writeSingleChunk(img, path, chunkIndex, chunkSize, fileDimensions, startInChunkCoords, startInImageCoords, chunkSize, compression);
								}
							}

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

			// TODO: This is not very efficient due to the copying of the block, and memory allocation, improve!
			template<typename pixel_t> void readFileIntoImageBlock(Image<pixel_t>& img, const string& filename, const Vec3c& imagePosition, const Vec3c& blockSize, NN5Compression compression, Image<pixel_t>& temp)
			{
				temp.ensureSize(blockSize);
				
				switch (compression)
				{
					case NN5Compression::Raw:
					{
						raw::read(temp, filename);
						break;
					}
					case NN5Compression::LZ4:
					{
						lz4::read(temp, filename);
						break;
					}
					default:
					{
						throw ITLException(string("Unsupported nn5 decompression algorithm: ") + toString(compression));
					}
				}

				copyValues(img, temp, imagePosition);
			}

			/**
			Reads single NN5 chunk file.
			*/
			template<typename pixel_t> void readSingleChunk(Image<pixel_t>& img, const std::string& path, const Vec3c& chunkIndex, const Vec3c& chunkStart, const Vec3c& chunkSize, NN5Compression compression, Image<pixel_t>& temp)
			{
				string dir = path;
				for (size_t n = 0; n < img.dimensionality(); n++)
					dir += string("/") + toString(chunkIndex[n]);

				Vec3c chunkEnd = chunkStart + chunkSize;
				for (size_t n = 0; n < chunkEnd.size(); n++)
				{
					if (chunkEnd[n] > img.dimension(n))
						chunkEnd[n] = img.dimension(n);
				}
				Vec3c realChunkSize = chunkEnd - chunkStart;

				// Search for files in the directory
				std::vector<string> files = getFileList(dir);

				if (files.size() <= 0)
				{
					// No file => all pixels in the block are zeroes.
					setValue(img, (pixel_t)0);
				}
				else if (files.size() == 1)
				{
					string filename = dir + "/" + files[0];
					readFileIntoImageBlock(img, filename, chunkStart, realChunkSize, compression, temp);
				}
				else
				{
					throw ITLException(string("Multiple image files found in block directory ") + dir);
				}
				
			}

			/**
			Reads NN5 chunk files.
			*/
			template<typename pixel_t> void readChunks(Image<pixel_t>& img, const std::string& path, const Vec3c& chunkSize, NN5Compression compression)
			{
				Image<pixel_t> temp;

				Vec3c chunkStart(0, 0, 0);
				Vec3c chunkIndex(0, 0, 0);
				while (chunkStart.z < img.depth())
				{
					while (chunkStart.y < img.height())
					{
						while (chunkStart.x < img.width())
						{
							readSingleChunk(img, path, chunkIndex, chunkStart, chunkSize, compression, temp);

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
			Reads NN5 chunk files.
			*/
			template<typename pixel_t> void readChunksInRange(Image<pixel_t>& img, const std::string& path, const Vec3c& chunkSize, NN5Compression compression,
				const Vec3c& datasetDimensions,
				const Vec3c& start, const Vec3c& end)
			{
				Image<pixel_t> temp;

				AABox<coord_t> imageBox(start, end);

				// This is a check-all-chunks algoritm. Alternatively, we could calculate the required chunk range.
				Vec3c chunkStart(0, 0, 0);
				Vec3c chunkIndex(0, 0, 0);
				while (chunkStart.z < datasetDimensions.z)
				{
					while (chunkStart.y < datasetDimensions.y)
					{
						while (chunkStart.x < datasetDimensions.x)
						{
							AABox<coord_t> currentChunk = AABox<coord_t>::fromPosSize(chunkStart, chunkSize);
							if(currentChunk.overlaps(imageBox))
								readSingleChunk(img, path, chunkIndex, chunkStart - start, chunkSize, compression, temp);

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
		}

		

	}

	template<>
	inline std::string toString(const nn5::NN5Compression& x)
	{
		switch (x)
		{
		case nn5::NN5Compression::Raw: return "Raw";
		case nn5::NN5Compression::LZ4: return "LZ4Raw";
		}
		throw ITLException("Invalid nn5 compression type.");
	}

	template<>
	inline nn5::NN5Compression fromString(const string& str0)
	{
		string str = str0;
		toLower(str);
		if (str == "raw")
			return nn5::NN5Compression::Raw;
		if (str == "lz4raw")
			return nn5::NN5Compression::LZ4;

		throw ITLException(string("Invalid nn5 compression type: ") + str);
	}

	namespace nn5
	{

		bool getInfo(const std::string& path, Vec3c& dimensions, bool& isNativeByteOrder, ImageDataType& dataType, Vec3c& chunkSize, NN5Compression& compression, std::string& reason);

		inline bool getInfo(const std::string& path, Vec3c& dimensions, ImageDataType& dataType, std::string& reason)
		{
			bool dummyIsNative;
			Vec3c dummyChunkSize;
			NN5Compression dummyCompression;
			return getInfo(path, dimensions, dummyIsNative, dataType, dummyChunkSize, dummyCompression, reason);
		}

		/**
		Write an image to an nn5 dataset.
		@param img Image to write.
		@param path Name of the top directory of the nn5 dataset.
		*/
		template<typename pixel_t> void write(const Image<pixel_t>& img, const std::string& path, const Vec3c& chunkSize, NN5Compression compression)
		{
			if (chunkSize.min() <= 0)
				throw ITLException(string("NN5 chunk size must be positive, but it is ") + toString(chunkSize));

			// Delete old dataset if it exists.
			if (fs::exists(path))
			{
				Vec3c dummyDimensions;
				ImageDataType dummyDatatype;
				string dummyReason;
				if (!io::getInfo(path, dummyDimensions, dummyDatatype, dummyReason))
					throw ITLException(string("Unable to write an NN5 as the dataset already exists but cannot be verified to be an image: ") + path + " Consider removing the existing dataset manually.");
				fs::remove_all(path);
			}

			fs::create_directories(path);

			// Write metadata
			internals::writeMetadata(path, img.dimensions(), img.dataType(), chunkSize, compression);
			
			// Write data
			internals::writeChunks(img, path, chunkSize, compression, img.dimensions());
		}


		/**
		Write an image to an nn5 dataset.
		@param img Image to write.
		@param path Name of the top directory of the nn5 dataset.
		*/
		template<typename pixel_t> void write(const Image<pixel_t>& img, const std::string& path, const Vec3c& chunkSize)
		{
			write(img, path, chunkSize, NN5Compression::LZ4);
		}

		/**
		Default chunk size for NN5 dataset.
		*/
		inline const Vec3c DEFAULT_CHUNK_SIZE = Vec3c(1536, 1536, 1536);

		/**
		Write an image to an nn5 dataset.
		@param img Image to write.
		@param path Name of the top directory of the nn5 dataset.
		*/
		template<typename pixel_t> void write(const Image<pixel_t>& img, const std::string& path)
		{
			write(img, path, DEFAULT_CHUNK_SIZE);
		}



		/**
		Reads an nn5 dataset file to the given image.
		@param img Image where the data is placed. The size of the image will be set based on the dataset contents.
		@param path Path to the root of the nn5 dataset.
		*/
		template<typename pixel_t> void read(Image<pixel_t>& img, const std::string& path)
		{
			bool isNativeByteOrder;
			Vec3c dimensions;
			ImageDataType dataType;
			Vec3c chunkSize;
			NN5Compression compression;
			string reason;
			if (!getInfo(path, dimensions, isNativeByteOrder, dataType, chunkSize, compression, reason))
				throw ITLException(string("Unable to read nn5 dataset: ") + reason);

			if (dataType != img.dataType())
				throw ITLException(string("Expected data type is ") + toString(img.dataType()) + " but the nn5 dataset contains data of type " + toString(dataType) + ".");

			img.ensureSize(dimensions);

			internals::readChunks(img, path, chunkSize, compression);

			if (!isNativeByteOrder)
				swapByteOrder(img);
		}

		/**
		Reads a part of a .nn5 dataset to the given image.
		NOTE: Does not support out of bounds start position.
		@param img Image where the data is placed. The size of the image defines the size of the block that is read.
		@param filename The name of the dataset to read.
		@param fileStart Start location of the read in the file. The size of the image defines the size of the block that is read.
		*/
		template<typename pixel_t> void readBlock(Image<pixel_t>& img, const std::string& path, const Vec3c& fileStart)
		{
			bool isNativeByteOrder;
			Vec3c fileDimensions;
			ImageDataType dataType;
			Vec3c chunkSize;
			NN5Compression compression;
			string reason;
			if (!getInfo(path, fileDimensions, isNativeByteOrder, dataType, chunkSize, compression, reason))
				throw ITLException(string("Unable to read nn5 dataset: ") + reason);

			if (dataType != img.dataType())
				throw ITLException(string("Expected data type is ") + toString(img.dataType()) + " but the nn5 dataset contains data of type " + toString(dataType) + ".");


			if (fileStart.x < 0 || fileStart.y < 0 || fileStart.z < 0 || fileStart.x >= fileDimensions.x || fileStart.y >= fileDimensions.y || fileStart.z >= fileDimensions.z)
				throw ITLException("Out of bounds start position in nn5::readBlock.");

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

			internals::readChunksInRange(img, path, chunkSize, compression, fileDimensions, cStart, cEnd);

			if (!isNativeByteOrder)
				swapByteOrder(img);
		}

		/**
		Writes a block of an image to the specified location in an .nn5 dataset.
		The output dataset is not truncated if it exists.
		If the output file does not exist, it is created.
		@param img Image to write.
		@param filename Name of file to write.
		@param filePosition Position in the file to write to.
		@param fileDimension Total dimensions of the entire output file.
		@param imagePosition Position in the image where the block to be written starts.
		@param blockDimensions Dimensions of the block of the source image to write.
		*/
		template<typename pixel_t> void writeBlock(Image<pixel_t>& img, const std::string& path, const Vec3c& chunkSize, NN5Compression compression,
			const Vec3c& filePosition, const Vec3c& fileDimensions, 
			const Vec3c& imagePosition,
			const Vec3c& blockDimensions)
		{
			fs::create_directories(path);

			// Write metadata
			internals::writeMetadata(path, fileDimensions, img.dataType(), chunkSize, compression);

			// Write data
			internals::writeChunksInRange(img, path, chunkSize, compression, filePosition, fileDimensions, imagePosition, blockDimensions);
		}


		//void startConcurrentWrite(const std::string& path, something that identifies reading and writing region for each process)
		//{

		//}

		//void endConcurrentWriteForBlock(const std::string& path, const Vec3c& blockIndex)
		//{
		//	endConcurrentWrite but for single block only
		//}

		//void endConcurrentWrite(const std::string& path)
		//{
		//	combine files in writes folders
		//}

		// TODO:
		// writeBlock (if writes folder -> write there; if not -> replace chunk)
		// startConcurrentWrite (create writes folders)
		// endConcurrentWrite (combine files in writes folders)
		// endConcurrentWriteForBlock (endConcurrentWrite but for single block only)

		namespace tests
		{
			void nn5Metadata();
			void nn5io();
			void nn5BlockIo();
		}
	}

}