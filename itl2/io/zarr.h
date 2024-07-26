#pragma once
#include <list>
#include <blosc.h>

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
			void writeMetadata(const std::string& path, const Vec3c& dimensions, ImageDataType dataType, const Vec3c& chunkSize, int fillValue, const std::list<ZarrCodec>& codecs);

			inline std::string chunkFile(const std::string& path, size_t dimensionality, const Vec3c& chunkIndex)
			{
				//TODO: use separator from zarr metadata
				string separator = string("/");
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

			/**
			Writes single Zarr chunk file.
			@param chunkIndex Index of the chunk to write. This is used to determine the correct output folder.
			@param chunkSize Chunk size of the dataset.
			@param startInChunkCoords Start position of the block to be written in the coordinates of the chunk.
			@param chunkStart Start position of the block to be written in the coordinates of image targetImg.
			@param writeSize Size of block to be written.
			*/
			template<typename pixel_t>
			void writeSingleChunk(const Image<pixel_t>& img, const std::string& path, const Vec3c& chunkIndex, const Vec3c& chunkShape, const Vec3c& datasetSize,
				Vec3c startInChunkCoords, const Vec3c& startInImageCoords, const Vec3c& writeSize, int fillValue, std::list<ZarrCodec>& codecs)
			{
				//todo: read this chunk and only update the requested region by const startInImageCoords,startInChunkCoords, writeSize
				Vec3c transposedChunkShape = chunkShape;
				Image<pixel_t> imgChunk(chunkShape);
				const Vec3c chunkStartPos = chunkIndex.componentwiseMultiply(chunkShape);
				//write all pixels of chunk in img to imgChunk
				forAllPixels(imgChunk, [&](coord_t x, coord_t y, coord_t z)
				{
				  Vec3c pos = Vec3c(x, y, z) + chunkStartPos;
				  imgChunk(x, y, z) = img(pos);
				});

				string filename = chunkFile(path, getDimensionality(datasetSize), chunkIndex);

				// Clamp write size to the size of the image.
				Vec3c imageChunkEnd = startInImageCoords + writeSize;
				for (size_t n = 0; n < imageChunkEnd.size(); n++)
				{
					if (imageChunkEnd[n] > img.dimension(n))
						imageChunkEnd[n] = img.dimension(n);
				}
				Vec3c realWriteSize = imageChunkEnd - startInImageCoords;
				Vec3c realChunkSize = clampedChunkSize(chunkIndex, chunkShape, datasetSize);
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

				std::list<ZarrCodec>::iterator codec = codecs.begin();
				for (; codec->type == ZarrCodecType::ArrayArrayCodec; ++codec)
				{
					if (codec->name == ZarrCodecName::Transpose)
					{
						Vec3c order = codec->transposeOrder();
						transposedChunkShape = transposedChunkShape.transposed(order);
						transpose(imgChunk, order, fillValue);
					}
					else throw ITLException("ArrayArrayCodec: " + toString(codec->name) + " not yet implemented: ");
				}
				assert(codec->type == ZarrCodecType::ArrayBytesCodec);
				std::vector<char> buffer;
				if (codec->name == ZarrCodecName::Bytes)
				{
					writeBytesCodec(imgChunk, buffer);
				}
				else throw ITLException("ArrayBytesCodec: " + toString(codec->name) + " not yet implemented: ");
				++codec;
				for (; codec != codecs.end(); ++codec)
				{
					assert(codec->type == ZarrCodecType::BytesBytesCodec);
					if (codec->name == ZarrCodecName::Blosc)
					{
						if (!codec->configuration.contains("typesize"))
							codec->configuration["typesize"] = sizeof(pixel_t);

						size_t destSize = buffer.size()+BLOSC_MIN_HEADER_LENGTH;
						size_t srcSize = buffer.size();
						cout << "compressing buffer.size()="<<buffer.size()<<" sizeof(pixel_t)= "<<sizeof(pixel_t)<<endl;
						std::vector<char> temp(destSize);
						string cname;
						int clevel;
						blosc::shuffle shuffle;
						size_t typesize;
						size_t blocksize;
						codec->getBloscConfiguration(cname, clevel, shuffle, typesize, blocksize);
						int rcode = blosc_set_compressor(cname.c_str());
						if (rcode < 0) {
						  throw ITLException("Blosc Error setting "+ cname + " compressor.  It really exists?");
						}
						cout << "Using blosc compressor" << cname << endl;
						blosc_set_blocksize(blocksize);
						//todo: use blosc_compress_ctx and no blosc_init for multithreaded
						size_t realDestSize = blosc_compress(clevel, (int)shuffle, typesize, srcSize, buffer.data(), temp.data(), destSize);

						if (realDestSize == 0) cout << "Buffer is incompressible.  Giving up." << endl;
						else if (realDestSize < 0) throw ITLException("Compression error.  Error code: " + toString(realDestSize));
						else cout << "Compression: " << srcSize << " -> " << realDestSize << " (" << static_cast<double>(srcSize) / realDestSize << "x)" << endl;

						buffer.resize(realDestSize);
						std::memcpy(buffer.data(), temp.data(), realDestSize);
					}
					else throw ITLException("BytesBytesCodec: " + toString(codec->name) + " not yet implemented: ");
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
							std::cout << "forAllChunks chunkIndex=" << chunkIndex << " chunkStart=" << chunkStart << std::endl;
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
			void writeChunksInRange(const Image<pixel_t>& img, const std::string& path, const Vec3c& chunkShape, int fillValue, std::list<ZarrCodec>& codecs,
				const Vec3c& filePosition, const Vec3c& fileDimensions,
				const Vec3c& imagePosition,
				const Vec3c& blockDimensions,
				bool showProgressInfo)
			{
				// Writes block of image defined by (imagePosition, blockDimensions) to file (defined by path),
				// to location defined by filePosition.

				bool is_blosc_init = false;
				for (auto codec : codecs)
				{
					if (codec.name == ZarrCodecName::Blosc && !is_blosc_init)
					{
						printf("Blosc version info: %s (%s)\n", BLOSC_VERSION_STRING, BLOSC_VERSION_DATE);
						blosc_init();
						is_blosc_init = true;
					}
				}

				AABoxc imageBlock = AABoxc::fromPosSize(imagePosition, blockDimensions);

				AABoxc fileTargetBlock = AABoxc::fromPosSize(filePosition, blockDimensions);

				forAllChunks(fileDimensions, chunkShape, showProgressInfo, [&](const Vec3c& chunkIndex, const Vec3c& chunkStart)
				{
				  // This is done for all chunks in the output file.
				  // We will need to update the chunk if the file block to be written (fileTargetBlock)
				  // overlaps the current chunk.
				  AABoxc currentChunk = AABoxc::fromPosSize(chunkStart, chunkShape);
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

					  // Write chunk data only if the update region is non-empty.
					  if (chunkUpdateRegion.size().min() > 0)
						  writeSingleChunk(img, path, chunkIndex, chunkShape, fileDimensions, chunkUpdateRegion.position(), imageUpdateRegion.position(), chunkUpdateRegion.size(), fillValue, codecs);
				  }
				});
				if (is_blosc_init)
				{
					blosc_destroy();
				}
			}

			/**
			Reads zarr chunk files.
			*/
			template<typename pixel_t>
			void readChunksInRange(Image<pixel_t>& img, const std::string& path,
				const Vec3c& chunkShape, int fillValue, std::list<ZarrCodec>& codecs,
				const Vec3c& datasetShape,
				const Vec3c& start, const Vec3c& end,
				bool showProgressInfo)
			{
				Vec3c transposedChunkShape = chunkShape;
				bool is_blosc_init = false;
				for (auto codec : codecs)
				{
					if (codec.name == ZarrCodecName::Blosc && !is_blosc_init)
					{
						printf("Blosc version info: %s (%s)\n", BLOSC_VERSION_STRING, BLOSC_VERSION_DATE);
						blosc_init();
						is_blosc_init = true;
					}
					if (codec.name == ZarrCodecName::Transpose)
					{
						transposedChunkShape = transposedChunkShape.transposed(codec.transposeOrder());
					}
				}
				cout << "transposedChunkShape= " << transposedChunkShape << endl;
				const AABoxc imageBox = AABoxc::fromMinMax(start, end);
				forAllChunks(datasetShape, chunkShape, showProgressInfo, [&](const Vec3c& chunkIndex, const Vec3c& chunkStart)
				{
				  AABox<coord_t> currentChunk = AABox<coord_t>::fromPosSize(chunkStart, chunkShape);
				  const Vec3c chunkStartInTarget = chunkStart - start;
				  const Vec3c readSize = currentChunk.intersection(imageBox).size();

				  if (currentChunk.overlapsExclusive(imageBox))
				  {
					  Vec3c currentChunkShape = transposedChunkShape;
					  Image<pixel_t> imgChunk(transposedChunkShape);
					  string filename = chunkFile(path, getDimensionality(datasetShape), chunkIndex);
					  if (fs::is_directory(filename))
					  {
						  throw ITLException(filename + string(" is a directory, but it should be a file."));
					  }
					  if (fs::is_regular_file(filename))
					  {
						  std::vector<char> buffer = readBytesOfFile(filename);
						  // iterate over codecs
						  std::list<ZarrCodec>::reverse_iterator codec = codecs.rbegin();
						  for (; codec->type == ZarrCodecType::BytesBytesCodec; ++codec)
						  {
							  if (codec->name == ZarrCodecName::Blosc)
							  {
								  size_t srcSize = buffer.size();
								  size_t destSize;
								  if (blosc_cbuffer_validate(buffer.data(), srcSize, &destSize) < 0){
									  throw ITLException("blosc_decompress error: \"Buffer does not contain valid blosc-encoded contents\"");
								  }
								  std::vector<char> temp(destSize);
								  size_t realDestSize = blosc_decompress(buffer.data(), temp.data(), destSize);
								  if (realDestSize < 0)
								  {
									  throw ITLException("blosc_decompress error.  Error code: " + toString(realDestSize));
								  }
								  cout << "buffer[0]= " << buffer[0] << "-> temp[0]= " << temp[0] << endl;
								  buffer.resize(realDestSize);
								  std::memcpy(buffer.data(), temp.data(), realDestSize);
							  }
							  else throw ITLException("BytesBytesCodec: " + toString(codec->name) + " not yet implemented: ");

						  }
						  assert(codec->type == ZarrCodecType::ArrayBytesCodec);
						  if (codec->name == ZarrCodecName::Bytes)
						  {
							  cout << "BytesCodec" << endl;
							  readBytesCodec(imgChunk, buffer);
						  }
						  else throw ITLException("ArrayBytesCodec: " + toString(codec->name) + " not yet implemented: ");
						  ++codec;
						  for (; codec != codecs.rend(); ++codec)
						  {
							  assert(codec->type == ZarrCodecType::ArrayArrayCodec);
							  if (codec->name == ZarrCodecName::Transpose)
							  {
								  if (!imgChunk.isInImage(currentChunkShape - Vec3c(1, 1, 1)))
								  {
									  throw ITLException("!imgChunk.isInImage " + toString(currentChunkShape - Vec3c(1, 1, 1)));
								  }
								  Vec3c order = codec->transposeOrder();
								  //todo: currentChunkShape not needed because equal to imgChunk.dims() ?
								  currentChunkShape = currentChunkShape.transposed(order.inverseOrder());
								  transpose(imgChunk, order.inverseOrder(), fillValue);
								  if (!imgChunk.sizeEquals(currentChunkShape))
								  {
									  throw ITLException("!imgChunk.sizeEquals" + toString(currentChunkShape));
								  }
							  }
							  else throw ITLException("ArrayArrayCodec: " + toString(codec->name) + " not yet implemented: ");
						  }
						  //write all pixels of chunk back to img
						  forAllPixels(imgChunk, [&](coord_t x, coord_t y, coord_t z)
						  {
							Vec3c pos = Vec3c(x, y, z) + chunkStartInTarget;
							//cout << "set img(" << pos << ") = imgChunk(" << Vec3c(x, y, z) << ")=" << imgChunk(x, y, z) << endl;
							img(pos) = imgChunk(x, y, z);
						  });
					  }
					  else
					  {
						  cout << "no file: " << filename << endl;
						  // No file => all pixels in the block are fillValue.
						  draw<pixel_t>(img, AABoxc::fromPosSize(chunkStartInTarget, readSize), (pixel_t)fillValue);
					  }
				  }
				});
				if (is_blosc_init)
				{
					blosc_destroy();
				}
			}

			/**
			Checks that provided zarr information is correct and writes metadata.
			@param deleteOldData Contents of the dataset are usually deleted when a write process is started. Set this to false to delete only if the path contains an incompatible dataset, but keep the dataset if it seems to be the same than the current image.
			*/
			void handleExisting(const Vec3c& imageDimensions,
				ImageDataType imageDataType,
				const std::string& path,
				const Vec3c& chunkSize,
				int fillValue,
				const std::list<ZarrCodec>& codecs,
				bool deleteOldData);

			template<typename pixel_t>
			void write(const Image<pixel_t>& img, const std::string& path, const Vec3c& chunkSize, int fillValue, std::list<ZarrCodec>& codecs, bool deleteOldData, bool showProgressInfo)
			{
				Vec3c dimensions = img.dimensions();
				internals::handleExisting(dimensions, img.dataType(), path, chunkSize, fillValue, codecs, deleteOldData);
				writeBlock(img, path, chunkSize, fillValue, codecs, Vec3c(0, 0, 0), dimensions, Vec3c(0, 0, 0), dimensions, showProgressInfo);
			}
		}

		bool getInfo(const std::string& path, Vec3c& shape, bool& isNativeByteOrder, ImageDataType& dataType, Vec3c& chunkSize, std::list<ZarrCodec>& codecs, int& fillValue, std::string& reason);

		inline bool getInfo(const std::string& path, Vec3c& shape, ImageDataType& dataType, std::string& reason)
		{
			bool dummyIsNative;
			Vec3c dummyChunkSize;
			int fillValue;
			std::list<ZarrCodec> codecs;
			return getInfo(path, shape, dummyIsNative, dataType, dummyChunkSize, codecs, fillValue, reason);
		}

		inline const Vec3c DEFAULT_CHUNK_SIZE = Vec3c(1536, 1536, 1536);
		inline const std::list<ZarrCodec> DEFAULT_CODECS = { ZarrCodec(ZarrCodecName::Bytes) };
		inline const nlohmann::json DEFAULT_CODECS_JSON = { ZarrCodec(ZarrCodecName::Bytes).toJSON() };
		/**
		Write an image to a zarr dataset.
		@param targetImg Image to write.
		@param path Name of the top directory of the nn5 dataset.
		*/
		template<typename pixel_t>
		void write(const Image<pixel_t>& img,
			const std::string& path,
			const Vec3c& chunkSize = DEFAULT_CHUNK_SIZE,
			std::list<ZarrCodec> codecs = DEFAULT_CODECS,
			int fillValue = 0,
			bool showProgressInfo = false)
		{
			Vec3c clampedChunkSize = chunkSize;
			clamp(clampedChunkSize, Vec3c(1, 1, 1), img.dimensions());
			internals::write(img, path, clampedChunkSize, fillValue, codecs, false, showProgressInfo);
		}

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
		void writeBlock(const Image<pixel_t>& img, const std::string& path, const Vec3c& chunkSize, int fillValue, std::list<ZarrCodec>& codecs,
			const Vec3c& filePosition, const Vec3c& fileDimensions,
			const Vec3c& imagePosition,
			const Vec3c& blockDimensions,
			bool showProgressInfo = false)
		{
			fs::create_directories(path);
			internals::writeMetadata(path, fileDimensions, img.dataType(), chunkSize, fillValue, codecs);
			internals::writeChunksInRange(img, path, chunkSize, fillValue, codecs, filePosition, fileDimensions, imagePosition, blockDimensions, showProgressInfo);
		}



		/**
		Reads a part of a .zarr dataset to the given image.
		NOTE: Does not support out of bounds start position.
		@param targetImg Image where the data is placed. The size of the image defines the size of the block that is read.
		@param filename The name of the dataset to read.
		@param fileStart Start location of the read in the file. The size of the image defines the size of the block that is read.
		*/
		//mark 5,6
		template<typename pixel_t>
		void readBlock(Image<pixel_t>& img, const std::string& path, const Vec3c& fileStart, bool showProgressInfo = false)
		{
			bool isNativeByteOrder;
			Vec3c datasetShape;
			ImageDataType dataType;
			Vec3c chunkSize;
			int fillValue;
			std::list<ZarrCodec> codecs;
			string reason;
			std::cout << "readBlock" << std::endl;
			if (!itl2::zarr::getInfo(path, datasetShape, isNativeByteOrder, dataType, chunkSize, codecs, fillValue, reason))
				throw ITLException(string("Unable to read zarr dataset: ") + reason);

			if (dataType != img.dataType())
				throw ITLException(string("Expected data type is ") + toString(img.dataType()) + " but the zarr dataset contains data of type " + toString(dataType) + ".");

			if (fileStart.x < 0 || fileStart.y < 0 || fileStart.z < 0 || fileStart.x >= datasetShape.x || fileStart.y >= datasetShape.y || fileStart.z >= datasetShape.z)
				throw ITLException("Out of bounds start position in readBlock.");

			Vec3c cStart = fileStart;
			clamp(cStart, Vec3c(0, 0, 0), datasetShape);
			Vec3c cEnd = fileStart + img.dimensions();
			clamp(cEnd, Vec3c(0, 0, 0), datasetShape);

			internals::readChunksInRange(img, path, chunkSize, fillValue, codecs, datasetShape, cStart, cEnd, showProgressInfo);

			if (!isNativeByteOrder)
				swapByteOrder(img);
		}
		/**
		Reads a zarr dataset file to the given image.
		@param targetImg Image where the data is placed. The size of the image will be set based on the dataset contents.
		@param path Path to the root of the nn5 dataset.
		*/
		//mark 5
		template<typename pixel_t>
		void read(Image<pixel_t>& img, const std::string& path, bool showProgressInfo = false)
		{
			bool isNativeByteOrder;
			Vec3c dimensions;
			ImageDataType dataType;
			Vec3c chunkSize;
			int fillValue;
			std::list<ZarrCodec> codecs;
			string reason;
			if (!itl2::zarr::getInfo(path, dimensions, isNativeByteOrder, dataType, chunkSize, codecs, fillValue, reason))
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
		}
	}
}