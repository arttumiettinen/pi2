
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

			void check(const Vec3c& chunkSize)
			{
				if (chunkSize.min() <= 0)
					throw ITLException(string("NN5 chunk size must be positive, but it is ") + toString(chunkSize));
			}

			void beginWrite(const Vec3c& imageDimensions, ImageDataType imageDataType, const std::string& path, const Vec3c& chunkSize, NN5Compression compression, bool deleteOldData)
			{
				check(chunkSize);

				// Delete old dataset if it exists.
				if (fs::exists(path))
				{
					bool oldIsNativeByteOrder;
					Vec3c oldDimensions;
					ImageDataType oldDataType;
					Vec3c oldChunkSize;
					NN5Compression oldCompression;
					string dummyReason;
					if (!nn5::getInfo(path, oldDimensions, oldIsNativeByteOrder, oldDataType, oldChunkSize, oldCompression, dummyReason))
					{
						// The path does not contain an NN5 dataset.
						// If it is no know image, do not delete it.
						if (!io::getInfo(path, oldDimensions, oldDataType, dummyReason))
							throw ITLException(string("Unable to write an NN5 as the output folder already exists but cannot be verified to be an image: ") + path + " Consider removing the existing dataset manually.");

						// Here the path does not contain NN5 but contains an image of known type.
						// Delete the old image.
						fs::remove_all(path);
					}
					else if (oldDimensions == imageDimensions &&
						oldIsNativeByteOrder == true &&
						oldDataType == imageDataType &&
						oldChunkSize == chunkSize &&
						oldCompression == compression)
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
				internals::writeMetadata(path, imageDimensions, imageDataType, chunkSize, compression);
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
		bool isChunkSafe(const AABoxc& box, const std::vector<NN5Process>& processes)
		{
			vector<size_t> readerIndices;
			vector<size_t> writerIndices;
			for (size_t n = 0; n < processes.size(); n++)
			{
				const auto& process = processes[n];
				if (process.readBlock.overlaps(box))
					readerIndices.push_back(n);
				if (process.writeBlock.overlaps(box))
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

		size_t startConcurrentWrite(const Vec3c& imageDimensions, ImageDataType imageDataType, const std::string& path, const Vec3c& chunkSize, NN5Compression compression, const std::vector<NN5Process>& processes)
		{
			// Find chunks that are
			// * written to by separate processes, or
			// * read from and written to by at least two separate processes,
			// and tag those unsafe by creating writes folder into the chunk folder.

			internals::beginWrite(imageDimensions, imageDataType, path, chunkSize, compression, false);

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
			NN5Compression compression;
			string reason;
			if (!getInfo(path, imageDimensions, isNativeByteOrder, dataType, chunkSize, compression, reason))
				throw ITLException(string("Unable to read nn5 dataset: ") + reason);

			return needsEndConcurrentWrite(path, getDimensionality(imageDimensions), chunkIndex);
		}

		vector<Vec3c> getChunksThatNeedEndConcurrentWrite(const std::string& path)
		{
			bool isNativeByteOrder;
			Vec3c imageDimensions;
			ImageDataType dataType;
			Vec3c chunkSize;
			NN5Compression compression;
			string reason;
			if (!getInfo(path, imageDimensions, isNativeByteOrder, dataType, chunkSize, compression, reason))
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

			template<typename pixel_t> void readAndAdd(Image<pixel_t>& img, const string& filename, NN5Compression compression)
			{
				Vec3c blockPos = parsePosition(filename);

				Image<pixel_t> block;
				readChunkFile(block, filename, compression);

				copyValues(img, block, blockPos);
			}

			template<typename pixel_t> struct CombineChunkWrites
			{
			public:
				static void run(const string& path, const Vec3c& datasetSize, const Vec3c& chunkSize, NN5Compression compression, const Vec3c& chunkIndex)
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
								readChunkFile(img, filename, compression);
							}
							else
							{
								throw ITLException(string("Multiple image files found in block directory ") + chunkFolder + " while combining chunk writes.");
							}

							// Modify data with the new writes.
							for (const string& file : writesFiles)
							{
								readAndAdd(img, file, compression);
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

			void endConcurrentWrite(const std::string& path, const Vec3c& imageDimensions, ImageDataType dataType, const Vec3c& chunkSize, NN5Compression compression, const Vec3c& chunkIndex)
			{
				pick<internals::CombineChunkWrites>(dataType, path, imageDimensions, chunkSize, compression, chunkIndex);
			}
		}

		

		void endConcurrentWrite(const std::string& path, const Vec3c& chunkIndex)
		{
			bool isNativeByteOrder;
			Vec3c imageDimensions;
			ImageDataType dataType;
			Vec3c chunkSize;
			NN5Compression compression;
			string reason;
			if (!getInfo(path, imageDimensions, isNativeByteOrder, dataType, chunkSize, compression, reason))
				throw ITLException(string("Unable to read nn5 dataset: ") + reason);

			internals::endConcurrentWrite(path, imageDimensions, dataType, chunkSize, compression, chunkIndex);
		}

		void endConcurrentWrite(const std::string& path, bool showProgressInfo)
		{
			bool isNativeByteOrder;
			Vec3c fileDimensions;
			ImageDataType dataType;
			Vec3c chunkSize;
			NN5Compression compression;
			string reason;
			if (!getInfo(path, fileDimensions, isNativeByteOrder, dataType, chunkSize, compression, reason))
				throw ITLException(string("Unable to read nn5 dataset: ") + reason);
			size_t dimensionality = getDimensionality(fileDimensions);

			internals::forAllChunks(fileDimensions, chunkSize, showProgressInfo, [&](const Vec3c& chunkIndex, const Vec3c& chunkStart)
				{
					if(needsEndConcurrentWrite(path, dimensionality, chunkIndex))
						internals::endConcurrentWrite(path, fileDimensions, dataType, chunkSize, compression, chunkIndex);
				});
			
			// Remove concurrent tag file after all blocks are processed such that if exception is thrown during processing,
			// the endConcurrentWrite can continue simply by re-running it.
			fs::remove_all(internals::concurrentTagFile(path));
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
				nn5::write(img, "./nn5_testimage", chunkSize, compression, true);

				// Read
				Image<uint16_t> read;
				nn5::read(read, "./nn5_testimage", true);

				raw::writed(read, "./nn5_results/read_from_disk");

				testAssert(equals(img, read), string("nn5 read image compared to written image, chunk size = ") + toString(chunkSize));
			}

			void nn5ioSimpleTest(const Vec3c& dimensions)
			{
				Image<uint16_t> img(dimensions);
				ramp3(img);
				add(img, 10);
				nn5::write(img, "./nn5_results/simpletest");

				Image<uint16_t> fromDisk;
				nn5::read(fromDisk, "./nn5_results/simpletest");

				testAssert(equals(img, fromDisk), string("NN5 simple test, dimensions = ") + toString(dimensions));
			}

			void nn5io()
			{
				// 0-dimensional image
				nn5ioSimpleTest(Vec3c(1, 1, 1));
				nn5ioSimpleTest(Vec3c(2, 1, 1));
				nn5ioSimpleTest(Vec3c(2, 2, 1));
				nn5ioSimpleTest(Vec3c(2, 2, 2));

				nn5ioSimpleTest(Vec3c(10, 1, 1));
				nn5ioSimpleTest(Vec3c(20, 1, 1));
				nn5ioSimpleTest(Vec3c(20, 20, 1));
				nn5ioSimpleTest(Vec3c(20, 20, 20));

				nn5ioSimpleTest(Vec3c(10, 10, 10));
				nn5ioSimpleTest(Vec3c(20, 10, 10));
				nn5ioSimpleTest(Vec3c(20, 20, 10));
				nn5ioSimpleTest(Vec3c(20, 20, 20));

				// Different chunk sizes
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
			}

			void nn5BlockIoOneTest(NN5Compression compression, const Vec3c& chunkSize,
				const Vec3c& blockStart,
				const Vec3c& blockSize)
			{
				cout << "Chunk size = " << chunkSize << endl;
				cout << "Compression = " << toString(compression) << endl;
				cout << "Block start = " << blockStart << endl;
				cout << "Block size = " << blockSize << endl;

				Vec3c dimensions(100, 200, 300);

				Image<uint16_t> img(dimensions);
				ramp3(img);

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
			
			void nn5BlockIoOneTest(const Vec3c& chunkSize,
				const Vec3c& blockStart,
				const Vec3c& blockSize)
			{
				nn5BlockIoOneTest(NN5Compression::Raw, chunkSize, blockStart, blockSize);
				nn5BlockIoOneTest(NN5Compression::LZ4, chunkSize, blockStart, blockSize);
			}

			void nn5BlockIoOneTest(
				const Vec3c& blockStart,
				const Vec3c& blockSize)
			{
				nn5BlockIoOneTest(Vec3c(50, 102, 99), blockStart, blockSize);
				nn5BlockIoOneTest(Vec3c(50, 40, 65), blockStart, blockSize);
				nn5BlockIoOneTest(Vec3c(40, 30, 20), blockStart, blockSize);
			}

			void nn5BlockIo()
			{
				//// Some specific edge cases
				//if(true)
				//{
				//	Image<uint16_t> img(256, 256, 129);
				//	ramp3(img);

				//	nn5::write(img, "./nn5_block_tests1/full_image", Vec3c(30, 32, 33), NN5Compression::LZ4);


				//	Image<uint16_t> block(256, 256, 64);
				//	nn5::readBlock(block, "./nn5_block_tests1/full_image", Vec3c(0, 0, 65));

				//	raw::writed(block, "./nn5_block_tests1/read_block_written_as_full_raw_file");
				//	raw::writeBlock(block, concatDimensions("./nn5_block_tests1/read_block_written_as_raw_block", img.dimensions()),
				//		Vec3c(0, 0, 65), img.dimensions(), Vec3c(0, 0, 0), Vec3c(256, 256, 64));
				//	//writerawblock("image_RJBTFVPJJR_1933", "../../testing/pi2py2/head_rotate_img_result_0.5233333333333333_[1, 0, 0]_[128, 128, 64]_[128, 128, 64]_distributed_256x256x129.raw", "[0, 0, 65]", "[256, 256, 129]", "[0, 0, 0]", "[256, 256, 64]")
				//}

				if (true)
				{
					nn5BlockIoOneTest(Vec3c(0, 0, 0), Vec3c(1, 200, 300));
					nn5BlockIoOneTest(Vec3c(0, 0, 0), Vec3c(100, 1, 300));
					nn5BlockIoOneTest(Vec3c(0, 0, 0), Vec3c(100, 200, 1));

					nn5BlockIoOneTest(Vec3c(99, 0, 0), Vec3c(1, 200, 300));
					nn5BlockIoOneTest(Vec3c(0, 199, 0), Vec3c(100, 1, 300));
					nn5BlockIoOneTest(Vec3c(0, 0, 299), Vec3c(100, 200, 1));

					nn5BlockIoOneTest(Vec3c(0, 0, 0), Vec3c(50, 102, 99));
					nn5BlockIoOneTest(Vec3c(1, 2, 3), Vec3c(50, 102, 99));
				}

				if(true)
				{
					Image<uint16_t> img(256, 256, 129);
					ramp3(img);

					//vector<NN5Process> processes;
					coord_t z = 0;
					while (z < img.depth())
					{
						coord_t y = 0;
						while (y < img.height())
						{
							//processes.push_back(NN5Process{ AABoxc::fromPosSize(Vec3c(), Vec3c()), AABoxc::fromPosSize(Vec3c(0, y, z), Vec3c(256, 128, 1)) });
							
							coord_t d = 2;
							if (z + d >= img.depth())
								d = 1;

							nn5::writeBlock(img, "./nn5_block_tests2/nn5", Vec3c(30, 32, 33), NN5Compression::LZ4,
								Vec3c(0, y, z), img.dimensions(), Vec3c(0, y, z), Vec3c(256, 128, d));

							y += 128;
						}
						z += 2;
					}
					//nn5::startConcurrentWrite(img, "./nn5_block_tests2/nn5", Vec3c(30, 32, 33), NN5Compression::LZ4, processes);
					//nn5::writeBlock(img, "./nn5_block_tests2/nn5", Vec3c(30, 32, 33), NN5Compression::LZ4,
					//	Vec3c(0,   0, 128), img.dimensions(), Vec3c(0, 0, 0), Vec3c(256, 128, 1));
					//nn5::writeBlock(img, "./nn5_block_tests2/nn5", Vec3c(30, 32, 33), NN5Compression::LZ4,
					//	Vec3c(0, 128, 128), img.dimensions(), Vec3c(0, 0, 0), Vec3c(256, 128, 1));
					//nn5::endConcurrentWrite("./nn5_block_tests2/nn5");



					Image<uint16_t> imgRaw, imgNN5;
					//raw::read(imgRaw, "./nn5_block_tests2/raw");
					nn5::read(imgNN5, "./nn5_block_tests2/nn5");
					raw::writed(imgNN5, "./nn5_block_tests2/nn5_to_raw");
					testAssert(equals(img, imgNN5), "orig vs NN5 writeBlock");
				}

			}


			void concurrencyOneTest(NN5Compression compression, const Vec3c& chunkSize)
			{
				cout << "Chunk size = " << chunkSize << ", compression = " << toString(compression) << endl;

				Vec3c dimensions(100, 200, 300);

				Image<uint16_t> img(dimensions);
				ramp3(img);

				Vec3c blockStart(10, 20, 30);
				Vec3c blockSize(50, 60, 70);

				string entireImageFile = "./nn5_concurrency/entire_image";

				// Write entire image and read.
				fs::remove_all(entireImageFile);
				{
					vector<nn5::NN5Process> processes;

					Vec3c processBlockSize(30, 30, 30);
					nn5::internals::forAllChunks(dimensions, processBlockSize, true, [&](const Vec3c& processBlockIndex, const Vec3c& processBlockStart)
						{
							processes.push_back(NN5Process{ AABoxc::fromPosSize(processBlockStart, processBlockSize + Vec3c(10, 10, 10)), AABoxc::fromPosSize(processBlockStart, processBlockSize) });
						});

					nn5::startConcurrentWrite(img, entireImageFile, chunkSize, compression, processes);
					nn5::write(img, entireImageFile, chunkSize, compression);
					nn5::endConcurrentWrite(entireImageFile);

					Image<uint16_t> entireFromDisk;
					nn5::read(entireFromDisk, entireImageFile);
					testAssert(equals(entireFromDisk, img), "NN5 entire image read/write cycle with concurrency enabled");
				}

				// Write in 2 blocks and read.
				fs::remove_all(entireImageFile);
				{
					Vec3c block1Start(0, 0, 0);
					Vec3c block1Size(dimensions.x / 2 + 3, dimensions.y, dimensions.z);
					Vec3c block2Start(block1Size.x, 0, 0);
					Vec3c block2Size(dimensions.x - block1Size.x, dimensions.y, dimensions.z);

					vector<nn5::NN5Process> processes;
					processes.push_back(NN5Process{ AABoxc::fromPosSize(block1Start, block1Size + Vec3c(10, 0, 0)), AABoxc::fromPosSize(block1Start, block1Size) });
					processes.push_back(NN5Process{ AABoxc::fromPosSize(block2Start, block2Size), AABoxc::fromPosSize(block2Start, block2Size) });

					nn5::startConcurrentWrite(img, entireImageFile, chunkSize, compression, processes);
					nn5::writeBlock(img, entireImageFile, chunkSize, compression, block1Start, dimensions, block1Start, block1Size);
					nn5::writeBlock(img, entireImageFile, chunkSize, compression, block2Start, dimensions, block2Start, block2Size);
					nn5::endConcurrentWrite(entireImageFile);

					Image<uint16_t> entireFromDisk;
					nn5::read(entireFromDisk, entireImageFile);
					testAssert(equals(entireFromDisk, img), "NN5 entire image read/write cycle in 2 pseudo-concurrent blocks");
				}

				// Write in multiple blocks and read.
				{
					vector<nn5::NN5Process> processes;

					Vec3c processBlockSize(30, 31, 32);
					nn5::internals::forAllChunks(dimensions, processBlockSize, true, [&](const Vec3c& processBlockIndex, const Vec3c& processBlockStart)
						{
							processes.push_back(NN5Process{ AABoxc::fromPosSize(processBlockStart, processBlockSize + Vec3c(10, 10, 10)), AABoxc::fromPosSize(processBlockStart, processBlockSize) });
						});

					nn5::startConcurrentWrite(img, entireImageFile, chunkSize, compression, processes);
					nn5::internals::forAllChunks(dimensions, processBlockSize, true, [&](const Vec3c& processBlockIndex, const Vec3c& processBlockStart)
						{
							nn5::writeBlock(img, entireImageFile, chunkSize, compression, processBlockStart, dimensions, processBlockStart, processBlockSize);
						});
					nn5::endConcurrentWrite(entireImageFile);

					Image<uint16_t> entireFromDisk;
					nn5::read(entireFromDisk, entireImageFile);
					testAssert(equals(entireFromDisk, img), "NN5 entire image read/write cycle with concurrency enabled");
				}
			}

			void concurrency()
			{
				concurrencyOneTest(NN5Compression::Raw, Vec3c(40, 200, 300));
				concurrencyOneTest(NN5Compression::LZ4, Vec3c(40, 200, 300));

				concurrencyOneTest(NN5Compression::Raw, Vec3c(40, 30, 20));
				concurrencyOneTest(NN5Compression::LZ4, Vec3c(40, 30, 20));
			}

			void concurrencyLong()
			{
				coord_t w = 100;
				coord_t h = 200;
				coord_t d = 300;

				srand(1212);

				for (coord_t n = 0; n < 100; n++)
				{
					Vec3c chunkSize(randc(10, w + 10), randc(10, h + 10), randc(10, d + 10));

					concurrencyOneTest(NN5Compression::Raw, chunkSize);
					concurrencyOneTest(NN5Compression::LZ4, chunkSize);
				}
			}
		}
	}
}