
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
		namespace internals
		{
			std::vector<char> readBytesOfFile(const std::string& filename)
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

		size_t startConcurrentWrite(const Vec3c& imageDimensions,
			ImageDataType imageDataType,
			const std::string& path,
			const Vec3c& chunkSize,
			const std::vector<io::DistributedImageProcess>& processes,
			const codecs::Pipeline& codecs)
		{
			ZarrMetadata metadata = DEFAULT_METADATA;
			metadata.chunkSize = chunkSize;
			metadata.dataType = imageDataType;
			metadata.codecs = codecs;
			internals::handleExisting(imageDimensions, path, metadata, false, true);
			fs::create_directories(path);
			internals::writeMetadata(path, imageDimensions, metadata);

			// Find chunks that are
			// * written to by separate processes, or
			// * read from and written to by at least two separate processes,
			// and tag those unsafe by creating writes folder into the chunk folder.

			// Tag the image as concurrently processed
			ofstream out(internals::concurrentTagFile(path), ios_base::out | ios_base::trunc | ios_base::binary);

			size_t unsafeChunkCount = 0;
			forAllChunks(imageDimensions, metadata.chunkSize, false, [&](const Vec3c& chunkIndex, const Vec3c& chunkStart)
				{
					string chunkFile = internals::chunkFile(path, getDimensionality(imageDimensions), chunkIndex, metadata.separator);
					//TODO: will that happen later? fs::create_directories(chunkFile);

					string writesFolder = internals::writesFolder(chunkFile);

					AABoxc chunkBox = AABoxc::fromPosSize(chunkStart, metadata.chunkSize);

					if (!io::isChunkSafe(chunkBox, processes))
					{
						// Mark the chunk as unsafe by creating writes folder.
						fs::create_directories(writesFolder);
						unsafeChunkCount++;
					}
				});
			return unsafeChunkCount;
		}


		template<typename pixel_t> struct CombineChunkWrites
		{
		public:
			static void run(const string& path, const Vec3c& datasetSize, const ZarrMetadata& metadata, const Vec3c& chunkIndex, const Vec3c& chunkStart)
			{
				string chunkFile = internals::chunkFile(path, getDimensionality(datasetSize), chunkIndex, metadata.separator);
				string writesFolder = internals::writesFolder(chunkFile);

				if (fs::exists(writesFolder))
				{
					// Find all files in the writes folder
					vector<string> sortedWritesFiles = buildFileList(writesFolder + "/");

					if (!sortedWritesFiles.empty())
					{
						Image<pixel_t> imgChunk(metadata.chunkSize);
						AABoxc chunkRegion = AABoxc::fromPosSize(chunkStart, metadata.chunkSize);
						internals::readChunkInBlock(imgChunk, chunkFile, chunkRegion, chunkStart,chunkStart, metadata);

						// Modify data with the new writes.
						for (const string& writesFile : sortedWritesFiles)
						{
							AABoxc updateRegion = internals::updateRegionOfWritesFile(writesFile);
							internals::readChunkInBlock(imgChunk, writesFile, updateRegion, chunkStart,chunkStart, metadata);
						}
						// Write back to disk.
						internals::writeSingleChunk(imgChunk, chunkFile, chunkRegion, metadata, true);
					}
					else cout << "writes folder empty: " << writesFolder << endl;

					// Remove the writes folder in order to mark this chunk processed.
					fs::remove_all(writesFolder);
				}
			}
		};

		void endConcurrentWrite(const std::string& path, bool showProgressInfo)
		{
			ZarrMetadata metadata;
			std::string reason;
			Vec3c dimensions;
			if(!internals::getInfo(path, dimensions, metadata, reason)){
				throw ITLException("Unable to read zarr dataset. path: " + path + " reason: " + reason);
			}

			forAllChunks(dimensions, metadata.chunkSize, showProgressInfo, [&](const Vec3c& chunkIndex, const Vec3c& chunkStart)
				{
				  pick<CombineChunkWrites>(metadata.dataType, path, dimensions, metadata, chunkIndex, chunkStart);
				});

			// Remove concurrent tag file after all blocks are processed such that if exception is thrown during processing,
			// the endConcurrentWrite can continue simply by re-running it.
			fs::remove_all(internals::concurrentTagFile(path));
		}

		namespace tests
		{
			void read()
			{
				Image<int32_t> fromDisk;
				zarr::read(fromDisk, "./testdata/zarrita.zarr");
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


			void writeBlockTest(Vec3c chunkSize, uint16_t changeValue)
			{
				string path = "./testoutput/writeBlock.zarr";
				Vec3c size = Vec3c(2, 3, 4);
				Vec3c startBlock(1, 2, 3);
				Vec3c blockSize(2, 2, 2);

				Image<uint16_t> img(size, 42);
				zarr::write(img, path, chunkSize);
				drawAll(img, changeValue);
				zarr::writeBlock(img, path, startBlock, blockSize, chunkSize);

				Image<uint16_t> fromDisk;
				zarr::read(fromDisk, path);

				Image<uint16_t> expected(size, 42);
				draw(expected, AABoxsc::fromPosSize(Vec3<int>(startBlock), Vec3<int>(blockSize)), changeValue);

				bool imgEquals = equals(expected, fromDisk);
				testAssert(imgEquals, string("zarr test writeBlock chunkSize=" + toString(chunkSize)));

				if(!imgEquals)
				{
					internals::printImg(expected);
					internals::printImg(fromDisk);
				}
			}

			void writeBlock()
			{
				writeBlockTest(Vec3c(1, 1, 1), 1);
				writeBlockTest(Vec3c(2, 2, 2), 1);
				writeBlockTest(DEFAULT_CHUNK_SIZE, 1);
				writeBlockTest(Vec3c(1, 1, 1), DEFAULT_FILLVALUE);
				writeBlockTest(Vec3c(2, 2, 2), DEFAULT_FILLVALUE);
				writeBlockTest(DEFAULT_CHUNK_SIZE, DEFAULT_FILLVALUE);
			}

			void readBlockTest(Vec3c chunkSize)
			{
				string path = "./testoutput/readBlock.zarr";
				Vec3c size = Vec3c(10, 10, 10);
				Vec3c startOnes(3, 4, 5);

				Vec3c startBlock(2, 2, 2);
				Vec3c blockSize(2, 4, 6);

				Image<uint16_t> img(size, 0);
				draw(img, AABoxsc::fromMinMax(Vec3<int>(startOnes), Vec3<int>(size)), (uint16_t)1);
				zarr::write(img, path, chunkSize);

				Image<uint16_t> fromDisk(blockSize, 0);
				zarr::readBlock(fromDisk, path, startBlock);

				Image<uint16_t> expected(blockSize, 0);
				draw(expected, AABoxsc::fromMinMax(Vec3<int>(startOnes - startBlock), Vec3<int>(blockSize)), (uint16_t)1);

				testAssert(equals(expected, fromDisk), string("zarr test readBlock chunkSize=" + toString(chunkSize)));
			}

			void readBlock()
			{
				readBlockTest(Vec3c(1, 1, 1));
				readBlockTest(Vec3c(2, 2, 2));
				readBlockTest(DEFAULT_CHUNK_SIZE);
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
			void zarrMetadataEquals(){
				ZarrMetadata metadata = {ImageDataType::Int32, Vec3c(1,1,1), { codecs::ZarrCodec(codecs::Name::Bytes) }, 0, "/"};
				ZarrMetadata equalMetadata = {ImageDataType::Int32, Vec3c(1,1,1), { codecs::ZarrCodec(codecs::Name::Bytes) }, 0, "/"};
				ZarrMetadata differentMetadata = {ImageDataType::Int32, Vec3c(1,1,1), { codecs::ZarrCodec(codecs::Name::Blosc) }, 0, "/"};
				testAssert(metadata == equalMetadata, string("zarr test equal metadata equals"));
				testAssert(!(metadata == differentMetadata), string("zarr test different metadata doesnt equal"));
			}

			void separatorTest(std::string separator, std::string label)
			{
				string path = "./testoutput/test_separator_" + label + ".zarr";

				Image<uint16_t> img(Vec3c(2, 3, 4));
				ramp3(img);
				add(img, 10);
				zarr::write(img,
					path,
					Vec3c(1, 1, 1),
					zarr::DEFAULT_CODECS,
					separator
				);
				Image<uint16_t> fromDisk;
				zarr::read(fromDisk, path);

				std::string filename = path + "/c" + separator + "1" + separator + "2" + separator + "3";
				testAssert(fs::exists(filename), string("zarr test read and write with separator " + separator + " file does not exist: " + filename));

				testAssert(equals(img, fromDisk), string("zarr test read and write with separator " + separator + " read failed"));
			}
			void separator(){
				separatorTest(".", "dot");
				separatorTest("/", "slash");
				separatorTest("-", "minus");
			}

			void shardingTest(std::string indexLocation, bool withBlosc, Vec3c imgShape = Vec3c(2, 5, 10), Vec3c shardShape = Vec3c(2, 5, 1), Vec3c chunkShape = Vec3c(2, 1, 1)){

				string path = "./testoutput/test_sharding_" + indexLocation;
				if(withBlosc) path+="_withBlosc";
				path+=".zarr";

				Image<uint16_t> img(imgShape);
				ramp(img, 0);
				add(img, 10);
				string bloscCodecConfig = R"({"cname": "lz4", "clevel": 1, "shuffle": "shuffle", "typesize": 4, "blocksize": 0})";

				nlohmann::json codecs = { codecs::ZarrCodec(codecs::Name::Bytes).toJSON()};
				if (withBlosc) codecs = { codecs::ZarrCodec(codecs::Name::Bytes).toJSON(), codecs::ZarrCodec(codecs::Name::Blosc, nlohmann::json::parse(bloscCodecConfig)).toJSON()};

				nlohmann::json shardingCodecConfigJSON = {
					{"chunk_shape", { chunkShape.x, chunkShape.y, chunkShape.z }},
					{"codecs", codecs},
					//{"index_codecs", { codecs::ZarrCodec(codecs::Name::Bytes).toJSON(), codecs::ZarrCodec(codecs::Name::Blosc, nlohmann::json::parse(bloscCodecConfig)).toJSON()}},
					{"index_codecs", { codecs::ZarrCodec(codecs::Name::Bytes).toJSON()}},
					{"index_location", indexLocation}
				};

				zarr::write(img,
					path,
					shardShape,
					{ codecs::ZarrCodec(codecs::Name::Sharding,  shardingCodecConfigJSON)});

				Image<uint16_t> fromDisk;
				zarr::read(fromDisk, path);

				testAssert(equals(img, fromDisk), string("zarr test write sharding with indexLocation=" + indexLocation + " withBlosc="+ toString(withBlosc)));
			}

			void shardingEmptyInnerChunksTest(string indexLocation)
			{
				Vec3c imgShape = Vec3c(10, 10,10);
				Vec3c shardShape = Vec3c(10, 10, 10);
				Vec3c chunkShape = Vec3c(5, 5, 5);

				string path = "./testoutput/test_sharding_empty_inner_chunks_" + indexLocation;
				Image<uint16_t> img(imgShape);
				add(img, DEFAULT_FILLVALUE);
				nlohmann::json shardingCodecConfigJSON = {
					{ "chunk_shape", { chunkShape.x, chunkShape.y, chunkShape.z }},
					{ "codecs", { codecs::ZarrCodec(codecs::Name::Bytes).toJSON() }},
					{ "index_codecs", { codecs::ZarrCodec(codecs::Name::Bytes).toJSON() }},
					{ "index_location", indexLocation }
				};

				img(0, 0, 0) = 1; //to not have complete shard empty
				zarr::write(img,
					path + "_empty",
					shardShape,
					{ codecs::ZarrCodec(codecs::Name::Sharding, shardingCodecConfigJSON) });

				ramp(img, 0);
				zarr::write(img,
					path + "_full",
					shardShape,
					{ codecs::ZarrCodec(codecs::Name::Sharding, shardingCodecConfigJSON) });

				auto sizeEmptyChunks = fileSize(internals::chunkFile(path + "_empty", 3, Vec3c(0, 0, 0), DEFAULT_SEPARATOR));
				auto sizeFullChunks = fileSize(internals::chunkFile(path + "_full", 3, Vec3c(0, 0, 0), DEFAULT_SEPARATOR));

				testAssert(sizeEmptyChunks < sizeFullChunks, "shardingEmptyInnerChunksTest" + indexLocation);
			}
			void sharding()
			{
				shardingTest("end", false);
				shardingTest("start", false);
				shardingTest("end", true);
				shardingTest("start", true);
				shardingTest("start", true, Vec3c(10, 10, 10), DEFAULT_CHUNK_SIZE, Vec3c(8,8,8));
				shardingEmptyInnerChunksTest("start");
				shardingEmptyInnerChunksTest("end");
			}
			void emptyChunks(){
				string path = "./testoutput/test_empty_chunks";
				Image<uint16_t> img(Vec3c(8, 8, 8));
				add(img, DEFAULT_FILLVALUE);

				bool expectedFiles[4][4][4] = {false};

				img(0, 0, 0) = 1;
				expectedFiles[0][0][0] = true;

				img(7, 7, 7) = 3;
				expectedFiles[3][3][3] = true;

				img(1, 2, 3) = 1;
				expectedFiles[0][1][1] = true;

				zarr::write(img,
					path,
					Vec3c(2,2,2));

				for (int i = 0; i < 4; ++i)
				{
					for (int j = 0; j < 4; ++j)
					{
						for (int k = 0; k < 4; ++k)
						{
							string filename = internals::chunkFile(path, 3, Vec3c(i,j,k), DEFAULT_SEPARATOR);
							testAssert(fs::exists(filename) == expectedFiles[i][j][k], "test emptyChunks at " + toString(i) + " " + toString(j) + " " + toString(k));

						}
					}
				}
			}
			void concurrencyOneTest(const codecs::Pipeline& codecs, const Vec3c& chunkSize)
			{
				cout << "Chunk size = " << chunkSize << endl;

				Vec3c dimensions(100, 200, 300);

				Image<uint16_t> img(dimensions);
				ramp3(img);

				string imageFilename = "./testoutput/zarr_concurrency/entire_image";
				ZarrMetadata metadata = { img.dataType(), chunkSize, codecs, DEFAULT_FILLVALUE, DEFAULT_SEPARATOR };

				// Write entire image and read.
				fs::remove_all(imageFilename);
				{
					vector<io::DistributedImageProcess> processes;

					Vec3c processBlockSize(30, 30, 30);
					forAllChunks(dimensions, processBlockSize, true, [&](const Vec3c& processBlockIndex, const Vec3c& processBlockStart)
					{
					  processes.push_back(io::DistributedImageProcess{ AABoxc::fromPosSize(processBlockStart, processBlockSize + Vec3c(4, 4, 4)),
																	   AABoxc::fromPosSize(processBlockStart, processBlockSize) });
					});

					zarr::startConcurrentWrite(img, imageFilename, chunkSize, processes, codecs);
					zarr::write(img, imageFilename, metadata);
					zarr::endConcurrentWrite(imageFilename);

					Image<uint16_t> entireFromDisk;
					zarr::read(entireFromDisk, imageFilename);
					testAssert(equals(entireFromDisk, img), "Zarr entire image read/write cycle with concurrency enabled");
				}
				// Write in 2 blocks and read.
				fs::remove_all(imageFilename);
				{
					Vec3c block1Start(0, 0, 0);
					Vec3c block1Size(dimensions.x / 2 + 3, dimensions.y, dimensions.z);
					Vec3c block2Start(block1Size.x, 0, 0);
					Vec3c block2Size(dimensions.x - block1Size.x, dimensions.y, dimensions.z);

					vector<io::DistributedImageProcess> processes;
					processes.push_back(io::DistributedImageProcess{ AABoxc::fromPosSize(block1Start, block1Size + Vec3c(10, 0, 0)), AABoxc::fromPosSize(block1Start, block1Size) });
					processes.push_back(io::DistributedImageProcess{ AABoxc::fromPosSize(block2Start, block2Size), AABoxc::fromPosSize(block2Start, block2Size) });

					zarr::startConcurrentWrite(img, imageFilename, chunkSize, processes, codecs);
					zarr::writeBlock(img, imageFilename, block1Start, block1Size, metadata);
					zarr::writeBlock(img, imageFilename, block2Start, block2Size, metadata);
					zarr::endConcurrentWrite(imageFilename);

					Image<uint16_t> entireFromDisk;
					zarr::read(entireFromDisk, imageFilename);
					testAssert(equals(entireFromDisk, img), "Zarr entire image read/write cycle in 2 pseudo-concurrent blocks");
				}

				// Write in multiple blocks and read.
				{
					vector<io::DistributedImageProcess> processes;

					Vec3c processBlockSize(30, 31, 32);
					forAllChunks(dimensions, processBlockSize, true, [&](const Vec3c& processBlockIndex, const Vec3c& processBlockStart)
						{
							processes.push_back(io::DistributedImageProcess{ AABoxc::fromPosSize(processBlockStart, processBlockSize + Vec3c(10, 10, 10)), AABoxc::fromPosSize(processBlockStart, processBlockSize) });
						});

					zarr::startConcurrentWrite(img, imageFilename, chunkSize, processes, codecs);
					forAllChunks(dimensions, processBlockSize, true, [&](const Vec3c& processBlockIndex, const Vec3c& processBlockStart)
						{
							zarr::writeBlock(img, imageFilename, processBlockStart, processBlockSize, metadata);
						});
					zarr::endConcurrentWrite(imageFilename);

					Image<uint16_t> entireFromDisk;
					zarr::read(entireFromDisk, imageFilename);
					testAssert(equals(entireFromDisk, img), "Zarr entire image read/write cycle with concurrency enabled");
				}
			}
			void concurrency()
			{
				nlohmann::json shardingCodecConfigJSON = {
					{ "chunk_shape", { 20, 10, 20 }},
					{ "codecs", { codecs::ZarrCodec(codecs::Name::Bytes).toJSON() }},
					{ "index_codecs", { codecs::ZarrCodec(codecs::Name::Bytes).toJSON() }},
					{ "index_location", "end" }
				};
				codecs::Pipeline codecsWithSharding = { codecs::ZarrCodec(codecs::Name::Sharding, shardingCodecConfigJSON) };

				concurrencyOneTest(DEFAULT_CODECS,   Vec3c(40, 200, 300));
				concurrencyOneTest(codecsWithSharding,Vec3c(40, 200, 300));

				concurrencyOneTest(DEFAULT_CODECS,   Vec3c(40, 30, 20));
				concurrencyOneTest(codecsWithSharding,Vec3c(40, 30, 20));
			}
		}
	}
}