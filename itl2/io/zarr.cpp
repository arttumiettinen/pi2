
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

				Image<uint16_t> img(size, 42);
				add(img, 10);
				zarr::writeBlock(img, path, Vec3c(0, 0, 0), size, startBlock, endBlock);

				Image<uint16_t> fromDisk;
				zarr::read(fromDisk, path);

				Image<uint16_t> expected(size, 42);
				draw(expected, AABoxsc::fromMinMax(Vec3<int>(startBlock), Vec3<int>(endBlock)), (uint16_t)52);

				testAssert(equals(img, fromDisk), string("zarr test writeBlock"));
			}

			void readBlockTest(Vec3c chunkSize)
			{
				cout << "readBlockTest chunkSize=" << chunkSize << endl;
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

				testAssert(equals(img, fromDisk), string("zarr test readBlock chunkSize=" + toString(chunkSize)));
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
				ZarrMetadata<uint16_t> metadata = {Vec3c(1,1,1), { codecs::ZarrCodec(codecs::Name::Bytes) }, 0, "/"};
				ZarrMetadata<uint16_t> equalMetadata = {Vec3c(1,1,1), { codecs::ZarrCodec(codecs::Name::Bytes) }, 0, "/"};
				ZarrMetadata<uint16_t> differentMetadata = {Vec3c(1,1,1), { codecs::ZarrCodec(codecs::Name::Blosc) }, 0, "/"};
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

			void shardingTest(std::string indexLocation, bool withBlosc){

				string path = "./testoutput/test_sharding_" + indexLocation;
				if(withBlosc) path+="_withBlosc";
				path+=".zarr";

				Image<uint16_t> img(Vec3c(2, 5, 10));
				ramp(img, 0);
				add(img, 10);
				string bloscCodecConfig = R"({"cname": "lz4", "clevel": 1, "shuffle": "shuffle", "typesize": 4, "blocksize": 0})";

				nlohmann::json codecs = { codecs::ZarrCodec(codecs::Name::Bytes).toJSON()};
				if (withBlosc) codecs = { codecs::ZarrCodec(codecs::Name::Bytes).toJSON(), codecs::ZarrCodec(codecs::Name::Blosc, nlohmann::json::parse(bloscCodecConfig)).toJSON()};

				nlohmann::json shardingCodecConfigJSON = {
					{"chunk_shape", { 2, 5, 1 }},
					{"codecs", codecs},
					//{"index_codecs", { codecs::ZarrCodec(codecs::Name::Bytes).toJSON(), codecs::ZarrCodec(codecs::Name::Blosc, nlohmann::json::parse(bloscCodecConfig)).toJSON()}},
					{"index_codecs", { codecs::ZarrCodec(codecs::Name::Bytes).toJSON()}},
					{"index_location", indexLocation}
				};

				zarr::write(img,
					path,
					zarr::DEFAULT_CHUNK_SIZE,
					{ codecs::ZarrCodec(codecs::Name::Sharding,  shardingCodecConfigJSON)});

				Image<uint16_t> fromDisk;
				zarr::read(fromDisk, path);

				testAssert(equals(img, fromDisk), string("zarr test write sharding with indexLocation=" + indexLocation + " withBlosc="+ toString(withBlosc)));
			}

			void shardingEmptyInnerChunksTest(string indexLocation)
			{
				string path = "./testoutput/test_sharding_empty_inner_chunks_" + indexLocation;
				Image<uint16_t> img(Vec3c(10, 10, 10));
				add(img, uint16_t());//TODO DEFAULT_FILLVALUE);
				nlohmann::json shardingCodecConfigJSON = {
					{ "chunk_shape", { 5, 5, 5 }},
					{ "codecs", { codecs::ZarrCodec(codecs::Name::Bytes).toJSON() }},
					{ "index_codecs", { codecs::ZarrCodec(codecs::Name::Bytes).toJSON() }},
					{ "index_location", indexLocation }
				};

				img(0, 0, 0) = 1; //to not have complete shard empty
				zarr::write(img,
					path + "_empty",
					zarr::DEFAULT_CHUNK_SIZE,
					{ codecs::ZarrCodec(codecs::Name::Sharding, shardingCodecConfigJSON) });

				ramp(img, 0);
				zarr::write(img,
					path + "_full",
					zarr::DEFAULT_CHUNK_SIZE,
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
				shardingEmptyInnerChunksTest("start");
				shardingEmptyInnerChunksTest("end");
			}
			void emptyChunks(){
				string path = "./testoutput/test_empty_chunks";
				Image<uint16_t> img(Vec3c(8, 8, 8));
				add(img, uint16_t());//TODO DEFAULT_FILLVALUE);

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
		}
	}
}