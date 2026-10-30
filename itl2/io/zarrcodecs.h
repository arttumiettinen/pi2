#pragma once

#include <string>
#include <sstream> // Ensure you include this for stringstream
#include <cstring>
#include <cstdint>
#include <blosc.h>
#include <zstd.h>

#include "json.h"
#include "utilities.h"
#include "iteration.h"

namespace itl2
{
	using std::cout, std::endl;

	namespace zarr
	{
		typedef coord_t fillValue_t; //TODO: allow other fillValue types
	}
	namespace zarr::codecs
	{

		enum class Type
		{
			ArrayArrayCodec,
			ArrayBytesCodec,
			BytesBytesCodec,
		};

		enum class Name
		{
			Bytes,
			Transpose,
			Blosc,
			Sharding,
			Zstd,
			Crc32c,
		};
		//TODO save codec in struct instead of json
		namespace blosc
		{
			enum class shuffle
			{
				noshuffle = 0,
				shuffle = 1,
				bitshuffle = 2,
			};
		}
		namespace sharding
		{
			enum class indexLocation
			{
				start,
				end,
			};
		}
	}

	template<>
	inline std::string toString(const zarr::codecs::Name& x)
	{
		switch (x)
		{
		case zarr::codecs::Name::Bytes:
			return "bytes";
		case zarr::codecs::Name::Transpose:
			return "transpose";
		case zarr::codecs::Name::Blosc:
			return "blosc";
		case zarr::codecs::Name::Sharding:
			return "sharding_indexed";
		case zarr::codecs::Name::Zstd:
			return "zstd";
		case zarr::codecs::Name::Crc32c:
			return "crc32c";
		}
		throw ITLException("Invalid zarr codec name.");
	}

	template<>
	inline zarr::codecs::Name fromString(const std::string& str0)
	{
		std::string str = str0;
		toLower(str);
		if (str == "bytes")
			return zarr::codecs::Name::Bytes;
		if (str == "transpose")
			return zarr::codecs::Name::Transpose;
		if (str == "blosc")
			return zarr::codecs::Name::Blosc;
		if (str == "sharding_indexed")
			return zarr::codecs::Name::Sharding;
		if (str == "zstd")
			return zarr::codecs::Name::Zstd;
		if (str == "crc32c")
			return zarr::codecs::Name::Crc32c;
		throw ITLException(std::string("Invalid zarr codec name: ") + str);
	}

	namespace zarr::codecs
	{
		class ZarrCodec
		{

		 public:
			typedef std::list<codecs::ZarrCodec> Pipeline;

			Type type;
			Name name;
			nlohmann::json configuration;

			ZarrCodec(const Name name, nlohmann::json config = nlohmann::json())
			{
				this->name = name;
				switch (name)
				{
				case Name::Bytes:
					this->type = Type::ArrayBytesCodec;
					parseBytesCodecConfig(config);
					break;
				case Name::Transpose:
					this->type = Type::ArrayArrayCodec;
					parseTransposeCodecConfig(config);
					break;
				case Name::Blosc:
					this->type = Type::BytesBytesCodec;
					parseBloscCodecConfig(config);
					break;
				case Name::Sharding:
					this->type = Type::ArrayBytesCodec;
					parseShardingCodecConfig(config);
					break;
				case Name::Zstd:
					this->type = Type::BytesBytesCodec;
					parseZstdCodecConfig(config);
					break;
				case Name::Crc32c:
					this->type = Type::BytesBytesCodec;
					// crc32c has no configuration.
					this->configuration = config;
					break;
				default:
					throw ITLException(std::string("Invalid zarr codec"));
				}
			}

			nlohmann::json toJSON() const
			{
				nlohmann::json j;
				j["name"] = toString(this->name);
				j["configuration"] = configuration; //todo: this does not return by reference to the configuration just a copy, right?
				return j;
			}

			bool operator==(const ZarrCodec& t) const
			{
				return toJSON() == t.toJSON();
			}

			void parseBloscCodecConfig(nlohmann::json config = nlohmann::json())
			{
				//TODO: validate
				this->configuration = config;
			}

			void parseTransposeCodecConfig(nlohmann::json config = nlohmann::json())
			{
				//TODO: validate
				this->configuration = config;
			}

			void parseZstdCodecConfig(nlohmann::json config = nlohmann::json())
			{
				//TODO: validate
				this->configuration = config;
			}

			void getZstdConfiguration(int& level, bool& checksum) const
			{
				if (this->name != Name::Zstd) throw ITLException("only zstd codec has zstd config");
				level = this->configuration.contains("level") ? this->configuration["level"].get<int>() : 0;
				checksum = this->configuration.contains("checksum") && this->configuration["checksum"].get<bool>();
			}

			void parseBytesCodecConfig(nlohmann::json config = nlohmann::json())
			{
				std::string endian = "little";
				for (auto it = config.begin(); it != config.end(); ++it)
				{
					if (it.key() == "endian")
					{
						endian = it.value();
						if (endian != "little" && endian != "big")
						{
							throw ITLException("Invalid endian in bytes codec config: " + endian);
						}
					}
					else
					{
						throw ITLException("Invalid key in bytes codec config: " + it.key());
					}
				}
				this->configuration = {
					{ "endian", endian }
				};
			}

			void getShardingConfiguration(Vec3c& chunkShape, Pipeline& codecs, Pipeline& indexCodecs, sharding::indexLocation& indexLocation) const;

			void parseShardingCodecConfig(nlohmann::json config = nlohmann::json())
			{
				//TODO: validate inner chunk size lies within shard size
				this->configuration = config;

				Vec3c dummyChunkShape;
				Pipeline dummyCodecs;
				Pipeline dummyIndexCodecs;
				sharding::indexLocation dummyIndexLocation;
				this->getShardingConfiguration(dummyChunkShape, dummyCodecs, dummyIndexCodecs, dummyIndexLocation);
			}

			void getBloscConfiguration(string& cname, int& clevel, blosc::shuffle& shuffle, size_t& typesize, size_t& blocksize) const
			{
				if (this->name != Name::Blosc) throw ITLException("only blosc codec has blosc config");
				try
				{
					cname = this->configuration["cname"];
					std::list<string> allowedCnames = { "blosclz", "lz4", "lz4hc", "zlib", "zstd" };
					if (!listContains<string>(allowedCnames, cname)) throw ITLException("invalid blosc cname: " + cname);
					clevel = this->configuration["clevel"];
					string shuffleName = this->configuration["shuffle"];
					if (shuffleName == "noshuffle") shuffle = blosc::shuffle::noshuffle;
					else if (shuffleName == "shuffle") shuffle = blosc::shuffle::shuffle;
					else if (shuffleName == "bitshuffle") shuffle = blosc::shuffle::bitshuffle;
					else throw ITLException("invalid blosc shuffle parameter: " + shuffleName);
					if (this->configuration.contains("typesize")) {
						typesize = this->configuration["typesize"];
						if (typesize < 1) typesize = 1;
					}
					else typesize = 1;
					blocksize = this->configuration["blocksize"];
				}
				catch (nlohmann::json::exception ex)
				{
					throw ITLException("error in reading blosc config: " + nlohmann::to_string(this->configuration) + " got exception: " + ex.what());
				}
			}

			void getTransposeConfiguration(Vec3c& order) const
			{
				try
				{
					if (this->name != Name::Transpose) throw ITLException("only transpose codec has transpose order");
					auto orderJSON = this->configuration["order"];
					order = Vec3c(0, 1, 2);
					order[0] = orderJSON[0].get<size_t>();
					if (orderJSON.size() >= 2)
						order[1] = orderJSON[1].get<size_t>();
					if (orderJSON.size() >= 3)
						order[2] = orderJSON[2].get<size_t>();
					if (!order.isPermutation()) throw ITLException("TransposeConfiguration error: invalid order: " + toString(order) + "expected a permutation of [0, 1, 2]");
				}
				catch (nlohmann::json::exception ex)
				{
					throw ITLException("error in reading transposeOrder of configuration: " + nlohmann::to_string(this->configuration) + " got exception: " + ex.what());
				}
			}

			friend std::ostream& operator<<(std::ostream& stream, const ZarrCodec& c)
			{
				stream << c.toJSON();
				return stream;
			}

			friend std::ostream& operator<<(std::ostream& stream, const Pipeline& p)
			{
				stream << "[";
				for (const ZarrCodec c : p)
				{
					stream << c;
					if (!(c == (*p.rbegin())))
						stream << ", ";
				}
				stream << "]";
				return stream;
			}

		};

		typedef ZarrCodec::Pipeline Pipeline;

		inline bool fromJSON(Pipeline& codecs, nlohmann::json codecsJSON, string& reason)
		{
			int numberArrayBytesCodecs = 0;
			for (auto& codec : codecsJSON)
			{
				if (!codec.contains("name"))
				{
					throw ITLException("codec name is missing in zarr metadata.");
				}

				nlohmann::json codecConfig = {};
				if (codec.contains("configuration"))
				{
					codecConfig = codec["configuration"];
				}
				Name zarrCodecName = fromString<Name>(codec["name"].get<string>());
				ZarrCodec zarrCodec = ZarrCodec(zarrCodecName, codecConfig);
				codecs.push_back(zarrCodec);
				switch (zarrCodec.type)
				{
				case Type::ArrayArrayCodec:
					if (numberArrayBytesCodecs > 0)
					{
						throw ITLException("ArrayArrayCodec cannot be used after ArrayBytesCodec.");
					}
					break;
				case Type::ArrayBytesCodec:
					numberArrayBytesCodecs++;
					break;
				case Type::BytesBytesCodec:
					if (numberArrayBytesCodecs < 1)
					{
						throw ITLException("ArrayBytesCodec must be used before BytesBytesCodec.");
					}
					break;
				default:
					reason = "Unknown codec type.";
					return false;
				}
			}

			if (numberArrayBytesCodecs != 1)
			{
				throw ITLException("Exactly one ArrayBytesCodec was expected in the codecs list, got " + std::to_string(numberArrayBytesCodecs) + ".");
			}
			return true;
		}

		// Rewrites a codecs JSON pipeline (in place) to drop the leading singleton axes
		// that were squeezed away when reading a >3-dimensional array as 3D. This keeps
		// the codec configurations consistent with the squeezed 3D shape:
		//  - transpose: remove the leading axes from the "order" and renumber the rest.
		//    Valid because size-1 axes do not affect the byte layout of the other axes.
		//  - sharding_indexed: squeeze the inner "chunk_shape" and recurse into its inner
		//    "codecs". The "index_codecs" operate on the shard index grid and are left
		//    unchanged.
		inline void squeezeLeadingSingletons(nlohmann::json& codecsJSON, size_t numLeadingSingletons)
		{
			if (numLeadingSingletons == 0)
				return;
			for (auto& codec : codecsJSON)
			{
				if (!codec.contains("name"))
					continue;
				string name = codec["name"].get<string>();
				toLower(name);
				if (name == "transpose"
					&& codec.contains("configuration")
					&& codec["configuration"].contains("order"))
				{
					nlohmann::json squeezed = nlohmann::json::array();
					for (auto& e : codec["configuration"]["order"])
					{
						size_t axis = e.get<size_t>();
						if (axis >= numLeadingSingletons)
							squeezed.push_back(axis - numLeadingSingletons);
					}
					codec["configuration"]["order"] = squeezed;
				}
				else if (name == "sharding_indexed" && codec.contains("configuration"))
				{
					auto& cfg = codec["configuration"];
					if (cfg.contains("chunk_shape"))
					{
						nlohmann::json squeezedShape = nlohmann::json::array();
						auto& cs = cfg["chunk_shape"];
						for (size_t axis = numLeadingSingletons; axis < cs.size(); axis++)
							squeezedShape.push_back(cs[axis]);
						cfg["chunk_shape"] = squeezedShape;
					}
					if (cfg.contains("codecs"))
						squeezeLeadingSingletons(cfg["codecs"], numLeadingSingletons);
				}
			}
		}

		template<typename pixel_t>
		void encodePipeline(const Pipeline& codecs, Image<pixel_t>& image, std::vector<char>& buffer, fillValue_t fillValue);

		template<typename pixel_t>
		void decodePipeline(const Pipeline& codecs, Image<pixel_t>& image, std::vector<char>& buffer, fillValue_t fillValue);

		template<typename pixel_t>
		void encodeTransposeCodec(const ZarrCodec& codec, Image<pixel_t>& image, fillValue_t fillValue)
		{
			Vec3c order;
			codec.getTransposeConfiguration(order);
			transpose(image, order, fillValue);
		}

		template<typename pixel_t>
		void decodeTransposeCodec(const ZarrCodec& codec, Image<pixel_t>& image, fillValue_t fillValue)
		{
			Vec3c order;
			codec.getTransposeConfiguration(order);
			transpose(image, order.inverseOrder(), fillValue);
		}

		inline void encodeBloscCodec(const ZarrCodec& codec, std::vector<char>& buffer)
		{
			size_t destSize = buffer.size() + BLOSC_MIN_HEADER_LENGTH;
			size_t srcSize = buffer.size();
			std::vector<char> temp(destSize);
			string cname;
			int clevel;
			codecs::blosc::shuffle shuffle;
			size_t typesize;
			size_t blocksize;
			codec.getBloscConfiguration(cname, clevel, shuffle, typesize, blocksize);
			int numinternalthreads = 1;

			size_t realDestSize = blosc_compress_ctx(clevel, (int)shuffle, typesize, srcSize, buffer.data(), temp.data(), destSize, cname.c_str(), blocksize, numinternalthreads);

			if (realDestSize == 0) cout << "Buffer is incompressible.  Giving up." << endl;
			else if (realDestSize < 0) throw ITLException("Compression error.  Error code: " + toString(realDestSize));

			buffer.resize(realDestSize);
			std::memcpy(buffer.data(), temp.data(), realDestSize);
		}

		inline void decodeBloscCodec(const ZarrCodec& codec, std::vector<char>& buffer)
		{
			size_t srcSize = buffer.size();
			size_t destSize;
			if (blosc_cbuffer_validate(buffer.data(), srcSize, &destSize) < 0)
			{
				throw ITLException("blosc_decompress error: \"Buffer does not contain valid blosc-encoded contents\"");
			}
			std::vector<char> temp(destSize);
			int numinternalthreads = 1;
			size_t realDestSize = blosc_decompress_ctx(buffer.data(), temp.data(), destSize, numinternalthreads);
			if (realDestSize < 0)
			{
				throw ITLException("blosc_decompress error.  Error code: " + toString(realDestSize));
			}
			buffer.resize(realDestSize);
			std::memcpy(buffer.data(), temp.data(), realDestSize);
		}

		inline void encodeZstdCodec(const ZarrCodec& codec, std::vector<char>& buffer)
		{
			int level;
			bool checksum;
			codec.getZstdConfiguration(level, checksum);
			if (level == 0)
				level = ZSTD_CLEVEL_DEFAULT;

			size_t destBound = ZSTD_compressBound(buffer.size());
			std::vector<char> temp(destBound);

			ZSTD_CCtx* cctx = ZSTD_createCCtx();
			if (cctx == nullptr)
				throw ITLException("zstd: could not create compression context.");
			ZSTD_CCtx_setParameter(cctx, ZSTD_c_compressionLevel, level);
			ZSTD_CCtx_setParameter(cctx, ZSTD_c_checksumFlag, checksum ? 1 : 0);
			size_t realDestSize = ZSTD_compress2(cctx, temp.data(), destBound, buffer.data(), buffer.size());
			ZSTD_freeCCtx(cctx);
			if (ZSTD_isError(realDestSize))
				throw ITLException(string("zstd compression error: ") + ZSTD_getErrorName(realDestSize));

			buffer.resize(realDestSize);
			std::memcpy(buffer.data(), temp.data(), realDestSize);
		}

		inline void decodeZstdCodec(const ZarrCodec& codec, std::vector<char>& buffer)
		{
			unsigned long long destSize = ZSTD_getFrameContentSize(buffer.data(), buffer.size());
			if (destSize == ZSTD_CONTENTSIZE_ERROR)
				throw ITLException("zstd: buffer does not contain a valid zstd frame.");
			if (destSize == ZSTD_CONTENTSIZE_UNKNOWN)
				throw ITLException("zstd: decompressed size is not stored in the frame header, which this implementation requires.");

			std::vector<char> temp(destSize);
			size_t realDestSize = ZSTD_decompress(temp.data(), destSize, buffer.data(), buffer.size());
			if (ZSTD_isError(realDestSize))
				throw ITLException(string("zstd decompression error: ") + ZSTD_getErrorName(realDestSize));

			buffer.resize(realDestSize);
			std::memcpy(buffer.data(), temp.data(), realDestSize);
		}

		// Software CRC32C (Castagnoli polynomial, reflected form 0x82F63B78).
		inline uint32_t crc32c(const uint8_t* data, size_t length)
		{
			uint32_t crc = 0xFFFFFFFFu;
			for (size_t i = 0; i < length; i++)
			{
				crc ^= data[i];
				for (int k = 0; k < 8; k++)
					crc = (crc >> 1) ^ (0x82F63B78u & (0u - (crc & 1u)));
			}
			return crc ^ 0xFFFFFFFFu;
		}

		// crc32c codec appends a 4-byte little-endian CRC32C of the preceding bytes.
		inline void encodeCrc32cCodec(std::vector<char>& buffer)
		{
			uint32_t crc = crc32c(reinterpret_cast<const uint8_t*>(buffer.data()), buffer.size());
			char crcBytes[4];
			for (int i = 0; i < 4; i++)
				crcBytes[i] = static_cast<char>((crc >> (8 * i)) & 0xFFu);
			buffer.insert(buffer.end(), crcBytes, crcBytes + 4);
		}

		inline void decodeCrc32cCodec(std::vector<char>& buffer)
		{
			if (buffer.size() < 4)
				throw ITLException("crc32c: buffer too small to contain a checksum.");
			size_t dataLen = buffer.size() - 4;
			uint32_t stored = 0;
			for (int i = 0; i < 4; i++)
				stored |= static_cast<uint32_t>(static_cast<uint8_t>(buffer[dataLen + i])) << (8 * i);
			uint32_t computed = crc32c(reinterpret_cast<const uint8_t*>(buffer.data()), dataLen);
			if (stored != computed)
				throw ITLException("crc32c checksum mismatch: the data is corrupt or the checksum is invalid.");
			buffer.resize(dataLen);
		}

		//todo: use swapByteOrder(img); depending on endian
		template<typename pixel_t>
		void decodeBytesCodec(Image<pixel_t>& image, std::vector<char>& buffer)
		{
			Vec3c shape = image.dimensions();
			if (shape.product() * sizeof(pixel_t) != buffer.size())
				throw ITLException("wrong buffer size for decoding bytes codec. buffersize=" + toString(buffer.size()) + " chunkSize=" + toString(shape.product()));

			pixel_t* buffer2 = (pixel_t*)buffer.data();
			size_t n = 0;
			for (coord_t x = 0; x < shape.x; x++)
			{
				for (coord_t y = 0; y < shape.y; y++)
				{
					for (coord_t z = 0; z < shape.z; z++)
					{
						//todo: might be faster to read data sequentially
						image(x, y, z) = buffer2[n++];
					}
				}
			}
		}

		template<typename pixel_t>
		void encodeBytesCodec(const Image<pixel_t>& image, std::vector<char>& buffer, size_t pixelSize = sizeof(pixel_t))
		{
			Vec3c shape = image.dimensions();
			size_t bufferSize = shape.product() * pixelSize;
			buffer = std::vector<char>(bufferSize);
			pixel_t* buffer2 = (pixel_t*)buffer.data();

			size_t n = 0;
			for (coord_t x = 0; x < shape.x; x++)
			{
				for (coord_t y = 0; y < shape.y; y++)
				{
					for (coord_t z = 0; z < shape.z; z++)
					{
						buffer2[n++] = image(x, y, z);
					}
				}
			}
		}

		//forward declaration (defined below)
		inline void decodeBytesBytesCodec(const ZarrCodec& codec, std::vector<char>& buffer);

		// Validates that a sharding index_codecs pipeline is one this implementation
		// supports: a single Bytes (ArrayBytes) codec optionally followed by BytesBytes
		// codecs (e.g. crc32c). Returns the number of extra bytes those BytesBytes
		// codecs append to the encoded index (e.g. 4 bytes for crc32c).
		inline coord_t validateShardingIndexCodecs(const Pipeline& indexCodecs)
		{
			if (indexCodecs.empty() || indexCodecs.begin()->name != Name::Bytes)
				throw ITLException("This zarr implementation requires the Bytes codec at the first position of the sharding index_codecs.");
			coord_t overhead = 0;
			for (auto it = std::next(indexCodecs.begin()); it != indexCodecs.end(); ++it)
			{
				if (it->name == Name::Crc32c)
					overhead += 4;
				else
					throw ITLException("This zarr implementation only supports the Bytes codec optionally followed by crc32c in the sharding index_codecs, but got: " + toString(it->name) + ".");
			}
			return overhead;
		}

		template<typename pixel_t>
		void decodeShardingCodec(const ZarrCodec& codec, Image<pixel_t>& shard, std::vector<char>& buffer, fillValue_t fillValue)
		{
			typedef uint64_t index_t;

			Vec3c innerChunkShape;
			Pipeline codecs;
			Pipeline indexCodecs;
			sharding::indexLocation indexLocation;
			codec.getShardingConfiguration(innerChunkShape, codecs, indexCodecs, indexLocation);

			Vec3c chunksPerShard = shard.dimensions().componentwiseDivide(innerChunkShape);
			coord_t chunkCount = chunksPerShard.product();
			coord_t indexCodecOverhead = validateShardingIndexCodecs(indexCodecs);
			coord_t rawIndexSize = 2 * sizeof(index_t) * chunkCount;
			coord_t indexSize = rawIndexSize + indexCodecOverhead;
			std::vector<char> indexBuffer;

			switch (indexLocation)
			{
			case sharding::indexLocation::start:
				indexBuffer.insert(indexBuffer.end(), buffer.begin(), buffer.begin() + indexSize);
				break;
			case sharding::indexLocation::end:
				indexBuffer.insert(indexBuffer.end(), buffer.end() - indexSize, buffer.end());
			}

			// Apply the BytesBytes index codecs in reverse to verify and strip any
			// trailing data (e.g. the crc32c checksum), leaving the raw index bytes.
			for (auto it = indexCodecs.rbegin(); it != indexCodecs.rend() && it->type == Type::BytesBytesCodec; ++it)
				decodeBytesBytesCodec(*it, indexBuffer);

			Image<index_t> shardIndexArrayOffsets(chunksPerShard);
			Image<index_t> shardIndexArrayNBytes(chunksPerShard);

			//decode indexArray buffer into shardIndexArrayOffsets and shardIndexArrayNBytes
			//TODO: extract to decodeBytesCodec
			std::vector<index_t> temp(chunkCount * 2);
			std::memcpy(temp.data(), indexBuffer.data(), rawIndexSize);

			size_t n = 0;
			for (coord_t x = 0; x < chunksPerShard.x; x++)
			{
				for (coord_t y = 0; y < chunksPerShard.y; y++)
				{
					for (coord_t z = 0; z < chunksPerShard.z; z++)
					{
						shardIndexArrayOffsets(x, y, z) = temp[n++];
						shardIndexArrayNBytes(x, y, z) = temp[n++];
					}
				}
			}
			forAllChunks(shard.dimensions(), innerChunkShape, [&](const Vec3c& chunkIndex, const Vec3c& chunkStart)
			{
			  index_t nBytes = shardIndexArrayNBytes(chunkIndex);
			  index_t offset = shardIndexArrayOffsets(chunkIndex);

			  //TODO: check both variables. for testing nBytes is sufficient as the library used to test against had partially wrong offset values
			  //if(nBytes==-1 && offset==-1){
			  if (nBytes == -1)
			  {
				  size_t ndrawn = draw(shard, AABoxc::fromMinMax(chunkStart, chunkStart + innerChunkShape), pixelRound<pixel_t>(fillValue));
#if defined(_DEBUG) || defined(BOUNDS_CHECK)
				  assert(ndrawn==innerChunkShape.product());
#endif

			  }
			  else
			  {
				  std::vector<char> chunkBuffer;
				  auto chunkBegin = buffer.begin() + offset;
				  auto chunkEnd = chunkBegin + nBytes;

				  chunkBuffer.insert(chunkBuffer.end(), chunkBegin, chunkEnd);

				  Image<pixel_t> innerChunk(innerChunkShape);
				  decodePipeline(codecs, innerChunk, chunkBuffer, fillValue);

				  forAllPixels(innerChunk, [&](coord_t x, coord_t y, coord_t z)
				  {
					Vec3c pos = Vec3c(x, y, z) + chunkStart;
					shard(pos) = innerChunk(x, y, z);
				  });
			  }
			});

		}

		//forward declaration
		inline void encodeBytesBytesCodec(const ZarrCodec& codec, std::vector<char>& buffer);

		template<typename pixel_t>
		void encodeShardingCodec(const ZarrCodec& codec, const Image<pixel_t>& shard, std::vector<char>& buffer, fillValue_t fillValue)
		{
			typedef uint64_t index_t;

			Vec3c innerChunkShape;
			Pipeline codecs;
			Pipeline indexCodecs;
			sharding::indexLocation indexLocation;
			codec.getShardingConfiguration(innerChunkShape, codecs, indexCodecs, indexLocation);

			coord_t indexCodecOverhead = validateShardingIndexCodecs(indexCodecs);

			if (!(shard.dimensions() >= innerChunkShape)) throw ITLException("inner chunk shape " + toString(innerChunkShape) + " does not fit into shard shape " + toString(shard.dimensions()));
			Vec3c chunksPerShard = shard.dimensions().componentwiseDivide(innerChunkShape);
			if (!(chunksPerShard.componentwiseMultiply(innerChunkShape) == shard.dimensions()))
				throw ITLException("inner chunk shape " + toString(innerChunkShape) + " does not evenly divide shard shape " + toString(shard.dimensions()));

			coord_t chunkCount = chunksPerShard.product();
			Image<index_t> shardIndexArrayOffsets(chunksPerShard);
			Image<index_t> shardIndexArrayNBytes(chunksPerShard);
			coord_t indexSize = 2 * sizeof(index_t) * chunkCount + indexCodecOverhead; //includes overhead of BytesBytes index codecs (e.g. crc32c)
			Image<std::vector<char>> chunkBytes(chunksPerShard);

			//TODO: concurrency needed? then working with fixed nbytes would be necessary
			//TODO: empty innerChunks? setting both values offset and nbytes to 2^64 - 1
			std::vector<char> shardBuffer;
			forAllChunks(shard.dimensions(), innerChunkShape, [&](const Vec3c& chunkIndex, const Vec3c& chunkStart)
			{
			  AABoxc currentInnerChunk = AABoxc::fromPosSize(chunkStart, innerChunkShape);
			  Image<pixel_t> innerChunk(innerChunkShape);
			  bool innerChunkEmpty = true;
			  //write pixels from shard to innerChunk
			  forAllPixels(innerChunk, [&](coord_t x, coord_t y, coord_t z)
			  {
				Vec3c pos = Vec3c(x, y, z) + chunkStart;
				innerChunk(x, y, z) = shard(pos);
				if (!(shard(pos) == pixelRound<pixel_t>(fillValue)))
				{
					innerChunkEmpty = false;
				}
			  });
			  if (innerChunkEmpty)
			  {
				  shardIndexArrayOffsets(chunkIndex) = -1;
				  shardIndexArrayNBytes(chunkIndex) = -1;
			  }
			  else
			  {
				  std::vector<char> chunkBuffer;
				  encodePipeline(codecs, innerChunk, chunkBuffer, fillValue);

				  shardIndexArrayOffsets(chunkIndex) = shardBuffer.size(); //position of shardBuffer.end() where we will write the data of chunkBuffer
				  if (indexLocation == sharding::indexLocation::start)
				  {
					  shardIndexArrayOffsets(chunkIndex) += indexSize;
				  }
				  shardIndexArrayNBytes(chunkIndex) = chunkBuffer.size();
				  shardBuffer.insert(shardBuffer.end(), chunkBuffer.begin(), chunkBuffer.end());
			  }
			});

			//apply encodePipeline for shardIndexArray, but we do not support 4d arrays
			if (std::find_if(indexCodecs.begin(), indexCodecs.end(), [](ZarrCodec codec)
			{ return codec.type == Type::ArrayArrayCodec; }) != indexCodecs.end())
				throw ITLException("This zarr implementation does not support ArrayArrayCodecs within the sharding index_codecs");
			auto indexCodec = indexCodecs.begin();
			if (indexCodec->name != Name::Bytes)
				throw ITLException("This zarr implementation only supports the Bytes Codec at the first position of the sharding index_codecs");

			//apply encodeBytesCodec for shardIndexArrayOffsets and shardIndexArrayNBytes combined
			//TODO: extract to encodeBytesCodec
			std::vector<char> indexBuffer;
			std::vector<index_t> temp(chunkCount * 2);
			size_t n = 0;
			for (coord_t x = 0; x < chunksPerShard.x; x++)
			{
				for (coord_t y = 0; y < chunksPerShard.y; y++)
				{
					for (coord_t z = 0; z < chunksPerShard.z; z++)
					{
						temp[n++] = shardIndexArrayOffsets(x, y, z);
						temp[n++] = shardIndexArrayNBytes(x, y, z);
					}
				}
			}
			size_t indexBufferSize = chunkCount * sizeof(index_t) * 2;
			indexBuffer = std::vector<char>(indexBufferSize);
			std::memcpy(indexBuffer.data(), temp.data(), indexBuffer.size());

			//encode BytesBytesCodecs
			++indexCodec;
			for (; indexCodec != indexCodecs.end(); ++indexCodec)
			{
				encodeBytesBytesCodec(*indexCodec, indexBuffer);
			}

			switch (indexLocation)
			{
			case sharding::indexLocation::start:
				buffer.insert(buffer.end(), indexBuffer.begin(), indexBuffer.end());
				buffer.insert(buffer.end(), shardBuffer.begin(), shardBuffer.end());
				break;
			case sharding::indexLocation::end:
				buffer.insert(buffer.end(), shardBuffer.begin(), shardBuffer.end());
				buffer.insert(buffer.end(), indexBuffer.begin(), indexBuffer.end());
			}
		}

		inline void decodeBytesBytesCodec(const ZarrCodec& codec, std::vector<char>& buffer)
		{
			assert(codec.type == codecs::Type::BytesBytesCodec);
			if (codec.name == codecs::Name::Blosc)
			{
				decodeBloscCodec(codec, buffer);
			}
			else if (codec.name == codecs::Name::Zstd)
			{
				decodeZstdCodec(codec, buffer);
			}
			else if (codec.name == codecs::Name::Crc32c)
			{
				decodeCrc32cCodec(buffer);
			}
			else throw ITLException("BytesBytesCodec: " + toString(codec.name) + " not yet implemented");
		}

		template<typename pixel_t>
		void decodeArrayBytesCodec(const ZarrCodec& codec, Image<pixel_t>& image, std::vector<char>& buffer, fillValue_t fillValue)
		{
			assert(codec.type == codecs::Type::ArrayBytesCodec);
			if (codec.name == codecs::Name::Bytes)
			{
				decodeBytesCodec(image, buffer);
			}
			else if (codec.name == codecs::Name::Sharding)
			{
				decodeShardingCodec(codec, image, buffer, fillValue);
			}
			else throw ITLException("ArrayBytesCodec: " + toString(codec.name) + " not yet implemented");
		}

		template<typename pixel_t>
		void decodeArrayArrayCodec(const ZarrCodec& codec, Image<pixel_t>& image, fillValue_t fillValue)
		{
			assert(codec.type == codecs::Type::ArrayArrayCodec);
			if (codec.name == codecs::Name::Transpose)
			{
				decodeTransposeCodec(codec, image, fillValue);
			}
			else throw ITLException("ArrayArrayCodec: " + toString(codec.name) + " not yet implemented");

		}

		inline void encodeBytesBytesCodec(const ZarrCodec& codec, std::vector<char>& buffer)
		{
			assert(codec.type == codecs::Type::BytesBytesCodec);
			if (codec.name == codecs::Name::Blosc)
			{
				encodeBloscCodec(codec, buffer);
			}
			else if (codec.name == codecs::Name::Zstd)
			{
				encodeZstdCodec(codec, buffer);
			}
			else if (codec.name == codecs::Name::Crc32c)
			{
				encodeCrc32cCodec(buffer);
			}
			else throw ITLException("BytesBytesCodec: " + toString(codec.name) + " not yet implemented");
		}

		template<typename pixel_t>
		void encodeArrayBytesCodec(const ZarrCodec& codec, Image<pixel_t>& image, std::vector<char>& buffer, fillValue_t fillValue)
		{
			assert(codec.type == codecs::Type::ArrayBytesCodec);
			if (codec.name == codecs::Name::Bytes)
			{
				encodeBytesCodec(image, buffer);
			}
			else if (codec.name == codecs::Name::Sharding)
			{
				encodeShardingCodec(codec, image, buffer, fillValue);
			}
			else throw ITLException("ArrayBytesCodec: " + toString(codec.name) + " not yet implemented");
		}

		template<typename pixel_t>
		void encodeArrayArrayCodec(const ZarrCodec& codec, Image<pixel_t>& image, fillValue_t fillValue)
		{
			assert(codec.type == codecs::Type::ArrayArrayCodec);
			if (codec.name == codecs::Name::Transpose)
			{
				encodeTransposeCodec(codec, image, fillValue);
			}
			else throw ITLException("ArrayArrayCodec: " + toString(codec.name) + " not yet implemented");
		}

		template<typename pixel_t>
		void encodePipeline(const Pipeline& codecs, Image<pixel_t>& image, std::vector<char>& buffer, fillValue_t fillValue)
		{
			codecs::Pipeline::const_iterator codec = codecs.begin();
			for (; codec->type == codecs::Type::ArrayArrayCodec; ++codec)
			{
				codecs::encodeArrayArrayCodec(*codec, image, fillValue);
			}
			codecs::encodeArrayBytesCodec(*codec, image, buffer, fillValue);
			++codec;
			for (; codec != codecs.end(); ++codec)
			{
				codecs::encodeBytesBytesCodec(*codec, buffer);
			}
		}

		template<typename pixel_t>
		void decodePipeline(const Pipeline& codecs, Image<pixel_t>& image, std::vector<char>& buffer, fillValue_t fillValue)
		{
			//todo: does this work with const codecs
			codecs::Pipeline::const_reverse_iterator codec = codecs.rbegin();
			for (; codec->type == codecs::Type::BytesBytesCodec; ++codec)
			{
				codecs::decodeBytesBytesCodec(*codec, buffer);
			}
			codecs::decodeArrayBytesCodec(*codec, image, buffer, fillValue);
			++codec;
			for (; codec != codecs.rend(); ++codec)
			{
				codecs::decodeArrayArrayCodec(*codec, image, fillValue);
			}
		}
	}
}

inline void itl2::zarr::codecs::ZarrCodec::getShardingConfiguration(Vec3c& chunkShape, Pipeline& codecs, Pipeline& indexCodecs, sharding::indexLocation& indexLocation) const
{
	try
	{
		if (this->name != Name::Sharding) throw ITLException("only sharding codec has sharding config");
		auto chunkShapeJSON = this->configuration["chunk_shape"];
		chunkShape = Vec3c(1, 1, 1);
		chunkShape[0] = chunkShapeJSON[0].get<size_t>();
		if (chunkShapeJSON.size() >= 2)
			chunkShape[1] = chunkShapeJSON[1].get<size_t>();
		if (chunkShapeJSON.size() >= 3)
			chunkShape[2] = chunkShapeJSON[2].get<size_t>();

		if (this->configuration.contains("codecs"))
		{
			string reason;
			if(!fromJSON(codecs, this->configuration["codecs"], reason))
				throw ITLException("Could not decode sharding codecs " + reason);
		}
		else
		{
			throw ITLException("No sharding codecs.");
		}

		if (this->configuration.contains("index_codecs"))
		{
			string reason;
			if(!fromJSON(indexCodecs, this->configuration["index_codecs"], reason))
				throw ITLException("could not decode sharding index_codecs " + reason);
		}
		else
		{
			throw ITLException("No sharding index_codecs.");
		}

		indexLocation = sharding::indexLocation::end;
		if (this->configuration.contains("index_location"))
		{
			string indexLocationString = this->configuration["index_location"];
			if (indexLocationString == "start")
				indexLocation = sharding::indexLocation::start;
			else if (indexLocationString == "end")
				indexLocation = sharding::indexLocation::end;
			else
				throw ITLException("Invalid sharding index_location: " + indexLocationString);
		}
	}
	catch (nlohmann::json::exception ex)
	{
		throw ITLException("error in reading sharding configuration: " + nlohmann::to_string(this->configuration) + " got exception: " + ex.what());
	}
}