#pragma once

#include <string>
#include <sstream> // Ensure you include this for stringstream
#include <blosc.h>

#include "json.h"
#include "utilities.h"
#include "iteration.h"

namespace itl2
{
	using std::cout, std::endl;
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
		namespace sharding{
			enum class indexLocation{
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
				default:
					throw ITLException(std::string("Invalid zarr codec"));
				}
				cout << "created codec: " << toString(this->name) << " " << (int)this->type << endl;
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

			void parseShardingCodecConfig(nlohmann::json config = nlohmann::json())
			{
				//TODO: validate
				this->configuration = config;
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
					if (shuffle != blosc::shuffle::noshuffle)
						typesize = this->configuration["typesize"];
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
					cout << "transposeOrder=" << order << endl;
					if (!order.isPermutation()) throw ITLException("TransposeConfiguration error: invalid order: " + toString(order) + "expected a permutation of [0, 1, 2]");
				}
				catch (nlohmann::json::exception ex)
				{
					throw ITLException("error in reading transposeOrder of configuration: " + nlohmann::to_string(this->configuration) + " got exception: " + ex.what());
				}
			}

			void getShardingConfiguration(Vec3c& chunkShape, Pipeline& codecs, Pipeline& indexCodecs, sharding::indexLocation& indexLocation) const;
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
		template<typename pixel_t>
		void encodePipeline(const Pipeline& codecs, Image <pixel_t>& image, std::vector<char>& buffer, int fillValue);

		template<typename pixel_t>
		void decodePipeline(const Pipeline& codecs, Image <pixel_t>& image, std::vector<char>& buffer, int fillValue);

		template<typename pixel_t>
		void encodeTransposeCodec(const ZarrCodec& codec, Image<pixel_t>& image, int fillValue)
		{
			Vec3c order;
			codec.getTransposeConfiguration(order);
			transpose(image, order, fillValue);
		}

		template<typename pixel_t>
		void decodeTransposeCodec(const ZarrCodec& codec, Image <pixel_t>& image, int fillValue)
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
			cout << "Using blosc compressor" << cname << endl;
			int numinternalthreads = 1;

			size_t realDestSize = blosc_compress_ctx(clevel, (int)shuffle, typesize, srcSize, buffer.data(), temp.data(), destSize, cname.c_str(), blocksize, numinternalthreads);

			if (realDestSize == 0) cout << "Buffer is incompressible.  Giving up." << endl;
			else if (realDestSize < 0) throw ITLException("Compression error.  Error code: " + toString(realDestSize));
			else cout << "Compression: " << srcSize << " -> " << realDestSize << " (" << static_cast<double>(srcSize) / realDestSize << "x)" << endl;

			buffer.resize(realDestSize);
			std::memcpy(buffer.data(), temp.data(), realDestSize);
		}

		//todo: use swapByteOrder(img); depending on endian
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

		template<typename pixel_t>
		void decodeBytesCodec(Image < pixel_t > &image, std::vector<char> & buffer)
		{
			Vec3c shape = image.dimensions();
			std::vector<pixel_t> temp(shape.product());
			std::memcpy(temp.data(), buffer.data(), buffer.size());

			assert(shape.product() * sizeof(pixel_t) == buffer.size());
			size_t n = 0;
			for (coord_t x = 0; x < shape.x; x++)
			{
				for (coord_t y = 0; y < shape.y; y++)
				{
					for (coord_t z = 0; z < shape.z; z++)
					{
						//todo: might be faster to read data sequentially
						image(x, y, z) = temp[n++];
						//cout << "set image(" << toString(Vec3c(x, y, z)) << ")=" << image(x, y, z) << endl;
					}
				}
			}
		}
		template<typename pixel_t>
		void encodeBytesCodec(const Image <pixel_t>& image, std::vector<char>& buffer)
		{
			Vec3c shape = image.dimensions();
			std::vector<pixel_t> temp(shape.product());

			size_t n = 0;
			for (coord_t x = 0; x < shape.x; x++)
			{
				for (coord_t y = 0; y < shape.y; y++)
				{
					for (coord_t z = 0; z < shape.z; z++)
					{
						temp[n++] = image(x, y, z);
					}
				}
			}

			size_t bufferSize = shape.product() * sizeof(pixel_t);
			buffer = std::vector<char>(bufferSize);
			std::memcpy(buffer.data(), temp.data(), buffer.size());
		}

		template<typename pixel_t>
		void decodeShardingCodec(const ZarrCodec& codec, Image < pixel_t > &shard, std::vector<char> & buffer, int fillValue)
		{
			Vec3c chunkShape;
			Pipeline codecs;
			Pipeline indexCodecs;
			sharding::indexLocation indexLocation;
			codec.getShardingConfiguration(chunkShape, codecs, indexCodecs, indexLocation);

		}

		template<typename pixel_t>
		void encodeShardingCodec(const ZarrCodec& codec, const Image <pixel_t>& shard, std::vector<char>& buffer, int fillValue)
		{
			Vec3c chunkShape;
			Pipeline codecs;
			Pipeline indexCodecs;
			sharding::indexLocation indexLocation;
			codec.getShardingConfiguration(chunkShape, codecs, indexCodecs, indexLocation);
			Vec3c chunksPerShard = shard.dimensions() / chunkShape;
			int chunkCount = chunksPerShard.product();
			std::vector<coord_t> shardIndexArray = {chunksPerShard.x, chunksPerShard.y, chunksPerShard.z, 2};
			std::vector<coord_t> chunkBytesList(chunkCount);
		}

		inline void decodeBytesBytesCodec(const ZarrCodec& codec, std::vector<char>& buffer)
		{
			assert(codec.type == codecs::Type::BytesBytesCodec);
			if (codec.name == codecs::Name::Blosc)
			{
				decodeBloscCodec(codec, buffer);
			}
			else throw ITLException("BytesBytesCodec: " + toString(codec.name) + " not yet implemented");
		}

		template<typename pixel_t>
		void decodeArrayBytesCodec(const ZarrCodec& codec, Image <pixel_t>& image, std::vector<char>& buffer, int fillValue)
		{
			assert(codec.type == codecs::Type::ArrayBytesCodec);
			if (codec.name == codecs::Name::Bytes)
			{
				decodeBytesCodec(image, buffer);
			}
			else if (codec.name == codecs::Name::Bytes)
			{
				decodeShardingCodec(codec, image, buffer, fillValue);
			}
			else throw ITLException("ArrayBytesCodec: " + toString(codec.name) + " not yet implemented");
		}

		template<typename pixel_t>
		void decodeArrayArrayCodec(const ZarrCodec& codec, Image <pixel_t>& image, int fillValue)
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
			else throw ITLException("BytesBytesCodec: " + toString(codec.name) + " not yet implemented");
		}

		template<typename pixel_t>
		void encodeArrayBytesCodec(const ZarrCodec& codec, Image <pixel_t>& image, std::vector<char>& buffer, int fillValue)
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
		void encodeArrayArrayCodec(const ZarrCodec& codec, Image <pixel_t>& image, int fillValue)
		{
			assert(codec.type == codecs::Type::ArrayArrayCodec);
			if (codec.name == codecs::Name::Transpose)
			{
				encodeTransposeCodec(codec, image, fillValue);
			}
			else throw ITLException("ArrayArrayCodec: " + toString(codec.name) + " not yet implemented");
		}

		template<typename pixel_t>
		void encodePipeline(const Pipeline& codecs, Image <pixel_t>& image, std::vector<char>& buffer, int fillValue){
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
		void decodePipeline(const Pipeline& codecs, Image <pixel_t>& image, std::vector<char>& buffer, int fillValue)
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
		cout << "chunkShape=" << chunkShape << endl;

		string reason;
		if (!fromJSON(codecs, this->configuration["codecs"], reason))
		{
			throw ITLException(reason);
		}
		if (!fromJSON(indexCodecs, this->configuration["index_codecs"], reason))
		{
			throw ITLException(reason);
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
		throw ITLException("error in reading transposeOrder of configuration: " + nlohmann::to_string(this->configuration) + " got exception: " + ex.what());
	}
}