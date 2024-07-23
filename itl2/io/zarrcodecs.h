#pragma once

#include <string>
#include <sstream> // Ensure you include this for stringstream

#include "json.h"
#include "utilities.h"
#include "zarrimagedatawrapper.h"

namespace itl2
{
	using std::cout, std::endl;
	namespace zarr
	{

		enum class ZarrCodecType
		{
			ArrayArrayCodec,
			ArrayBytesCodec,
			BytesBytesCodec,
		};

		enum class ZarrCodecName
		{
			Bytes,
			Transpose,
			Blosc,
		};
		//TODO save codec in struct instead of json
	}

	template<>
	inline std::string toString(const zarr::ZarrCodecName& x)
	{
		switch (x)
		{
		case zarr::ZarrCodecName::Bytes:
			return "bytes";
		case zarr::ZarrCodecName::Transpose:
			return "transpose";
		case zarr::ZarrCodecName::Blosc:
			return "blosc";
		}
		throw ITLException("Invalid ZarrCodecName.");
	}

	template<>
	inline zarr::ZarrCodecName fromString(const std::string& str0)
	{
		std::string str = str0;
		toLower(str);
		if (str == "bytes")
			return zarr::ZarrCodecName::Bytes;
		if (str == "transpose")
			return zarr::ZarrCodecName::Transpose;
		if(str == "blosc")
			return zarr::ZarrCodecName::Blosc;

		throw ITLException(std::string("Invalid zarr codec name: ") + str);
	}

	namespace zarr
	{
		class ZarrCodec
		{
		 public:

			ZarrCodecType type;
			ZarrCodecName name;
			nlohmann::json configuration;

			ZarrCodec(const ZarrCodecName name, nlohmann::json config = nlohmann::json())
			{
				this->name = name;
				switch (name)
				{
				case ZarrCodecName::Bytes:
					this->type = ZarrCodecType::ArrayBytesCodec;
					parseBytesCodecConfig(config);
					break;
				case ZarrCodecName::Transpose:
					this->type = ZarrCodecType::ArrayArrayCodec;
					parseTransposeCodecConfig(config);
					break;
				case ZarrCodecName::Blosc:
					this->type = ZarrCodecType::BytesBytesCodec;
					parseBloscCodecConfig(config);
					break;
				default:
					throw ITLException(std::string("Invalid zarr codec"));
				}
				cout << "created codec: " << toString(this->name) << " " << (int)this->type << endl;
			}

			void parseBloscCodecConfig(nlohmann::json config = nlohmann::json())
			{
				this->configuration = config;
			}


			void parseTransposeCodecConfig(nlohmann::json config = nlohmann::json())
			{
				this->configuration = config;
			}

			nlohmann::json toJSON() const
			{
				nlohmann::json j;
				j["name"] = toString(this->name);
				j["configuration"] = configuration; //this does not return by reference to the configuration just a copy, right?
				return j;
			}

			bool operator==(const ZarrCodec& t) const
			{
				return toJSON() == t.toJSON();
			}

			// BytesCodec
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

			Vec3c transposeOrder(){
				if (this->name!=ZarrCodecName::Transpose) throw ITLException("only transpose codec has transpose order");
				auto orderJSON = this->configuration["order"];
				//TODO check if valid order
				Vec3c order = Vec3c(0, 1, 2);
				order[0] = orderJSON[0].get<size_t>();
				if (orderJSON.size() >= 2)
					order[1] = orderJSON[1].get<size_t>();
				if (orderJSON.size() >= 3)
					order[2] = orderJSON[2].get<size_t>();
				cout << "transposeOrder=" << order << endl;
				return order;
			}
		};

		namespace internals
		{
			template<typename pixel_t>
			std::vector<pixel_t> readAllBytes(std::string filename)
			{
				std::ifstream ifs(filename, std::ios_base::binary|std::ios::ate);
				if (!ifs)
				{
					throw ITLException(std::string("Unable to open ") + filename + std::string(", ") + getStreamErrorMessage());
				}
				std::ifstream::pos_type pos = ifs.tellg();

				std::vector<pixel_t>  result(pos);

				ifs.seekg(0, std::ios::beg);
				ifs.read((char*)&result[0], pos);

				return result;
			}

			template<typename pixel_t>
			void readBytesCodec(Image<pixel_t>& image, std::vector<pixel_t>& buffer)
			{
				Vec3c shape = image.dimensions();
				pixel_t data[buffer.size()];
				std::copy(buffer.begin(), buffer.end(), data);

				size_t n = 0;
				for (coord_t x = 0; x < shape.x; x++)
				{
					for (coord_t y = 0; y < shape.y; y++)
					{
						for (coord_t z = 0; z < shape.z; z++)
						{
						  	//todo: might be faster to read data sequentially
							image(x, y, z) = data[n++];
							cout << "set image(" << toString(Vec3c(x, y, z)) << ")=" <<image(x, y, z)<<endl;
						}
					}
				}
			}
			template<typename pixel_t>
			void writeBytesCodecBlock(const Image<pixel_t>& img, const std::string& filename,
				const Vec3c& filePosition, const Vec3c& fileDimensions,
				const Vec3c& imagePosition,
				const Vec3c& blockDimensions,
				bool showProgressInfo = false)
			{
				Vec3c fileStartPos = filePosition;
				clamp(fileStartPos, Vec3c(0, 0, 0), fileDimensions);
				Vec3c fileEndPos = filePosition + blockDimensions;
				clamp(fileEndPos, Vec3c(0, 0, 0), fileDimensions);
				if (!img.isInImage(imagePosition))
					throw ITLException("Block start position must be inside the image.");
				if (!img.isInImage(imagePosition + blockDimensions - Vec3c(1, 1, 1)))
					throw ITLException("Block end position must be inside the image.");

				createFoldersFor(filename);

				// Create file if it does not exist, otherwise set file size to the correct value.
				setFileSize(filename, fileDimensions.x * fileDimensions.y * fileDimensions.z * sizeof(pixel_t));

				std::ofstream out(filename.c_str(), std::ios_base::in | std::ios_base::out | std::ios_base::binary);

				if (!out)
					throw ITLException(std::string("Unable to open ") + filename + std::string(", ") + getStreamErrorMessage());

				const pixel_t* pBuffer = img.getData();

				{
					ProgressIndicator prog(fileEndPos.z - fileStartPos.z, showProgressInfo);

					for (coord_t x = fileStartPos.x; x < fileEndPos.x; x++)
					{
						for (coord_t y = fileStartPos.y; y < fileEndPos.y; y++)
						{
							for (coord_t z = fileStartPos.z; z < fileEndPos.z; z++)
							{
								Vec3c imgPos = Vec3c(x, y, z) - fileStartPos + imagePosition;
								size_t linearIndex = img.getLinearIndex(imgPos);

								//size_t filePos = ((z * fileDimensions.x * fileDimensions.y) + (y * fileDimensions.x) + x) * sizeof(pixel_t);
								//out.seekp(filePos);

								if (!out)
									throw ITLException(std::string("Seek failed for file ") + filename + std::string(", ") + getStreamErrorMessage());

								//todo: might be faster to read pBuffer sequentially
								out.write((char*)&pBuffer[linearIndex], sizeof(pixel_t));

								if (!out)
									throw ITLException(std::string("Unable to write to ") + filename + std::string(", ") + getStreamErrorMessage());
							}
						}
						prog.step();
					}
				}
			}
		}

	}
}