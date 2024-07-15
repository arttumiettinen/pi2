#pragma once

#include <string>
#include <sstream> // Ensure you include this for stringstream

#include "json.h"
#include "utilities.h"

namespace itl2
{
	namespace zarr
	{

		enum class ZarrCodecType
		{
			None,
			ArrayArrayCodec,
			ArrayBytesCodec,
			BytesBytesCodec,
		};

		enum class ZarrCodecName
		{
			None,
			Bytes,
			Transpose,
		};
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
				default:
					throw ITLException(std::string("Invalid zarr codec"));
				}
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
		};
		namespace internals
		{
			template<typename pixel_t, typename ReadPixel = decltype(raw::readPixel<pixel_t>)>
			void readBytesCodec(Image<pixel_t>& img, std::string filename, size_t bytesToSkip = 0, ReadPixel readPixel = raw::readPixel<pixel_t>)
			{
				std::ifstream in(filename.c_str(), std::ios_base::in | std::ios_base::binary);

				if (!in)
				{
					throw ITLException(std::string("Unable to open ") + filename + std::string(", ") + getStreamErrorMessage());
				}
				in.seekg(bytesToSkip, std::ios::beg);

				for (coord_t x = 0; x < img.width(); x++)
				{
					for (coord_t y = 0; y < img.height(); y++)
					{
						for (coord_t z = 0; z < img.depth(); z++)
						{
							readPixel(in, img(x, y, z));
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
				std::cout << "writeBlock fileStartPos=" << fileStartPos << " fileEndPos=" << fileEndPos << std::endl;
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
								std::cout << "writeBlock imgPos=" << imgPos << " x=" << x << " y=" << y << " z=" << z << std::endl;
								size_t linearIndex = img.getLinearIndex(imgPos);

								std::cout << "writeBlock linearIndex=" << linearIndex << " x=" << x << " y=" << y << " z=" << z << std::endl;
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