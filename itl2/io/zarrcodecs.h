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
			Bytes
		};
	}

	template<>
	inline std::string toString(const zarr::ZarrCodecName& x)
	{
		switch (x)
		{
		case zarr::ZarrCodecName::Bytes:
			return "bytes";
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

			ZarrCodec(const ZarrCodecName name)
			{
				std::cout << "creating new zarrCodec" << std::endl;
				this->name = name;
				switch (name)
				{
				case ZarrCodecName::Bytes:
					this->type = ZarrCodecType::ArrayBytesCodec;
					readBytesCodecConfig();
					break;
				default:
					throw ITLException(std::string("Invalid zarr codec"));
				}
			}

			void readConfig(nlohmann::json config)
			{
				switch (this->name)
				{
				case ZarrCodecName::Bytes:
					readBytesCodecConfig(config);
					break;
				default:
					throw ITLException(std::string("Invalid zarr codec"));
				}
			}

			void readBytesCodecConfig(nlohmann::json config = nlohmann::json())
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
		};
		/**
			 * copy from raw::read
			Reads a .raw file to the given image, initializes the image to correct size read from file name.
			The file name must be in format image_name_100x200x300.raw or image_name_100x200.raw.
			@param img Image where the data is placed. The size of the image will be set based on the .raw file name.
			@param filename The name of the file to read.
			@param bytesToSkip Skip this many bytes from the beginning of the file.
			@param readPixel Function that reads one pixel. Relevant only for non-trivially copyable pixel data types.
			*/
		template<typename pixel_t, typename ReadPixel = decltype(raw::readPixel<pixel_t>)> void readBytesCodec(Image<pixel_t>& img, std::string filename, size_t bytesToSkip = 0, ReadPixel readPixel = raw::readPixel<pixel_t>)
		{
			raw::readNoParse<pixel_t, ReadPixel>(img, filename, bytesToSkip, readPixel);
		}
	}

	template<>
	inline std::string toString(const zarr::ZarrCodec& x)
	{
		return toString(x.name);
	}

	template<>
	inline zarr::ZarrCodec fromString(const std::string& str0)
	{
		std::string str = str0;
		toLower(str);
		zarr::ZarrCodecName name = fromString<zarr::ZarrCodecName>(str);
		return *new zarr::ZarrCodec(name);

	}
}