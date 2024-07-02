#pragma once

#include <string>

#include "json.h"
#include "utilities.h"

namespace itl2
{
	namespace zarr
	{
		/**
		Enumerates compression methods supported by the NN5 file format.
		*/
		enum class ZarrCodecType{
            ArrayArrayCodec,
            ArrayBytesCodec,
            BytesBytesCodec
        };

		class ZarrCodec
		{
        public:
            static const ZarrCodecType type;
            static const std::string name;
            void readConfig(nlohmann::json config);
            nlohmann::json toJSON() const;
            bool operator==(const ZarrCodec &t) const
            {
                return toJSON() == t.toJSON();
            }
        };

        class ZarrBytesCodec: public ZarrCodec
        {
        public:
            ZarrCodecType type = ZarrCodecType::ArrayBytesCodec;
            std::string name = "bytes";
            std::string endian;
            ZarrBytesCodec()
            {
                endian = "little";
            }
            void readConfig(nlohmann::json config)
            {
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

            }
            nlohmann::json toJSON()
            {
                nlohmann::json j;
                j["name"] = name;
                j["configuration"]["endian"] = endian;
                return j;
            }

            bool operator==(ZarrCodec &t)
            {
                if (t.name != name)
                {
                    return false;
                }
                ZarrBytesCodec t1 = (ZarrBytesCodec&)t;
                return t1.endian == endian;
            }
        };
	}

	template<>
	inline std::string toString(const zarr::ZarrCodec& codec)
	{
		if (codec.name == "") throw ITLException("Invalid zarr codec.");
        return codec.name;
	}

	template<>
	inline zarr::ZarrCodec fromString(const string& str0)
	{
		string str = str0;
		toLower(str);
		if (str == "bytes")
			return *new zarr::ZarrBytesCodec();

		throw ITLException(string("Invalid zarr codec: ") + str);
	}
}