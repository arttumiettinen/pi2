#pragma once

#include <string>
#include <sstream> // Ensure you include this for stringstream

#include "json.h"
#include "utilities.h"

namespace itl2 {
    namespace zarr {

        enum class ZarrCodecType {
            None,
            ArrayArrayCodec,
            ArrayBytesCodec,
            BytesBytesCodec,
        };

        enum class ZarrCodecName {
            None,
            Bytes
        };
    }

    template<>
    inline std::string toString(const zarr::ZarrCodecName &x) {
        switch (x) {
            case zarr::ZarrCodecName::Bytes:
                return "bytes";
        }
        throw ITLException("Invalid ZarrCodecName.");
    }

    template<>
    inline zarr::ZarrCodecName fromString(const std::string &str0) {
        std::string str = str0;
        toLower(str);
        if (str == "bytes")
            return zarr::ZarrCodecName::Bytes;

        throw ITLException(std::string("Invalid zarr codec name: ") + str);
    }

    namespace zarr {

        class ZarrCodec {
        public:
            ZarrCodecType type;
            ZarrCodecName name;
            nlohmann::json configuration;

            ZarrCodec(const ZarrCodecName name, const nlohmann::json configuration) {
                this->name = name;
                this->configuration = configuration;
                switch (name) {
                    case ZarrCodecName::Bytes:
                        this->type = ZarrCodecType::ArrayBytesCodec;

                        break;
                    default:
                        throw ITLException(std::string("Invalid zarr codec"));

                }
            }

            void readBytesCodecConfig(nlohmann::json config) {
                std::string endian = "little";
                for (auto it = config.begin(); it != config.end(); ++it) {
                    if (it.key() == "endian") {
                        endian = it.value();
                        if (endian != "little" && endian != "big") {
                            throw ITLException("Invalid endian in bytes codec config: " + endian);
                        }
                    } else {
                        throw ITLException("Invalid key in bytes codec config: " + it.key());
                    }
                }
                //TODO
            }

            nlohmann::json toJSON() const {
                nlohmann::json j;
                j["name"] = toString(this->name);
                j["configuration"] = configuration; //this does not return the reference to the configuration just a copy, right?
                return j;
            }

            bool operator==(const ZarrCodec &t) const {
                return toJSON() == t.toJSON();
            }
        };
    }

    template<>
    inline std::string toString(const zarr::ZarrCodec &x) {
        return toString(x.name);
    }

    template<>
    inline zarr::ZarrCodec fromString(const std::string &str0) {
        std::string str = str0;
        toLower(str);
        if (str == "bytes")
            return *new zarr::ZarrCodec(zarr::ZarrCodecName::Bytes, "");

        throw ITLException(std::string("Invalid zarr codec: ") + str);
    }
}