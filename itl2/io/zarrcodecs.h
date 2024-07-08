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
    inline std::string toString(const zarr::ZarrCodecName& x) {
        switch (x) {
            case zarr::ZarrCodecName::Bytes: return "bytes";
        }
        throw ITLException("Invalid ZarrCodecName.");
    }

    template<>
    inline zarr::ZarrCodecName fromString(const std::string& str0) {
        std::string str = str0;
        toLower(str);
        if (str == "bytes")
            return zarr::ZarrCodecName::Bytes;

        throw ITLException(std::string("Invalid zarr codec name: ") + str);
    }

    namespace zarr {

        class ZarrCodec {
            // TODO: make this an abstract class
        public:
            static const ZarrCodecType type = ZarrCodecType::None;
            static const ZarrCodecName name = ZarrCodecName::None;

            void readConfig(nlohmann::json config) {}

            nlohmann::json toJSON() const {
                nlohmann::json j;
                j["name"] = toString(this->name);
                return j;
            }

            bool operator==(const ZarrCodec &t) const {
                return toJSON() == t.toJSON();
            }
        };

        class ZarrBytesCodec: public ZarrCodec {
        public:
            static const ZarrCodecType type = ZarrCodecType::ArrayBytesCodec;
            static const ZarrCodecName name = ZarrCodecName::Bytes;
            std::string endian;

            ZarrBytesCodec() {
                endian = "little";
            }

            void readConfig(nlohmann::json config) {
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
            }

            nlohmann::json toJSON() {
                nlohmann::json j = ZarrCodec::toJSON();
                j["configuration"]["endian"] = endian;
                return j;
            }
        };
    }

    template<>
    inline std::string toString(const zarr::ZarrCodec& x) {
        return toString(x.name);
    }

    template<>
    inline zarr::ZarrCodec fromString(const std::string& str0) {
        std::string str = str0;
        toLower(str);
        if (str == "bytes")
            return *new zarr::ZarrBytesCodec();

        throw ITLException(std::string("Invalid zarr codec: ") + str);
    }
}