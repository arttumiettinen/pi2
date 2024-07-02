#pragma once

#include <string>

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
            ZarrCodecType type;
            std::string name = "invalid";
		};
        class ZarrBytesCodec: public ZarrCodec
        {
        public:
            ZarrCodecType type = ZarrCodecType::ArrayBytesCodec;
            std::string name = "bytes";
        };

	}

	template<>
	inline std::string toString(const zarr::ZarrCodec& codec)
	{
		if (codec.name == "invalid") throw ITLException("Invalid zarr codec.");
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