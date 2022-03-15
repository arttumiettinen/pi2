#pragma once

#include <string>

#include "utilities.h"

namespace itl2
{
	namespace nn5
	{
		/**
		Enumerates compression methods supported by the NN5 file format.
		*/
		enum class NN5Compression
		{
			/**
			Uncompressed .raw format.
			*/
			Raw,
			/**
			LZ4 compressed format.
			*/
			LZ4
		};
	}

	template<>
	inline std::string toString(const nn5::NN5Compression& x)
	{
		switch (x)
		{
		case nn5::NN5Compression::Raw: return "Raw";
		case nn5::NN5Compression::LZ4: return "LZ4Raw";
		}
		throw ITLException("Invalid nn5 compression type.");
	}

	template<>
	inline nn5::NN5Compression fromString(const string& str0)
	{
		string str = str0;
		toLower(str);
		if (str == "raw")
			return nn5::NN5Compression::Raw;
		if (str == "lz4raw")
			return nn5::NN5Compression::LZ4;

		throw ITLException(string("Invalid nn5 compression type: ") + str);
	}
}