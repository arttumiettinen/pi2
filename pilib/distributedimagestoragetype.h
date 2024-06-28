#pragma once

#include "itlexception.h"

namespace pilib
{
	enum class DistributedImageStorageType
	{
		NN5,
		Raw,
		Sequence,
		Zarr
	};

	inline string toString(const DistributedImageStorageType type)
	{
		switch (type)
		{
		case DistributedImageStorageType::NN5: return "NN5";
		case DistributedImageStorageType::Zarr: return "Zarr";
		case DistributedImageStorageType::Raw: return "Raw";
		case DistributedImageStorageType::Sequence: return "Sequence";
		default: throw itl2::ITLException("Unknown distributed image storage type.");
		}
	}
}