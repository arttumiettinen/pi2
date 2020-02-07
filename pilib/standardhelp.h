#pragma once

#include <string>

namespace itl2
{

	inline std::string connectivityHelp()
	{
		return "Can be Nearest for connectivity to nearest neighbours only, or All for connectivity to all neighbours.";
	}

	inline std::string interpolationHelp()
	{
		return "Can be Nearest for nearest neighbour interpolation, Linear for linear interpolation, or Cubic for cubic interpolation.";
	}

	inline std::string neighbourhoodTypeHelp()
	{
		return "Can be Ellipsoidal for ellipsoidal or spherical neighbourhood; or Rectangular for rectangular neighbourhood.";
	}

	inline std::string boundaryConditionHelp()
	{
		return "Zero indicates that values outside of image bounds are taken to be zero. Nearest indicates that the nearest value inside the image is to be used in place of values outside of image bounds.";
	}

	inline std::string rawFilenameFormatHelp()
	{
		return "Dimensions of the .raw file do not need to be specified if the file name is in format name_WxHxD.raw, where [W, H, D] are the dimensions of the image. "
			"The system tries to guess the pixel data type, too, based on the physical size and dimensions of the file as follows. "
			"If pixel size in bytes is 1, the system sets the pixel type to uint8. If pixel size in bytes is 2, the system sets the pixel type to uint16. "
			"If pixel size in bytes is 4, float32 pixel data is assumed (instead of e.g. int32 or uint32). "
			"If pixel size in bytes is 8, pixels are assumed to be of type uint64 (instead of e.g. int64 or complex32). "
			"If the guess is wrong, the pixel data type must be explicitly specified using the corresponding argument. "
			"Even in this case the dimensions can be read from the name of the file if the file name contains the dimensions.";
	}
}