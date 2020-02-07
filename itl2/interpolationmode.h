#pragma once

#include "utilities.h"

namespace itl2
{

	/*
	Enumerates supported interpolation modes.
	*/
	enum class InterpolationMode
	{
		/**
		Nearest neighbour interpolation.
		*/
		Nearest,
		/**
		Linear interpolation.
		*/
		Linear,
		/**
		Cubic interpolation
		*/
		Cubic
	};

	template<>
	inline std::string toString(const InterpolationMode& x)
	{
		switch (x)
		{
		case InterpolationMode::Nearest: return "Nearest neighbour";
		case InterpolationMode::Linear: return "Linear";
		case InterpolationMode::Cubic: return "Cubic";
		}
		throw ITLException("Invalid boundary condition.");
	}

	template<>
	inline InterpolationMode fromString(const string& dt)
	{
		string str = dt;
		trim(str);
		toLower(str);
		if (startsWith(str, "nearest")
			|| str == "0" || str == "no" || str == "none" || str == "off")
			return InterpolationMode::Nearest;

		if (startsWith(str, "linear")
			|| str == "1")
			return InterpolationMode::Linear;

		if (startsWith(str, "cubic")
			|| str == "2")
			return InterpolationMode::Cubic;

		throw ITLException("Invalid connectivity: " + dt);
	}

}