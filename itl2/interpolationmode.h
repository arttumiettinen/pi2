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

	

}