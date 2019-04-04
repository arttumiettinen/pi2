#pragma once

#include "utilities.h"

namespace itl2
{
	/*
	Enumerates supported image boundary conditions.
	NOTE: If boundary conditions are added, at least these places must be checked:
	getNeighbourhood
	separable filtering functions
	interpolation
	*/
	enum class BoundaryCondition
	{
		/**
		Boundary condition where zero is returned for locations outside the image.
		*/
		Zero,
		/**
		Boundary condition where the nearest value inside the image is used instead of out-of-bounds value.
		*/
		Nearest
	};

	template<>
	inline string toString(const BoundaryCondition& x)
	{
		switch (x)
		{
		case BoundaryCondition::Zero: return "BoundaryCondition::Zero";
		case BoundaryCondition::Nearest: return "Nearest";
		}
		throw ITLException("Invalid boundary condition.");
	}
}
