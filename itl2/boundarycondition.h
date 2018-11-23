#pragma once

namespace itl2
{
	/*
	Enumerates supported image boundary conditions.
	NOTE: If boundary conditions are added, at least these places must be checked:
	getNeighbourhood
	separable filtering functions
	interpolation
	*/
	enum BoundaryCondition
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
}
