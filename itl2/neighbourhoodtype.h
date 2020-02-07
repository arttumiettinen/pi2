#pragma once

#include "utilities.h"

namespace itl2
{
	/*
	Enumerates possible neighbourhood types.
	*/
	enum class NeighbourhoodType
	{
		/**
		Indicates rectangular neighbourhood with different radius in all dimensions.
		The size (pixel count) of the rectangular neighbourhood is 2.0 * radius + 1.0.
		*/
		Rectangular,
		/**
		Indicates ellipsoidal neighbourhood with different radius in all dimensions.
		*/
		Ellipsoidal
	};

	template<>
	inline string toString(const NeighbourhoodType& x)
	{
		switch (x)
		{
		case NeighbourhoodType::Rectangular: return "Rectangular";
		case NeighbourhoodType::Ellipsoidal: return "Ellipsoidal";
		}
		throw ITLException("Invalid neighbourhood type value.");
	}

	template<>
	inline NeighbourhoodType fromString(const string& dt)
	{
		string dt2 = dt;
		trim(dt2);
		toLower(dt2);
		if (dt2 == "rectangular" || dt2 == "rect" || dt2 == "box")
			return NeighbourhoodType::Rectangular;

		if (dt2 == "ell" || dt2 == "ellipsoidal" || dt2 == "ellipsoid" || dt2 == "spherical" || dt2 == "sphere" || dt2 == "sph")
			return NeighbourhoodType::Ellipsoidal;

		throw ITLException("Invalid neighbourhood type string: " + dt);
	}
}
