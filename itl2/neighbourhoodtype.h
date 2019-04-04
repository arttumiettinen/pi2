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
}
