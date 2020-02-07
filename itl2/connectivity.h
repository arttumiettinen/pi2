#pragma once

#include "utilities.h"

namespace itl2
{
	/*
	Defines possible connectivity modes (etc. for flood fill).
	*/
	enum class Connectivity
	{
		/**
		Declares connectivity where only nearest (in 2D 4, in 3D 6) neighbours are connected.
		*/
		NearestNeighbours,

		/**
		Declares connectivity where all neighbours are connected (in 2D 8, in 3D 27).
		*/
		AllNeighbours
	};

	template<>
	inline std::string toString(const Connectivity& x)
	{
		switch (x)
		{
		case Connectivity::NearestNeighbours: return "Nearest";
		case Connectivity::AllNeighbours: return "All";
		}
		throw ITLException("Invalid connectivity.");
	}

	template<>
	inline Connectivity fromString(const string& dt)
	{
		string str = dt;
		trim(str);
		toLower(str);
		if (str == "nearest" || str == "nearestneigbours" || str == "nearest_neighbours" ||
			str == "nearestneigbors" || str == "nearest_neighbors" ||
			str == "4" || str == "6" || str == "0")
			return Connectivity::NearestNeighbours;

		if (str == "all" || str == "allneigbours" || str == "all_neighbours" ||
			str == "allneigbors" || str == "all_neighbors" ||
			str == "8" || str == "27" || str == "1")
			return Connectivity::AllNeighbours;

		throw ITLException("Invalid connectivity: " + dt);
	}
}