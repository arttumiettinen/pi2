#pragma once

#include "line.h"

namespace itl2
{
	/**
	Represents capsule spanning between two points.
	*/
	template<typename T>
	class Capsule : public Line<T>
	{
	public:
		T radius;

		Capsule(const Vec3<T>& start, const Vec3<T>& end, T radius) :
			Line<T>(start, end),
			radius(radius)
		{
		}

		/**
		Tests if this capsule contains the given point.
		Assumes that all points x that satisfy |centerline - x| < r belong to the capsule.
		*/
		bool contains(const Vec3<T>& x) const
		{
			return Line<T>::distance(x) < radius;
		}
	};

	typedef Capsule<double> Capsuled;
	typedef Capsule<float32_t> Capsulef;
}
