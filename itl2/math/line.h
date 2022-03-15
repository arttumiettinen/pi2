#pragma once

#include "math/vec3.h"
#include "math/geometry.h"

namespace itl2
{
	/**
	Represents line between two points.
	*/
	template<typename T>
	class Line
	{
	public:
		Vec3<T> start;
		Vec3<T> end;

		Line(const Vec3<T>& start, const Vec3<T>& end) :
			start(start),
			end(end)
		{
		}

		/**
		Calculates distance between this line and the given point.
		*/
		T distance(const Vec3<T>& x) const
		{
			return (T)pointLineSegmentDistance(start, end, x);
		}
	};
}