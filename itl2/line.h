#pragma once

#include "math/vec3.h"

using math::Vec3;

namespace itl2
{
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
	};
}