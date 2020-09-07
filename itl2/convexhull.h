#pragma once

#include <vector>

namespace itl2
{

	/**
	Returns true if 'c' is left of line 'a'-'b'.
	*/
	template<typename POINT> int left(const POINT& a, const POINT& b, const POINT& c)
	{
		double value = (c.x - a.x) * (b.y - a.y) - (c.y - a.y) * (b.x - a.x);
		if (value < -1e-10)
			return -1;
		else if (value > 1e-10)
			return 1;
		else
			return 0;
	}

	/**
	Returns true if 'a' is more lower-left than 'b'.
	*/
	template<typename POINT> bool cmppt(const POINT& a, const POINT& b)
	{

		if (a[0] < b[0])
			return true;
		if (a[0] > b[0])
			return false;
		if (a[1] < b[1])
			return true;
		if (a[1] > b[1])
			return false;
		return false;

		// 1|0  
		//--b--
		// 1|1  
		//if(a[0] >= b[0] && a[1] >= b[1])
		//	return false;
		return true;
	}

	/**
	Calculates convex hull of the inputs points and outputs its vertices to out vector.
	Uses gift-wrapping algorithm.
	*/
	template<typename POINT> void convexHull2D(const std::vector<POINT>& points, std::vector<POINT>& out)
	{
		out.clear();

		// Two points form convex hull in all cases. (except if they are the same)
		if (points.size() <= 2)
		{
			for (size_t i = 0; i < points.size(); ++i)
				out.push_back(points[i]);
			return;
		}

		// Find lower-leftmost point.
		size_t pointOnHull = 0;
		for (size_t i = 1; i < points.size(); ++i)
		{
			if (cmppt(points[i], points[pointOnHull]))
			{
				pointOnHull = i;
			}
		}

		// Gift wrap
		size_t out0 = pointOnHull;
		do
		{
			out.push_back(points[pointOnHull]);

			size_t endpointCandidate = 0;
			for (size_t j = 1; j < points.size(); ++j)
			{
				if (pointOnHull != endpointCandidate) // If the endPointCandidate is the last point on the hull, there's no line that point j could be compared against.
				{
					if (j != pointOnHull)
					{
						int pos = left(points[pointOnHull], points[endpointCandidate], points[j]);
						if (pos < 0)
						{
							endpointCandidate = j;
						}
						else if (pos == 0)
						{
							if ((points[pointOnHull] - points[endpointCandidate]).normSquared() < (points[pointOnHull] - points[j]).normSquared())
								endpointCandidate = j;
						}
					}
				}
				else
				{
					endpointCandidate = j;
				}
			}

			pointOnHull = endpointCandidate;
		} while (pointOnHull != out0);

	}


}

