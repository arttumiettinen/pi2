#pragma once

#include "math/vec3.h"
#include "math/numberutils.h"
#include "test.h"

namespace itl2
{
	/**
	Find nearest point of some point p on line segment from l1 to l2.
	Assumes l1 != l2.
	*/
	template<typename T> Vec3<T> nearestPointOnLineSegment(const Vec3<T>& l1, const Vec3<T>& l2, const Vec3<T> &p)
	{
		const Vec3<T> d = l2 - l1;
		T t = (p - l1).dot(d) / d.dot(d);
		if (t < 0)
			return l1;
		if (t > 1)
			return l2;
		return l1 + d * t;
	}

	/**
	Calculates distance between p and line segment from l1 to l2.
	*/
	template<typename T> double pointLineSegmentDistance(const Vec3<T>& l1, const Vec3<T>& l2, const Vec3<T> &p)
	{
		return sqrt(pointLineSegmentDistanceSquared(l1, l2, p));
	}

	/**
	Calculates squared distance between p and line segment from l1 to l2.
	*/
	template<typename T> double pointLineSegmentDistanceSquared(const Vec3<T>& l1, const Vec3<T>& l2, const Vec3<T> &p)
	{
		return (p - nearestPointOnLineSegment(l1, l2, p)).normSquared();
	}

	/**
	Find nearest point of some point p on infinite line that goes through points l1 and l2.
	Assumes l1 != l2.
	*/
	template<typename T> Vec3<T> nearestPointOnLine(const Vec3<T>& l1, const Vec3<T>& l2, const Vec3<T> &p)
	{
		const Vec3<T> d = l1 - l2;
		return l1 + d * (p - l1).dot(d) / d.dot(d);
	}

	/**
	Calculates distance between p and infinite line that goes through points l1 and l2.
	*/
	template<typename T> double pointLineDistance(const Vec3<T>& l1, const Vec3<T>& l2, const Vec3<T> &p)
	{
		return sqrt(pointLineDistanceSquared(l1, l2, p));
	}

	/**
	Calculates squared distance between p and infinite line that goes through points l1 and l2.
	*/
	template<typename T> double pointLineDistanceSquared(const Vec3<T>& l1, const Vec3<T>& l2, const Vec3<T> &p)
	{
		return (p - nearestPointOnLine(l1, l2, p)).normSquared();
	}

	/// <summary>
	/// Perform one ray slab test step.
	/// </summary>
	/// <param name="start">Starting coordinate for a ray.</param>
	/// <param name="vel">Velocity of ray. (Direction component normally)</param>
	/// <param name="min">Minimum coordinate for the slab interval.</param>
	/// <param name="max">Maximum coordinate for the slab interval.</param>
	/// <param name="tmin">Minimal time value from previous tests. This is replaced with new minimum value.</param>
	/// <param name="tmax">Maximal time value from previous tests. This is replaced with new maximum value.</param>
	/// <returns>true if ray intersects the slab interval; false otherwise. (false=no ray-box intersection.)</returns>
	template<typename real_t> bool raySlab(const real_t start, const real_t vel, const real_t min, const real_t max, real_t& tmin, real_t& tmax)
	{
		// See ClipLine below.

		// If velocity is zero, test if the starting point is in the interval.
		if (NumberUtils<real_t>::equals(vel, 0))
		{
			return min <= start && start <= max;
		}

		real_t newtmin = (min - start) / vel;
		real_t newtmax = (max - start) / vel;

		if (newtmin > newtmax)
			std::swap(newtmin, newtmax);

		if (newtmax < tmin || newtmin > tmax)
			return false;

		if (newtmin > tmin)
			tmin = newtmin;
		if (newtmax < tmax)
			tmax = newtmax;

		return true;
	}

	/// <summary>
	/// Clip 3D line.
	/// Uses ray slab algorithm. Discussed in
	/// http://www.gamedev.net/community/forums/topic.asp?topic_id=429443 and
	/// http://www.siggraph.org/education/materials/HyperGraph/raytrace/rtinter3.htm.
	/// </summary>
	/// <param name="lineStart">Starting point of the line.</param>
	/// <param name="lineEnd">Ending point of the line.</param>
	/// <param name="clipMin">Minimum x- and y-values for the clip rectangle.</param>
	/// <param name="clipMax">Maximum x- and y-values for the clip rectangle.</param>
	/// <returns>true if the line intersects the clipping rectangle, false otherwise.</returns>
	template<typename real_t> bool clipLine(Vec3<real_t>& lineStart, Vec3<real_t>& lineEnd, const Vec3<real_t>& clipMin, const Vec3<real_t>& clipMax)
	{
		Vec3<real_t> p0 = lineStart;
		Vec3<real_t> dir = (lineEnd - lineStart).normalized();

		Vec3<real_t> halfDimensions = (clipMax - clipMin) / 2.0;

		Vec3<real_t> b0 = (clipMin + clipMax) / 2.0;

		Vec3<real_t> b1(1, 0, 0);

		real_t tmin = 0.0;
		real_t tmax = (lineEnd - lineStart).norm();

		if (!raySlab(p0.dot(b1), dir.dot(b1), b0.dot(b1) - halfDimensions.x, b0.dot(b1) + halfDimensions.x, tmin, tmax))
			return false;

		Vec3<real_t> b2(0, 1, 0);

		if (!raySlab(p0.dot(b2), dir.dot(b2), b0.dot(b2) - halfDimensions.y, b0.dot(b2) + halfDimensions.y, tmin, tmax))
			return false;

		Vec3<real_t> b3(0, 0, 1);

		if (!raySlab(p0.dot(b3), dir.dot(b3), b0.dot(b3) - halfDimensions.z, b0.dot(b3) + halfDimensions.z, tmin, tmax))
			return false;

		//distance = tmin;
		lineStart = p0 + tmin * dir;
		lineEnd = p0 + tmax * dir;
		return true;
	}

	namespace tests
	{
		inline void geometry()
		{
			Vec3d start, end;

			start = Vec3d(-1, -1, -1);
			end = Vec3d(2, 2, 2);
			itl2::testAssert(clipLine(start, end, Vec3d(0, 0, 0), Vec3d(1, 1, 1)) == true, "clipLine return value");
			itl2::testAssert(start == Vec3d(0, 0, 0), "start point");
			itl2::testAssert(end == Vec3d(1, 1, 1), "end point");

			double d1 = pointLineSegmentDistance(start, end, start);
			double d2 = pointLineSegmentDistance(start, end, end);
			double d3 = pointLineSegmentDistance(start, end, (start + end) / 2);
			double d4 = pointLineSegmentDistance(start, end, Vec3d(-2, -2, -2));
			double d4_true = (start - Vec3d(-2, -2, -2)).norm();
			itl2::testAssert(d1 == 0, "point line segment distance 1");
			itl2::testAssert(d2 == 0, "point line segment distance 2");
			itl2::testAssert(d3 == 0, "point line segment distance 3");
			itl2::testAssert(NumberUtils<double>::equals(d4, d4_true), "point line segment distance 4");

		}
	}
}
