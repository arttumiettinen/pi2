#pragma once

#include "math/vec3.h"
#include "math/matrix3x3.h"
#include "aabox.h"

namespace itl2
{
	/**
	Box in generic orientation (not axis-aligned).
	See also class AABox that represents an axis-aligned box.
	*/
	class Box
	{
	private:
		/**
		Position of the center of the box.
		*/
		Vec3d center;

		/**
		Half-width, height and depth of the box.
		*/
		Vec3d radius;

		/**
		Rotation matrix that rotates the box so that its axes are aligned with the coordinate axes.
		*/
		Matrix3x3d Rinv;

	public:

		/**
		Constructs box located at origin and whose semi-axis lengths are zeroes.
		*/
		Box() :
			center(0, 0, 0),
			radius(0, 0, 0),
			Rinv(Matrix3x3d::identity())
		{
		}

		/**
		Constructs axis-aligned box with given center points and radius.
		*/
		Box(const Vec3d& center, const Vec3d& radius) :
			center(center),
			radius(radius),
			Rinv(Matrix3x3d::identity())
		{
		}

		/**
		Constructs box located at the given center point, with given radius and
		orientation vectors.
		*/
		Box(const Vec3d& center, const Vec3d& radius,
			const Vec3d& u1, const Vec3d& u2) :
			center(center),
			radius(radius)
		{
			Rinv = Matrix3x3d::rotationMatrix(u1, u2);
			Rinv.transpose(); // Invert the rotation matrix
		}

		/**
		Tests if given point is inside this box.
		The test is exclusive, i.e. returns false at the edge of the box.
		*/
		bool contains(const Vec3d& p) const
		{
			// Transform to coordinates where the box is centered at the origin
			Vec3d pdot = p - center;

			// Rotate pdot such that box axes are aligned with coordinate axes
			pdot = Rinv * pdot;

			// See if pdot is inside box radius
			return NumberUtils<double>::lessThan(abs(pdot.x), radius.x) &&
				   NumberUtils<double>::lessThan(abs(pdot.y), radius.y) &&
				   NumberUtils<double>::lessThan(abs(pdot.z), radius.z);
		}

		/**
		Calculates the axis-aligned bounding box of this box.
		*/
		AABox<double> boundingBox() const
		{
			// Corners of the box in coordinates where the box is axis-aligned at the origin.
			Vec3d corners[] = {
				Vec3d(-radius.x, -radius.y, -radius.z),
				Vec3d( radius.x, -radius.y, -radius.z),
				Vec3d(-radius.x,  radius.y, -radius.z),
				Vec3d( radius.x,  radius.y, -radius.z),
				Vec3d(-radius.x, -radius.y,  radius.z),
				Vec3d( radius.x, -radius.y,  radius.z),
				Vec3d(-radius.x,  radius.y,  radius.z),
				Vec3d( radius.x,  radius.y,  radius.z)
			};

			// Rotate and translate the corners of the axis-aligned box to the correct orientation
			// and location
			Matrix3x3d R = Rinv;
			R.transpose();
			for (size_t n = 0; n < 8; n++)
				corners[n] = R * corners[n] + center;

			// Find componentwise bounds
			Vec3d minc = corners[0];
			Vec3d maxc = corners[0];
			for (size_t n = 1; n < 8; n++)
			{
				minc = min(minc, corners[n]);
				maxc = max(maxc, corners[n]);
			}

			return AABox<double>(minc, maxc);
		}
	};
}