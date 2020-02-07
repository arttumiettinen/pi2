#pragma once

#include "math/vec3.h"
#include "math/mathutils.h"
#include "math/matrix3x3.h"
#include "aabox.h"

namespace itl2
{

	/**
	Gets value of left side of ellipsoid equation at p for ellipsoid located at c that has
	semi-axis lengths l1, l2 and l3	and whose semi-axis orientations are given by (phiN, thetaN), N=1..3.
	Points inside the ellipsoid are characterized by return value <= 1.
	*/
	inline double getEllipsoidFunctionValue(const Vec3d& p,
								const Vec3d& c,
								double l1, double l2, double l3,
								double phi1, double theta1,
								double phi2, double theta2,
								double phi3, double theta3)
	{
		// Make the ellipsoid centered to origin
		Vec3d pdot = p - c;

		// Rotate pdot such that ellipsoid axes are aligned with coordinate axes
		Vec3d u1 = toCartesian(1.0, phi1, theta1);
		Vec3d u2 = toCartesian(1.0, phi2, theta2);
		Vec3d u3 = toCartesian(1.0, phi3, theta3);

		Matrix3x3d R(u1.x, u2.x, u3.x,
					 u1.y, u2.y, u3.y,
					 u1.z, u2.z, u3.z);

		R.transpose();
		pdot = R * pdot;

		//Matrix3x3d Ri;
		//R.inverse(Ri);
		//pdot = Ri * pdot;

		// Use ellipsoid equation
		double f = (pdot.x * pdot.x) / (l1 * l1) + (pdot.y * pdot.y) / (l2 * l2) + (pdot.z * pdot.z) / (l3 * l3);

		return f;
	}

	/**
	Test if point p is in ellipsoid located at c that has semi-axis lengths l1, l2 and l3
	and whose semi-axis orientations are given by (phiN, thetaN), N=1..3.
	*/
	inline bool isInEllipsoid(const Vec3d& p,
								const Vec3d& c,
								double l1, double l2, double l3,
								double phi1, double theta1,
								double phi2, double theta2,
								double phi3, double theta3)
	{
		double f = getEllipsoidFunctionValue(p, c, l1, l2, l3, phi1, theta1, phi2, theta2, phi3, theta3);
		return NumberUtils<double>::lessThan(f, 1);
	}

	/**
	Test if point p is in an ellipsoid located at c that has semi-axis lengths l1, l2 and l3
	and that can be transformed to axis-aligned ellipsoid by rotating with rotation matrix Rinv
	The test is exclusive.
	*/
	inline bool isInEllipsoid(const Vec3d& p,
								const Vec3d& c,
								double l1, double l2, double l3,
								const Matrix3x3d& Rinv)
	{
		// Make the ellipsoid centered to origin
		Vec3d pdot = p - c;

		// Rotate pdot such that ellipsoid axes are aligned with coordinate axes
		pdot = Rinv * pdot;

		// Use ellipsoid equation
		double f = (pdot.x * pdot.x) / (l1 * l1) + (pdot.y * pdot.y) / (l2 * l2) + (pdot.z * pdot.z) / (l3 * l3);

		return NumberUtils<double>::lessThan(f, 1);
	}


	class Ellipsoid
	{
	private:
		/**
		Position of the center of the ellipsoid.
		*/
		Vec3d center;

		/**
		Lengths of the three semi-axes.
		*/
		Vec3d semiAxisLengths;
		
		/**
		Rotation matrix that rotates the ellipsoid so that its axes are aligned with the coordinate axes.
		*/
		Matrix3x3d Rinv;

	public:

		/**
		Constructs ellipsoid located at origin and whose semi-axis lengths are zeroes.
		*/
		Ellipsoid() :
			center(0, 0, 0),
			semiAxisLengths(0, 0, 0),
			Rinv(Matrix3x3d::identity())
		{
		}

		/**
		Constructs axis-aligned ellipsoid with given center points and semi-axis lengths.
		*/
		Ellipsoid(const Vec3d& center, const Vec3d& axes) :
			center(center),
			semiAxisLengths(axes),
			Rinv(Matrix3x3d::identity())
		{
		}

		/**
		Constructs ellipsoid located at the given center point, with given semi-axis lengths and
		orientation vectors.
		*/
		Ellipsoid(const Vec3d& center, const Vec3d& axes,
			const Vec3d& u1, const Vec3d& u2) :
			center(center),
			semiAxisLengths(axes)
		{
			Rinv = Matrix3x3d::rotationMatrix(u1, u2);
			Rinv.transpose(); // Invert the rotation matrix
		}

		/**
		Constructs ellipsoid located at origin and whose semi-axis lengths are zeroes, and orientations given in spherical
		coordinates by (phiN, thetaN) pair for two first semi-axes. The orientation of the third axis is given by t1 x t3
		where t1 and t2 are the orientations of the first two axes.
		*/
		Ellipsoid(const Vec3d& center, const Vec3d& axes,
			double phi1, double theta1,
			double phi2, double theta2) :
			Ellipsoid(center, axes,
				toCartesian(1.0, phi1, theta1),
				toCartesian(1.0, phi2, theta2))
		{
		}

		/**
		Tests if given point is inside this ellipsoid.
		The test is exclusive, i.e. returns false at the edge of the ellipsoid.
		*/
		bool contains(const Vec3d& p) const
		{
			return isInEllipsoid(p, center, semiAxisLengths.x, semiAxisLengths.y, semiAxisLengths.z, Rinv);
		}

		/**
		Calculates bounding box of this ellipsoid.
		*/
		AABox<double> boundingBox() const
		{
			double m11 = Rinv.a00 * semiAxisLengths.x;
			double m12 = Rinv.a10 * semiAxisLengths.y;
			double m13 = Rinv.a20 * semiAxisLengths.z;
			double dx = sqrt(m11 * m11 + m12 * m12 + m13 * m13);


			double m21 = Rinv.a01 * semiAxisLengths.x;
			double m22 = Rinv.a11 * semiAxisLengths.y;
			double m23 = Rinv.a21 * semiAxisLengths.z;
			double dy = sqrt(m21 * m21 + m22 * m22 + m23 * m23);

			double m31 = Rinv.a02 * semiAxisLengths.x;
			double m32 = Rinv.a12 * semiAxisLengths.y;
			double m33 = Rinv.a22 * semiAxisLengths.z;
			double dz = sqrt(m31 * m31 + m32 * m32 + m33 * m33);

			Vec3d D(dx, dy, dz);
			return AABox<double>(center - D, center + D);
		}

	};

	namespace tests
	{
		void ellipsoid();
	}

}