#pragma once

#include "math/vec3.h"
#include "math/matrix3x3.h"
#include "math/numberutils.h"

#include <fstream>
#include <vector>

namespace itl2
{
	
	template<typename T> class Sphere
	{
	private:

		static double radiusEpsilon()
		{
			return 1.0e-6;
		}

		void init(const math::Vec3<T>& p1, const math::Vec3<T>& p2)
		{
			center = (p1 + p2) / 2;
			radius = (p1 - p2).norm() / 2 + radiusEpsilon();
		}



	public:
		
		/**
		Center of the sphere.
		*/
		math::Vec3<T> center;
		
		/**
		Radius of the sphere.
		*/
		T radius;


		/**
		Constructs sphere whose position = 0 and radius = 0.
		*/
		Sphere() :
			center(0, 0, 0),
			radius(0)
		{
		}

		/**
		Copy constructor
		*/
		Sphere(const Sphere& s) :
			center(s.center),
			radius(s.radius)
		{
		}

		/**
		Constructs a sphere whose position is the given point and radius = 0.
		*/
		Sphere(const math::Vec3<T>& center) :
			center(center),
			radius(0 + radiusEpsilon())
		{
		}

		/**
		Constructs a sphere located at the given point and having the given radius.
		*/
		Sphere(const math::Vec3<T>& center, T radius) :
			center(center),
			radius(radius)
		{
		}

		/**
		Constructs a sphere through the two given points.
		*/
		Sphere(const math::Vec3<T>& p1, const math::Vec3<T>& p2)
		{
			init(p1, p2);
		}

		/**
		Constructs a sphere through the three given points.
		*/
		Sphere(const math::Vec3<T>& p1, const math::Vec3<T>& p2, const math::Vec3<T>& p3)
		{
			math::Vec3<T> a = p2 - p1;
			math::Vec3<T> b = p3 - p1;

			T denom = 2 * (a.cross(b)).dot(a.cross(b));

			if(!math::NumberUtils<T>::equals(denom, 0))
			{
				math::Vec3<T> o = ((b.normSquared()) * ((a.cross(b)).cross(a)) +
							(a.normSquared()) * (b.cross(a.cross(b)))) / denom;

				radius = o.norm() + radiusEpsilon();
				center = p1 + o;
			}
			else
			{
				// The three points are on the same line. Take the pair of points with largest distance
				// between them and construct the sphere based on the pair.
				
				T p1p2 = (p1 - p2).norm();
				T p2p3 = (p2 - p3).norm();
				T p1p3 = (p1 - p3).norm();

				if(p1p2 >= p1p3 && p1p2 >= p2p3)
				{
					init(p1, p2);
				}
				else if(p2p3 >= p1p2 && p2p3 >= p1p3)
				{
					init(p2, p3);
				}
				else
				{
					init(p1, p3);
				}
			}
		}

		/**
		Constructs a sphere through four given points.
		*/
		Sphere(const math::Vec3<T>& p1, const math::Vec3<T>& p2, const math::Vec3<T>& p3, const math::Vec3<T>& p4)
		{
			math::Vec3<T> a = p2 - p1;
			math::Vec3<T> b = p3 - p1;
			math::Vec3<T> c = p4 - p1;

			math::Matrix3x3<T> A(a.x, a.y, a.z,
						   b.x, b.y, b.z,
						   c.x, c.y, c.z);

			T denom = 2.0f * A.det();

			if(!math::NumberUtils<T>::equals(denom, 0.0))
			{
				math::Vec3<T> o = ((c.normSquared()) * (a.cross(b)) +
							(b.normSquared()) * (c.cross(a)) +
							(a.normSquared()) * (b.cross(c))) / denom;

				radius = o.norm() + radiusEpsilon();
				center = p1 + o;
			}
			else
			{
				// The four points are in the same plane or three of them lie on the same line.
				// Try all spheres through the three points.
				Sphere s1(p1, p2, p3);
				Sphere s2(p2, p3, p4);
				Sphere s3(p3, p4, p1);
				Sphere s4(p4, p1, p2);

				if(s1.radius <= s2.radius && s1.radius <= s3.radius && s1.radius <= s4.radius &&
					s1.containsInclusive(p4))
				{
					*this = s1;
				}
				else if(s2.radius <= s1.radius && s2.radius <= s3.radius && s2.radius <= s4.radius &&
					s2.containsInclusive(p1))
				{
					*this = s2;
				}
				else if(s3.radius <= s1.radius && s3.radius <= s2.radius && s3.radius <= s4.radius &&
					s3.containsInclusive(p2))
				{
					*this = s3;
				}
				else if(s4.radius <= s1.radius && s4.radius <= s2.radius && s4.radius <= s3.radius &&
					s4.containsInclusive(p3))
				{
					*this = s4;
				}
				else
				{
					// All four points are on the edge of a circle in the plane where the points lie.
					// TODO: Does this ever happen?
					center = (p1 + p2 + p3 + p4) / 4;
					radius = math::max(math::max((p1 - center).norm(), (p2 - center).norm()), math::max((p3 - center).norm(), (p4 - center).norm()));
				}
			}
		}

		/**
		Tests if this sphere contains the given point, assumes that all points x that satisfy (center - x) <= r belong to the sphere.
		*/
		bool containsInclusive(const math::Vec3<T>& x) const
		{
			return (x - center).norm() <= radius;
		}

		/**
		Tests if this sphere contains the given point, assumes that all points x that satisfy (center - x) < r belong to the sphere.
		*/
		bool contains(const math::Vec3<T>& x) const
		{
			return (x - center).norm() < radius;
		}

		/**
        Converts this object to string.
        */
        friend std::ostream& operator<<(std::ostream& stream, const Sphere<T>& v)
        {
            stream << "[Sphere center = " << v.center << ", radius = " << v.radius << "]";
            return stream;
        }

		/**
		Calculates signed distance of the given point and the surface of this sphere.
		The result is negative if the point is inside the sphere, zero if it is on the surface of the sphere
		and positive if it is outside of the sphere.
		*/
		double distance(const math::Vec3<T>& p) const
		{
			return (p - center).norm() - radius;
		} 



	private:
		/**
		Helper method for miniball method.
		*/
		template<typename POINTT> static Sphere<T> recurseMini(math::Vec3<POINTT> const* P[], size_t p, size_t b = 0)
		{
			Sphere<T> MB;

			// TODO: Remove weird indexing.
			switch(b)
			{
			case 0:
				MB = Sphere<T>();
				break;
			case 1:
				MB = Sphere<T>(*P[-1]);
				break;
			case 2:
				MB = Sphere<T>(*P[-1], *P[-2]);
				break;
			case 3:
				MB = Sphere<T>(*P[-1], *P[-2], *P[-3]);
				break;
			case 4:
				MB = Sphere<T>(*P[-1], *P[-2], *P[-3], *P[-4]);
				return MB;
			}

			for(size_t i = 0; i < p; i++)
			{
				if(MB.distance(*P[i]) > 0)
				{
					for(size_t j = i; j > 0; j--)
					{
						math::Vec3<POINTT> const* temp = P[j];
						P[j] = P[j - 1];
						P[j - 1] = temp;
					}

					MB = recurseMini(P + 1, i, b + 1);
				}
			}

			return MB;
		}


	public:
		/**
		Finds minimum ball over set of points.
		See Welzl, E. - Smallest enclosing disks (balls and ellipsoids)
		and
		http://www.flipcode.com/archives/Smallest_Enclosing_Spheres.shtml
		TODO: Remove unsafe memory allocation.
		*/
		template<typename POINTT> static Sphere<T> miniball(const std::vector<math::Vec3<POINTT>>& P)
		{
			std::vector<math::Vec3<T> > conv;
			conv.reserve(P.size());
			for (size_t n = 0; n < P.size(); n++)
				conv.push_back(math::Vec3<T>(P[n]));

			math::Vec3<T> const* *L = new math::Vec3<T> const*[P.size()];

			for (size_t i = 0; i < P.size(); i++)
				L[i] = &conv[i];
				//L[i] = &P[i];

			Sphere<T> MB = recurseMini(L, P.size());

			delete[] L;

			return MB;
		}

	};

}
