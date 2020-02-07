#pragma once

#include "math/vec3.h"
#include "math/numberutils.h"

#include <vector>
#include <limits>

namespace itl2
{

	/**
	Represents rectangular region between points minc and maxc, i.e. axis aligned box.
	*/
    template<typename T> class AABox
    {
	private:
		/**
		Returns the maximum of the given three numbers.
		*/
		static T max3(T a, T b, T c)
		{
			if (a >= b && a >= c)
				return a;
			if (b >= a && b >= c)
				return b;
			return c;
		}

		/**
		Calculates squared distance between box [xmin, xmax] x [ymin, ymax] x [zmin, zmax] and point (x, y, z).
		If (x, y, z) is inside the box, returns 0.
		The algorithm is from https://stackoverflow.com/questions/5254838/calculating-distance-between-a-point-and-a-rectangular-box-nearest-point
		*/
		static T boxDistanceSquared(T xmin, T xmax, T ymin, T ymax, T zmin, T zmax, T x, T y, T z)
		{
			T dx = max3(xmin - x, 0, x - xmax);
			T dy = max3(ymin - y, 0, y - ymax);
			T dz = max3(zmin - z, 0, z - zmax);
			return dx * dx + dy * dy + dz * dz;
		}

    public:
        /**
        Minimum and maximum coordinates of the box.
        */
        Vec3<T> minc, maxc;
        
		/**
		Constructor.
		Creates empty box at origin.
		*/
		AABox()
		{

		}

		/**
		Constructor.
		Creates box from minimum and maximum coordinates.
		*/
        AABox(const Vec3<T>& minc, const Vec3<T>& maxc) :
            minc(minc),
            maxc(maxc)
        {
        }

		/**
		Constructor.
		Creates box from box of other type by rounding and casting the coordinates.
		*/
		template<typename other_t> AABox(const AABox<other_t>& r) :
			minc(r.minc),
			maxc(r.maxc)
		{
		}

		/**
		Constructs box from position and size.
		*/
		static AABox<T> fromPosSize(const Vec3<T>& pos, const Vec3<T>& size)
		{
			return AABox(pos, pos + size);
		}

		/**
		Constructs box from center point and radius.
		*/
		static AABox<T> fromCenterRadius(const Vec3<T>& center, const Vec3<T>& radius)
		{
			return AABox(center - radius, center + radius);
		}
        
		/**
		Constructs box from center point and radius.
		*/
		static AABox<T> fromCenterRadius(const Vec3<T>& center, T radius)
		{
			Vec3<T> r(radius, radius, radius);
			return AABox<T>::fromCenterRadius(center, r);
		}

		/**
		Calculates squared distance between this box and point p.
		Returns 0 if the point is inside the box.
		*/
		T distanceSquared(const Vec3<T>& p) const
		{
			return boxDistanceSquared(minc.x, maxc.x, minc.y, maxc.y, minc.z, maxc.z, p.x, p.y, p.z);
		}

		/**
		Calculates distance between point p and this box.
		Returns 0 if the point is inside the box or on the edge of it.
		*/
		typename NumberUtils<T>::RealFloatType distance(const Vec3<T>& p) const
		{
			return sqrt(distanceSquared(p));
		}

		/**
		Tests if p.x_i is in range [box.min_i, box.max_i[ for each dimension i.
		*/
		bool contains(const Vec3<T>& p) const
		{
			return minc.x <= p.x && p.x < maxc.x &&
				minc.y <= p.y && p.y < maxc.y &&
				minc.z <= p.z && p.z < maxc.z;
		}

        /**
        Enlarges the box in all directions by given amount.
        */
        void inflate(T amount)
        {
            inflate(Vec3<T>(amount, amount, amount));
        }
        
        /**
        Enlarges the box in all directions by given amount.
        */
        void inflate(const Vec3<T>& amount)
        {
            minc -= amount;
            maxc += amount;
        }
        
        /**
        Tests if this box and the given box overlap.
        */
        bool overlaps(const AABox<T>& r) const
        {
            return maxc.x >= r.minc.x && r.maxc.x >= minc.x &&
                    maxc.y >= r.minc.y && r.maxc.y >= minc.y &&
                    maxc.z >= r.minc.z && r.maxc.z >= minc.z;
        }

		/**
		Calculates intersection of this box and the given box.
		*/
		AABox intersection(const AABox<T>& r) const
		{
			return AABox(max(minc, r.minc), min(maxc, r.maxc));
		}

		/**
		Calculates width of the box.
		*/
		T width() const
		{
			return maxc.x - minc.x;
		}

		/**
		Calculates height of the box.
		*/
		T height() const
		{
			return maxc.y - minc.y;
		}

		/**
		Calculates depth of the box.
		*/
		T depth() const
		{
			return maxc.z - minc.z;
		}

		/**
		Calculates volume of the box.
		*/
		T volume() const
		{
			return (T)abs(width()) * (T)abs(height()) * (T)abs(depth());
		}

		/**
		Calculates the size of the box.
		*/
		Vec3<T> size() const
		{
			return maxc - minc;
		}

		/**
		Gets position of the box.
		*/
		Vec3<T> position() const
		{
			return minc;
		}
        
        /**
        Calculates bounding box of points in the given list.
        */
        static AABox<T> boundingBox(const std::vector<Vec3<T> >& points)
        {
            Vec3<T> minc(std::numeric_limits<T>::max(), std::numeric_limits<T>::max(), std::numeric_limits<T>::max());
            Vec3<T> maxc(std::numeric_limits<T>::lowest(), std::numeric_limits<T>::lowest(), std::numeric_limits<T>::lowest());
            
            for(const Vec3<T>& p : points)
            {
                minc = min(minc, p);
                maxc = max(maxc, p);
            }
            
            return AABox(minc, maxc);
        }

    };

}
