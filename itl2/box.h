#pragma once

#include "math/vec3.h"

#include <vector>
#include <limits>

namespace itl2
{

    template<typename T>
    class Box
    {
    public:
        /**
        Minimum and maximum coordinates of the box.
        */
        math::Vec3<T> minc, maxc;
        
        Box(const math::Vec3<T>& minc, const math::Vec3<T>& maxc) :
            minc(minc),
            maxc(maxc)
        {
        }
        
        /**
        Enlarges the box in all directions by given amount.
        */
        void inflate(T amount)
        {
            inflate(math::Vec3<T>(amount, amount, amount));
        }
        
        /**
        Enlarges the box in all directions by given amount.
        */
        void inflate(const math::Vec3<T>& amount)
        {
            minc -= amount;
            maxc += amount;
        }
        
        /**
        Tests if this box and the given box overlap.
        */
        bool overlaps(const Box<T>& r) const
        {
            return maxc.x >= r.minc.x && r.maxc.x >= minc.x &&
                    maxc.y >= r.minc.y && r.maxc.y >= minc.y &&
                    maxc.z >= r.minc.z && r.maxc.z >= minc.z;
        }

		/**
		Calculates intersection of this box and the given box.
		*/
		Box intersection(const Box<T>& r) const
		{
			return Box(componentwiseMax(minc, r.minc), componentwiseMin(maxc, r.maxc));
		}

		/**
		Calculates width of the box.
		*/
		coord_t width() const
		{
			return maxc.x - minc.x;
		}

		/**
		Calculates height of the box.
		*/
		coord_t height() const
		{
			return maxc.y - minc.y;
		}

		/**
		Calculates depth of the box.
		*/
		coord_t depth() const
		{
			return maxc.z - minc.z;
		}

		/**
		Calculates volume of the box.
		*/
		size_t volume() const
		{
			return (size_t)abs(width()) * (size_t)abs(height()) * (size_t)abs(depth());
		}
        
        /**
        Calculates bounding box of points in the given list.
        */
        static Box<T> boundingBox(const std::vector<math::Vec3<T> >& points)
        {
            math::Vec3<T> minc(std::numeric_limits<T>::max(), std::numeric_limits<T>::max(), std::numeric_limits<T>::max());
            math::Vec3<T> maxc(std::numeric_limits<T>::lowest(), std::numeric_limits<T>::lowest(), std::numeric_limits<T>::lowest());
            
            for(const math::Vec3<T>& p : points)
            {
                minc = min(minc, p);
                maxc = max(maxc, p);
            }
            
            return Box(minc, maxc);
        }
    };

}
