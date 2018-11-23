#pragma once

#include "math/vec3.h"

using math::Vec3;

namespace itl2
{

    template<typename T>
    class Box
    {
    public:
        /**
        Minimum and maximum coordinates of the box.
        */
        Vec3<T> minc, maxc;
        
        Box(const Vec3<T>& minc, const Vec3<T>& maxc) :
            minc(minc),
            maxc(maxc)
        {
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
        bool overlaps(const Box<T>& r) const
        {
            return maxc.x >= r.minc.x && r.maxc.x >= minc.x &&
                    maxc.y >= r.minc.y && r.maxc.y >= minc.y &&
                    maxc.z >= r.minc.z && r.maxc.z >= minc.z;
        }
        
        /**
        Calculates bounding box of points in the given list.
        */
        static Box<T> boundingBox(const vector<Vec3<T> >& points)
        {
            Vec3<T> minc(numeric_limits<T>::max(), numeric_limits<T>::max(), numeric_limits<T>::max());
            Vec3<T> maxc(numeric_limits<T>::lowest(), numeric_limits<T>::lowest(), numeric_limits<T>::lowest());
            
            for(const Vec3<T>& p : points)
            {
                minc = min(minc, p);
                maxc = max(maxc, p);
            }
            
            return Box(minc, maxc);
        }
    };

}
