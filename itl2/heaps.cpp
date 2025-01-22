
#include <iostream>

#include "heaps.h"
#include "math/vec3.h"
#include "test.h"

namespace itl2
{
	namespace tests
	{
		inline void hheap()
		{
			//HHeap<Vec3c, int> points(0, 10, 10);
			//GHeap<Vec3c, int> points;
			BucketMap<Vec3c, int> points;


			points.push(1, Vec3c(1, 0, 0));
			points.push(1, Vec3c(2, 0, 0));
			points.push(1, Vec3c(3, 0, 0));

			points.push(0, Vec3c(4, 0, 0));
			points.push(0, Vec3c(5, 0, 0));
			points.push(0, Vec3c(6, 0, 0));

			points.push(7, Vec3c(-2, 0, 0));
			points.push(7, Vec3c(-1, 0, 0));
			points.push(7, Vec3c(0, 0, 0));


			points.push(-1, Vec3c(7, 0, 0));

			points.push(15, Vec3c(-3, 0, 0));

			int topPriority = points.topPriority();
			while (!points.empty())
			{
				testAssert(points.topPriority() <= topPriority, "Invalid retrieval order");

				topPriority = points.topPriority();
				std::cout << points.topItem() << ", priority = " << points.topPriority() << std::endl;

				points.pop();
			}

			return;
		}
	}
}