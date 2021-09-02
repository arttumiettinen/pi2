
#include "sdmap.h"
#include "generation.h"
#include "pointprocess.h"

namespace itl2
{
	namespace tests
	{
		void seededDMap()
		{
			Image<float32_t> distance;

			Image<uint8_t> geometry(50, 50, 50);
			draw(geometry, Sphere(Vec3f(15, 15, 25), 10.0f), (uint8_t)255);
			draw(geometry, Sphere(Vec3f(15, 15 + 18, 25), 10.0f), (uint8_t)255);
			draw(geometry, Sphere(Vec3f(15 + 18, 15 + 18, 25), 10.0f), (uint8_t)255);
			draw(geometry, Sphere(Vec3f(15 + 22, 15, 25), 10.0f), (uint8_t)255);

			raw::writed(geometry, "./sdmap/geometry");

			Image<uint8_t> seeds(geometry.dimensions());
			seeds(24, 19, 25) = 255;
			raw::writed(seeds, "./sdmap/seeds");

			seededDistanceMap(seeds, geometry, distance);
			raw::writed(distance, "./sdmap/seeded_dmap");
		}
	}
}