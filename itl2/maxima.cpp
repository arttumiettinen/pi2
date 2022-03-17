
#include "maxima.h"
#include "io/raw.h"
#include "dmap.h"
#include "thickmap.h"

namespace itl2
{
	namespace tests
	{
		void localMaxima()
		{
			// Make an image
			Image<uint32_t> geom(100, 100, 100);
			draw(geom, Sphere(Vec3c(50, 13, 50), (coord_t)18), (uint32_t)255);
			draw(geom, AABoxc::fromMinMax(Vec3c(20, 30, 30), Vec3c(80, 70, 70)), (uint32_t)255);

			raw::writed(geom, "maxima/geometry");

			// Calculate its distance map
			Image<float32_t> dmap;
			distanceTransform(geom, dmap);

			// Find and label local maxima
			auto results = findLocalMaxima(dmap);

			Image<uint8_t> vis;
			setValue(vis, geom);
			draw(vis, results);
			raw::writed(vis, "maxima/dmap_local_maxima");

			// Create thickness map. NOTE: Thickness map has only 1 local maxima (inside the rectangle) as the sphere has smaller diameter.
			Image<float32_t> tmap;
			thicknessMap(geom, tmap);

			// Find and label local maxima
			results = findLocalMaxima(tmap);

			setValue(vis, geom);
			draw(vis, results);
			raw::writed(vis, "maxima/tmap_local_maxima");
		}
	}
}