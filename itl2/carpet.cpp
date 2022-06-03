
#include "carpet.h"
#include "io/raw.h"
#include "transform.h"

namespace itl2
{
	namespace tests
	{
		void carpet()
		{
			Image<uint16_t> orig;
			raw::read(orig, "../test_input_data/t1-head");

			// Crop smaller piece so that the carpet does not fall through the head as there's nothing in the edges
			// of the image.
			Image<uint16_t> img(85, 120, orig.depth());
			crop(orig, img, Vec3c(80, 90, 0));

			Image<float32_t> hmap;
			Image<uint16_t> vis;

			findSurface(img, hmap, 100, Direction::Down, 1.0, 150, &vis, img.height() / 2);

			raw::writed(vis, "./carpet/vis");
			raw::writed(hmap, "./carpet/hmap");

			Image<uint16_t> vis1;
			setValue(vis1, img);
			drawHeightMap(vis1, hmap, (uint16_t)800);
			raw::writed(vis1, "./carpet/full_carpet_vis");


			Image<uint16_t> vis2;
			setValue(vis2, img);
			setBeforeHeightMap(vis2, hmap, (uint16_t)800);
			raw::writed(vis2, "./carpet/full_carpet_set_before");


			Image<uint16_t> vis3;
			setValue(vis3, img);
			setAfterHeightMap(vis3, hmap, (uint16_t)800);
			raw::writed(vis3, "./carpet/full_carpet_set_after");
		}
	}
}