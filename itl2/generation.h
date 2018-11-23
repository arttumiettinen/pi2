#pragma once

#include "image.h"

namespace itl2
{

	/**
	Fills given image with ramp in given dimension.
	*/
	template<typename pixel_t> void ramp(Image<pixel_t>& img, size_t dimension)
	{
		if (dimension > 2)
			throw ITLException("Invalid ramp dimension.");

		#pragma omp parallel for if(!omp_in_parallel())
		for (coord_t z = 0; z < img.depth(); z++)
		{
			for (coord_t y = 0; y < img.height(); y++)
			{
				for (coord_t x = 0; x < img.width(); x++)
				{
					Vec3c pos(x, y, z);
					img(x, y, z) = pixelRound<pixel_t>(pos[dimension]);
				}
			}
		}
	}

}