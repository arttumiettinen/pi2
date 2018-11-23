#pragma once

#include "image.h"
#include "math/mathutils.h"

namespace itl2
{
	/**
	Copies pixels of input image to output image, and converts data type in the process.
	Resizes output image to the size of the input image if necessary.
	@param in Input image.
	@param out Output image.
	*/
	template<typename pixel_t, typename newPixel_t> void convert(const Image<pixel_t>& in, Image<newPixel_t>& out)
	{
		out.ensureSize(in);

		size_t counter = 0;
		#pragma omp parallel for if(in.pixelCount() > PARALLELIZATION_THRESHOLD)
		for (coord_t n = 0; n < in.pixelCount(); n++)
		{
			out(n) = math::pixelRound<newPixel_t>(in(n));

			// Showing progress info here would induce more processing than is done in the whole loop.
		}
	}

}
