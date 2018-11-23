#pragma once

#include "image.h"
#include "projections.h"

namespace itl2
{
	/**
	Sets edge pixels to given value.
	@param img Image whose edges are to be set.
	@param value Value that is set to the edges of the image.
	@param r Thickness of strip of values at image edges to set.
	*/
	template<typename pixel_t, typename value_t> void setEdges(Image<pixel_t>& img, value_t value, Vec3c r)
	{
		pixel_t val = pixelRound<pixel_t, value_t>(value);
		
		r = componentwiseMax(r, Vec3c(1, 1, 1));

		for (size_t skip = 0; skip < img.dimensionality(); skip++)
		{
			Vec3c reducedsize = img.dimensions();
			reducedsize[skip] = 1;

			for (coord_t z = 0; z < reducedsize.z; z++)
			{
				for (coord_t y = 0; y < reducedsize.y; y++)
				{
					for (coord_t x = 0; x < reducedsize.x; x++)
					{
						Vec3c coords(x, y, z);

						coord_t s = img.dimension(skip) - 1;
						for (coord_t n = 0; n < math::min(r[skip], s); n++)
						{
							coords[skip] = n;
							img(coords) = val;
							coords[skip] = s - n;
							img(coords) = val;
						}
					}
				}
			}
		}
	}

	/**
	Sets edge pixels to given value.
	@param img Image whose edges are to be set.
	@param value Value that is set to the edges of the image.
	@param r Thickness of strip of values at image edges to set.
	*/
	template<typename pixel_t, typename value_t> void setEdges(Image<pixel_t>& img, value_t value = 0, coord_t r = 1)
	{
		setEdges(img, value, Vec3c(r, r, r));
	}


	/**
	Ensures that all xy-slices have the same mean.
	@param img Image to process.
	@param globalMean The desired mean value. NaN corresponds to the global mean of the image.
	*/
	template<typename pixel_t> void normalizeZ(Image<pixel_t>& img, float32_t globalMean = numeric_limits<float32_t>::signaling_NaN())
	{
		Image<float32_t> tmp, zmean, allmean;
		mean(img, 0, tmp);
		mean(tmp, 1, zmean);
		if(math::isnan(globalMean))
			globalMean = (float32_t)mean(zmean);

		#pragma omp parallel for if(img.pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		for (coord_t z = 0; z < img.depth(); z++)
		{
			// NOTE: The order of indices in zmean is designed for visual inspection and from this point of view the order is wrong...
			float32_t shift = globalMean - zmean(img.depth() - 1 - z);

			for (coord_t y = 0; y < img.height(); y++)
			{
				for (coord_t x = 0; x < img.width(); x++)
				{
					img(x, y, z) = pixelRound<pixel_t>((float32_t)img(x, y, z) + shift);
				}
			}
		}

	}


	namespace tests
	{
		void edges();
		void normalizeZ();
	}
}
