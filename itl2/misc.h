#pragma once

#include "image.h"
#include "projections.h"
#include "pointprocess.h"

#include <set>
#include <type_traits>
#include <limits>

namespace itl2
{
	/**
	Finds all unique colors in the given image.
	*/
	template<typename pixel_t> std::set<pixel_t> colors(const Image<pixel_t>& img)
	{
		std::set<pixel_t> c;
		for (coord_t i = 0; i < img.pixelCount(); i++)
		{
			c.insert(img(i));
		}
		return c;
	}

	/**
	Finds the highest possible value that is lower than the given value.
	*/
	template<typename T, typename std::enable_if_t<std::is_integral_v<T>, int> = 0> T decrement(T value)
	{
		return value - 1;
	}

	/**
	Finds the highest possible value that is lower than the given value.
	*/
	template<typename T, typename std::enable_if_t<std::is_floating_point_v<T>, int> = 0> T decrement(T value)
	{
		return nextafter(value, std::numeric_limits<T>::lowest());
	}

	/**
	Finds value that is not used by any of the pixels in the given image.
	*/
	template<typename pixel_t> pixel_t findUnusedValue(const Image<pixel_t>& image)
	{
		std::set<pixel_t> unique = colors(image);
		pixel_t val = std::numeric_limits<pixel_t>::max();
		while (unique.find(val) != unique.end())
		{
			if (val <= std::numeric_limits<pixel_t>::lowest())
				throw ITLException("The image does not contain any free color to be used as a temporary color.");

			decrement(val);
		}
		return val;
	}


	/**
	Dual thresholding with tracking.
	Results in an image where all regions that have value over upperThreshold are white, and additionally
	those regions that have value over lowerThreshold and are connected to a region with value over upperThreshold.
	*/
	template<typename pixel_t> void dualThreshold(Image<pixel_t>& img, pixel_t lowerThreshold, pixel_t upperThreshold)
	{
		// Multi-threshold to two classes.
		std::vector<pixel_t> th = { lowerThreshold, upperThreshold };
		multiThreshold(img, th);

		// Convert all those structures to "sure" that touch a "sure" structure.
		grow(img, pixelRound<pixel_t>(2), pixelRound<pixel_t>(1));

		// Threshold so that only "sure" structures are left.
		threshold<pixel_t>(img, 1);
	}


	/**
	Sets edge pixels to given value.
	@param img Image whose edges are to be set.
	@param value Value that is set to the edges of the image.
	@param r Thickness of strip of values at image edges to set.
	*/
	template<typename pixel_t, typename value_t> void setEdges(Image<pixel_t>& img, value_t value, Vec3c r)
	{
		pixel_t val = pixelRound<pixel_t, value_t>(value);
		
		r = max(r, Vec3c(1, 1, 1));

		// NOTE: This algorithm sets some pixels multiple times, beware if you use it for something else than setting pixel values!
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
						for (coord_t n = 0; n < std::min(r[skip], s); n++)
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
	template<typename pixel_t> void normalizeZ(Image<pixel_t>& img, float32_t globalMean = std::numeric_limits<float32_t>::signaling_NaN())
	{
		Image<float32_t> tmp, zmean, allmean;
		mean(img, 0, tmp);
		mean(tmp, 1, zmean);
		if(std::isnan(globalMean))
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


	///**
	//Scales image such that its mean and standard deviation have the given values.
	//*/
	//template<typename pixel_t> void normalize(Image<pixel_t>& img, pixel_t targetMean, pixel_t targetStdDev)
	//{

	//}


	namespace tests
	{
		void edges();
		void normalizeZ();
	}
}
