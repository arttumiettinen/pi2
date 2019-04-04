#pragma once

#include "image.h"
#include "box.h"
#include "sphere.h"
#include "line.h"
#include "raytrace.h"
#include "math/vec3.h"

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
					math::Vec3c pos(x, y, z);
					img(x, y, z) = pixelRound<pixel_t>(pos[dimension]);
				}
			}
		}
	}

	/**
	Draws a box into the given image with given color.
	Clips the box to the image coordinates.
	@return Count of filled pixels.
	*/
	template<typename pixel_t> size_t draw(Image<pixel_t>& img, const Box<coord_t>& box, pixel_t value)
	{

		Box region = box.intersection(Box(math::Vec3c(), img.dimensions()));

		#pragma omp parallel for if(!omp_in_parallel() && region.volume() > PARALLELIZATION_THRESHOLD)
		for (coord_t z = region.minc.z; z < region.maxc.z; z++)
		{
			for (coord_t y = region.minc.y; y < region.maxc.y; y++)
			{
				for (coord_t x = region.minc.x; x < region.maxc.x; x++)
				{
					img(x, y, z) = value;
				}
			}
		}

		return box.volume();
	}

	/**
	Draws a sphere.
	@return Count of filled pixels.
	*/
	template<typename pixel_t, typename real_t> size_t draw(Image<pixel_t>& image, const Sphere<real_t>& sphere, pixel_t color)
	{
		math::Vec3c minPos = round(sphere.center - sphere.radius * math::Vec3<real_t>(1, 1, 1));
		math::Vec3c maxPos = round(sphere.center + sphere.radius * math::Vec3<real_t>(1, 1, 1));

		clamp(minPos, math::Vec3c(0, 0, 0), image.dimensions() - math::Vec3c(1, 1, 1));
		clamp(maxPos, math::Vec3c(0, 0, 0), image.dimensions() - math::Vec3c(1, 1, 1));

		size_t filledCount = 0;
		#pragma omp parallel for if(!omp_in_parallel() && Box(minPos, maxPos).volume() > PARALLELIZATION_THRESHOLD) reduction(+:filledCount)
		for(coord_t z = minPos.z; z <= maxPos.z; z++)
		{
			for(coord_t y = minPos.y; y <= maxPos.y; y++)
			{
				for(coord_t x = minPos.x; x <= maxPos.x; x++)
				{
					if(sphere.contains(math::Vec3<real_t>((real_t)x, (real_t)y, (real_t)z)))
					{
						image(x, y, z) = color;
						filledCount++;
					}
				}
			}
		}

		return filledCount;
	}

	/**
	Draws a line.
	*/
	template<typename pixel_t, typename real_t> void draw(Image<pixel_t>& image, const Line<real_t>& line, pixel_t color)
	{
		Vec3<real_t> start = line.start;
		Vec3<real_t> end = line.end;
		LineColorPlotter<pixel_t, real_t> plotter(image, color);
		siddonLineClip(start, end, plotter, Vec3<real_t>(image.dimensions()));
	}
}
