#pragma once

// NOTE: This file contains experiments on "new-style" iteration that is to replace most loops over the image,
// and in particular replace pointprocess.h.

#include "image.h"
#include "aabox.h"

namespace itl2
{

	/**
	Call lambda(x, y, z) for all (x, y, z) in range [block.minc, block.maxc[.
	*/
	template<typename F>
	void forAllInBox(const AABox<coord_t>& block, F&& lambda)
	{
		//#pragma omp parallel for if(block.volume() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		//for (coord_t z = block.minc.z; z < block.maxc.z; z++)
		//{
		//	for (coord_t y = block.minc.y; y < block.maxc.y; y++)
		//	{
		//		for (coord_t x = block.minc.x; x < block.maxc.x; x++)
		//		{
		//			lambda(x, y, z);
		//		}
		//	}
		//}

		//coord_t w = block.width();
		//coord_t h = block.height();
		//coord_t d = block.depth();
		//coord_t dh = d * h;
		//#pragma omp parallel for if(block.volume() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		//for (coord_t k = 0; k < dh; k++)
		//{
		//	coord_t z = k / h + block.minc.z;
		//	coord_t y = k % h + block.minc.y;
		//	for (coord_t x = block.minc.x; x < block.maxc.x; x++)
		//	{
		//		lambda(x, y, z);
		//	}
		//}

		coord_t w = block.width();
		coord_t h = block.height();
		coord_t d = block.depth();

		if (d <= 1)
		{
			// 2D image
			#pragma omp parallel for if(block.volume() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
			for (coord_t y = block.minc.y; y < block.maxc.y; y++)
			{
				for (coord_t x = block.minc.x; x < block.maxc.x; x++)
				{
					lambda(x, y, block.minc.z);
				}
			}
		}
		else
		{
			// 3D image
			#pragma omp parallel for if(block.volume() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
			for (coord_t z = block.minc.z; z < block.maxc.z; z++)
			{
				for (coord_t y = block.minc.y; y < block.maxc.y; y++)
				{
					for (coord_t x = block.minc.x; x < block.maxc.x; x++)
					{
						lambda(x, y, z);
					}
				}
			}
		}
	}

	/**
	Call lambda(coord_t x, coord_t y, coord_t z) for all pixels in the image.
	*/
	template<typename F>
	void forAllPixels(const ImageBase& img, F&& lambda)
	{
		forAllInBox(AABox<coord_t>(Vec3c(), img.dimensions()), lambda);
	}

	/**
	Perform img(x, y, z) = pixelRound<pixel_t>(lambda(img(x, y, z))) for all (x, y, z) in the image.
	*/
	template<typename pixel_t, typename F>
	void forAll(Image<pixel_t>& img, F&& lambda)
	{
		forAllPixels(img, [&](coord_t x, coord_t y, coord_t z)
			{
				img(x, y, z) = pixelRound<pixel_t>(lambda(img(x, y, z)));
			});
	}

	/**
	Perform l(x, y, z) = pixelRound<pixel_t>(lambda(l(x, y, z), r(x, y, z))) for all (x, y, z) in the image.
	*/
	template<typename pixel1_t, typename pixel2_t, typename F>
	void forAll(Image<pixel1_t>& l, const Image<pixel2_t>& r, F&& lambda)
	{
		l.ensureSize(r);
		forAllPixels(l, [&](coord_t x, coord_t y, coord_t z)
			{
				l(x, y, z) = pixelRound<pixel1_t>(lambda(l(x, y, z), r(x, y, z)));
			});
	}
}