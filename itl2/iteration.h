#pragma once

// NOTE: This file contains experiments on "new-style" iteration that is to replace most loops over the image,
// and in particular replace pointprocess.h.

#include "image.h"
#include "math/aabox.h"
#include "progress.h"

namespace itl2
{

	/**
	Call lambda(x, y, z) for all (x, y, z) in range [block.minc, block.maxc[.
	*/
	template<typename F>
	void forAllInBox(const AABox<coord_t>& block, F&& lambda, bool showProgressIndicator = false)
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

		if (d <= 0)
		{
			return;
		}
		else if (d <= 1)
		{
			// 2D image
			ProgressIndicator progress(block.maxc.y - block.minc.y, showProgressIndicator);
			#pragma omp parallel for if(block.volume() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
			for (coord_t y = block.minc.y; y < block.maxc.y; y++)
			{
				for (coord_t x = block.minc.x; x < block.maxc.x; x++)
				{
					lambda(x, y, block.minc.z);
				}
				progress.step();
			}
		}
		else
		{
			// 3D image
			ProgressIndicator progress(block.maxc.z - block.minc.z, showProgressIndicator);
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
				progress.step();
			}
		}
	}

	/**
	Call lambda(coord_t x, coord_t y, coord_t z) for all pixels in the image.
	*/
	template<typename F>
	void forAllPixels(const ImageBase& img, F&& lambda, bool showProgressIndicator = false)
	{
		forAllInBox(AABox<coord_t>(Vec3c(), img.dimensions()), lambda, showProgressIndicator);
	}

	/**
	Perform img(x, y, z) = pixelRound<pixel_t>(lambda(img(x, y, z))) for all (x, y, z) in the image.
	*/
	template<typename pixel_t, typename F>
	void forAll(Image<pixel_t>& img, F&& lambda, bool showProgressIndicator = false)
	{
		forAllPixels(img, [&](coord_t x, coord_t y, coord_t z)
			{
				img(x, y, z) = pixelRound<pixel_t>(lambda(img(x, y, z)));
			},
			showProgressIndicator);
	}

	/**
	Perform l(x, y, z) = pixelRound<pixel_t>(lambda(l(x, y, z), r(x, y, z))) for all (x, y, z) in the image.
	*/
	template<typename pixel1_t, typename pixel2_t, typename F>
	void forAll(Image<pixel1_t>& l, const Image<pixel2_t>& r, F&& lambda, bool showProgressIndicator = false)
	{
		l.ensureSize(r);
		forAllPixels(l, [&](coord_t x, coord_t y, coord_t z)
			{
				l(x, y, z) = pixelRound<pixel1_t>(lambda(l(x, y, z), r(x, y, z)));
			},
			showProgressIndicator);
	}


	namespace internals
	{

		/**
		Helper for forEdges function.
		Processes x and y edges for given z. The z coordinate must be in [block.minc.z+1, block.maxc.z-2].
		This code is from StackOverflow: https://stackoverflow.com/questions/53432767/how-to-iterate-over-pixels-on-edge-of-a-square-in-1-iteration
		*/
		template<typename F>
		void forEdgesHelperXY(const AABox<coord_t>& block, F&& lambda, coord_t z, size_t boxDimensionality)
		{
			// y = block.minc.y
			if (boxDimensionality >= 2)
			{
				for (coord_t x = block.minc.x; x < block.maxc.x; x++)
				{
					lambda(x, block.minc.y, z);
				}

				// y in block
				for (coord_t y = block.minc.y + 1; y < block.maxc.y - 1; y++)
				{
					lambda(block.minc.x, y, z);
					lambda(block.maxc.x - 1, y, z);
				}

				// y = block.maxc.y - 1
				for (coord_t x = block.minc.x; x < block.maxc.x; x++)
				{
					lambda(x, block.maxc.y - 1, z);
				}
			}
			else
			{
				for (coord_t y = block.minc.y; y < block.maxc.y; y++)
				{
					lambda(block.minc.x, y, z);
					lambda(block.maxc.x - 1, y, z);
				}
			}
		}
	}

	/**
	Call lambda(x, y, z) for all (x, y, z) in the edges of the box that spans [block.minc, block.maxc[.
	This function is not multi-threaded and processes the edge points in the same order.
	Each edge point is processed once.
	@param boxDimensionality If 1, fills only x-directional edges. If 2, fills x- and y-directional edges. If 3, fills x-, y-, and z-directional edges.
	*/
	template<typename F>
	void forEdges(const AABox<coord_t>& block, size_t boxDimensionality, F&& lambda)
	{
		if (boxDimensionality >= 3)
		{
			// z = block.minc.z
			coord_t z = block.minc.z;
			for (coord_t y = block.minc.y; y < block.maxc.y; y++)
			{
				for (coord_t x = block.minc.x; x < block.maxc.x; x++)
				{
					lambda(x, y, z);
				}
			}

			// z between block.minc.z + 1 and block.maxc.z - 1
			for (coord_t z = block.minc.z + 1; z < block.maxc.z - 1; z++)
			{
				internals::forEdgesHelperXY(block, lambda, z, boxDimensionality);
			}

			// z = block.maxc.z - 1
			z = block.maxc.z - 1;
			for (coord_t y = block.minc.y; y < block.maxc.y; y++)
			{
				for (coord_t x = block.minc.x; x < block.maxc.x; x++)
				{
					lambda(x, y, z);
				}
			}
		}
		else
		{
			for (coord_t z = block.minc.z; z < block.maxc.z; z++)
			{
				internals::forEdgesHelperXY(block, lambda, z, boxDimensionality);
			}
		}
	}
	
}