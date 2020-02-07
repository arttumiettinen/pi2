#pragma once

#include <omp.h>

#include "image.h"
#include "math/mathutils.h"
#include "boundarycondition.h"
#include "neighbourhoodtype.h"

namespace itl2
{

	/**
	Tests whether a position relative to the center of the neighbourhood is inside the neighbourhood.
	Zero elements in nbRadius are not tested at all.
	*/
	inline bool isInNeighbourhood(const Vec3c& relPos, NeighbourhoodType nbType, const Vec3c& nbRadius)
	{
		if (nbType == NeighbourhoodType::Rectangular)
		{
			size_t size = nbRadius.size();
			for (size_t n = 0; n < size; n++)
			{
				if (fabs(relPos[n]) > nbRadius[n] && nbRadius[n] > 0)
					return false;
			}
			return true;
		}
		else if (nbType == NeighbourhoodType::Ellipsoidal)
		{
			double sum = 0;
			size_t size = nbRadius.size();
			for (size_t n = 0; n < size; n++)
			{
				double ri = (double)relPos[n];
				double ai = (double)nbRadius[n];
				if (ai > 0)
					
				sum += (ri * ri) / (ai * ai);
			}

			return sum <= 1;
		}

		throw ITLException("Invalid or unsupported neighbourhood type.");
	}

	/**
	Creates a mask image whose size is the size of the neighbourhood and that contains 1 in pixels that are inside the neighbourhood
	and 0 everywhere else.
	*/
	template<typename pixel_t> void createNeighbourhoodMask(NeighbourhoodType nbType, coord_t radius, Image<pixel_t>& mask)
	{
		createNeighbourhoodMask(nbType, Vec3c(radius, radius, radius), mask);
	}

	/**
	Creates a mask image whose size is the size of the neighbourhood and that contains 1 in pixels that are inside the neighbourhood
	and 0 everywhere else.
	*/
	template<typename pixel_t> void createNeighbourhoodMask(NeighbourhoodType nbType, const Vec3c& nbRadius, Image<pixel_t>& mask)
	{
		mask.ensureSize(2 * nbRadius + Vec3c(1, 1, 1));

		#pragma omp parallel for if(mask.pixelCount() > PARALLELIZATION_THRESHOLD)
		for (coord_t z = 0; z < mask.depth(); z++)
		{
			for (coord_t y = 0; y < mask.height(); y++)
			{
				for (coord_t x = 0; x < mask.width(); x++)
				{
					Vec3c relPos = Vec3c(x, y, z) - nbRadius;
					mask(x, y, z) = isInNeighbourhood(relPos, nbType, nbRadius) ? (pixel_t)1 : (pixel_t)0;
				}
			}
		}
	}



	namespace internals
	{
		/**
		Extracts pixels in neighbourhood from a big image.
		Sets those neighbourhood pixels to zero that are outside the image.
		@param nb Neighbourhood pixels are assigned to this image. The size of this image is not checked but it must be 2 * nbRadius + 1.
		*/
		template<typename pixel_t, typename out_t = pixel_t> void getNeighbourhoodZero(const Image<pixel_t>& img, const Vec3c& nbCenter, const Vec3c& nbRadius, Image<out_t>& nb)
		{
			Vec3c start = nbCenter - nbRadius;
			Vec3c start0 = start;
			Vec3c end = nbCenter + nbRadius;

			// Make sure start point is in the image.
			// If start is > image dimensions, checks made for end will prevent the loops from running,
			// so there's no need to check that start < image dimensions.
			if (start.x < 0)
				start.x = 0;
			if (start.y < 0)
				start.y = 0;
			if (start.z < 0)
				start.z = 0;

			coord_t w = img.width();
			coord_t h = img.height();
			coord_t d = img.depth();

			// Make sure end point is in the image.
			if (end.x > w - 1)
				end.x = w - 1;
			if (end.y > h - 1)
				end.y = h - 1;
			if (end.z > d - 1)
				end.z = d - 1;

			// Zero neighbourhood if there's not enough data to fill it.
			if (start.x == 0 || start.y == 0 || start.z == 0 ||
				end.x == w - 1 || end.y == h - 1 || end.z == d - 1)
			{
//				#pragma omp parallel for if(nb.pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
				for (coord_t n = 0; n < nb.pixelCount(); n++)
					nb(n) = 0;
			}

			// Copy data to neighbourhood image.
//			#pragma omp parallel for if(nb.pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
			for (coord_t z = start.z; z <= end.z; z++)
			{
				for (coord_t y = start.y; y <= end.y; y++)
				{
					for (coord_t x = start.x; x <= end.x; x++)
					{
						pixel_t val = img(x, y, z);
						nb(x - start0.x, y - start0.y, z - start0.z) = pixelRound<out_t>(val);
					}
				}
			}
		}

		/**
		Extracts pixels in neighbourhood from a big image.
		Neighbourhood pixel outside of the image is set to the value of the nearest edge pixel.
		@param nb Neighbourhood pixels are assigned to this image. The size of this image is not checked but it must be 2 * nbRadius + 1.
		*/
		template<typename pixel_t, typename out_t = pixel_t> void getNeighbourhoodClamp(const Image<pixel_t>& img, const Vec3c& nbCenter, const Vec3c& nbRadius, Image<out_t>& nb)
		{
			Vec3c start = nbCenter - nbRadius;
			Vec3c end = nbCenter + nbRadius;

			// Copy data to neighbourhood image.
//			#pragma omp parallel for if(nb.pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
			for (coord_t z = start.z; z <= end.z; z++)
			{
				coord_t imgz = z;
				clamp<coord_t>(imgz, 0, img.depth() - 1);

				for (coord_t y = start.y; y <= end.y; y++)
				{
					coord_t imgy = y;
					clamp<coord_t>(imgy, 0, img.height() - 1);

					for (coord_t x = start.x; x <= end.x; x++)
					{
						coord_t imgx = x;
						clamp<coord_t>(imgx, 0, img.width() - 1);

						pixel_t val = img(imgx, imgy, imgz);
						nb(x - start.x, y - start.y, z - start.z) = pixelRound<out_t>(val);
					}
				}
			}
		}
	}

	/**
	Extracts pixels in a small neighbourhood from a larger image.
	Sets those neighbourhood pixels to zero that are outside the image.
	@param img The larger image.
	@param nbCenter Location of the center of the neighbourhood in the larger image.
	@param nbRadius Radius of the neighbourhood.
	@param nb Neighbourhood pixels are assigned to this image. The size of this image is not checked but it must be 2 * nbRadius + 1.
	@param bc Boundary condition.
	*/
	template<typename pixel_t, typename out_t = pixel_t> void getNeighbourhood(const Image<pixel_t>& img, const Vec3c& nbCenter, const Vec3c& nbRadius, Image<out_t>& nb, BoundaryCondition bc)
	{
		if (bc == BoundaryCondition::Zero)
			internals::getNeighbourhoodZero<pixel_t, out_t>(img, nbCenter, nbRadius, nb);
		else // if bc == Nearest
			internals::getNeighbourhoodClamp<pixel_t, out_t>(img, nbCenter, nbRadius, nb);
	}

	namespace tests
	{
		void neighbourhoodTools();
	}
}
