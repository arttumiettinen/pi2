#pragma once

/*
This file contains code that has been compared to code in surfaceskeleton.h
Both pieces of code seem to give the same result, but this one is more clear so
I have left it available here if it is needed in the future.
*/

#include "image.h"
#include "floodfill.h"

#include "surfaceskeleton.h"

namespace itl2
{
	namespace experimental
	{
		namespace internals
		{
			/**
			Tests if the pixel in the center of the given 3x3x3 neighbourhood has background 6-neighbour.

			This test equals internals:N, internals:S, etc. calls in itl2::surfaceSkeleton function,
			but for 6-subiteration algorithm to work one must test only single direction per sub-iteration.

			*/
			//template<typename pixel_t> bool hasBackground6Neighbour(const Image<pixel_t>& nb)
			//{
			//	constexpr coord_t c = 1;
			//	return nb(c - 1, c, c) == 0 ||
			//		nb(c + 1, c, c) == 0 ||
			//		nb(c, c - 1, c) == 0 ||
			//		nb(c, c + 1, c) == 0 ||
			//		nb(c, c, c - 1) == 0 ||
			//		nb(c, c, c + 1) == 0;
			//}

			/**
			Tests if the pixel in the center of the given 3x3x3 neighbourhood has foreground 26-neighbour.
			*/
			template<typename pixel_t> bool hasForeground26Neighbour(const Image<pixel_t>& nb)
			{
				coord_t fgCount = 0;
				for (coord_t n = 0; n < nb.pixelCount(); n++)
					if (nb(n) != 0)
						fgCount++;

				return fgCount >= 2; // Center pixel is always foreground, so if count is 2 there is at least one foreground neighbour.
			}

			/**
			Tests if the set of foreground 26-neighbors of center point in the given 3x3x3 neighbourhood is 26-connected.

			This function PROBABLY gives the same result than itl2::isSimplePointHybrid, but I have not tested that very well.

			*/
			template<typename pixel_t> bool isForeground26NeighbourSet26Connected(const Image<pixel_t>& nb)
			{
				Image<pixel_t> nb2;
				setValue(nb2, nb);

				nb2(1, 1, 1) = 0;

				// Flood fill from any foreground pixel
				for (coord_t z = 0; z < nb2.depth(); z++)
				{
					for (coord_t y = 0; y < nb2.height(); y++)
					{
						for (coord_t x = 0; x < nb2.width(); x++)
						{
							Vec3c p(x, y, z);
							if (nb2(p) != 0)
							{
								slowFloodfill(nb2, p, (pixel_t)0, (pixel_t)0, Connectivity::AllNeighbours);

								// Break all loops
								x = nb2.width() + 1;
								y = nb2.height() + 1;
								z = nb2.depth() + 1;
								break;
							}
						}
					}
				}

				// Test if the result is zero. If yes, all the 26-neighbours are 26-connected as the 26-connected fill filled all of them when started from one of the points.
				return max(nb2) == 0;
			}

			/**
			Tests if the set of background 6-neighbors of center point of given neighbourhood is 6-connected in the set of background 18-neighbors.
			*/
			template<typename pixel_t> bool areBackground6Neighbours6ConnectedInBackground18Neighbours(const Image<pixel_t>& nb)
			{
				// Create neighbourhood that contains background 6-neighbours and background 18-neighbours as foreground points (= everything except corners)
				Image<pixel_t> nb2(3, 3, 3);

				constexpr coord_t c = 1;
				//nb2(c - 1, c - 1, c - 1) = nb(c - 1, c - 1, c - 1) == 0 ? 1 : 0;
				nb2(c, c - 1, c - 1) = nb(c, c - 1, c - 1) == 0 ? 1 : 0;
				//nb2(c + 1, c - 1, c - 1) = nb(c + 1, c - 1, c - 1) == 0 ? 1 : 0;
				nb2(c - 1, c, c - 1) = nb(c - 1, c, c - 1) == 0 ? 1 : 0;
				nb2(c, c, c - 1) = nb(c, c, c - 1) == 0 ? 1 : 0;
				nb2(c + 1, c, c - 1) = nb(c + 1, c, c - 1) == 0 ? 1 : 0;
				//nb2(c - 1, c + 1, c - 1) = nb(c - 1, c + 1, c - 1) == 0 ? 1 : 0;
				nb2(c, c + 1, c - 1) = nb(c, c + 1, c - 1) == 0 ? 1 : 0;
				//nb2(c + 1, c + 1, c - 1) = nb(c + 1, c + 1, c - 1) == 0 ? 1 : 0;

				nb2(c - 1, c - 1, c) = nb(c - 1, c - 1, c) == 0 ? 1 : 0;
				nb2(c, c - 1, c) = nb(c, c - 1, c) == 0 ? 1 : 0;
				nb2(c + 1, c - 1, c) = nb(c + 1, c - 1, c) == 0 ? 1 : 0;
				nb2(c - 1, c, c) = nb(c - 1, c, c) == 0 ? 1 : 0;
				//nb2(c    , c    , c    ) = nb(c    , c    , c    ) == 0 ? 1 : 0;
				nb2(c + 1, c, c) = nb(c + 1, c, c) == 0 ? 1 : 0;
				nb2(c - 1, c + 1, c) = nb(c - 1, c + 1, c) == 0 ? 1 : 0;
				nb2(c, c + 1, c) = nb(c, c + 1, c) == 0 ? 1 : 0;
				nb2(c + 1, c + 1, c) = nb(c + 1, c + 1, c) == 0 ? 1 : 0;

				//nb2(c - 1, c - 1, c + 1) = nb(c - 1, c - 1, c + 1) == 0 ? 1 : 0;
				nb2(c, c - 1, c + 1) = nb(c, c - 1, c + 1) == 0 ? 1 : 0;
				//nb2(c + 1, c - 1, c + 1) = nb(c + 1, c - 1, c + 1) == 0 ? 1 : 0;
				nb2(c - 1, c, c + 1) = nb(c - 1, c, c + 1) == 0 ? 1 : 0;
				nb2(c, c, c + 1) = nb(c, c, c + 1) == 0 ? 1 : 0;
				nb2(c + 1, c, c + 1) = nb(c + 1, c, c + 1) == 0 ? 1 : 0;
				//nb2(c - 1, c + 1, c + 1) = nb(c - 1, c + 1, c + 1) == 0 ? 1 : 0;
				nb2(c, c + 1, c + 1) = nb(c, c + 1, c + 1) == 0 ? 1 : 0;
				//nb2(c + 1, c + 1, c + 1) = nb(c + 1, c + 1, c + 1) == 0 ? 1 : 0;

				// Fill with zero from any 6-neighbour
				if (nb2(c - 1, c, c) != 0)
					slowFloodfill(nb2, Vec3c(c - 1, c, c), (pixel_t)0, (pixel_t)0, Connectivity::NearestNeighbours);
				else if (nb2(c + 1, c, c) != 0)
					slowFloodfill(nb2, Vec3c(c + 1, c, c), (pixel_t)0, (pixel_t)0, Connectivity::NearestNeighbours);
				else if (nb2(c, c - 1, c) != 0)
					slowFloodfill(nb2, Vec3c(c, c - 1, c), (pixel_t)0, (pixel_t)0, Connectivity::NearestNeighbours);
				else if (nb2(c, c + 1, c) != 0)
					slowFloodfill(nb2, Vec3c(c, c + 1, c), (pixel_t)0, (pixel_t)0, Connectivity::NearestNeighbours);
				else if (nb2(c, c, c - 1) != 0)
					slowFloodfill(nb2, Vec3c(c, c, c - 1), (pixel_t)0, (pixel_t)0, Connectivity::NearestNeighbours);
				else if (nb2(c, c, c + 1) != 0)
					slowFloodfill(nb2, Vec3c(c, c, c + 1), (pixel_t)0, (pixel_t)0, Connectivity::NearestNeighbours);

				// Test if all 6-neighbours are zero
				// If yes, all of them could be reached from one of them through background 18-neighbours with 6-connected route.
				return nb2(c - 1, c, c) == 0 &&
					nb2(c + 1, c, c) == 0 &&
					nb2(c, c - 1, c) == 0 &&
					nb2(c, c + 1, c) == 0 &&
					nb2(c, c, c - 1) == 0 &&
					nb2(c, c, c + 1) == 0;
			}
		}

		/*
		Performs one thinning iteration and returns count of pixels removed from the image.
		Repeating this thinning process until the image does not change results in hybrid skeleton.

		This method gives the same results as itl2::surfaceThin.

		@param img Image that should be thinned.
		*/
		template<typename pixel_t> size_t surfaceThin2(Image<pixel_t>& img)
		{
			itl2::internals::createAllowedHashList();

			size_t changed = 0;
			ProgressIndicator progress(6 * img.depth());
			std::vector<Vec3c> pointsToRemove;
			pointsToRemove.reserve(10000);

			for (int currentBorder = 1; currentBorder <= 6; currentBorder++)
			{
				pointsToRemove.clear();

				#pragma omp parallel
				{
					Image<uint8_t> nbPrivate(3, 3, 3);

					#pragma omp for
					for (coord_t z = 0; z < img.depth(); z++)
					{
						for (coord_t y = 0; y < img.height(); y++)
						{
							for (coord_t x = 0; x < img.width(); x++)
							{
								pixel_t c = img(x, y, z);
								if (c != (pixel_t)0)
								{

									// check 6-neighbors to see if the point is a border point of type currentBorder
									bool isBorderPoint = false;
									switch (currentBorder)
									{
									case 1: isBorderPoint = itl2::internals::N(img, x, y, z); break;
									case 2: isBorderPoint = itl2::internals::S(img, x, y, z); break;
									case 3: isBorderPoint = itl2::internals::E(img, x, y, z); break;
									case 4: isBorderPoint = itl2::internals::W(img, x, y, z); break;
									case 5: isBorderPoint = itl2::internals::U(img, x, y, z); break;
									case 6: isBorderPoint = itl2::internals::B(img, x, y, z); break;
									}

									if (isBorderPoint)
									{
										getNeighbourhood(img, Vec3c(x, y, z), Vec3c(1, 1, 1), nbPrivate, BoundaryCondition::Zero);

										if (//internals::hasBackground6Neighbour(nbPrivate) && // This condition is replaced by isBorderPoint test above.
											internals::hasForeground26Neighbour(nbPrivate) &&
											internals::isForeground26NeighbourSet26Connected(nbPrivate) &&
											internals::areBackground6Neighbours6ConnectedInBackground18Neighbours(nbPrivate))
										{
											#pragma omp critical(surfacethin2_insert)
												pointsToRemove.push_back(Vec3c(x, y, z));
										}
									}
								}
							}

						}

						progress.step();
					}
				}

				// This is required to make the result exactly the same than in Fiji (non-threaded version)
				// Otherwise, the point removal order may change the skeleton points.
				std::sort(pointsToRemove.begin(), pointsToRemove.end(), vecComparer<coord_t>);

				// Re-check while removing points so that connectivity is preserved.
				Image<uint8_t> nb(3, 3, 3);
				for (size_t n = 0; n < pointsToRemove.size(); n++)
				{
					Vec3c p = pointsToRemove[n];

					getNeighbourhood(img, p, Vec3c(1, 1, 1), nb, BoundaryCondition::Zero);
					//nb(1, 1, 1) = 0; // This must not be done. Function isSimplePointHybrid does not consider the center pixel, but it is needed for isSurfacePoint.

					if (internals::isForeground26NeighbourSet26Connected(nb) &&
						!itl2::internals::isSurfacePoint(nb))
					{
						img(p) = 0;
						changed++;
					}
				}
			}

			return changed;
		}


		/**
		Calculate skeleton of objects in the given binary image.
		Uses simple point conditions in algorithm in Theorem 2 of
		Saha - A survey on skeletonization algorithms and their applications.
		The thinning process is implemented as 6-directional subiteration.
		Planar points are detected using idea from
		Lee et al. Building skeleton models via 3-D medial surface/axis thinning algorithms. Computer Vision, Graphics, and Image Processing, 56(6):462–478, 1994.
		but additional conditions are added as those in the paper are not sufficient in many practical cases.
		See internals::isSurfacePoint and internals::createAllowedHashList for details.

		This method gives the same results as itl2::surfaceSkeleton, but is slower (although probably more understandable).
		*/
		template<typename pixel_t> void surfaceSkeleton2(Image<pixel_t>& img, size_t maxIterations = std::numeric_limits<size_t>::max())
		{

			size_t changes;
			size_t it = 0;
			do
			{
				changes = surfaceThin2(img);

				std::cout << std::endl << changes << " pixels removed." << std::endl;

				it++;
				//raw::writed(img, string("./skeleton/surface_skeleton2_sequence") + toString(it));
				if (it >= maxIterations)
					break;

			} while (changes > 0);

		}

		namespace tests
		{
			void surfaceSkeleton2();
		}
	}

}
