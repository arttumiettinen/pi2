#pragma once

#include "math/vec3.h"
#include "image.h"
#include "connectivity.h"
#include "floodfill.h"
#include "generation.h"
#include "misc.h"
#include "indexforest.h"

#include <vector>
#include <set>

namespace itl2
{
	/**
	Finds local maxima in a 3D image.
	@param orig Maxima of this image are to be found.
	@return List of pixel locations for each maximum.
	*/
	template<typename pixel_t> std::vector<std::vector<Vec3sc> > findLocalMaxima(Image<pixel_t>& orig, Connectivity connectivity = Connectivity::AllNeighbours)
	{
		pixel_t tempColor = findUnusedValue(orig);

		Image<uint8_t> processed(orig.dimensions());

		std::vector<std::vector<Vec3sc> > results;
		std::vector<Vec3sc> filledPoints;
		std::set<pixel_t> neighbourValues;

		ProgressIndicator progress(orig.depth());
		for (coord_t z = 0; z < orig.depth(); z++)
		{
			for (coord_t y = 0; y < orig.height(); y++)
			{
				for (coord_t x = 0; x < orig.width(); x++)
				{
					Vec3c p(x, y, z);
					pixel_t v = orig(p);

					// Only continue if the current pixel is in a non-processed & non-zero region
					if (v != 0 && processed(p) == 0)
					{
						filledPoints.clear();
						neighbourValues.clear();

						// Find all points in the current region
						floodfillSingleThreaded(orig, p, tempColor, tempColor, connectivity, nullptr, &filledPoints, 0, &neighbourValues);

						// Put the correct value back
						draw(orig, filledPoints, v);

						// Mark the current region as processed
						draw(processed, filledPoints, (uint8_t)1);

						// Check if this region is a local maximum
						bool isMax = true;
						for (pixel_t val : neighbourValues)
						{
							if (val != 0 &&		// Not background
								val > v)		// > region value
							{
								// Region is not a local maximum as it has a
								// non-background neighbour with larger value
								isMax = false;
								break;
							}
						}

						if (isMax)
							results.push_back(filledPoints);
					}
				}
			}
			
			progress.step();
		}

		return results;
	}


	/**
	Removes all maxima that are smaller in radius than neighbouring maximum.
	Maximum is neighbour to another maximum if distance between them is less than radius of the larger maximum multiplied by radiusMultiplier.
	Removes all maxima m that satisfy distance(m, n) < radiusMultiplier * radius(n) and radius(n) > radius(m) for some n.
	The distance is measured between centroids of the maxima.
	The maxima are removed by combining them to the larger maxima.
	@param maximaList List of local maxima, see findLocalMaxima.
	@param orig Original image from which the maxima were found.
	@param radiusMultiplier Multiplier to apply to the radius of the larger maximum.
	*/
	template<typename pixel_t> void removeMaximaInsideLargerOnes(std::vector<std::vector<Vec3sc> >& maximaList, const Image<pixel_t>& orig, double radiusMultiplier = 1)
	{
		std::vector<Vec3d> centroids(maximaList.size(), Vec3d());
		std::vector<double> radii(maximaList.size(), 0.0);
		std::vector<bool> removalFlags(maximaList.size(), false);
		IndexForest sets(maximaList.size());

		// Calculate centroid and radius for each maximum
		{
			ProgressIndicator progress(maximaList.size());
#pragma omp parallel for if(maximaList.size() > PARALLELIZATION_THRESHOLD)
			for (int n = 0; n < maximaList.size(); n++)
			{
				const auto& points = maximaList[n];
				double radius = (double)orig(points[0]);
				radii[n] = radius;

				Vec3d centroid;
				for (size_t m = 0; m < points.size(); m++)
				{
					centroid += Vec3d(points[m]);
				}
				centroid /= (double)points.size();
				centroids[n] = centroid;

				// removalFlags[n] is already false due to the initialization above.

				progress.step();
			}
		}

		double radiusMultiplier2 = radiusMultiplier * radiusMultiplier;

		// Calculate removal flag for each maximum
		{
			ProgressIndicator progress(maximaList.size());
#pragma omp parallel for if(maximaList.size() > PARALLELIZATION_THRESHOLD)
			for (int n = 0; n < maximaList.size(); n++)
			{
				// Find all maxima inside radius of this maximum and mark those to be removed if their radius is less than radius of the current maximum.
				Vec3d cn = centroids[n];
				double rn = radii[n];
				double rn2 = rn * rn;

				for (size_t m = n + 1; m < maximaList.size(); m++)
				{
					// Version 2, remove small maxima and combine ones with equal radius.
					double distance2 = (centroids[m] - cn).normSquared();
					double rm = radii[m];
					

					// n -> m comparison
					if (distance2 < radiusMultiplier2 * rn2)
					{
						// At this point maximum m is near maximum n.

						if (NumberUtils<double>::equals(rm, rn))
						{
							// Maxima m and n are equal in radius. Combine them.
#pragma omp critical(cleanmaxima_union)
							{
								sets.union_sets(n, m);
							}
						}
						else if(rm < rn)
						{
							// Maximum m is smaller than maximum n => remove maximum m.
							removalFlags[m] = true;
						}
					}

					double rm2 = rm * rm;

					// m -> n comparison (for symmetry due to multiplier * radii[m] term)
					if (distance2 < radiusMultiplier2 * rm2)
					{
						// At this point maximum n is near maximum m.

						if (NumberUtils<double>::equals(rn, rm))
						{
							// Maxima n and m are equal in radius. Combine them.
#pragma omp critical(cleanmaxima_union)
							{
								sets.union_sets(m, n);
							}
						}
						else if(rn < rm)
						{
							// Maximum n is smaller than maximum m => remove maximum n.
							removalFlags[n] = true;
						}
					}
				}

				progress.step();
			}
		}

		// Version 2: No maxima points are removed, but they are combined instead.
		{
			ProgressIndicator progress(maximaList.size());
			for (size_t n = 0; n < maximaList.size(); n++)
			{
				size_t root = sets.find_set(n);
				if (root != n)
				{
					// NOTE: Enable this if to make the entire maximum be removed if one of its members is marked to be removed.
					// If not enabled, the result depends on the order of maxima, e.g., in case, where max A erases max B,
					// and max B is combined to max C. Here, if root is max B, the entire combined max will be removed,
					// but if root is max C, it will remain. This is usually not good.
					if (removalFlags[n])
						removalFlags[root] = true;

					removalFlags[n] = true;
					for (size_t m = 0; m < maximaList[n].size(); m++)
						maximaList[root].push_back(maximaList[n][m]);
				}

				progress.step();
			}
		}

		// Remove all flagged maxima
		std::vector<std::vector<Vec3sc> > newMaximaList;
		for (size_t n = 0; n < maximaList.size(); n++)
		{
			if (!removalFlags[n])
				newMaximaList.push_back(maximaList[n]);
		}

		maximaList = newMaximaList;
	}


	namespace tests
	{
		void localMaxima();
	}
}
