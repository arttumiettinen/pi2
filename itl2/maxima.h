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
						floodfill(orig, p, tempColor, tempColor, connectivity, nullptr, &filledPoints, 0, &neighbourValues);

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

			showProgress(z, orig.depth());
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
		std::vector<Vec3d> centroids;
		centroids.reserve(maximaList.size());

		std::vector<double> radii;
		radii.reserve(maximaList.size());

		std::vector<bool> removalFlags;
		removalFlags.reserve(maximaList.size());

		IndexForest sets(maximaList.size());

		// Calculate centroid and radius for each maximum
		for (size_t n = 0; n < maximaList.size(); n++)
		{
			const auto& points = maximaList[n];
			double radius = (double)orig(points[0]);
			radii.push_back(radius);

			Vec3d centroid;
			for (size_t m = 0; m < points.size(); m++)
			{
				centroid += Vec3d(points[m]);
			}
			centroid /= (double)points.size();
			centroids.push_back(centroid);

			removalFlags.push_back(false);
		}

		// Calculate removal flag for each maximum
		for (size_t n = 0; n < maximaList.size(); n++)
		{
			// Find all maxima inside radius of this maximum and mark those to be removed if their radius is less than radius of the current maximum.

			Vec3d p0 = centroids[n];
			double radius = radii[n];

			for (size_t m = 0; m < maximaList.size(); m++)
			{
				/*
				Version 1, only remove smaller maxima but do not combine.
				if(m != n && radii[m] < radius)
				{
					if((centroids[m] - p0).norm() < radius)
					{
						removalFlags[m] = true;
					}
				}
				*/
				// Version 2, remove small maxima and combine ones with equal radius.
				if (m != n)
				{
					double distance = (centroids[m] - p0).norm();
					if (distance < radiusMultiplier * radius)
					{
						// At this point maximum m is near maximum n.

						if (NumberUtils<double>::equals(radii[m], radius))
						{
							// Maxima m and n are equal in radius. Combine them.
							sets.union_sets(n, m);
						}
						else if (radii[m] < radius)
						{
							// Maximum m is smaller than maximum n => remove maximum m.
							removalFlags[m] = true;
						}
					}
				}
			}
		}

		/*
		// Version 1: combined maximum = single point at common centroid of all combined maxima.

		// Remove non-root maxima
		for(size_t n = 0; n < maximaList.size(); n++)
		{
			size_t root = sets.find_set(n);
			if(root != n)
				removalFlags[n] = true;
		}

		// Calculate new centroid for each maximum
		std::vector<Tvecd> newCentroids;
		std::vector<double> counts;

		newCentroids.reserve(maximaList.size());
		counts.reserve(maximaList.size());

		for(size_t n = 0; n < maximaList.size(); n++)
		{
			newCentroids.push_back(Tvecd());
			counts.push_back(0);
		}

		for(size_t n = 0; n < maximaList.size(); n++)
		{
			size_t root = sets.find_set(n);
			newCentroids[root] += centroids[n];
			counts[root]++;
		}

		for(size_t n = 0; n < maximaList.size(); n++)
		{
			newCentroids[n] /= counts[n];
		}

		// If new centroid is different from old one, replace points of the maximum by one point corresponding
		// to the new centroid
		for(size_t n = 0; n < maximaList.size(); n++)
		{
			if(!centroids[n].equals(newCentroids[n]))
			{
				// Replace points of the maximum
				maximaList[n].clear();
				maximaList[n].push_back(round(newCentroids[n]));
				centroids[n] = newCentroids[n];
			}
		}
		*/

		// Version 2: No maxima points are removed, but they are combined instead.
		for (size_t n = 0; n < maximaList.size(); n++)
		{
			size_t root = sets.find_set(n);
			if (root != n)
			{
				removalFlags[n] = true;
				for (size_t m = 0; m < maximaList[n].size(); m++)
					maximaList[root].push_back(maximaList[n][m]);
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
