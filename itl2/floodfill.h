#pragma once

#include <set>
#include <vector>
#include <queue>
#include <tuple>
#include <iostream>
#include <limits>

#include "image.h"
#include "math/vec3.h"
#include "connectivity.h"
#include "progress.h"
#include "heaps.h"

namespace itl2
{

	namespace internals
	{
		template<typename pixel_t> bool processNeighbours(coord_t x, coord_t y, coord_t z, std::queue<Vec3sc>& points, std::vector<std::tuple<coord_t, coord_t, bool> >& nbs, const Image<pixel_t>& image, pixel_t fillColor, pixel_t origColor, pixel_t stopColor, std::set<pixel_t>* pNeighbouringColors)
		{
			for (auto& nb : nbs)
			{
				coord_t deltay = std::get<0>(nb);
				coord_t deltaz = std::get<1>(nb);
				bool& active = std::get<2>(nb);

				Vec3c xn(x, y + deltay, z + deltaz);
				if (image.isInImage(xn))
				{
					pixel_t p = image(xn);
					if (!active && p == origColor)
					{
						points.push(Vec3sc(xn));
						active = true;
					}
					else if (active && p != origColor)
					{
						active = false;
					}

					if (p != origColor)
					{
						if (fillColor != stopColor && p == stopColor)
							return false;

						if (pNeighbouringColors != 0 && p != origColor && p != fillColor)
							pNeighbouringColors->insert(p);
					}
				}
			}

			return true;
		}
	}

	/**
	Flood fill beginning from the given seed points.
	@param origColor Original color that we are filling. (the color of the region where the fill is allowed to proceed)
	@param fillColor Fill color. The filled pixels will be colored with this color.
	@param stopColor Set to value different from fillColor to stop filling when a pixel of this color is encountered. This argument is used for efficient implementation of small region removal.
	@param pNeighbouringColors Pointer to a set that will contain colors neighbouring the filled region. Set to null not to collect this information. Values of seed points that are not origColor, fillColor, or stopColor are added to the set, too.
	@return True if the fill was terminated naturally; false if the fill was terminated by reaching fillLimit in filled pixel count; by encountering pixel with stopColor value; or if the origColor is fillColor.
	*/
	template<typename pixel_t> bool floodfillSingleThreaded(Image<pixel_t>& image, const std::vector<Vec3sc>& seeds, pixel_t origColor, pixel_t fillColor, pixel_t stopColor, Connectivity connectivity = Connectivity::NearestNeighbours, size_t* pFilledPointCount = nullptr, std::vector<Vec3sc>* pFilledPoints = nullptr, size_t fillLimit = 0, std::set<pixel_t>* pNeighbouringColors = nullptr)
	{
		if (pFilledPointCount)
			*pFilledPointCount = 0;

		if (pFilledPoints)
			pFilledPoints->clear();

		if (pNeighbouringColors)
			pNeighbouringColors->clear();

		if (origColor == fillColor)
			return false;

		if (fillLimit <= 0)
			fillLimit = std::numeric_limits<size_t>::max();

		// Contains {deltay, deltaz, active} for all neighbouring scanlines.
		std::vector<std::tuple<coord_t, coord_t, bool> > nbs;
		if (connectivity == Connectivity::NearestNeighbours)
		{
			nbs = { {1, 0, false}, {-1, 0, false}, {0, 1, false,}, {0, -1, false} };
		}
		else if (connectivity == Connectivity::AllNeighbours)
		{
			nbs = { {1, 0, false}, {-1, 0, false}, {0, 1, false,}, {0, -1, false}, {1, 1, false}, {1, -1, false}, {-1, 1, false}, {-1, -1, false} };
		}
		else
		{
			throw ITLException("Unsupported connectivity value.");
		}

		std::queue<Vec3sc> points;
		for (const Vec3sc& v : seeds)
		{
			if (image.isInImage(v))
			{
				pixel_t p = image(v);
				if (fillColor != stopColor && p == stopColor)
					return false;

				if (pNeighbouringColors != 0 && p != origColor && p != fillColor)
					pNeighbouringColors->insert(p);

				points.push(v);
			}
		}

		size_t tmp = 0;
		size_t* pCount = &tmp;
		if (pFilledPointCount)
			pCount = pFilledPointCount;

		ProgressIndicator progress("");
		while (!points.empty())
		{
			const Vec3c p = Vec3c(points.front());
			
			// Check that this point has not been filled before (there might be multiple routes to the same location).
			if (image(p) == origColor)
			{

				coord_t xl = p.x;
				coord_t y = p.y;
				coord_t z = p.z;

				while (xl >= 0 && image(xl, y, z) == origColor)
					xl--;

				// Stop color check and neighbouring point set update
				if (xl >= 0)
				{
					pixel_t p = image(xl, y, z);
					if (fillColor != stopColor && p == stopColor)
						return false;

					if (pNeighbouringColors != 0 && p != origColor && p != fillColor)
						pNeighbouringColors->insert(p);
				}

				xl++;

				// Set Active flags to zero for all neighbour directions.
				for (auto& nb : nbs)
					std::get<2>(nb) = false;

				// Fill neighbouring rows (don't fill the first and last pixels at xl-1 and end pos+1 if doing filling with All connectivity.
				if (connectivity == Connectivity::AllNeighbours && xl > 0)
				{
					if (!internals::processNeighbours(xl - 1, y, z, points, nbs, image, fillColor, origColor, stopColor, pNeighbouringColors))
						return false;
				}

				while (xl < image.width() && image(xl, y, z) == origColor)
				{
					image(xl, y, z) = fillColor;
					(*pCount)++;
					if (pFilledPoints)
						pFilledPoints->push_back(Vec3sc((int32_t)xl, (int32_t)y, (int32_t)z));

					// Fill volume limit check
					if (*pCount >= fillLimit)
						return false;

					if(!internals::processNeighbours(xl, y, z, points, nbs, image, fillColor, origColor, stopColor, pNeighbouringColors))
						return false;

					xl++;
				}

				if (connectivity == Connectivity::AllNeighbours && xl < image.width())
				{
					if (!internals::processNeighbours(xl, y, z, points, nbs, image, fillColor, origColor, stopColor, pNeighbouringColors))
						return false;
				}


				// Stop color check and neighbouring point set update
				if (xl < image.width())
				{
					pixel_t p = image(xl, y, z);
					if (fillColor != stopColor && p == stopColor)
						return false;

					if (pNeighbouringColors != 0 && p != origColor && p != fillColor)
						pNeighbouringColors->insert(p);
				}
			}

			points.pop();

			// Progress report for large fills
			size_t s = points.size();
			if (s > 0 && s % 50000 == 0)
			{
				progress.showMessage(toString(s) + " seeds");
			}
		}

		return true;
	}

	/**
	Perform flood fill.
	Uses single-threaded scanline fill algorithm.
	@param image Image containing the geometry to be filled.
	@param start Starting position.
	@param fillColor Fill color. The filled pixels will be colored with this color.
	@param stopColor Set to value different from fillColor to stop filling when a pixel of this color is encountered.
	@param connectivity Connectivity of the fill.
	@param pFilledPoints Pointer to vector that will receive the coordinates of the filled points. Set to zero if this information is not required.
	@param fillLimit Set to value to limit count of filled points to that value. Used in regionremoval code. Set to std::numeric_limits<size_t>::max() to allow (practically) any number of particles.
	@param pNeighbouringColors Pointer to set that will contain the colors of non-filled points neighbouring the filled region.
	@return True if the fill was terminated naturally; false if the fill was terminated by reaching fillLimit in filled pixel count; by encountering pixel with stopColor value; or if the origColor is fillColor.
	*/
	template<typename pixel_t> bool floodfillSingleThreaded(Image<pixel_t>& image, const Vec3c& start, pixel_t fillColor, pixel_t stopColor, Connectivity connectivity = Connectivity::NearestNeighbours, size_t* pFilledPointCount = nullptr, std::vector<Vec3sc>* pFilledPoints = nullptr, size_t fillLimit = 0, std::set<pixel_t>* pNeighbouringColors = nullptr)
	{
		if (!image.isInImage(start))
			return true;

		pixel_t origColor = image(start);
		std::vector<Vec3sc> seeds;
		seeds.push_back(Vec3sc(start));

		return floodfillSingleThreaded(image, seeds, origColor, fillColor, stopColor, connectivity, pFilledPointCount, pFilledPoints, fillLimit, pNeighbouringColors);
	}

	/**
	Perform flood fill using a multi-threaded algorithm.
	@param image Image containing the geometry to be filled.
	@param seeds Initial seed points.
	@param origColor The color to fill.
	@param fillColor Fill color. The filled pixels will be colored with this color.
	@param stopColor Set to value different from fillColor to stop filling when a pixel of this color is encountered.
	@param connectivity Connectivity of the fill.
	@param pFilledPoints Pointer to vector that will receive the coordinates of the filled points. Set to zero if this information is not required.
	@param fillLimit Set to value to limit count of filled points to that value. Used in regionremoval code. Set to std::numeric_limits<size_t>::max() to allow (practically) any number of particles.
	@param pNeighbouringColors Pointer to set that will contain the colors of non-filled points neighbouring the filled region.
	@param minBlockSize Minimum z-directional size of a block assigned to a single thread.
	@return True if the fill was terminated naturally; false if the fill was terminated by reaching fillLimit in filled pixel count; by encountering pixel with stopColor value; or if the origColor is fillColor.
	*/
	template<typename pixel_t> bool floodfillBlocks(Image<pixel_t>& image, const std::vector<Vec3sc>& seeds, pixel_t origColor, pixel_t fillColor, pixel_t stopColor, Connectivity connectivity = Connectivity::NearestNeighbours, size_t* pFilledPointCount = nullptr, std::vector<Vec3sc>* pFilledPoints = nullptr, size_t fillLimit = 0, std::set<pixel_t>* pNeighbouringColors = nullptr,
		coord_t minBlockSize = 700)
	{
		// Algorithm:
		// divide image into blocks
		// process each block that contains seed points
		//	flood fill
		// 	store points that are on the edge of the block and that changed as new seed points
		// remove seed points on the edge of the full image
		// repeat fill for the new set of seed points

		if (pFilledPointCount)
			*pFilledPointCount = 0;

		if (pFilledPoints)
			pFilledPoints->clear();

		if (pNeighbouringColors)
			pNeighbouringColors->clear();

		if (origColor == fillColor)
			return false;

		if (fillLimit <= 0)
			fillLimit = std::numeric_limits<size_t>::max();

		// The threads will use prevSeeds array as seed points
		std::vector<Vec3sc> prevSeeds;
		prevSeeds.insert(prevSeeds.end(), seeds.begin(), seeds.end());

		// The threads fill this array with the new seeds for the next filling round
		std::vector<Vec3sc> newSeeds;

		bool totalResult = true;

		ProgressIndicator progress("");

		size_t round = 0;
		do
		{
			round++;

			#pragma omp parallel
			{
				int idx = omp_get_thread_num();
				int count = omp_get_num_threads();

				if (idx == 0)
				{
					progress.showMessage(string("Flood-filling in maximum of ") + toString(count) + " blocks, iteration " + toString(round));
				}

				// This barrier is required such that it is guaranteed that outputs above are shown before thread-specific
				// outputs.
				#pragma omp barrier

				// Calculate amount of slices single thread should process.
				// The value is capped so that we don't divide the work unnecessarily too much if the image is small.
				coord_t size = image.depth() / count;
				if (size < minBlockSize)
					size = minBlockSize;

				// If there's nothing to do for all the threads, the excess threads will just skip processing.
				coord_t minZ = idx * size;
				if (minZ < 0)
					minZ = 0;
				if (minZ < image.depth())
				{

					coord_t maxZ = minZ + size - 1;
					if (maxZ >= image.depth())
						maxZ = image.depth() - 1;

					if (idx == count - 1)
					{
						// The last thread processes possible "rounding error" slices
						maxZ = image.depth() - 1;
					}

					//#pragma omp critical
					//					{
					//						std::cout << "Thread " << idx << " analyzes z-range [" << minZ << ", " << maxZ << "[" << std::endl;
					//					}

					// Get view of a part of the image
					Image<pixel_t> block(image, minZ, maxZ);

					// Find seeds for this block, and convert to block coordinates.
					std::vector<Vec3sc> blockSeeds;
					for (Vec3sc v : prevSeeds)
					{
						v.z -= (int32_t)minZ;
						if (block.isInImage(v))
							blockSeeds.push_back(v);
					}

					if (blockSeeds.size() > 0)
					{
						progress.showMessage(string("Thread ") + toString(idx) + " has " + toString(blockSeeds.size()) + " seeds to process.");
					}

					if (blockSeeds.size() > 0)
					{

						// Iterate block edge values and save to a list.
						// NOTE: It is necessary to evaluate only the z-faces as the block division is done in z only.
						std::vector<pixel_t> edgeValues;
						//edgeValues.reserve(block.edgePixelCount());
						//forEdges(block.bounds(), image.dimensionality(), [&](coord_t x, coord_t y, coord_t z)
						edgeValues.reserve(block.width() * block.height() * 2);
						for (coord_t z = 0; z < block.depth(); z += block.depth() - 1)
						{
							for (coord_t y = 0; y < block.height(); y++)
							{
								for (coord_t x = 0; x < block.width(); x++)
								{
									edgeValues.push_back(block(x, y, z));
								}
							}
						}

						// Flood fill in block
						size_t blockFilledPointCount;
						size_t* pBlockFilledPointCount = pFilledPointCount ? &blockFilledPointCount : nullptr;
						std::vector<Vec3sc> blockFilledPoints;
						std::vector<Vec3sc>* pBlockFilledPoints = pFilledPoints ? &blockFilledPoints : nullptr;
						std::set<pixel_t> blockNeighbouringColors;
						std::set<pixel_t>* pBlockNeighbouringColors = pNeighbouringColors ? &blockNeighbouringColors : nullptr;
						bool blockResult = floodfillSingleThreaded(block, blockSeeds, origColor, fillColor, stopColor, connectivity, pBlockFilledPointCount, pBlockFilledPoints, fillLimit, pBlockNeighbouringColors);

						// Iterate block edge values and compare to saved, make a list of points that changed.
						// Neighbours of the changed points that are not in the current block are new seeds for the neighbouring block.
						std::vector<Vec3sc> blockNewSeeds;
						size_t n = 0;
						//forEdges(block.bounds(), image.dimensionality(), [&](coord_t x, coord_t y, coord_t z)
						for (coord_t z = 0; z < block.depth(); z += block.depth() - 1)
						{
							int32_t dz = 0;
							if (z <= 0)
							{
								dz = -1;
							}
							else if (z >= block.depth() - 1)
							{
								dz = 1;
							}
							else
							{
								// This happens when iterating over non-z faces
							}

							int32_t newZ = (int32_t)z + (int32_t)minZ + dz;
							
							for (coord_t y = 0; y < block.height(); y++)
							{
								for (coord_t x = 0; x < block.width(); x++)
								{
									if (dz != 0 && newZ >= 0 && newZ < image.depth())
									{
										if (edgeValues[n] != block(x, y, z))
										{
											Vec3sc p((int32_t)x, (int32_t)y, newZ);
											if (connectivity == Connectivity::NearestNeighbours)
											{
												blockNewSeeds.push_back(p + Vec3sc(0, 0, 0));
											}
											else if (connectivity == Connectivity::AllNeighbours)
											{
												blockNewSeeds.push_back(p + Vec3sc(-1, -1, 0));
												blockNewSeeds.push_back(p + Vec3sc(0, -1, 0));
												blockNewSeeds.push_back(p + Vec3sc(1, -1, 0));
												blockNewSeeds.push_back(p + Vec3sc(-1, 0, 0));
												blockNewSeeds.push_back(p + Vec3sc(0, 0, 0));
												blockNewSeeds.push_back(p + Vec3sc(1, 0, 0));
												blockNewSeeds.push_back(p + Vec3sc(-1, 1, 0));
												blockNewSeeds.push_back(p + Vec3sc(0, 1, 0));
												blockNewSeeds.push_back(p + Vec3sc(1, 1, 0));
											}
											else
											{
												throw std::logic_error("Unsupported connectivity value in flood fill.");
											}
										}
									}
									n++;
								}
							}
						}


						#pragma omp critical(floodfill_blocks)
						{
							newSeeds.insert(newSeeds.end(), blockNewSeeds.begin(), blockNewSeeds.end());
							if (pFilledPointCount)
								*pFilledPointCount += blockFilledPointCount;
							if (pFilledPoints)
							{
								for (Vec3sc v : blockFilledPoints)
									pFilledPoints->push_back(v + Vec3sc(0, 0, (int32_t)minZ));
							}
							if (pNeighbouringColors)
							{
								pNeighbouringColors->merge(blockNeighbouringColors);
							}
							if (blockResult == false && totalResult == true)
								totalResult = false;
						}
					}

				}
			}

			std::swap(prevSeeds, newSeeds);
			newSeeds.clear();
		} while (totalResult == true && prevSeeds.size() > 0);

		return totalResult;
	}

	/**
	Perform flood fill using a multi-threaded algorithm.
	@param image Image containing the geometry to be filled.
	@param start Starting position.
	@param fillColor Fill color. The filled pixels will be colored with this color.
	@param stopColor Set to value different from fillColor to stop filling when a pixel of this color is encountered.
	@param connectivity Connectivity of the fill.
	@param pFilledPoints Pointer to vector that will receive the coordinates of the filled points. Set to zero if this information is not required.
	@param fillLimit Set to value to limit count of filled points to that value. Used in regionremoval code. Set to std::numeric_limits<size_t>::max() to allow (practically) any number of particles.
	@param pNeighbouringColors Pointer to set that will contain the colors of non-filled points neighbouring the filled region.
	@param minBlockSize Minimum z-directional size of a block assigned to a single thread.
	@return True if the fill was terminated naturally; false if the fill was terminated by reaching fillLimit in filled pixel count; by encountering pixel with stopColor value; or if the origColor is fillColor.
	*/
	template<typename pixel_t> bool floodfillBlocks(Image<pixel_t>& image, const Vec3c& start, pixel_t fillColor, pixel_t stopColor, Connectivity connectivity = Connectivity::NearestNeighbours, size_t* pFilledPointCount = nullptr, std::vector<Vec3sc>* pFilledPoints = nullptr, size_t fillLimit = 0, std::set<pixel_t>* pNeighbouringColors = nullptr,
		bool showProgressInfo = true, coord_t minBlockSize = 700)
	{
		if (!image.isInImage(start))
			return true;

		pixel_t origColor = image(start);
		std::vector<Vec3sc> seeds;
		seeds.push_back(Vec3sc(start));

		return floodfillBlocks(image, seeds, origColor, fillColor, stopColor, connectivity, pFilledPointCount, pFilledPoints, fillLimit, pNeighbouringColors, minBlockSize);
	}

	/**
	Perform flood fill.
	Chooses between single- and multi-threaded implementation based on image size.
	@param image Image containing the geometry to be filled.
	@param start Starting position.
	@param fillColor Fill color. The filled pixels will be colored with this color.
	@param stopColor Set to value different from fillColor to stop filling when a pixel of this color is encountered.
	@param connectivity Connectivity of the fill.
	@param pFilledPoints Pointer to vector that will receive the coordinates of the filled points. Set to zero if this information is not required.
	@param fillLimit Set to value to limit count of filled points to that value. Used in regionremoval code. Set to std::numeric_limits<size_t>::max() to allow (practically) any number of particles.
	@param pNeighbouringColors Pointer to set that will contain the colors of non-filled points neighbouring the filled region.
	@return True if the fill was terminated naturally; false if the fill was terminated by reaching fillLimit in filled pixel count; by encountering pixel with stopColor value; or if the origColor is fillColor.
	*/
	template<typename pixel_t> bool floodfill(Image<pixel_t>& image, const Vec3c& start, pixel_t fillColor, pixel_t stopColor, Connectivity connectivity = Connectivity::NearestNeighbours, size_t* pFilledPointCount = nullptr, std::vector<Vec3sc>* pFilledPoints = nullptr, size_t fillLimit = 0, std::set<pixel_t>* pNeighbouringColors = nullptr)
	{
		if (image.pixelCount() < 700 * 700 * 700 || image.dimensionality() < 3)
			return floodfillSingleThreaded(image, start, fillColor, stopColor, connectivity, pFilledPointCount, pFilledPoints, fillLimit, pNeighbouringColors);
		else
			return floodfillBlocks(image, start, fillColor, stopColor, connectivity, pFilledPointCount, pFilledPoints, fillLimit, pNeighbouringColors);
	}

	/**
	Perform flood fill.
	Chooses between single- and multi-threaded implementation based on image size.
	@param image Image containing the geometry to be filled.
	@param seeds Initial seed points.
	@param fillColor Fill color. The filled pixels will be colored with this color.
	@param stopColor Set to value different from fillColor to stop filling when a pixel of this color is encountered.
	@param connectivity Connectivity of the fill.
	@param pFilledPoints Pointer to vector that will receive the coordinates of the filled points. Set to zero if this information is not required.
	@param fillLimit Set to value to limit count of filled points to that value. Used in regionremoval code. Set to std::numeric_limits<size_t>::max() to allow (practically) any number of particles.
	@param pNeighbouringColors Pointer to set that will contain the colors of non-filled points neighbouring the filled region.
	@return True if the fill was terminated naturally; false if the fill was terminated by reaching fillLimit in filled pixel count; by encountering pixel with stopColor value; or if the origColor is fillColor.
	*/
	template<typename pixel_t> bool floodfill(Image<pixel_t>& image, const std::vector<Vec3sc>& seeds, pixel_t origColor, pixel_t fillColor, pixel_t stopColor, Connectivity connectivity = Connectivity::NearestNeighbours, size_t* pFilledPointCount = nullptr, std::vector<Vec3sc>* pFilledPoints = nullptr, size_t fillLimit = 0, std::set<pixel_t>* pNeighbouringColors = nullptr)
	{
		if (image.pixelCount() < 700 * 700 * 700 || image.dimensionality() < 3)
			return floodfillSingleThreaded(image, seeds, origColor, fillColor, stopColor, connectivity, pFilledPointCount, pFilledPoints, fillLimit, pNeighbouringColors);
		else
			return floodfillBlocks(image, seeds, origColor, fillColor, stopColor, connectivity, pFilledPointCount, pFilledPoints, fillLimit, pNeighbouringColors);
	}

	/**
	Flood fill beginning from the given seed points.
	Does not use scan line fill algorithm, may be faster than scan line fill for very small images.
	@param origColor Original color that we are filling. (the color of the region where the fill is allowed to proceed)
	@param fillColor Fill color. The filled pixels will be colored with this color.
	@param stopColor Set to value different from fillColor to stop filling when a pixel of this color is encountered.
	@return Count of pixels filled if the fill was terminated naturally; -1 times count of pixels filled if the fill was terminated by reaching fillLimit in filled pixel count or by encountering pixel with stopColor value.
	*/
	template<typename pixel_t> coord_t slowFloodfill(Image<pixel_t>& image, const std::vector<Vec3sc>& seeds, pixel_t origColor, pixel_t fillColor, pixel_t stopColor, Connectivity connectivity = Connectivity::NearestNeighbours, std::vector<Vec3sc>* pFilledPoints = nullptr, size_t fillLimit = std::numeric_limits<size_t>::max(), std::set<pixel_t>* pNeighbouringColors = nullptr)
	{
		std::queue<Vec3sc> points;
		for (const Vec3sc& p : seeds)
			points.push(p);

		size_t lastPrinted = 0;
		coord_t filledPoints = 0;
		while (!points.empty())
		{
			const Vec3sc p = points.front();

			pixel_t pixelValue = image(p);
			if (pixelValue == origColor)
			{
				filledPoints++;
				image(p) = fillColor;
				if (pFilledPoints)
					pFilledPoints->push_back(p);

				// Add items to queue
				if (connectivity == Connectivity::NearestNeighbours)
				{
					for (size_t n = 0; n < 3; n++)
					{
						if (p[n] > 0)
						{
							Vec3sc np = p;
							np[n]--;
							points.push(np);
						}

						if (p[n] < image.dimension(n) - 1)
						{
							Vec3sc np = p;
							np[n]++;
							points.push(np);
						}
					}
				}
				else
				{
					// All neighbours
					for (int32_t dx = -1; dx <= 1; dx++)
					{
						for (int32_t dy = -1; dy <= 1; dy++)
						{
							for (int32_t dz = -1; dz <= 1; dz++)
							{
								Vec3sc np(p.x + dx, p.y + dy, p.z + dz);
								if (np.x >= 0 && np.y >= 0 && np.z >= 0 &&
									np.x < image.dimension(0) &&
									np.y < image.dimension(1) &&
									np.z < image.dimension(2)
									)
									points.push(np);
							}
						}
					}
				}
			}
			else if (fillColor != stopColor && pixelValue == stopColor)
			{
				// Stop color has been encountered.
				return -filledPoints;
			}
			else
			{
				if (pNeighbouringColors != 0)
					pNeighbouringColors->insert(pixelValue);
			}

			points.pop();

			if ((size_t)filledPoints > fillLimit)
				return -filledPoints; // Fill volume limit has been encountered.

			// Progress report for large fills
			size_t s = points.size();
			if (s > 0 && s % 50000 == 0 && lastPrinted != s)
			{
				lastPrinted = s;
				std::cout << s << " seed points...\r" << std::flush;
			}
		}

		if (lastPrinted != 0)
			std::cout << std::endl;

		return filledPoints;
	}

	template<typename pixel_t> coord_t slowFloodfill(Image<pixel_t>& image, const Vec3c& start, pixel_t fillColor, pixel_t stopColor, Connectivity connectivity = Connectivity::NearestNeighbours, std::vector<Vec3sc>* pFilledPoints = nullptr, size_t fillLimit = std::numeric_limits<size_t>::max(), std::set<pixel_t>* pNeighbouringColors = nullptr)
	{
		pixel_t origColor = image(start);
		if (origColor == fillColor)
			return true;

		std::vector<Vec3sc> seeds;
		seeds.push_back(Vec3sc(start));

		return slowFloodfill(image, seeds, origColor, fillColor, stopColor, connectivity, pFilledPoints, fillLimit, pNeighbouringColors);
	}

	/**
	Grows regions colored with sourceColor into regions colored with allowedColor.
	@return Number of pixels whose color changed.
	*/
	template<typename pixel_t> coord_t grow(Image<pixel_t>& image, const pixel_t sourceColor, const pixel_t allowedColor, Connectivity connectivity = Connectivity::NearestNeighbours)
	{
		// First find seed points
		std::vector<Vec3sc> seeds;
		for (coord_t z = 0; z < image.depth(); z++)
		{
			for (coord_t y = 0; y < image.height(); y++)
			{
				for (coord_t x = 0; x < image.width(); x++)
				{
					if (image(x, y, z) == sourceColor)
					{
						// Only add a seed if it has neighbour having value 'targetColor'

						Vec3sc p((int32_t)x, (int32_t)y, (int32_t)z);
						if (connectivity == Connectivity::NearestNeighbours)
						{
							for (size_t n = 0; n < 3; n++)
							{
								if (p[n] > 0)
								{
									Vec3c np(p);
									np[n]--;
									if (image.isInImage(np) && image(np) == allowedColor)
									{
										seeds.push_back(p);

										// Skip further processing of the loops
										break;
									}
								}

								if (p[n] < image.dimension(n) - 1)
								{
									Vec3c np(p);
									np[n]++;
									if (image.isInImage(np) && image(np) == allowedColor)
									{
										seeds.push_back(p);

										// Skip further processing of the loops
										break;
									}
								}
							}
						}
						else
						{
							for (coord_t dz = -1; dz <= 1; dz++)
							{
								for (coord_t dy = -1; dy <= 1; dy++)
								{
									for (coord_t dx = -1; dx <= 1; dx++)
									{
										Vec3c np(x + dx, y + dy, z + dz);
										if (image.isInImage(np) && image(np) == allowedColor)
										{
											seeds.push_back(p);

											// Skip further processing of the loops
											dz = 2;
											dy = 2;
											dx = 2;
											break;
										}
									}
								}
							}
						}

					}
				}
			}
		}

		// Set all seed points to allowedColor so that the flood fill does not stop immediately because of wrongly colored seeds.
		size_t seedCount = seeds.size();
		if (seedCount <= 0)
			return 0;

		for (const Vec3sc& p : seeds)
			image(p) = allowedColor;

		// Flood fill from all the seed points at once.
		size_t filled;
		floodfill(image, seeds, allowedColor, sourceColor, sourceColor, connectivity, &filled);
		return filled - seedCount;
	}


	namespace internals
	{
		/**
		Storage class for fill point priority queue in Meyer's algorithm.
		*/
		template<typename label_t, typename weight_t> class MeyerSeed
		{
		private:
			Vec3sc pos;

			label_t targetLabel;
			weight_t myWeight;
			size_t birthday;

		public:

			/**
			Constructor
			@param p The point.
			@param label Label color of the point.
			@param w Weight of the point. Used to prioritize points with larger weight before points with smaller weight.
			@param birthday The filling round number. Used to prioritize older points before newer points so that points near seeds are filled first.
			*/
			MeyerSeed(const Vec3sc& p, label_t label, weight_t w, size_t birthday) :
				pos(p),
				targetLabel(label),
				myWeight(w),
				birthday(birthday)
			{
			}

			/**
			Gets position.
			*/
			const Vec3sc& position() const
			{
				return pos;
			}

			/**
			Compares weights.
			*/
			bool operator < (const MeyerSeed& right) const
			{
				if (myWeight != right.myWeight)
					return myWeight < right.myWeight;
				else
					return birthday > right.birthday;
			}

			/**
			Gets label value.
			*/
			const label_t label() const
			{
				return targetLabel;
			}

			const weight_t weight() const
			{
				return myWeight;
			}

		};
	}


	/**
Region grow segmentation.
The argument images must be of the same size.
Uses Meyer's flooding algorithm, implementation inspired by OpenCV - see below.
This function is a single-threaded implementation, and will give different results than the multi-threaded version,
as this function is not invariant to image rotations and therefore not invariant to order of seeds.
@param labels Image containing the labels of distinct areas. At input, the image must contain the seed points as nonzero pixels and background as zero pixels; after the algorithm finishes, the image will contain the segmented regions corresponding to the seed points. Multiple seeds may have the same value.
@param weights Image containing the filling priority of each pixel. This image is not modified. If priority is zero or negative, the pixel is never filled.



This function is inspired by OpenCV implementation of watershed, found in file \opencv\modules\imgproc\src\segmentation.cpp.
This was chosen as it is quite fast according to Kornilov - An Overview of Watershed Algorithm Implementations in Open Source Libraries.
OpenCV implementation is licensed with the following text:

							   License Agreement
				For Open Source Computer Vision Library

 Copyright (C) 2000, Intel Corporation, all rights reserved.
 Copyright (C) 2013, OpenCV Foundation, all rights reserved.
 Third party copyrights are property of their respective owners.

 Redistribution and use in source and binary forms, with or without modification,
 are permitted provided that the following conditions are met:

   * Redistribution's of source code must retain the above copyright notice,
	 this list of conditions and the following disclaimer.

   * Redistribution's in binary form must reproduce the above copyright notice,
	 this list of conditions and the following disclaimer in the documentation
	 and/or other materials provided with the distribution.

   * The name of the copyright holders may not be used to endorse or promote products
	 derived from this software without specific prior written permission.

 This software is provided by the copyright holders and contributors "as is" and
 any express or implied warranties, including, but not limited to, the implied
 warranties of merchantability and fitness for a particular purpose are disclaimed.
 In no event shall the Intel Corporation or contributors be liable for any direct,
 indirect, incidental, special, exemplary, or consequential damages
 (including, but not limited to, procurement of substitute goods or services;
 loss of use, data, or profits; or business interruption) however caused
 and on any theory of liability, whether in contract, strict liability,
 or tort (including negligence or otherwise) arising in any way out of
 the use of this software, even if advised of the possibility of such damage.
*/
	template<typename label_t, typename weight_t> void growSingleThreaded(Image<label_t>& labels, const Image<weight_t>& weights, label_t inQueue = std::numeric_limits<label_t>::max())
	{
		// NOTE: The result depends on the order of the seed points_shared, i.e., on the orientation of the image!

		weights.checkSize(labels);
		weights.mustNotBe(labels);

		// Note: This function works with HHeap and GHeap, too, but BucketMap seems to be the fastest choice.
		//HHeap<Vec3c, weight_t> points_shared(min(weights), max(weights), 100);
		//GHeap<Vec3c, weight_t> points_shared;
		BucketMap<Vec3c, weight_t> points;

		// Add all neighbours of seed points_shared to the priority queue
		for (coord_t z = 0; z < labels.depth(); z++)
			//for(coord_t z = labels.depth() - 1; z >= 0; z--)
		{
			for (coord_t y = 0; y < labels.height(); y++)
			{
				for (coord_t x = 0; x < labels.width(); x++)
				{
					Vec3 p(x, y, z);

					if (labels(p) == inQueue)
						throw ITLException(string("Value ") + toString(inQueue) + " must not be used in the labels image, it is needed as a temporary value.");

					if (labels(p) == 0) // Only background points_shared can be neighbours of seeds
					{
						bool hasSeedNeighbour = false;

						for (size_t n = 0; n < labels.dimensionality(); n++)
						{
							for (int delta = -1; delta <= 1; delta += 2)
							{
								Vec3c np = p;
								np[n] += delta;

								if (np[n] >= 0 && np[n] < labels.dimension(n))
								{
									label_t lbl = labels(np);
									weight_t w = weights(np);
									if (w > 0 && lbl != 0 && lbl != inQueue)
									{
										hasSeedNeighbour = true;
									}
								}
							}
						}

						if (hasSeedNeighbour)
						{
							points.push(weights(p), p);
							labels(p) = inQueue;
						}
					}
				}
			}
		}

		size_t round = 0;

		// Grow from the point p to all directions if they are not filled yet.
		while (!points.empty())
		{
			Vec3c p = points.topItem();
			points.pop();

			// Find label for this pixel from the surrounding pixels.
			// Add adjacent unlabeled pixels to the queue.
			label_t currPixelLabel = 0;
			weight_t currPixelLabelWeight = 0;
			for (size_t n = 0; n < labels.dimensionality(); n++)
			{
				for (int delta = -1; delta <= 1; delta += 2)
				{
					Vec3c np = p;
					np[n] += delta;

					if (np[n] >= 0 && np[n] < labels.dimension(n))
					{
						label_t& lbl = labels(np);

						if (lbl == 0)
						{
							// This is a non-labeled neighbour. It goes to the queue.
							weight_t w = weights(np);
							if (w > 0)
							{
								points.push(w, np);

								// Mark the pixel to avoid re-entering it into the queue.
								lbl = inQueue;
							}
						}
						else
						{
							// This is a labeled neighbour, find if the label of this point should be the label of the neighbour.
							if (lbl != inQueue)
							{
								weight_t w = weights(np);
								if (currPixelLabelWeight == 0 || w > currPixelLabelWeight) // If there are multiple possibilities, choose the one with the highest weight.
								{
									currPixelLabel = lbl;
									currPixelLabelWeight = w;
								}
							}
						}
					}
				}
			}

			// Fill current pixel
			labels(p) = currPixelLabel;

			round++;
		}
	}

	/**
	Region grow segmentation.
	The argument images must be of the same size.
	Uses Meyer's flooding algorithm, implementation inspired by OpenCV - see below.
	@param labels Image containing the labels of distinct areas. At input, the image must contain the seed points as nonzero pixels and background as zero pixels; after the algorithm finishes, the image will contain the segmented regions corresponding to the seed points. Multiple seeds may have the same value.
	@param weights Image containing the filling priority of each pixel. This image is not modified. If priority is zero or negative, the pixel is never filled.



	This function is inspired by OpenCV implementation of watershed, found in file \opencv\modules\imgproc\src\segmentation.cpp.
	This was chosen as it is quite fast according to Kornilov - An Overview of Watershed Algorithm Implementations in Open Source Libraries.
	OpenCV implementation is licensed with the following text:

								   License Agreement
					For Open Source Computer Vision Library

	 Copyright (C) 2000, Intel Corporation, all rights reserved.
	 Copyright (C) 2013, OpenCV Foundation, all rights reserved.
	 Third party copyrights are property of their respective owners.

	 Redistribution and use in source and binary forms, with or without modification,
	 are permitted provided that the following conditions are met:

	   * Redistribution's of source code must retain the above copyright notice,
		 this list of conditions and the following disclaimer.

	   * Redistribution's in binary form must reproduce the above copyright notice,
		 this list of conditions and the following disclaimer in the documentation
		 and/or other materials provided with the distribution.

	   * The name of the copyright holders may not be used to endorse or promote products
		 derived from this software without specific prior written permission.

	 This software is provided by the copyright holders and contributors "as is" and
	 any express or implied warranties, including, but not limited to, the implied
	 warranties of merchantability and fitness for a particular purpose are disclaimed.
	 In no event shall the Intel Corporation or contributors be liable for any direct,
	 indirect, incidental, special, exemplary, or consequential damages
	 (including, but not limited to, procurement of substitute goods or services;
	 loss of use, data, or profits; or business interruption) however caused
	 and on any theory of liability, whether in contract, strict liability,
	 or tort (including negligence or otherwise) arising in any way out of
	 the use of this software, even if advised of the possibility of such damage.
	*/
	template<typename label_t, typename weight_t> void grow(Image<label_t>& labels, const Image<weight_t>& weights, label_t inQueue = std::numeric_limits<label_t>::max())
	{
		// NOTE: This function is invariant to the ordering of the seed points_shared & image orientation - compare
		// to the singlethreaded one that is not.

		weights.checkSize(labels);
		weights.mustNotBe(labels);

		// Note: This function works with HHeap and GHeap, too, but BucketMap seems to be the fastest choice.
		//HHeap<Vec3sc, weight_t> points_shared(min(weights), max(weights), 100);
		//GHeap<Vec3sc, weight_t> points_shared;
		BucketMap<Vec3sc, weight_t> points_shared;


		// Add all neighbours of seed points_shared to the priority queue
		forAllPixels(labels, [&](coord_t x, coord_t y, coord_t z)
			{
				Vec3sc p((int32_t)x, (int32_t)y, (int32_t)z);

				if (labels(p) == inQueue)
					throw ITLException(string("Value ") + toString(inQueue) + " must not be used in the labels image, it is needed as a temporary value.");

				if (labels(p) == 0) // Only background points_shared can be neighbours of seeds
				{
					bool hasSeedNeighbour = false;

					for (size_t n = 0; n < labels.dimensionality(); n++)
					{
						for (int delta = -1; delta <= 1; delta += 2)
						{
							Vec3sc np = p;
							np[n] += delta;

							if (np[n] >= 0 && np[n] < labels.dimension(n))
							{
								label_t lbl = labels(np);
								weight_t w = weights(np);
								if (w > 0 && lbl != 0 && lbl != inQueue)
								{
									hasSeedNeighbour = true;
								}
							}
						}
					}

					if (hasSeedNeighbour)
					{
#pragma omp critical(growInit)
						{
							points_shared.push(weights(p), p);
						}
						labels(p) = inQueue;
					}
				}
			});


		//vector<vector<Vec3sc>*> newPointsLists;

		ProgressIndicator progress("");
		size_t prevSize = 0;

		// Grow from the point p to all directions if they are not filled yet.
#pragma omp parallel
		{
			std::vector<Vec3sc> newPoints;
			newPoints.reserve(250);

			//#pragma omp critical
			//			{
			//				newPointsLists.push_back(&newPoints);
			//			}

			std::vector<label_t> newLabels;
			newLabels.reserve(250);

			coord_t idx = omp_get_thread_num();
			coord_t threadCount = omp_get_num_threads();

			while (true)//(!points_shared.empty())
			{
				// NOTE: Loop start and end indices cannot be changed inside the omp parallel region.
				// NOTE: For reason or another, while(!points_shared.empty()) does not seem to work correctly.
				//			Separate if(points_shared.empty()) break; works as expected.


				if (points_shared.empty())
					break;

				newPoints.clear();

				// This block ensures that topItems_shared is not used after it will be destroyed in popTopItems.
				{
					const std::deque<Vec3sc>& topItems_shared = points_shared.topItems();
					coord_t N = topItems_shared.size();

					const coord_t minPointsPerThread = 1;
					coord_t countPerThread = std::max(minPointsPerThread, N / threadCount + 1);
					coord_t starti = idx * countPerThread;
					coord_t endi = starti + countPerThread;
					if (endi > N)
						endi = N;

					if (starti < N)
					{
						newLabels.clear();
						for (coord_t i = starti; i < endi; i++)
							newLabels.push_back(0);

						for (coord_t i = starti; i < endi; i++)
						{
							Vec3sc p = topItems_shared[i];

							label_t currPixelLabel = 0;
							weight_t currPixelLabelWeight = 0;

							// Find label for this pixel from the surrounding pixels.
							// Add adjacent unlabeled pixels to the queue.
							for (size_t n = 0; n < labels.dimensionality(); n++)
							{
								for (int delta = -1; delta <= 1; delta += 2)
								{
									Vec3sc np = p;
									np[n] += delta;

									if (np[n] >= 0 && np[n] < labels.dimension(n))
									{
										label_t& lbl = labels(np);

										if (lbl == 0)
										{
											// This is a non-labeled neighbour. It goes to the queue.
											weight_t w = weights(np);
											if (w > 0)
											{
												newPoints.push_back(np);
											}
										}
										else
										{
											// This is a labeled neighbour, find if the label of this point should be the label of the neighbour.
											if (lbl != inQueue)
											{
												weight_t w = weights(np);
												if (currPixelLabelWeight == 0 || w > currPixelLabelWeight) // If there are multiple possibilities, choose the one with the highest weight.
												{
													currPixelLabel = lbl;
													currPixelLabelWeight = w;
												}
											}
										}
									}
								}
							}

							newLabels[i - starti] = currPixelLabel;
							//newLabels[i] = currPixelLabel;						
						}
					}

#pragma omp barrier

					// Set new labels to the image
					if (starti < N)
					{
						for (coord_t n = starti; n < endi; n++)
						{
							Vec3sc p = topItems_shared[n];

							// Fill current pixel
							labels(p) = newLabels[n - starti];
							//labels(p) = newLabels[n];
						}
					}
				}

#pragma omp barrier

#pragma omp master
				{
					points_shared.popTopItems();


					size_t s = points_shared.size();
					if (std::abs((coord_t)s - (coord_t)prevSize) > 500000)
					{
						progress.showMessage(toString(s) + " seeds");
						prevSize = s;
					}
				}

#pragma omp barrier


				// Alternative 1: Process points in main thread
				//				// Combine all new_points lists, store minimum w for each point.
				//				// Then insert to points_shared.
				//
				//#pragma omp master
				//				{
				//					for (const vector<Vec3sc>* pList : newPointsLists)
				//					{
				//						for (size_t n = 0; n < pList->size(); n++)
				//						{
				//							Vec3sc np = (*pList)[n];
				//
				//							weight_t w = weights(np);
				//
				//							if (labels(np) != inQueue)
				//							{
				//								points_shared.push(w, np);
				//								labels(np) = inQueue;
				//							}
				//						}
				//					}
				//				}


				// Alternative 2: Process points using a critical section. This seems to be a tiny bit faster.
								// Add new points to the points_shared queue.
#pragma omp critical(grow_collect)
				{
					for (size_t n = 0; n < newPoints.size(); n++)
					{
						Vec3sc np = newPoints[n];

						weight_t w = weights(np);

						if (labels(np) != inQueue)
						{
							points_shared.push(w, np);
							labels(np) = inQueue;
						}
					}
				}

				// Do not allow the next round to start until all the points are pushed into the points_shared queue.
#pragma omp barrier

			}
		}


	}

	/**
	Region grow segmentation.
	The argument images must be of the same size.
	Uses Meyer's flooding algorithm.
	This is the "old" implementation from early versions of pi2.
	@param labels Image containing the labels of distinct areas. At input, the image must contain the seed points as nonzero pixels and background as zero pixels; after the algorithm finishes, the image will contain the segmented regions corresponding to the seed points. Multiple seeds may have the same value.
	@param weights Image containing the filling priority of each pixel. This image is not modified. If priority is zero or negative, the pixel is never filled.
	*/
	template<typename label_t, typename weight_t> void growOld(Image<label_t>& labels, const Image<weight_t>& weights)
	{
		weights.checkSize(labels);
		weights.mustNotBe(labels);

		std::priority_queue<internals::MeyerSeed<label_t, weight_t> > points;

		// Add all seed points to the priority queue

		for (coord_t z = 0; z < labels.depth(); z++)
		{
			for (coord_t y = 0; y < labels.height(); y++)
			{
				for (coord_t x = 0; x < labels.width(); x++)
				{
					label_t& p = labels(x, y, z);
					if (p != 0)
					{
						points.push(internals::MeyerSeed<label_t, weight_t>(Vec3sc((int32_t)x, (int32_t)y, (int32_t)z), p, std::numeric_limits<weight_t>::max(), 0));
						//points.push(internals::MeyerSeed<label_t, weight_t>(Vec3sc(x, y, z), p, weights(x, y, z), 0));
						p = 0;
					}
				}
			}
		}

		size_t lastPrinted = 0;
		long round = 0;

		//size_t maxQueueDepth = points.size();
		// Grow from the point p to all directions if they are not filled yet.
		while (!points.empty())
		{
			//maxQueueDepth = std::max(maxQueueDepth, points.size());

			round++;

			const internals::MeyerSeed<label_t, weight_t>& obj = points.top();
			Vec3sc p = obj.position();
			label_t targetLabel = obj.label();

			points.pop();

			// Only proceed if the label of this pixel has not been set yet.
			if (labels(p) == 0)
			{

				labels(p) = targetLabel;

				// Insert neighbours into the priority queue.
				for (size_t n = 0; n < p.size(); n++)
				{
					if (p[n] > 0)
					{
						Vec3sc np = p;
						np[n]--;

						label_t lbl = labels(np);
						if (lbl == 0)
						{
							weight_t w = weights(np);
							if (w > 0)
							{
								points.push(internals::MeyerSeed<label_t, weight_t>(np, targetLabel, w, round));
							}
						}
					}

					if (p[n] < labels.dimension(n) - 1)
					{
						Vec3sc np = p;
						np[n]++;

						label_t lbl = labels(np);
						if (lbl == 0)
						{
							weight_t w = weights(np);
							if (w > 0)
							{
								points.push(internals::MeyerSeed<label_t, weight_t>(np, targetLabel, w, round));
							}
						}
					}
				}


			}


			// Progress report for large fills
			size_t s = points.size();
			if (s > 0 && s % 50000 == 0 && lastPrinted != s)
			{
				lastPrinted = s;
				std::cout << s << " seeds...\r" << std::flush;
			}
		}

		if (lastPrinted != 0)
			std::cout << std::endl;

		//cout << "Max queue depth = " << maxQueueDepth << std::endl;
	}

	



	namespace internals
	{
		/**
		Storage class for fill point priority queue in Meyer's algorithm.
		*/
		template<typename label_t> class MeyerSeedNoWeight
		{
		private:
			Vec3sc pos;

			label_t targetLabel;
			size_t birthday;

		public:

			/**
			Constructor
			@param p The point.
			@param label Label color of the point.
			@param birthday The filling round number. Used to prioritize older points before newer points so that points near seeds are filled first.
			*/
			MeyerSeedNoWeight(const Vec3sc& p, label_t label, size_t birthday) :
				pos(p),
				targetLabel(label),
				birthday(birthday)
			{
			}

			/**
			Gets position.
			*/
			const Vec3sc& position() const
			{
				return pos;
			}

			/**
			Compares weights.
			*/
			bool operator < (const MeyerSeedNoWeight& right) const
			{
				return birthday > right.birthday;
			}

			/**
			Gets label value.
			*/
			const label_t label() const
			{
				return targetLabel;
			}

		};
	}



	/**
	Region grow segmentation.
	Grows all colored regions towards specific color.
	@param labels Image containing the labels of distinct areas to be grown and allowed regions marked with allowedColor.
	@param allowedColor Labels will be grown only to pixels that have this color.
	@param backgroundColor No pixels having this color will be filled. Set to allowedColor to fill to all pixels.
	*/
	template<typename label_t> size_t growAll(Image<label_t>& labels, label_t allowedColor, label_t backgroundColor, Connectivity connectivity = Connectivity::NearestNeighbours)
	{
		// Find all colors
		std::set<label_t> values;
		unique(labels, values);

		size_t totalChanged = 0;

		// Grow each of them
		ProgressIndicator progress(values.size());
		for (label_t srcColor : values)
		{
			if(srcColor != backgroundColor && srcColor != allowedColor)
				totalChanged += grow(labels, srcColor, allowedColor, connectivity);

			progress.step();
		}

		return totalChanged;

		/*
		std::priority_queue<internals::MeyerSeedNoWeight<label_t> > points;

		// Add all seed points to the priority queue

		for (coord_t z = 0; z < labels.depth(); z++)
		{
			for (coord_t y = 0; y < labels.height(); y++)
			{
				for (coord_t x = 0; x < labels.width(); x++)
				{
					label_t& p = labels(x, y, z);
					if (p != allowedColor && p != backgroundColor)
					{
						points.push(internals::MeyerSeedNoWeight<label_t>(Vec3sc((int32_t)x, (int32_t)y, (int32_t)z), p, 0));
						p = allowedColor;
					}
				}
			}
		}

		size_t lastPrinted = 0;
		long round = 0;
		size_t changed = 0;
		
		size_t seedCount = points.size();

		// Grow from the point p to all directions if they are not filled yet.
		while (!points.empty())
		{
			round++;

			const internals::MeyerSeedNoWeight<label_t>& obj = points.top();
			Vec3sc p = obj.position();
			label_t targetLabel = obj.label();

			points.pop();

			// Only proceed if the label of this pixel has not been set yet.
			if (labels(p) == allowedColor)
			{
				labels(p) = targetLabel;
				changed++;

				// Insert neighbours into the priority queue.
				if (connectivity == Connectivity::NearestNeighbours)
				{
					for (size_t n = 0; n < 3; n++)
					{
						if (p[n] > 0)
						{
							Vec3sc np(p);
							np[n]--;

							label_t lbl = labels(np);
							if (lbl == allowedColor)
							{
								points.push(internals::MeyerSeedNoWeight<label_t>(np, targetLabel, round));
							}
						}

						if (p[n] < labels.dimension(n) - 1)
						{
							Vec3sc np(p);
							np[n]++;

							label_t lbl = labels(np);
							if (lbl == allowedColor)
							{
								points.push(internals::MeyerSeedNoWeight<label_t>(np, targetLabel, round));
							}
						}
					}
				}
				else
				{
					for (int32_t dz = -1; dz <= 1; dz++)
					{
						for (int32_t dy = -1; dy <= 1; dy++)
						{
							for (int32_t dx = -1; dx <= 1; dx++)
							{
								Vec3sc np(p.x + dx, p.y + dy, p.z + dz);

								
								if (labels.isInImage(np) && labels(np) == allowedColor)
								{
									points.push(internals::MeyerSeedNoWeight<label_t>(np, targetLabel, round));
								}
							}
						}
					}
				}


			}


			// Progress report for large fills
			size_t s = points.size();
			if (s > 0 && s % 50000 == 0 && lastPrinted != s)
			{
				lastPrinted = s;
				std::cout << s << " seeds...\r" << std::flush;
				//std::cout << s << " seeds = " << bytesToString((double)(s * sizeof(internals::MeyerSeedNoWeight<label_t>))) << "...\r" << std::flush;
			}
		}

		if (lastPrinted != 0)
			std::cout << std::endl;

		return changed - seedCount;
		*/
	}

	namespace tests
	{
		void floodfillSanityChecks();
		void floodfill();
		void floodfillLeaks();
		void floodfillThreading();
		void growPriority();
		void growAll();
		void growComparison();
		void growCustomHeap();
	}

}
