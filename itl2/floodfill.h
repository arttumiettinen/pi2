#pragma once

#include <set>
#include <vector>
#include <queue>
#include <tuple>
#include <iostream>

#include "image.h"
#include "math/vec3.h"
#include "connectivity.h"

using math::Vec3c;
using math::Vec3sc;

namespace itl2
{

	/**
	Perform flood fill.
	@param image Image containing the geometry to be filled.
	@param start Starting position.
	@param fillColor Fill color. The filled pixels will be colored with this color.
	@param stopColor Set to value different from fillColor to stop filling when a pixel of this color is encountered.
	@param connectivity Connectivity of the fill.
	@param pFilledPoints Pointer to vector that will receive the coordinates of the filled points. Set to zero if this information is not required.
	@param fillLimit Set to value to limit count of filled points to that value. Used in regionremoval code. Set to numeric_limits<size_t>::max() to allow (practically) any number of particles.
	@param pNeighbouringColors Pointer to set that will contain the colors of non-filled points neighbouring the filled region.
	@return True if the fill was terminated naturally; false if the fill was terminated by reaching fillLimit in filled pixel count; by encountering pixel with stopColor value; or if the origColor is fillColor.
	*/
	template<typename pixel_t> bool floodfill(Image<pixel_t>& image, const math::Vec3c& start, pixel_t fillColor, pixel_t stopColor, Connectivity connectivity = Connectivity::NearestNeighbours, size_t* pFilledPointCount = 0, std::vector<math::Vec3sc>* pFilledPoints = 0, size_t fillLimit = 0, std::set<pixel_t>* pNeighbouringColors = 0)
	{
		if (!image.isInImage(start))
			return true;

		pixel_t origColor = image(start);
		std::vector<Vec3sc> seeds;
		seeds.push_back(Vec3sc(start));

		return floodfill(image, seeds, origColor, fillColor, stopColor, connectivity, pFilledPointCount, pFilledPoints, fillLimit, pNeighbouringColors);
	}

	namespace internals
	{
		template<typename pixel_t> bool processNeighbours(coord_t x, coord_t y, coord_t z, queue<Vec3sc>& points, std::vector<tuple<coord_t, coord_t, bool> >& nbs, const Image<pixel_t>& image, pixel_t fillColor, pixel_t origColor, pixel_t stopColor, std::set<pixel_t>* pNeighbouringColors)
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
	@return True if the fill was terminated naturally; false if the fill was terminated by reaching fillLimit in filled pixel count; by encountering pixel with stopColor value; or if the origColor is fillColor.
	*/
	template<typename pixel_t> bool floodfill(Image<pixel_t>& image, const std::vector<Vec3sc>& seeds, pixel_t origColor, pixel_t fillColor, pixel_t stopColor, Connectivity connectivity = Connectivity::NearestNeighbours, size_t* pFilledPointCount = 0, std::vector<math::Vec3sc>* pFilledPoints = 0, size_t fillLimit = 0, std::set<pixel_t>* pNeighbouringColors = 0)
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
			fillLimit = numeric_limits<size_t>::max();

		// Contains {deltay, deltaz, active} for all neighbouring scanlines.
		std::vector<tuple<coord_t, coord_t, bool> > nbs;
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
		for (const Vec3sc& p : seeds)
		{
			if(image.isInImage(p))
				points.push(p);
		}

		size_t lastPrinted = 0;
		size_t tmp = 0;
		size_t* pCount = &tmp;
		if (pFilledPointCount)
			pCount = pFilledPointCount;

		while (!points.empty())
		{
			const math::Vec3c p = Vec3c(points.front());
			
			// Check that this point has not been filled before (there might be multiple routes to the same location).
			if (image(p) == origColor) // Can we get rid of this check?
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
					get<2>(nb) = false;

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
			if (s > 0 && s % 50000 == 0 && lastPrinted != s)
			{
				lastPrinted = s;
				std::cout << s << " seeds...\r" << std::flush;
			}
		}

		return true;
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

	namespace tests
	{
		void floodfillSanityChecks();
		void floodfill();
		void floodfillLeaks();
	}

}
