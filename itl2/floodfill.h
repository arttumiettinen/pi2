#pragma once

#include <set>
#include <vector>
#include <queue>

#include "image.h"
#include "math/vec3.h"

using math::Vec3c;

namespace itl2
{
	/*
	* Defines possible connectivity modes (etc. for flood fill).
	*/
	enum Connectivity
	{
		/**
		Declares connectivity where only nearest (in 2D 4, in 3D 6) neighbours are connected.
		*/
		NearestNeighbours,

		/**
		Declares connectivity where all neighbours are connected (in 2D 8, in 3D 27).
		*/
		AllNeighbours
	};

	/**
	Perform flood fill.
	@param image Image containing the geometry to be filled.
	@param start Starting position.
	@param fillColor Fill color.
	@param stopColor Set to value different from fillColor to stop filling when a pixel of this color is encountered.
	@param connectivity Connectivity of the fill. (NearestNeighbours of AllNeighbours).
	@param pFilledPoints Pointer to vector that will receive the coordinates of the filled points. Set to zero if this information is not required.
	@param fillLimit Set to value to limit count of filled points to that value. Used in regionremoval code. Set to numeric_limits<size_t>::max() to allow (practically) any number of particles.
	@param pNeighbouringColors Pointer to set that will contain the colors of non-filled points neighbouring the filled region.
	@return true if the fill was terminated naturally; false if the fill was terminated by reaching fillLimit in filled pixel count or by encountering pixel with stopColor value.
	*/
	template<typename pixel_t> bool floodfill(Image<pixel_t>& image, const math::Vec3c& start, pixel_t fillColor, pixel_t stopColor, Connectivity connectivity = NearestNeighbours, std::vector<math::Vec3c>* pFilledPoints = 0, size_t fillLimit = numeric_limits<size_t>::max(), std::set<pixel_t>* pNeighbouringColors = 0)
	{
		pixel_t origColor = image(start);
		if (origColor == fillColor)
			return true;

		std::queue<Vec3c> points;
		points.push(start);

		size_t filledPoints = 0;
		while (!points.empty())
		{
			const math::Vec3c p = points.front();

			pixel_t pixelValue = image(p);
			if (pixelValue == origColor)
			{
				filledPoints++;
				image(p) = fillColor;
				if (pFilledPoints)
					pFilledPoints->push_back(p);

				// Add items to queue
				if (connectivity == NearestNeighbours)
				{
					for (size_t n = 0; n < 3; n++)
					{
						if (p[n] > 0)
						{
							Vec3c np = p;
							np[n]--;
							points.push(np);
						}

						if (p[n] < image.dimension(n) - 1)
						{
							Vec3c np = p;
							np[n]++;
							points.push(np);
						}
					}
				}
				else
				{
					// All neighbours
					for (coord_t dx = -1; dx <= 1; dx++)
					{
						for (coord_t dy = -1; dy <= 1; dy++)
						{
							for (coord_t dz = -1; dz <= 1; dz++)
							{
								Vec3c np(p.x + dx, p.y + dy, p.z + dz);
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
				return false;
			}
			else
			{
				if (pNeighbouringColors != 0)
					pNeighbouringColors->insert(pixelValue);
			}

			points.pop();

			if (filledPoints > fillLimit)
				return false; // Fill volume limit has been encountered.
		}

		return true;
	}

	namespace tests
	{
		void floodfill();
	}

}
