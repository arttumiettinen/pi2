#pragma once

#include "image.h"
#include "filters.h"
#include "projections.h"
#include "utilities.h"
#include "interpolation.h"

namespace itl2
{
	namespace internals
	{
		template<typename pixel_t> void visualizeLayer(const Image<pixel_t>& geometry, Image<float32_t>& heightMap, Image<pixel_t>* pVisualization, coord_t visualizeY, pixel_t color, coord_t iteration)
		{
			if (pVisualization != nullptr && visualizeY >= 0 && visualizeY < geometry.height())
			{
				// Fill the slice from the original
				for (coord_t z = 0; z < geometry.depth(); z++)
				{
					for (coord_t x = 0; x < geometry.width(); x++)
					{
						(*pVisualization)(x, z, iteration) = geometry(x, visualizeY, z);
					}
				}

				// Mark position of the surface
				for (coord_t x = 0; x < heightMap.width(); x++)
				{
					coord_t z = itl2::pixelRound<coord_t>(heightMap(x, visualizeY));
					Vec3c p(x, z, iteration);
					if(pVisualization->isInImage(p))
						(*pVisualization)(p) = color;
				}
			}
		}
	}


	/**
	Direction of carpet movement.
	*/
	enum class Direction
	{
		/**
		Towards larger z.
		*/
		Down,
		/**
		Towards smaller z.
		*/
		Up
	};

	template<>
	inline Direction fromString(const string& dt)
	{
		string dt2 = dt;
		trim(dt2);
		toLower(dt2);

		if (dt2 == "up")
			return Direction::Up;

		return Direction::Down;
	}

	inline string toString(Direction dt)
	{
		switch (dt)
		{
		case Direction::Up: return "Up";
		case Direction::Down: return "Down";
		default: return "Unknown";
		}
	}

	/**
	Surface recognition algorithm 'Carpet' according to Turpeinen - Interface Detection Using a Quenched-Noise Version of the Edwards–Wilkinson Equation. 
	The algorithm places a surface above (alternatively below) the image, and moves it towards larger (alternatively smaller) $z$ values according to the E-W equation. 
	The movement of the surface stops when it encounters enough pixels with value above specific stopping gray level. 
	The surface does not move through small holes in the object as it has controllable amount of surface tension.
	@param geometry Geometry image.
	@param heightMap At input, this image should be empty or contain initial height map of the surface. The size of the height map must be w x h x 1 where w and h are width and height of the geometry image.
	@param stoppingValue Stopping gray-value.
	@param direction Direction where the surface moves. Down corresponds to direction towards larger z values, and Up corresponds to direction towards smaller z values.
	@param surfaceTension Value that indicates how smooth the surface will be.
	@param iterations Count of iterations to perform.
	@param pVisualization Pointer to an image where a visualization of the evolution of the surface will be made. The dimensions of the visualization will be set to w x d x iterations, where w and d are width and depth of the geometry image, and iterations is the count of iterations to perform. Set to nullptr to skip creation of the visualization.
	@param visualizeY If pVisualization is not nullptr, this value indicates the y-coordinate of the xz-slice that will be visualized.
	@param visColor Color of the surface in the visualization. If set to zero, the color will be set to one above the maximum in geometry image.
	*/
	template<typename pixel_t> void findSurface(const Image<pixel_t>& geometry, Image<float32_t>& heightMap, double stoppingValue, Direction direction = Direction::Down, double surfaceTension = 1.0, size_t iterations = 150, Image<pixel_t>* pVisualization = nullptr, coord_t visualizeY = 0, pixel_t visColor = 0, float32_t maxMove = 1.0)
	{
		float32_t upForceFactor = (float32_t)(-1.0 / stoppingValue);
		float32_t downForce = 1.0;
		if (direction == Direction::Up)
		{
			downForce *= -1;
			upForceFactor *= -1;
		}

		heightMap.ensureSize(geometry.width(), geometry.height(), 1);

		if (pVisualization != nullptr && visualizeY >= 0 && visualizeY < geometry.height())
		{
			pVisualization->ensureSize(geometry.width(), geometry.depth(), iterations);

			if (visColor == 0)
				visColor = NumberUtils<pixel_t>::saturatingAdd(max(geometry), 1);
		}

		// Linear interpolation may make the carpet a bit more stable, Zero boundary condition reduces bulging out of image.
		LinearInterpolator<float32_t, pixel_t> interp(BoundaryCondition::Zero);

		ProgressIndicator progress(iterations);
		for (size_t n = 0; n < iterations; n++)
		{
			// Move
			for (coord_t y = 0; y < heightMap.height(); y++)
			{
				for (coord_t x = 0; x < heightMap.width(); x++)
				{
					float32_t imgValue = interp(geometry, (float32_t)x, (float32_t)y, heightMap(x, y));

					float32_t newDepth = heightMap(x, y) + downForce + upForceFactor * imgValue;

					// Constraint amount of movement in either direction.
					if (newDepth > heightMap(x, y) + maxMove)
						newDepth = heightMap(x, y) + maxMove;
					else if(newDepth < heightMap(x, y) - maxMove)
						newDepth = heightMap(x, y) - maxMove;

					// Ensure that the carpet does not bulge out of the image.
					if (newDepth < 0)
						newDepth = 0;
					else if (newDepth > geometry.depth() - 1)
						newDepth = (float32_t)geometry.depth() - 1;
					
					heightMap(x, y) = newDepth;
				}
			}

			// Apply surface tension
			gaussFilter(heightMap, surfaceTension, BoundaryCondition::Nearest);

			internals::visualizeLayer(geometry, heightMap, pVisualization, visualizeY, visColor, (coord_t)n);

			progress.step();
		}
	}

	/**
	Draw height map to given image.
	@param image Image to draw to. Image width and height must match those of the height map.
	@param heightMap Height map to draw.
	@param color Color that is used to draw the height map.
	*/
	template<typename pixel_t> void drawHeightMap(Image<pixel_t>& image, const Image<float32_t>& heightMap, pixel_t color)
	{
		if (image.width() != heightMap.width() || image.height() != heightMap.height())
			throw ITLException("Image width and height must match those of the height map.");

		for (coord_t y = 0; y < heightMap.height(); y++)
		{
			for (coord_t x = 0; x < heightMap.width(); x++)
			{
				coord_t z = itl2::pixelRound<coord_t>(heightMap(x, y));
				Vec3c p(x, y, z);
				if (image.isInImage(p))
					image(p) = color;
			}
		}
	}

	/**
	Finds pixels of given image located before the height map and sets them to the given value.
	@param image Image whose pixels will be set. The width and the height of the image must be the same than those of the height map.
	@param heightMap Height map.
	@param color Color to which pixels are set.
	*/
	template<typename pixel_t> void setBeforeHeightMap(Image<pixel_t>& image, const Image<float32_t>& heightMap, pixel_t color)
	{
		if (image.width() != heightMap.width() || image.height() != heightMap.height())
			throw ITLException("Image width and height must match those of the height map.");

		coord_t maxZ = itl2::ceil(max(heightMap));
		clamp(maxZ, (coord_t)0, image.depth());

		for (coord_t z = 0; z < maxZ; z++)
		{
			for (coord_t y = 0; y < heightMap.height(); y++)
			{
				for (coord_t x = 0; x < heightMap.width(); x++)
				{
					float32_t zh = heightMap(x, y);
					if (z < zh)
						image(x, y, z) = color;
				}
			}
		}
	}


	/**
	Finds pixels of given image located after the height map and sets them to the given value.
	@param image Image whose pixels will be set. The width and the height of the image must be the same than those of the height map.
	@param heightMap Height map.
	@param color Color to which pixels are set.
	*/
	template<typename pixel_t> void setAfterHeightMap(Image<pixel_t>& image, const Image<float32_t>& heightMap, pixel_t color)
	{
		if (image.width() != heightMap.width() || image.height() != heightMap.height())
			throw ITLException("Image width and height must match those of the height map.");

		coord_t minZ = itl2::floor(min(heightMap));
		clamp(minZ, (coord_t)0, image.depth());

		for (coord_t z = minZ; z < image.depth(); z++)
		{
			for (coord_t y = 0; y < heightMap.height(); y++)
			{
				for (coord_t x = 0; x < heightMap.width(); x++)
				{
					float32_t zh = heightMap(x, y);
					if (z > zh)
						image(x, y, z) = color;
				}
			}
		}
	}


	/**
	Shift each z-directional column of input image by amount given in the shift map.
	@param img Image to shift.
	@param shiftMap Image that contains the amount of shift for each (x, y)-position of the input image. The width and height of this image and the input image must be equal.
	@param subtractMean Set to true to automatically subtract the average value of the shiftMap from each shift. This is useful if the shift map is negation of a surface map found using the findsurface command, and the intent is to make the surface straight.
	@param interpolate Interpolator
	*/
	template<typename pixel_t> void shiftZ(Image<pixel_t>& img, const Image<float32_t>& shiftMap, bool subtractMean = false, const Interpolator<pixel_t, pixel_t>& interpolate = LinearInterpolator<pixel_t, pixel_t>(BoundaryCondition::Zero))
	{
		if(img.width() != shiftMap.width() || img.height() != shiftMap.height())
			throw ITLException("Shift map width and height do not correspond to the geometry image width and height.");

		float32_t extraShift = 0;
		if (subtractMean)
		{
			extraShift = -(float32_t)mean(shiftMap);
		}

		#pragma omp parallel if(!omp_in_parallel())
		{
			Image<pixel_t> buffer(img.depth());

			#pragma omp for
			for (coord_t y = 0; y < img.height(); y++)
			{
				for (coord_t x = 0; x < img.width(); x++)
				{
					float32_t shift = shiftMap(x, y) + extraShift;

					// Get shifted values into a buffer
					for (coord_t z = 0; z < img.depth(); z++)
					{
						buffer(z) = interpolate(img, (float32_t)x, (float32_t)y, (float32_t)z - shift);
					}

					// Write back to image
					for (coord_t z = 0; z < img.depth(); z++)
					{
						img(x, y, z) = buffer(z);
					}
				}
			}
		}
	}


	namespace tests
	{
		void carpet();
	}

}
