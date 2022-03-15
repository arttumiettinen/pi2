#pragma once

#include "image.h"
#include "math/aabox.h"
#include "math/box.h"
#include "sphere.h"
#include "math/ellipsoid.h"
#include "math/line.h"
#include "math/capsule.h"
#include "raytrace.h"
#include "math/vec3.h"
#include "network.h"
#include "pointprocess.h"
#include "projections.h"

#include <iostream>
#include <vector>

namespace itl2
{

	/**
	Fills given image with ramp in given dimension.
	*/
	template<typename pixel_t> void ramp(Image<pixel_t>& img, size_t dimension, coord_t shift = 0)
	{
		if (dimension > 2)
			throw ITLException("Invalid ramp dimension.");

		#pragma omp parallel for if(!omp_in_parallel())
		for (coord_t z = 0; z < img.depth(); z++)
		{
			for (coord_t y = 0; y < img.height(); y++)
			{
				for (coord_t x = 0; x < img.width(); x++)
				{
					Vec3c pos(x, y, z);
					img(x, y, z) = pixelRound<pixel_t>(pos[dimension] + shift);
				}
			}
		}
	}

	/**
	Fills given image with ramp in all dimensions
	*/
	template<typename pixel_t> void ramp3(Image<pixel_t>& img, coord_t shift = 0)
	{
		#pragma omp parallel for if(!omp_in_parallel())
		for (coord_t z = 0; z < img.depth(); z++)
		{
			for (coord_t y = 0; y < img.height(); y++)
			{
				for (coord_t x = 0; x < img.width(); x++)
				{
					img(x, y, z) = pixelRound<pixel_t>(x + y + z + shift);
				}
			}
		}
	}


	/**
	Sets pixel values in an image to a constant.
	@param img Image.
	@param locations List of locations to set.
	@param value Drawing color.
	*/
	template<typename pixel_t> void draw(Image<pixel_t>& img, const std::vector<Vec3sc>& locations, pixel_t value)
	{
		for (const auto& p : locations)
			img(p) = value;
	}

	/**
	Draws each region in the given list of regions with different color.
	Color of the first region is 1, the second region 2, etc.
	Throws exception if there are more regions than colors.
	*/
	template<typename pixel_t> void draw(Image<pixel_t>& img, const std::vector<std::vector<Vec3sc> >& regions)
	{
		if (regions.size() > (size_t)std::numeric_limits<pixel_t>::max() + 1)
			throw ITLException("Image data type cannot hold enough values to store all the regions with different color.");

		for (size_t n = 0; n < regions.size(); n++)
		{
			draw(img, regions[n], pixelRound<pixel_t>(n + 1));
		}
	}

	/**
	Draws a box into the given image with given color.
	Clips the box to the image coordinates.
	@return Count of filled pixels.
	*/
	template<typename pixel_t, typename box_t> size_t draw(Image<pixel_t>& img, const AABox<box_t>& box, pixel_t value)
	{
		AABox<coord_t> rounded(box);
		AABox<coord_t> region = rounded.intersection(AABox(Vec3c(), img.dimensions()));

		#pragma omp parallel for if(!omp_in_parallel() && region.volume() > PARALLELIZATION_THRESHOLD)
		for (coord_t z = region.minc.z; z < region.maxc.z; z++)
		{
			for (coord_t y = region.minc.y; y < region.maxc.y; y++)
			{
				for (coord_t x = region.minc.x; x < region.maxc.x; x++)
				{
					img(x, y, z) = value;
				}
			}
		}

		return (size_t)region.volume();
	}

	/**
	Draws a sphere.
	@return Count of filled pixels.
	*/
	template<typename pixel_t, typename real_t> size_t draw(Image<pixel_t>& image, const Sphere<real_t>& sphere, pixel_t color)
	{
		Vec3c minPos = round(sphere.center - sphere.radius * Vec3<real_t>(1, 1, 1));
		Vec3c maxPos = round(sphere.center + sphere.radius * Vec3<real_t>(1, 1, 1));

		clamp(minPos, Vec3c(0, 0, 0), image.dimensions() - Vec3c(1, 1, 1));
		clamp(maxPos, Vec3c(0, 0, 0), image.dimensions() - Vec3c(1, 1, 1));

		size_t filledCount = 0;
		#pragma omp parallel for if(!omp_in_parallel() && AABox(minPos, maxPos).volume() > PARALLELIZATION_THRESHOLD) reduction(+:filledCount)
		for(coord_t z = minPos.z; z <= maxPos.z; z++)
		{
			for(coord_t y = minPos.y; y <= maxPos.y; y++)
			{
				for(coord_t x = minPos.x; x <= maxPos.x; x++)
				{
					if(sphere.contains(Vec3<real_t>((real_t)x, (real_t)y, (real_t)z)))
					{
						image(x, y, z) = color;
						filledCount++;
					}
				}
			}
		}

		return filledCount;
	}

	/**
	Draws a sphere based on center point and squared radius.
	@return Count of filled pixels.
	*/
	template<typename pixel_t> size_t draw(Image<pixel_t>& image, const Sphere2& sphere, pixel_t color)
	{
		const Vec3sc& center = sphere.pos;
		int32_t r2 = sphere.r2;

		int32_t ri = (int32_t)ceil(sqrt(r2));
		Vec3sc minPos = center - ri * Vec3sc(1, 1, 1);
		Vec3sc maxPos = center + ri * Vec3sc(1, 1, 1);

		clamp(minPos, Vec3sc(0, 0, 0), Vec3sc(image.dimensions()) - Vec3sc(1, 1, 1));
		clamp(maxPos, Vec3sc(0, 0, 0), Vec3sc(image.dimensions()) - Vec3sc(1, 1, 1));

		size_t filledCount = 0;
		#pragma omp parallel for if(!omp_in_parallel() && AABox(minPos, maxPos).volume() > PARALLELIZATION_THRESHOLD) reduction(+:filledCount)
		for (int32_t z = minPos.z; z <= maxPos.z; z++)
		{
			for (int32_t y = minPos.y; y <= maxPos.y; y++)
			{
				for (int32_t x = minPos.x; x <= maxPos.x; x++)
				{
					Vec3sc d = center - Vec3sc(x, y, z);

					if (d.x * d.x + d.y * d.y + d.z * d.z < r2)
					{
						image(x, y, z) = color;
						filledCount++;
					}
				}
			}
		}

		return filledCount;
	}

	/**
	Draws an ellipsoid.
	@return Count of filled pixels.
	*/
	template<typename pixel_t> size_t draw(Image<pixel_t>& image, const Ellipsoid& ellipsoid, pixel_t color)
	{
		AABox<double> bounds = ellipsoid.boundingBox();
		Vec3sc minPos(itl2::floor(bounds.minc));
		Vec3sc maxPos(itl2::ceil(bounds.maxc));

		minPos -= Vec3sc(1, 1, 1);
		maxPos += Vec3sc(1, 1, 1);

		clamp(minPos, Vec3sc(0, 0, 0), Vec3sc(image.dimensions()) - Vec3sc(1, 1, 1));
		clamp(maxPos, Vec3sc(0, 0, 0), Vec3sc(image.dimensions()) - Vec3sc(1, 1, 1));

		size_t filledCount = 0;
		#pragma omp parallel for if(!omp_in_parallel() && AABox(minPos, maxPos).volume() > PARALLELIZATION_THRESHOLD) reduction(+:filledCount)
		for (int32_t z = minPos.z; z <= maxPos.z; z++)
		{
			for (int32_t y = minPos.y; y <= maxPos.y; y++)
			{
				for (int32_t x = minPos.x; x <= maxPos.x; x++)
				{
					if(ellipsoid.contains(Vec3d(x, y, z)))
					{
						image(x, y, z) = color;
						filledCount++;
					}
				}
			}
		}

		return filledCount;
	}

	/**
	Draws a non-axis-aligned box.
	@return Count of filled pixels.
	*/
	template<typename pixel_t> size_t draw(Image<pixel_t>& image, const Box& box, pixel_t color)
	{
		AABox<double> bounds = box.boundingBox();
		Vec3sc minPos(itl2::floor(bounds.minc));
		Vec3sc maxPos(itl2::ceil(bounds.maxc));

		minPos -= Vec3sc(1, 1, 1);
		maxPos += Vec3sc(1, 1, 1);

		clamp(minPos, Vec3sc(0, 0, 0), Vec3sc(image.dimensions()) - Vec3sc(1, 1, 1));
		clamp(maxPos, Vec3sc(0, 0, 0), Vec3sc(image.dimensions()) - Vec3sc(1, 1, 1));

		size_t filledCount = 0;
#pragma omp parallel for if(!omp_in_parallel() && AABox(minPos, maxPos).volume() > PARALLELIZATION_THRESHOLD) reduction(+:filledCount)
		for (int32_t z = minPos.z; z <= maxPos.z; z++)
		{
			for (int32_t y = minPos.y; y <= maxPos.y; y++)
			{
				for (int32_t x = minPos.x; x <= maxPos.x; x++)
				{
					if (box.contains(Vec3d(x, y, z)))
					{
						image(x, y, z) = color;
						filledCount++;
					}
				}
			}
		}

		return filledCount;
	}

	/**
	Draws a line.
	*/
	template<typename pixel_t, typename real_t> void draw(Image<pixel_t>& image, const Line<real_t>& line, pixel_t color)
	{
		Vec3<real_t> start = line.start;
		Vec3<real_t> end = line.end;
		LineColorPlotter<pixel_t, real_t> plotter(image, color);
		siddonLineClip(start, end, plotter, Vec3<real_t>(image.dimensions()));
	}

	/**
	Draws a capsule.
	*/
	template<typename pixel_t, typename real_t> size_t draw(Image<pixel_t>& image, const Capsule<real_t>& capsule, pixel_t color)
	{
		Vec3c minPosStart = round(capsule.start - capsule.radius * Vec3<real_t>(1, 1, 1));
		Vec3c maxPosStart = round(capsule.start + capsule.radius * Vec3<real_t>(1, 1, 1));
		Vec3c minPosEnd = round(capsule.end - capsule.radius * Vec3<real_t>(1, 1, 1));
		Vec3c maxPosEnd = round(capsule.end + capsule.radius * Vec3<real_t>(1, 1, 1));

		Vec3c minPos = min(minPosStart, minPosEnd);
		Vec3c maxPos = max(maxPosStart, maxPosEnd);

		clamp(minPos, Vec3c(0, 0, 0), image.dimensions() - Vec3c(1, 1, 1));
		clamp(maxPos, Vec3c(0, 0, 0), image.dimensions() - Vec3c(1, 1, 1));

		size_t filledCount = 0;
		#pragma omp parallel for if(!omp_in_parallel() && AABox(minPos, maxPos).volume() > PARALLELIZATION_THRESHOLD) reduction(+:filledCount)
		for (coord_t z = minPos.z; z <= maxPos.z; z++)
		{
			for (coord_t y = minPos.y; y <= maxPos.y; y++)
			{
				for (coord_t x = minPos.x; x <= maxPos.x; x++)
				{
					if (capsule.contains(Vec3<real_t>((real_t)x, (real_t)y, (real_t)z)))
					{
						image(x, y, z) = color;
						filledCount++;
					}
				}
			}
		}

		return filledCount;
	}

	/**
	Draws a network to given image.
	Vertices are drawn as spheres and edges are drawn as lines connecting the spheres.
	@param useMeasuredEdgeRadius Set to true to draw edges as capsules with cross-sectional area read from network measurements. Edges with no nan area are drawn as lines. Set to false to draw the edges as lines.
	*/
	template<typename pixel_t> void draw(Image<pixel_t>& image, const Network& network, bool drawVerts = true, float32_t vertexRadius = 2, pixel_t vertexColor = std::numeric_limits<pixel_t>::max(), bool drawEdges = true, bool useMeasuredEdgeRadius = false, pixel_t edgeColor = std::numeric_limits<pixel_t>::max())
	{
		// Draw edges
		if (drawEdges)
		{
			for (const auto& edge : network.edges)
			{
				if (edge.verts[0] >= 0 && (size_t)edge.verts[0] < network.vertices.size() &&
					edge.verts[1] >= 0 && (size_t)edge.verts[1] < network.vertices.size())
				{
					float32_t r = std::sqrt(edge.properties.area / PIf);
					if (!useMeasuredEdgeRadius || std::isnan(r) || r < 0.75)
					{
						Line l(network.vertices[edge.verts[0]], network.vertices[edge.verts[1]]);
						draw(image, l, edgeColor);
					}
					else
					{
						Capsule c(network.vertices[edge.verts[0]], network.vertices[edge.verts[1]], r);
						draw(image, c, edgeColor);
					}
				}
			}
		}

		// Draw vertices
		if (drawVerts)
		{
			for (const auto& vertex : network.vertices)
			{
				Sphere<float32_t> s(vertex, vertexRadius);
				draw(image, s, vertexColor);
			}
		}
	}


	/**
	Draws boxes and spheres to the image.
	*/
	inline void generateSimpleGeometry(Image<uint8_t>& geom, unsigned int seed)
	{
		setValue(geom, 0);

		srand(seed);

		// Extract a few random numbers from the generator as otherwise sphere count is almost always less than box count.
		frand(0, 1000);
		frand(0, 1000);
		frand(0, 1000);

		coord_t sphereCount = itl2::round(frand(0, 1000) / (200.0 * 200.0 * 200.0) * geom.pixelCount());
		coord_t boxCount = itl2::round(frand(0, 1000) / (200.0 * 200.0 * 200.0) * geom.pixelCount());

		std::cout << "Generating " << sphereCount << " spheres..." << std::endl;
		for (coord_t n = 0; n < sphereCount; n++)
		{
			double r = frand(1, 20);
			double x = frand(-r, geom.width() + r);
			double y = frand(-r, geom.height() + r);
			double z = frand(-r, geom.depth() + r);

			draw(geom, Sphere(Vec3d(x, y, z), r), (uint8_t)1);
		}

		std::cout << "Generating " << boxCount << " boxes..." << std::endl;
		for (coord_t n = 0; n < boxCount; n++)
		{
			coord_t rx = randc(1, 20);
			coord_t ry = randc(1, 20);
			coord_t rz = randc(1, 20);
			coord_t x = randc(-rx, geom.width() + rx);
			coord_t y = randc(-ry, geom.height() + ry);
			coord_t z = randc(-rz, geom.depth() + rz);

			Vec3c c(x, y, z);
			Vec3c r(rx, ry, rz);
			draw(geom, AABox(c - r, c + r), (uint8_t)1);
		}

	}

	///**
	//Draws spheres and boxes to the image until volume fraction of filled pixels is at least targetVolumeFraction.
	//*/
	//inline void generateGeometry(Image<uint8_t>& geom, unsigned int seed, double targetVolumeFraction)
	//{

	//	setValue(geom, 0);

	//	srand(seed);

	//	// Extract a few random numbers from the generator to reduce bias.
	//	frand(0, 1000);
	//	frand(0, 1000);
	//	frand(0, 1000);

	//	size_t n = 0;
	//	double filledPixels;
	//	do
	//	{
	//		if (n % 2 == 0)
	//		{
	//			double r = frand(1, 50);
	//			double x = frand(-r, geom.width() + r);
	//			double y = frand(-r, geom.height() + r);
	//			double z = frand(-r, geom.depth() + r);

	//			draw(geom, Sphere(Vec3d(x, y, z), r), (uint8_t)1);
	//		}
	//		else
	//		{
	//			coord_t rx = randc(1, 30);
	//			coord_t ry = randc(1, 30);
	//			coord_t rz = randc(1, 30);
	//			coord_t x = randc(-rx, geom.width() + rx);
	//			coord_t y = randc(-ry, geom.height() + ry);
	//			coord_t z = randc(-rz, geom.depth() + rz);

	//			Vec3c c(x, y, z);
	//			Vec3c r(rx, ry, rz);
	//			draw(geom, AABox(c - r, c + r), (uint8_t)1);
	//		}

	//		n++;

	//		filledPixels = sum(geom);
	//	} while (filledPixels < geom.pixelCount() * targetVolumeFraction); // Generate until specified volume fraction
	//}

	/**
	Draws spheres and boxes to the image until volume fraction of filled pixels is at least targetVolumeFraction.
	*/
	inline void generateGeometry(Image<uint8_t>& geom, unsigned int seed, double targetVolumeFraction, bool indicateProgress = false)
	{

		setValue(geom, 0);

		srand(seed);

		// Extract a few random numbers from the generator to reduce bias.
		frand(0, 1000);
		frand(0, 1000);
		frand(0, 1000);

		size_t n = 0;
		double filledPixels;
		do
		{
			if (n % 2 == 0)
			{
				double r = frand(1, (double)geom.width() / 2);
				double x = frand(-r, geom.width() + r);
				double y = frand(-r, geom.height() + r);
				double z = frand(-r, geom.depth() + r);

				draw(geom, Sphere(Vec3d(x, y, z), r), (uint8_t)1);
			}
			else
			{
				coord_t rx = randc(1, geom.width() / 2);
				coord_t ry = randc(1, geom.width() / 2);
				coord_t rz = randc(1, geom.width() / 2);
				coord_t x = randc(-rx, geom.width() + rx);
				coord_t y = randc(-ry, geom.height() + ry);
				coord_t z = randc(-rz, geom.depth() + rz);

				Vec3c c(x, y, z);
				Vec3c r(rx, ry, rz);
				draw(geom, AABox(c - r, c + r), (uint8_t)1);
			}

			n++;

			filledPixels = sum(geom);

			std::cout << ((double)filledPixels / geom.pixelCount()) << " / " << (targetVolumeFraction) << "\r" << std::flush;
		} while (filledPixels < geom.pixelCount() * targetVolumeFraction); // Generate until specified volume fraction
	}
}
