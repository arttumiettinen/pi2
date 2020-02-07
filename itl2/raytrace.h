#pragma once

#include "image.h"
#include "math/vec3.h"
#include "math/numberutils.h"
#include "io/raw.h"
#include "math/geometry.h"

namespace itl2
{

	/**
	Calculate ray trace through image using Siddon method.
	Assumes pixel at location (x, y, z) spans region [x, x+1[ x [y, y+1[ x [z, z+1].
	Calls processing function/functor 'operation' for each pixel on the path.
	Call is like operation(const Vec3c& pixel_coordinates, const real_t length_in_that_pixel, const real_t total_line_length)
	*/
	template<typename real_t, typename op> void siddonLine(const Vec3<real_t>& start, const Vec3<real_t>& end, op& operation)
	{
		Vec3<real_t> ray = end - start;
		real_t l = ray.norm();
		Vec3<real_t> dir = ray / l;

		// Divide length of ray to steps in each direction
		Vec3<real_t> deltat(l / fabs(ray.x), l / fabs(ray.y), l / fabs(ray.z));

		// Calculate first intersection with coordinate plane in each direction
		Vec3<real_t> nextt(
			(ray.x >= 0 ? (std::ceil(start.x) - start.x) : (start.x - std::floor(start.x))) * deltat.x,
			(ray.y >= 0 ? (std::ceil(start.y) - start.y) : (start.y - std::floor(start.y))) * deltat.y,
			(ray.z >= 0 ? (std::ceil(start.z) - start.z) : (start.z - std::floor(start.z))) * deltat.z);

		// If there's no movement in some direction, make sure that the next t in that direction is infinity.
		if (std::isinf(deltat.x))
			nextt.x = std::numeric_limits<real_t>::infinity();
		if (std::isinf(deltat.y))
			nextt.y = std::numeric_limits<real_t>::infinity();
		if (std::isinf(deltat.z))
			nextt.z = std::numeric_limits<real_t>::infinity();

		real_t t = 0;
		real_t prevt = 0;
		while (NumberUtils<real_t>::lessThan(t, l))
		{
			// The next intersection point with pixel edge
			// is the one that is nearest to the previous one along the ray.
			if (nextt.x <= nextt.y && nextt.x <= nextt.z)
			{
				t = nextt.x;
				nextt.x += deltat.x;
			}
			else if (nextt.y <= nextt.x && nextt.y <= nextt.z)
			{
				t = nextt.y;
				nextt.y += deltat.y;
			}
			else
			{
				t = nextt.z;
				nextt.z += deltat.z;
			}

			// Make sure that we don't pass the end of the line
			if (NumberUtils<real_t>::greaterThanOrEqual(t, l))
				t = l;

			// Length of intersection between line and the current pixel
			real_t inpixel = t - prevt;

			// Calculate location of the current pixel
			Vec3<real_t> tmp = start + (t + prevt) / 2 * dir;
			Vec3c pixelpos((coord_t)tmp.x, (coord_t)tmp.y, (coord_t)tmp.z);

			// Process current pixel if there's really nonzero length of the line in the current pixel.
			// This takes care of not processing same pixel multiple times.
			if (NumberUtils<real_t>::greaterThan(inpixel, 0, (real_t)1e-3))
				operation(pixelpos, inpixel, l);

			prevt = t;
		}
	}

	/**
	Similar to siddonLine method above but clips the line to the [0, 0, 0] x dimensions.
	The clipped line is returned in start and end vectors. If the line does not intersect the region, start and end will be undefined and
	the line is not drawn at all.
	*/
	template<typename real_t, typename op> void siddonLineClip(Vec3<real_t>& start, Vec3<real_t>& end, op& operation, const Vec3<real_t>& dimensions)
	{
		if (clipLine(start, end, Vec3<real_t>(0, 0, 0), Vec3<real_t>(dimensions)))// - Vec3<real_t>(1, 1, 1)))
			siddonLine(start, end, operation);
	}

	/**
	Plots constant color value to image
	*/
	template<typename pixel_t, typename real_t> struct LineColorPlotter
	{
	private:
		Image<pixel_t> & img;
		pixel_t color;

	public:
		LineColorPlotter(Image<pixel_t>& img, pixel_t color) : img(img), color(color)
		{

		}

		void operator()(const Vec3c& pos, const real_t length, const real_t totalLength)
		{
			if (img.isInImage(pos))
				img(pos) = color;
		}
	};

	/**
	Plots in-pixel line lengths to given image.
	*/
	template<typename pixel_t, typename real_t> struct LinePlotter
	{
	private:
		Image<pixel_t> & img;

	public:
		LinePlotter(Image<pixel_t>& img) : img(img)
		{

		}

		void operator()(const Vec3c& pos, const real_t length, const real_t totalLength)
		{
			if (img.isInImage(pos))
				img(pos) = pixelRound<pixel_t>(length);
		}
	};

	template<typename pixel_t, typename real_t> struct LineProjector
	{
	private:
		real_t projectionValue;
	public:
		const Image<pixel_t> & img;

		LineProjector(const Image<pixel_t>& img) : img(img), projectionValue(0)
		{

		}

		real_t getValue() const
		{
			return projectionValue;
		}

		void reset()
		{
			projectionValue = 0;
		}

		void operator()(const Vec3c& pos, const real_t length, const real_t totalLength)
		{
			if (img.isInImage(pos))
				projectionValue += img(pos) * length;
		}
	};

	template<typename pixel_t, typename real_t> struct LineBackProjector
	{
	private:
		real_t projectionValue;
	public:
		Image<pixel_t> & img;

		LineBackProjector(Image<pixel_t>& img, const real_t projectionValue) : img(img), projectionValue(projectionValue)
		{

		}

		void setProjectionValue(real_t val)
		{
			projectionValue = val;
		}

		void operator()(const Vec3c& pos, const real_t length, const real_t totalLength)
		{
			if (img.isInImage(pos))
				img(pos) += projectionValue * length / totalLength;
		}
	};

	namespace tests
	{

		inline void siddonProject()
		{

			Image<float32_t> img1(30, 30);
			LinePlotter<float32_t, double> plotter1(img1);
			siddonLine<double, LinePlotter<float32_t, double> >(Vec3d(10, 10, 0), Vec3d(20, 19, 0), plotter1);

			siddonLine<double, LinePlotter<float32_t, double> >(Vec3d(25, 3, 0), Vec3d(25, 15, 0), plotter1);

			siddonLine<double, LinePlotter<float32_t, double> >(Vec3d(10.1, 15, 0), Vec3d(15, 25.2, 0), plotter1);

			Image<float32_t> img2(30, 30);
			LinePlotter<float32_t, double> plotter2(img2);
			siddonLine<double, LinePlotter<float32_t, double> >(Vec3d(20, 19, 0), Vec3d(10, 10, 0), plotter2);

			// Compare img1 and img2, they should be the same.
			// Sum of image should be the length of the line.

			raw::writed(img1, "./siddon/siddon_line_1");
			raw::writed(img2, "./siddon/siddon_line_2");

			Image<float32_t> img3(30, 30);
			LinePlotter<float32_t, double> plotter3(img3);
			siddonLine<double, LinePlotter<float32_t, double> >(Vec3d(10.5, 10, 0), Vec3d(20.5, 10, 0), plotter3);
			Vec3d start(-10.5, -10, 0);
			Vec3d end(20.5, 35, 0);
			siddonLineClip<double, LinePlotter<float32_t, double> >(start, end, plotter3, Vec3d((double)img3.width(), (double)img3.height(), 0));
			raw::writed(img3, "./siddon/siddon_line_3");

			Image<float32_t> img4(100, 100, 100);
			LinePlotter<float32_t, double> plotter4(img4);
			siddonLine<double, LinePlotter<float32_t, double> >(Vec3d(10.2, 15.5, 5.7), Vec3d(80.3, 75, 90.27), plotter4);
			raw::writed(img4, "./siddon/siddon_line_4");
		}
	}
}
