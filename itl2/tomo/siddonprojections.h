#pragma once

#include "datatypes.h"
#include "math/vec3.h"
#include "math/vec2.h"
#include "image.h"
#include "io/raw.h"
#include "pointprocess.h"
#include "raytrace.h"
#include "generation.h"

namespace itl2
{
	/**
	Rotates vector v given source rotation and tilt angles.
	*/
	template<typename real_t> void rotate(Vec3<real_t>& v, real_t rotationAngle, real_t tiltAngle)
	{
		// Rotate
		v = v.rotate(Vec3<real_t>(0, 0, 1), rotationAngle);

		// Tilt
		Vec3<real_t> rotatedTiltAxis(0, 1, 0);
		rotatedTiltAxis = rotatedTiltAxis.rotate(Vec3<real_t>(0, 0, 1), rotationAngle);
		v = v.rotate(rotatedTiltAxis, tiltAngle);
	}

	/**
	Calculates projections or backprojections as determined by operation argument.
	Calls operation(projectionx, projectiony, projectionindex, startpoint, endpoint) for all lines
	from source to center of every pixel on every projection.
	*/
	template<typename pixel_t, typename real_t, typename op>
	void project(const Vec3<coord_t>& imgDimensions,
		const Vec3<coord_t>& projectionDimensions,
		const real_t sourceDistance,
		const Vec2<real_t>& angularRange,
		const Vec2<real_t>& tiltRange,
		op& operation)
	{

		coord_t projectionWidth = projectionDimensions.x;
		coord_t projectionHeight = projectionDimensions.y;
		coord_t projectionCount = projectionDimensions.z;

		for (coord_t i = 0; i < projectionCount; i++)
		{
			real_t rotationAngle = projectionCount > 1 ? angularRange.x + (angularRange.y - angularRange.x) / (projectionCount - 1) * i : angularRange.x;
			real_t tiltAngle = projectionCount > 1 ? tiltRange.x + (tiltRange.y - tiltRange.x) / (projectionCount - 1) * i : tiltRange.x;

			Vec3<real_t> srcPos(-sourceDistance, 0, 0);
			Vec3<real_t> cameraX(0, 1, 0);
			Vec3<real_t> cameraY(0, 0, 1);

			rotate(srcPos, rotationAngle, tiltAngle);
			rotate(cameraX, rotationAngle, tiltAngle);
			rotate(cameraY, rotationAngle, tiltAngle);

			Vec3<real_t> center = ((Vec3<real_t>)imgDimensions + Vec3<real_t>(1, 1, 1)) / 2;
			srcPos += center;

			for (coord_t y = 0; y < projectionHeight; y++)
			{
				for (coord_t x = 0; x < projectionWidth; x++)
				{
					// Camera is placed at the center of the image perpendicular to non-tilted optical axis.
					// Pixel size of projection and original image are equal.
					Vec3<real_t> destPixelPos = center + (x - projectionWidth / (real_t)2.0) * cameraX + (y - projectionHeight / (real_t)2.0) * cameraY;
					Vec3<real_t> start = srcPos;
					Vec3<real_t> end = srcPos + 5 * (destPixelPos - srcPos);

					operation(x, y, i, start, end);
				}
			}
		}

	}

	/**
	Siddon ray-tracing based forward projection operation for project method.
	*/
	template<typename pixel_t, typename real_t> struct ForwardProjectionSiddon
	{
	private:
		LineProjector<pixel_t, real_t> projector;
		Image<real_t>& projections;
	public:

		ForwardProjectionSiddon(const Image<pixel_t>& img, Image<real_t>& projections) :
			projector(img),
			projections(projections)
		{

		}

		void operator()(coord_t x, coord_t y, coord_t i, Vec3<real_t>& start, Vec3<real_t>& end)
		{
			projector.reset();
			siddonLineClip<real_t, LineProjector<pixel_t, real_t> >(start, end, projector, (Vec3<real_t>)projector.img.dimensions());
			projections(x, y, i) = projector.getValue();
		}
	};

	/**
	Siddon ray-tracing based backprojection operator for project method.
	*/
	template<typename pixel_t, typename real_t> struct BackProjectionSiddon
	{
	private:
		LineBackProjector<pixel_t, real_t> projector;
		const Image<real_t>& projections;
	public:

		BackProjectionSiddon(Image<pixel_t>& img, const Image<real_t>& projections) :
			projector(img, 0),
			projections(projections)
		{

		}

		void operator()(coord_t x, coord_t y, coord_t i, Vec3<real_t>& start, Vec3<real_t>& end)
		{
			projector.setProjectionValue(projections(x, y, i));
			siddonLineClip<real_t, LineBackProjector<pixel_t, real_t> >(start, end, projector, (Vec3<real_t>)projector.img.dimensions());
		}
	};

	/**
	Creates forward projections of image img to image projections with Siddon algorithm.
	*/
	template<typename pixel_t, typename real_t> void createForwardProjectionSiddon(
		const Image<pixel_t>& img,
		Image<real_t>& projections,
		const real_t sourceDistance,
		const Vec2<real_t>& angularRange,
		const Vec2<real_t>& tiltRange)
	{
		ForwardProjectionSiddon<pixel_t, real_t> forwardProjection(img, projections);
		project<pixel_t, real_t, ForwardProjectionSiddon<pixel_t, real_t> >(
			img.dimensions(),
			projections.dimensions(),
			sourceDistance,
			angularRange,
			tiltRange,
			forwardProjection);
	}

	/**
	Creates backprojection of projections to img with Siddon algorithm.
	Produces aliasing artefacts.
	*/
	template<typename real_t> void createBackProjectionSiddon(
		Image<real_t>& img,
		const Image<real_t>& projections,
		const real_t sourceDistance,
		const Vec2<real_t>& angularRange,
		const Vec2<real_t>& tiltRange)
	{
		BackProjectionSiddon<real_t, real_t> backProjection(img, projections);
		project<real_t, real_t, BackProjectionSiddon<real_t, real_t> >(
			img.dimensions(),
			projections.dimensions(),
			sourceDistance,
			angularRange,
			tiltRange,
			backProjection);
	}

	namespace tests
	{
		inline void projectionConsistency()
		{
			string filename = "input_data/plates";
			Image<float32_t> img(100, 100, 100);
			raw::read(img, filename);

			Image<float32_t> projections(100, 100, 1);
			createForwardProjectionSiddon<float32_t, float32_t>(
				img, projections,
				300,
				Vec2f(0, 0),
				Vec2f(0, 0));

			raw::writed(projections, "siddon_consistency/projections");

			// Projection of backprojection should equal projection
			Image<float32_t> backprojection(img.dimensions());
			createBackProjectionSiddon<float32_t>(
				backprojection, projections,
				300,
				Vec2f(0, 0),
				Vec2f(0, 0));

			Image<float32_t> projectionsb(100, 100, 1);
			createForwardProjectionSiddon<float32_t, float32_t>(
				backprojection, projectionsb,
				300,
				Vec2f(0, 0),
				Vec2f(0, 0));

			raw::writed(projectionsb, "siddon_consistency/projections_of_backprojection");

			subtract(projections, projectionsb);
			raw::writed(projections, "siddon_consistency/difference");
		}

		inline void createPlates()
		{
			Image<float32_t> plates(100, 100, 100);
			int r = 5;
			int step = 3 * r;
			int rxy = 30;
			double alpha = 0;
			for (int z = r + 3; z < plates.depth() - r - 3; z += step)
			{
				double c = cos(alpha / 180.0 * 3.1415);
				double s = sin(alpha / 180.0 * 3.1415);
				draw(plates, Box(Vec3d(plates.width() / 2.0, plates.height() / 2.0, z), Vec3d(rxy, rxy, r), Vec3d(c, s, 0), Vec3d(-s, c, 0)), 1.0f);
				alpha += 45;
			}

			raw::writed(plates, "projections/plates");
		}

		inline void createProjection()
		{
			string filename = "projections/plates";
			Image<float32_t> img(100, 100, 100);
			raw::read(img, filename);

			Image<float32_t> projections(100, 100, 1);
			createForwardProjectionSiddon<float32_t, float32_t>(
				img, projections,
				300,
				Vec2f(0, 0),
				Vec2f(0, 0));

			raw::writed(projections, "projections/single_projection");
		}

		inline void createProjections()
		{
			string filename = "projections/plates";
			Vec2f angularRange(-180, 180);
			Vec2f tiltRange(0, 0);
			float32_t sourceDistance = 300;

			Image<float32_t> img(100, 100, 100);
			raw::read(img, filename);

			Image<float32_t> projections(100, 100, 180);
			createForwardProjectionSiddon<float32_t, float32_t>(
				img, projections,
				sourceDistance,
				angularRange / 180 * PIf,
				tiltRange / 180 * PIf);

			raw::writed(projections, "projections/projected");
		}

		inline void create10Projections()
		{
			string filename = "projections/plates";
			Vec2f angularRange(-180, 180);
			Vec2f tiltRange(0, 0);
			float32_t sourceDistance = 300;

			Image<float32_t> img(100, 100, 100);
			raw::read(img, filename);

			Image<float32_t> projections(100, 100, 10);
			createForwardProjectionSiddon<float32_t, float32_t>(
				img, projections,
				sourceDistance,
				angularRange / 180 * PIf,
				tiltRange / 180 * PIf);

			raw::writed(projections, "projections/projected10");
		}

		inline void createMoreProjections()
		{
			string filename = "projections/plates";
			Vec2f angularRange(-180, 180);
			Vec2f tiltRange(0, 0);
			float32_t sourceDistance = 300;

			Image<float32_t> img(100, 100, 100);
			raw::read(img, filename);

			Image<float32_t> projections(100, 100, 120);
			createForwardProjectionSiddon<float32_t, float32_t>(
				img, projections,
				sourceDistance,
				angularRange / 180 * PIf,
				tiltRange / 180 * PIf);

			raw::writed(projections, "projections/projected120");
		}

		inline void create36ProjectionsSheppLogan()
		{
			string filename = "input_data/shepp-logan_179x179x179.raw";
			Vec2f angularRange(-180, 180);
			Vec2f tiltRange(0, 0);
			float32_t sourceDistance = 300;

			Image<float32_t> img(179, 179, 179);
			raw::read(img, filename);

			Image<float32_t> projections(179, 179, 36);
			createForwardProjectionSiddon<float32_t, float32_t>(
				img, projections,
				sourceDistance,
				angularRange / 180 * PIf,
				tiltRange / 180 * PIf);

			raw::writed(projections, "shepp-logan/projected36sh");
		}

		inline void createBackprojection()
		{
			// This file is generated by create*Projections methods
			string filename = "projections/projected_100x100x180.raw";
			Vec2f angularRange(-180, 180);
			Vec2f tiltRange(0, 0);
			float32_t sourceDistance = 300;

			Image<float32_t> img(100, 100, 100);
			Image<float32_t> projections(100, 100, 180);

			// First create line length projections
			Image<float32_t> lineLength(img.dimensions());

			setValue(projections, 1);

			createBackProjectionSiddon<float32_t>(
				lineLength, projections,
				sourceDistance,
				angularRange / 180 * PIf,
				tiltRange / 180 * PIf);

			raw::writed(lineLength, "projections/line_length");

			// Now read projections and back-project
			raw::read(projections, filename);

			Image<float32_t> backprojection(img.dimensions());
			createBackProjectionSiddon<float32_t>(
				backprojection, projections,
				sourceDistance,
				angularRange / 180 * PIf,
				tiltRange / 180 * PIf);

			raw::writed(backprojection, "projections/backprojection");
		}
	}
}
