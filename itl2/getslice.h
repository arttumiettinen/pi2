#pragma once

#include "image.h"
#include "math/matrix3x3.h"
#include "interpolation.h"

namespace itl2
{
	namespace internals
	{
		inline Matrix3x3d sliceRotationMatrix(Vec3d dir)
		{
			// Create orthogonal base by choosing two more directions
			dir = dir.normalized();
			Vec3d up(0, 0, 1);
			Vec3d right = dir.cross(up);
			if (right.norm() < 0.0001)
			{
				up = Vec3d(1, 0, 0);
				right = dir.cross(up);
			}
			right = right.normalized();
			up = right.cross(dir);

			// Build rotation matrix
			Matrix3x3d rot(up.x, right.x, dir.x,
				up.y, right.y, dir.y,
				up.z, right.z, dir.z);

			if (!NumberUtils<double>::equals(rot.det(), 1.0, 100 * NumberUtils<double>::tolerance()))
				throw ITLException("Invalid rotation matrix.");

			return rot;
		}
	}

	/**
	Extracts 2D slice from 3D image.
	@param img Image where the pixel data is extracted.
	@param pos Position of the center point of the slice.
	@param dir Direction of slice normal.
	@param slice Slice pixels will be set to this image. The size of the image must be the required size of the slice.
	@param pTouchesOriginalEdge If not 0, a value indicating whether the extracted slice touches the edge of the original is stored here.
	*/
	template<typename pixel_t> void getSlice(const Image<pixel_t>& img, const Vec3d& pos, Vec3d dir, Image<pixel_t>& slice, bool* pTouchesOriginalEdge = 0, const Interpolator<pixel_t, pixel_t, double>& interpolator = LinearInterpolator<pixel_t, pixel_t, double, double>(BoundaryCondition::Zero))
	{
		Matrix3x3d rot = internals::sliceRotationMatrix(dir);

		Vec2d sliceRadius = Vec2d((double)slice.width() - 1, (double)slice.height() - 1) / 2.0;

		for (coord_t yi = 0; yi < (coord_t)slice.height(); yi++)
		{
			for (coord_t xi = 0; xi < (coord_t)slice.width(); xi++)
			{
				Vec3d p(xi - (double)sliceRadius.x, (double)yi - (double)sliceRadius.y, 0);
				Vec3d x = rot * p + pos;
				Vec3c xc = round(x);

				if (img.isInImage(xc))
				{
					pixel_t pixel = interpolator(img, x);
					slice(xi, yi) = pixel;
				}
				else
				{
					slice(xi, yi) = 0;
					if (pTouchesOriginalEdge)
						*pTouchesOriginalEdge = true;
				}
			}
		}
	}

	/**
	Draws the slice that getSlice function would return, into the original image.
	*/
	template<typename pixel_t> void drawSlice(Image<pixel_t>& img, const Vec3d& pos, Vec3d dir, const Vec2d& sliceDimensions, pixel_t color)
	{
		Matrix3x3d rot = internals::sliceRotationMatrix(dir);

		Vec2d sliceRadius = Vec2d((double)sliceDimensions.x - 1, (double)sliceDimensions.y - 1) / 2.0;

		for (coord_t yi = 0; yi < (coord_t)sliceDimensions.y; yi++)
		{
			for (coord_t xi = 0; xi < (coord_t)sliceDimensions.x; xi++)
			{
				Vec3d p(xi - (double)sliceRadius.x, (double)yi - (double)sliceRadius.y, 0);
				Vec3d x = rot * p + pos;
				Vec3c xc = round(x);

				if (img.isInImage(xc))
					img(xc) = color;
			}
		}
	}
}