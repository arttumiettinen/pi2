#pragma once

#include "image.h"
#include "transform.h"

namespace itl2
{
	/**
	Creates a 2D montage from a 3D image.
	@param in Input image where the montage is created from. Usually a 3D stack.
	@param out The montage will be stored here.
	@param columns, rows Width and height of the montage. Unit = image count.
	@param firstSlice, lastSlice The first and the last slice to include in the montage.
	@param step Step between slices to include in the montage.
	@param borderWidth Width of border between slices in the montage.
	@param boderColor Color of border between slices in the montage.
	*/
	template<typename pixel_t> void montage(const Image<pixel_t>& in, Image<pixel_t>& out,
		size_t columns, size_t rows,
		double scale = 1.0,
		size_t firstSlice = 0, size_t lastSlice = std::numeric_limits<size_t>::max(),
		size_t step = 0,
		size_t borderWidth = 0,
		pixel_t borderColor = 0)
	{
		if (firstSlice > (size_t)in.depth())
			return;

		if (lastSlice > (size_t)in.depth() - 1)
			lastSlice = (size_t)in.depth() - 1;

		// Set step such that the whole stack is visible in the output
		if (step <= 0)
		{
			step = (lastSlice - firstSlice) / (rows * columns);
			if (step < 1)
				step = 1;
		}

		Vec3c scaledDimensions(itl2::round(in.width() * scale), itl2::round(in.height() * scale), 1);
		Vec3c outDimensions((scaledDimensions.x + borderWidth) * columns - borderWidth, (scaledDimensions.y + borderWidth) * rows - borderWidth, 1);

		out.ensureSize(outDimensions);
		if(borderColor != (pixel_t)0)
			setValue(out, borderColor);

		Image<pixel_t> slice(in.width(), in.height());
		Image<pixel_t> scaled(scaledDimensions);

		size_t row = 0;
		size_t column = 0;
		for (size_t n = firstSlice; n <= lastSlice; n += step)
		{
			crop(in, slice, Vec3c(0, 0, n));
			itl2::scale(slice, scaled);
			copyValues(out, scaled, Vec3c((scaledDimensions.x + borderWidth) * column, (scaledDimensions.y + borderWidth) * row, 0));
			column++;
			if (column >= columns)
			{
				column = 0;
				row++;
				if (row >= rows)
					return;
			}
		}
	}

	namespace tests
	{
		void montage();
	}
}