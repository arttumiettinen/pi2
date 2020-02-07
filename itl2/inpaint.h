#pragma once

#include <cmath>
#include <omp.h>
#include <limits>

#include "image.h"
#include "dmap.h"
#include "math/mathutils.h"
#include "math/vec3.h"
#include "fft.h"
#include "conversions.h"
#include "pointprocess.h"

namespace itl2
{

	template<typename pixel_t> bool isFlag(pixel_t pix, pixel_t val)
	{
		return pix == val;
	}

	template<> inline bool isFlag(float32_t pix, float32_t val)
	{
		if (std::isnan(val))
		{
			return std::isnan(pix);
		}

		if (std::isinf(val))
		{
			return std::isinf(pix);
		}

		return pix == val;
	}

	template<> inline bool isFlag(double pix, double val)
	{
		if (std::isnan(val))
		{
			return std::isnan(pix);
		}

		if (std::isinf(val))
		{
			return std::isinf(pix);
		}

		return pix == val;
	}

	/**
	Replaces val in the img by nearest non-val value.
	@param image Image to process.
	@param val Value that marks missing pixels.
	*/
	template<typename pixel_t> void inpaintNearest(Image<pixel_t>& image, pixel_t val = 0)
	{
		bool process = false;
		for (coord_t n = 0; n < image.pixelCount(); n++)
		{
			if (isFlag(image(n), val))
			{
				process = true;
				break;
			}
		}

		if (process)
		{

			Image<float32_t> distance(image.dimensions());
			Image<Vec3c> coords(image.dimensions());

			for (coord_t n = 0; n < image.pixelCount(); n++)
			{
				if (isFlag(image(n), val))
				{
					distance(n) = std::numeric_limits<float32_t>::max();
				}
				else
				{
					distance(n) = 0;
				}
			}

			distanceTransform(distance, &coords);

			#pragma omp parallel for if(image.pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
			for (coord_t n = 0; n < image.pixelCount(); n++)
			{
				Vec3c p = coords(n);
				if (image.isInImage(p))
					image(n) = image(p);
			}
		}
	}

	/**
	Replaces value val in the img by a value interpolated from nearby non-val pixels.

	The inpainting algorithm is from
    Garcia D, Robust smoothing of gridded data in one and higher dimensions with missing values. Computational Statistics & Data Analysis, 2010;54:1167-1178
    and
    Wang G, Garcia D et al. A three-dimensional gap filling method for large geophysical datasets: Application to global satellite soil moisture observations. Environ Modell Softw, 2012;30:139-142.

    The implementation is based on the reference implementation in Matlab.

	@param x Image that is to be inpainted.
	@param val Value that marks missing pixels.
	@param indicateProgress Set to true to show progress indicator.
	@param n Number of iterations.
	@param RF Relaxation factor.
	@param m Some undocumented parameter.
	*/
	template<typename pixel_t> void inpaintGarcia(Image<pixel_t>& img, pixel_t val = 0, bool indicateProgress = false, float32_t tolerance = 0, int n = 100, float32_t RF = 2, float32_t m = 2)
	{
		// Skip processing if there are no zeroes in the image
		bool process = false;
		for (coord_t m = 0; m < img.pixelCount(); m++)
		{
			if (isFlag(img(m), val))
			{
				process = true;
				break;
			}
		}

		if (!process)
			return;

		// Lambda matrix
		Image<float32_t> lambda(img.dimensions());

		#pragma omp parallel for if(img.pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		for (coord_t z = 0; z < img.depth(); z++)
		{
			float32_t lz = (float32_t)std::cos(PI * z / img.depth());
			if (img.depth() <= 1)
				lz = 0;
			for (coord_t y = 0; y < img.height(); y++)
			{
				float32_t ly = (float32_t)std::cos(PI * y / img.height());
				if (img.height() <= 1)
					lz = 0;
				for (coord_t x = 0; x < img.width(); x++)
				{
					float32_t lx = (float32_t)std::cos(PI * x / img.width());

					lambda(x, y, z) = (float32_t)::pow(2 * (img.dimensionality() - (lz + ly + lx)), m);
				}
			}
		}

		// Initial guess
		Image<float32_t> y(img.dimensions());
		convert(img, y);
		inpaintNearest<float32_t>(y, pixelRound<float32_t>(val)); // TODO: This might not work for all data types and val values.

		// Inpainting iterations
		Image<float32_t> tmp(img.dimensions());
		for (coord_t i = 0; i < n; i++)
		{
			// Smoothness parameter range
			constexpr float32_t start = 3;
			constexpr float32_t end = -6;
			float32_t si = (float32_t)::pow(10.0, start + (end - start) / (n - 1) * i);

			#pragma omp parallel for if(img.pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
			for (coord_t m = 0; m < img.pixelCount(); m++)
			{
				tmp(m) = y(m);
				if (!isFlag(img(m), val))
				{
					tmp(m) += (float32_t)(img(m) - y(m));
				}
			}

			dct(tmp);

			#pragma omp parallel for if(img.pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
			for (coord_t m = 0; m < img.pixelCount(); m++)
			{
				float32_t gamma = 1 / (1 + si * lambda(m));
				tmp(m) *= gamma;
			}

			idct(tmp);


			// Calculate
			// y_new = y * (1 - RF) + tmp * RF
			// with simple algorithm
			//multiply(tmp, RF);
			//multiply(y, 1 - RF);
			//add(y, tmp);

			// This calculation should reduce overhead and provide maximum absolute change in y.
			float32_t maxDiff = 0;
			#pragma omp parallel if(img.pixelCount() > PARALLELIZATION_THRESHOLD)
			{
				float32_t res_private = 0;
				#pragma omp for nowait
				for (coord_t m = 0; m < img.pixelCount(); m++)
				{
					float32_t y_new = y(m) * (1 - RF) + tmp(m) * RF;
					float32_t diff = ::abs(y(m) - y_new);
					y(m) = y_new;

					if (isFlag(img(m), val) && diff > res_private)
						res_private = diff;
				}

				#pragma omp critical(inpaint_maximum_reduction)
				{
					if (res_private > maxDiff)
						maxDiff = res_private;
				}
			}

			//cout << maxDiff << endl;

			if (maxDiff <= tolerance)
			{ 
				if (indicateProgress)
					std::cout << "\r" << std::flush;
				break;
			}

			if (indicateProgress)
				showProgress(i, n);
		}

		// Pick output pixels
		#pragma omp parallel for if(img.pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		for (coord_t m = 0; m < img.pixelCount(); m++)
		{
			if (isFlag(img(m), val))
				img(m) = pixelRound<pixel_t>(y(m));
		}
	}

	namespace tests
	{
		void inpaintNearest();
		void inpaintGarcia();
		void inpaintGarcia2();
	}
}

