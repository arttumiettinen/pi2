#pragma once

#include <iostream>
#include <cmath>
#include <array>
#include <tuple>
#include <numeric>

#include "image.h"
#include "utilities.h"
#include "math/vec2.h"
#include "math/vec3.h"
#include "pointprocess.h"
#include "progress.h"

namespace itl2
{
	/**
	Selects suitable type for intermediate sum image in histogram calculation:
	hist_t or weight_t are floating point -> double
	hist_t and weight_t are integers -> uint64_t
	*/
	template<class hist_t, class weight_t> struct histogram_intermediate_type {
		using type = typename std::conditional <
			std::is_floating_point_v<hist_t> || std::is_floating_point_v<weight_t>,
			double,
			uint64_t
		>::type;
	};

	/**
	Calculates unweighted or weighted histogram of input image.
	@param img Image whose histogram is calculated.
	@param histogram Image that will store the histogram. Bin count of the histogram is pixel count of this image. Make sure that the pixel data type of the image can hold big enough values; uint64 or float32 image is recommended.
	@param range Gray value range for the histogram. Pixels out of range are counted in first or last bins.
	@param edgeSkip This many pixels at the image edge are not considered in the histogram.
	@param pWeight Pointer to image that stores weight of each pixel.
	*/
	template<typename pixel_t, typename hist_t, typename weight_t = hist_t> void histogram(const Image<pixel_t>& img, Image<hist_t>& histogram, const Vec2d& range, coord_t edgeSkip = 0, const Image<weight_t>* pWeight = nullptr)
	{
		if (pWeight)
			pWeight->checkSize(img);

		coord_t dim = histogram.pixelCount();
		size_t counter = 0;

		using sum_t = typename histogram_intermediate_type<hist_t, weight_t>::type;
		Image<sum_t> sums(histogram.dimensions());

		ProgressIndicator progress(img.depth());
		#pragma omp parallel if(img.pixelCount() > PARALLELIZATION_THRESHOLD)
		{
			Image<sum_t> privateHist(sums.dimensions());
			#pragma omp for nowait
			for (coord_t z = 0; z < img.depth(); z++)
			{
				for (coord_t y = 0; y < img.height(); y++)
				{
					for (coord_t x = 0; x < img.width(); x++)
					{
						if (img.edgeDistance(Vec3c(x, y, z)) >= edgeSkip)
						{
							pixel_t pix = img(x, y, z);
							coord_t bin = itl2::floor(((pix - range.x) / (range.y - range.x)) * (double)dim);

							// Problem with above expression is that the terms inside floor() may give, e.g. 0.2899999998 for pixel
							// that should go to bin 290-300. The floor makes it end in bin 280-290.
							// Check that the pixel really belongs to the bin determined using above expression, and adjust if necessary.
							// TODO: There is probably some numerically stable algorithm that does not need this check.
							pixel_t binMin = pixelRound<pixel_t>(range.x + (double)bin / (double)dim * (range.y - range.x));
							pixel_t binMax = pixelRound<pixel_t>(range.x + (double)(bin + 1) / (double)dim * (range.y - range.x));
							if (pix < binMin)
								bin--;
							else if (pix >= binMax)
								bin++;

							if (bin < 0)
								bin = 0;
							else if (bin >= dim)
								bin = dim - 1;

							if (!pWeight)
								privateHist(bin)++;
							else
								privateHist(bin) += (sum_t)(*pWeight)(x, y, z);
						}
					}
				}

				progress.step();
			}

			#pragma omp critical(histogram_reduction)
			{
				for (coord_t n = 0; n < privateHist.pixelCount(); n++)
				{
					sums(n) += privateHist(n);
				}
			}
		}

		setValue(histogram, sums);

		//coord_t dim = histogram.pixelCount();
		//size_t counter = 0;
		//for (coord_t z = 0; z < img.depth(); z++)
		//{
		//	for (coord_t y = 0; y < img.height(); y++)
		//	{
		//		for (coord_t x = 0; x < img.width(); x++)
		//		{
		//			if (edgeDistance(Vec3c(x, y, z), img.dimensions()) >= edgeSkip)
		//			{
		//				pixel_t pix = img(x, y, z);
		//				coord_t bin = (coord_t)floor(((pix - range.min) / (range.max - range.min)) * (double)dim);

		//				// Problem with above expression is that the terms inside floor() may give, e.g. 0.2899999998 for pixel
		//				// that should go to bin 290-300. The floor makes it end in bin 280-290.
		//				// Check that the pixel really belongs to the bin determined using above expression, and adjust if necessary.
		//				// TODO: There is probably some numerically stable algorithm that does not need this check.
		//				pixel_t binMin = pixelRound<pixel_t>(range.min + (double)bin / (double)dim * (range.max - range.min));
		//				pixel_t binMax = pixelRound<pixel_t>(range.min + (double)(bin + 1) / (double)dim * (range.max - range.min));
		//				if (pix < binMin)
		//					bin--;
		//				else if (pix >= binMax)
		//					bin++;

		//				if (bin < 0)
		//					bin = 0;
		//				else if (bin >= dim)
		//					bin = dim - 1;

		//				histogram(bin)++;
		//			}
		//		}
		//	}
		//}
	}

	/**
	Encapsulates image and corresponding histogram range and bin count.
	*/
	template<typename pixel_t> struct ImageAndRange
	{
		const Image<pixel_t>& image;
		const Vec2d range;
		size_t binCount;

		ImageAndRange(const Image<pixel_t>& img, const Vec2d& range, size_t binCount) :
			image(img),
			range(range),
			binCount(binCount)
		{

		}
	};

	namespace internals
	{
		template <std::size_t From, size_t... indices, typename T1, typename T2, typename Func>
		void transform(T1&& s, T2& t, Func f, std::index_sequence<indices...>)
		{
			(void)std::initializer_list<int>{
				(std::get<indices + From>(t) = f(std::get<indices>(std::forward<T1>(s))), 0)...};
		}


		/**
		Transform range from index 'From' to and including 'To' of tuple or array s to tuple or array t.
		*/
		template <std::size_t From, std::size_t To, typename T1, typename T2, typename Func>
		void transform(T1&& s, T2& t, Func f)
		{
			transform<From>(std::forward<T1>(s), t, f, std::make_index_sequence<To - From + 1>());
		}
	}


	/**
	Calculates unweighted or weighted multivariate histogram of input images.
	TODO: As Image class supports only up to 3-dimensional images, this method is limited to 3-dimensional histograms.
	@param histogram Image that will store the histogram. Dimensionality of this image must equal count of images in imgs list. Size in each dimension is bin count in that dimension. Make sure that the pixel data type of the image can hold big enough values; floating point image is recommended.
	@param edgeSkip This many pixels at the image edge are not considered in the histogram.
	@param pWeight Pointer to image where weight for each pixel will be read. Pass nullptr to use uniform weight for each pixel. The size of this image must equal the size of the input images.
	@param imgs Images whose joint histogram is to be calculated, together with histogram range and bin count for each image. Each image will become one dimension in the output histogram. Pixels out of range are counted in first or last bins.
	*/
	template<typename hist_t, typename weight_t, typename ... pixel_t> void multiHistogram(Image<hist_t>& histogram, coord_t edgeSkip, const Image<weight_t>* pWeight = nullptr, const ImageAndRange<pixel_t>&... imgs)
	{
		constexpr coord_t N = sizeof...(pixel_t);
		if constexpr (N <= 0)
			return;

		auto&& t = std::forward_as_tuple(imgs...);
		
		// Check that all input images have the same dimensions.
		//std::array<int, N> dummy;
		//internals::transform<0, N - 1>(t, dummy, [&](auto i) { i.image.checkSize(dims); return 0; });
		//Vec3c dims = get<0>(t).image.dimensions();

		// Calculate minimum dimensions where we have data available.
		std::array<Vec3c, N> allDimensions;
		internals::transform<0, N - 1>(t, allDimensions, [&](auto i) { return i.image.dimensions(); });
		//Vec3c minDims = componentwiseMin(allDimensions);
		Vec3c minDims = std::accumulate<typename std::array<Vec3c, N>::const_iterator, Vec3c, Vec3c (const Vec3c&, const Vec3c&)>(allDimensions.begin(), allDimensions.end(), Vec3c(std::numeric_limits<coord_t>::max(), std::numeric_limits<coord_t>::max(), std::numeric_limits<coord_t>::max()), min<coord_t>);

		if (pWeight)
			pWeight->checkSize(minDims);

		if (edgeSkip > 2 * minDims.min())
			throw ITLException("Too large edge skip value.");

		// Initialize output image
		std::array<size_t, N> binCounts;
		internals::transform<0, N - 1>(t, binCounts, [&](auto i) { return i.binCount; });

		if constexpr (N == 1)
			histogram.ensureSize(binCounts[0]);
		else if constexpr (N == 2)
			histogram.ensureSize(binCounts[0], binCounts[1]);
		else if constexpr (N == 3)
			histogram.ensureSize(binCounts[0], binCounts[1], binCounts[2]);
		else
			static_assert("Only 1-, 2-, and 3-dimensional histograms are supported.");

		using sum_t = typename histogram_intermediate_type<hist_t, weight_t>::type;
		Image<sum_t> sums(histogram.dimensions());

		ProgressIndicator progress(std::max((coord_t)0, minDims.z - 2 * edgeSkip));
		#pragma omp parallel if(minDims.x * minDims.y * minDims.z > PARALLELIZATION_THRESHOLD)
		{
			Image<sum_t> privateHist(sums.dimensions());
			#pragma omp for nowait
			for (coord_t z = edgeSkip; z < minDims.z - edgeSkip; z++)
			{
				for (coord_t y = edgeSkip; y < minDims.y - edgeSkip; y++)
				{
					for (coord_t x = edgeSkip; x < minDims.x - edgeSkip; x++)
					{
						// Calculate bin for each input image
						std::array<coord_t, N> ndbin;
						internals::transform<0, N - 1>(t, ndbin, [&](auto i) {
							auto& img = i.image;
							Vec2d range = i.range;
							double dim = (double)i.binCount;

							auto pix = img(x, y, z);
							coord_t bin = (coord_t)floor((((double)pix - range.x) / (range.y - range.x)) * dim);

							// Problem with above expression is that the terms inside floor() may give, e.g. 0.2899999998 for pixel
							// that should go to bin 290-300. The floor makes it end in bin 280-290.
							// Check that the pixel really belongs to the bin determined using above expression, and adjust if necessary.
							// TODO: There is probably some numerically stable algorithm that does not need this check.
							auto binMin = pixelRound<decltype(pix)>(range.x + (double)bin / dim * (range.y - range.x));
							auto binMax = pixelRound<decltype(pix)>(range.x + (double)(bin + 1) / dim * (range.y - range.x));
							if (pix < binMin)
								bin--;
							else if (pix >= binMax)
								bin++;

							if (bin < 0)
								bin = 0;
							else if ((size_t)bin >= i.binCount)
								bin = i.binCount - 1;

							return bin;
						});

						Vec3c ndbinv;
						for (size_t n = 0; n < N; n++)
							ndbinv[n] = ndbin[n];

						if (!pWeight)
							privateHist(ndbinv)++;
						else
							privateHist(ndbinv) += (sum_t)(*pWeight)(x, y, z);
					}
				}

				progress.step();
			}

			#pragma omp critical(ndhistogram_reduction)
			{
				for (coord_t n = 0; n < privateHist.pixelCount(); n++)
				{
					sums(n) += privateHist(n);
				}
			}
		}

		setValue(histogram, sums);
	}

	template<typename hist_t, typename ... pixel_t> void multiHistogram(Image<hist_t>& histogram, coord_t edgeSkip, const ImageAndRange<pixel_t>&&... imgs)
	{
		multiHistogram<hist_t, hist_t, pixel_t...>(histogram, edgeSkip, nullptr, imgs...);
	}


	namespace tests
	{
		void histogram();
		void histogram2d();
		void histogramIntermediateType();
	}

}
