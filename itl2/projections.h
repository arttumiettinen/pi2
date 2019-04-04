#pragma once

#include "image.h"
#include "utilities.h"
#include "math/mathutils.h"
#include "pointprocess.h"

#include <set>

using namespace math;

namespace itl2
{

	/**
	Project one dimension of the image.
	Use to create x, y and z projections.
	*/
	template<typename pixel_t, typename out_t, void process(pixel_t, double&)> void projectDimension(const Image<pixel_t>& img, size_t dimension, Image<out_t>& out, double initialValue, bool showProgressInfo = true)
	{
		//// This version works but is relatively complicated compared to the simple version below.
		//// Z project -> zstep = inf, ystep = 1, xstep = 1, zstep2 = 1, ystep2 = inf, xstep2 = inf; resultsize = [xdim, ydim]
		//// Y project -> zstep = 1, ystep = inf, xstep = 1, zstep2 = inf, ystep2 = 1, xstep2 = inf; resultsize = [xdim, zdim]
		//// X project -> zstep = 1, ystep = 1, xstep = inf, zstep2 = inf, ystep2 = inf, xstep2 = 1; resultsize = [ydim, zdim]

		//// Step sizes for outer loops
		//coord_t steps[3];
		//steps[0] = 1;
		//steps[1] = 1;
		//steps[2] = 1;
		//steps[dimension] = img.dimension(dimension);

		//// Step sizes for inner loops
		//coord_t steps2[3];
		//steps2[0] = img.dimension(0);
		//steps2[1] = img.dimension(1);
		//steps2[2] = img.dimension(2);
		//steps2[dimension] = 1;

		//// Size of result image
		//vector<coord_t> resultsize;
		//resultsize.push_back(img.dimension(0));
		//resultsize.push_back(img.dimension(1));
		//resultsize.push_back(img.dimension(2));
		//resultsize[dimension] = 1;

		//out.init(resultsize[0], resultsize[1], resultsize[2]);

		//for (coord_t z = 0; z < img.depth(); z += steps[2])
		//{
		//	for (coord_t y = 0; y < img.height(); y += steps[1])
		//	{
		//		for (coord_t x = 0; x < img.width(); x += steps[0])
		//		{
		//			double res = 0;

		//			for (coord_t z2 = steps[2] == 1 ? z : 0; z2 < img.depth(); z2 += steps2[2])
		//			{
		//				for (coord_t y2 = steps[1] == 1 ? y : 0; y2 < img.height(); y2 += steps2[1])
		//				{
		//					for (coord_t x2 = steps[0] == 1 ? x : 0; x2 < img.width(); x2 += steps2[0])
		//					{
		//						res += img(x2, y2, z2);
		//					}
		//				}
		//			}

		//			out(x, y, z) = pixelRound<out_t>(res);
		//		}
		//	}
		//}
		//out.squeeze();



		// Simpler version
		size_t counter = 0;
		if (dimension == 2)
		{
			// Z project
			out.init(img.width(), img.height());

			#pragma omp parallel for if(img.pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
			for (coord_t y = 0; y < img.height(); y++)
			{
				for (coord_t x = 0; x < img.width(); x++)
				{
					double res = initialValue;
					for (coord_t z = 0; z < img.depth(); z++)
					{
						process(img(x, y, z), res);
					}
					out(x, y) = pixelRound<out_t>(res);
				}

				showThreadProgress(counter, img.height(), showProgressInfo);
			}
		}
		else if (dimension == 1)
		{

			// Y project
			out.init(img.width(), img.depth());

			#pragma omp parallel for if(img.pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
			for (coord_t z = 0; z < img.depth(); z++)
			{
				for (coord_t x = 0; x < img.width(); x++)
				{
					double res = initialValue;
					for (coord_t y = 0; y < img.height(); y++)
					{
						process(img(x, y, z), res);
					}
					out(x, z) = pixelRound<out_t>(res);
				}

				showThreadProgress(counter, img.depth(), showProgressInfo);
			}
		}
		else if (dimension == 0)
		{
			// X project. Swap z and y in output to make the image more logical (in some sense...)
			out.init(img.depth(), img.height());

			#pragma omp parallel for if(img.pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
			for (coord_t y = 0; y < img.height(); y++)
			{
				for (coord_t z = 0; z < img.depth(); z++)
				{
					double res = initialValue;
					for (coord_t x = 0; x < img.width(); x++)
					{
						process(img(x, y, z), res);
					}
					out(out.width() - z - 1, y) = pixelRound<out_t>(res);
				}

				showThreadProgress(counter, img.height(), showProgressInfo);
			}
		}
		else
		{
			throw ITLException("Invalid dimension.");
		}

	}


	namespace internals
	{
		template<typename pixel_t> void sumProjectionOp(pixel_t p, double& result)
		{
			result += (double)p;
		}

		template<typename pixel_t> void minProjectionOp(pixel_t p, double& result)
		{
			if (p < result)
				result = (double)p;
		}

		template<typename pixel_t> void maxProjectionOp(pixel_t p, double& result)
		{
			if (p > result)
				result = (double)p;
		}

		template<typename pixel_t> void squareSumProjectionOp(pixel_t val, double& result)
		{
			result = result + (double)val * (double)val;
		}

		template<typename pixel_t> void crossSumProjectionOp(pixel_t val1, pixel_t val2, double& result)
		{
			result = result + (double)val1 * (double)val2;
		}

		template<typename pixel_t> void absDifferenceProjectionOp(pixel_t val1, pixel_t val2, double& result)
		{
			result = result + std::abs((double)val1 - (double)val2);
		}
	}


	/**
	Determines if pixel values in the two images are equal.
	@param a, b Images to compare
	*/
	template<typename pixel1_t, typename pixel2_t> bool equals(const Image<pixel1_t>& a, const Image<pixel2_t>& b)
	{
		a.checkSize(b);

		bool eq = true;
		#pragma omp parallel for if(a.pixelCount() > PARALLELIZATION_THRESHOLD)
		for (coord_t n = 0; n < a.pixelCount(); n++)
		{
			if (eq && a(n) != b(n))
				eq = false;

			// Showing progress info here would induce more processing than is done in the whole loop.
		}

		return eq;
	}

	/**
	Determines if pixel values in the two images are not equal.
	@param a, b Images to compare
	*/
	template<typename pixel1_t, typename pixel2_t> bool differs(const Image<pixel1_t>& a, const Image<pixel2_t>& b)
	{
		return !equals(a, b);
	}

	/**
	Finds all unique pixel values in the image and adds them to the set.
	*/
	template<typename pixel_t> void unique(const Image<pixel_t>& img, std::set<pixel_t>& values)
	{
		// TODO: Parallelize
		for (coord_t n = 0; n < img.pixelCount(); n++)
		{
			values.emplace(img(n));
		}
	}


	/**
	Calculates sum of all pixels in the image.
	@param img Image to process.
	@return Sum of all pixels.
	*/
	template<typename pixel_t, typename out_t = double> out_t sum(const Image<pixel_t>& img)
	{
		out_t res = out_t();
		#pragma omp parallel if(img.pixelCount() > PARALLELIZATION_THRESHOLD)
		{
			out_t res_private = out_t();
			#pragma omp for nowait
			for (coord_t n = 0; n < img.pixelCount(); n++)
			{
				res_private += img(n);

				// Showing progress info here would induce more processing than is done in the whole loop.
			}

			#pragma omp critical(sum_reduction)
			{
				res += res_private;
			}
		}

		return res;

		// OpenMP reduction does not support non-scalar data types
		//out_t res = 0;
		//#pragma omp parallel for if(img.pixelCount() > PARALLELIZATION_THRESHOLD) reduction(+:res)
		//for (coord_t n = 0; n < img.pixelCount(); n++)
		//{
		//	res += img(n);

		//	// Showing progress info here would induce more processing than is done in the whole loop.
		//}
		//return res;

	}

	/*
	Calculates sum of all pixels except those whose value is ignoreValue.
	Counts pixels that contributes to the sum and places the count to count argument.
	@param img Image to process.
	@param ignoreValue Value that should be ignored while calculating sum.
	@param count Count of pixels that do not have ignoreValue.
	@return Sum of all pixels except those whose value is ignoreValue.
	*/
	template<typename pixel_t, typename out_t = double> out_t maskedsum(const Image<pixel_t>& img, pixel_t ignoreValue, out_t& count)
	{
		if (!NumberUtils<pixel_t>::isnan(ignoreValue))
		{
			out_t res = out_t();
			count = 0;
#pragma omp parallel if(img.pixelCount() > PARALLELIZATION_THRESHOLD)
			{
				out_t res_private = out_t();
				out_t count_private = out_t();
#pragma omp for nowait
				for (coord_t n = 0; n < img.pixelCount(); n++)
				{
					pixel_t p = img(n);
					if (p != ignoreValue)
					{
						res_private += p;
						count_private++;
					}

					// Showing progress info here would induce more processing than is done in the whole loop.
				}

#pragma omp critical(maskedsum_reduction)
				{
					res += res_private;
					count += count_private;
				}
			}

			return res;
		}
		else
		{
			// ignoreValue is NaN, comparison logic must be different
			out_t res = out_t();
			count = 0;
#pragma omp parallel if(img.pixelCount() > PARALLELIZATION_THRESHOLD)
			{
				out_t res_private = out_t();
				out_t count_private = out_t();
#pragma omp for nowait
				for (coord_t n = 0; n < img.pixelCount(); n++)
				{
					pixel_t p = img(n);
					if (!NumberUtils<pixel_t>::isnan(p))
					{
						res_private += p;
						count_private++;
					}

					// Showing progress info here would induce more processing than is done in the whole loop.
				}

#pragma omp critical(maskedsum_reduction)
				{
					res += res_private;
					count += count_private;
				}
			}

			return res;
		}

		// OpenMP reduction does not support non-scalar data types
		//out_t res = 0;
		//#pragma omp parallel for if(img.pixelCount() > PARALLELIZATION_THRESHOLD) reduction(+:res)
		//for (coord_t n = 0; n < img.pixelCount(); n++)
		//{
		//	res += img(n);

		//	// Showing progress info here would induce more processing than is done in the whole loop.
		//}
		//return res;

	}

	/**
	Calculates minimum of all pixels in the image.
	@param img Image to process.
	@return Minimum of all pixels.
	*/
	template<typename pixel_t> pixel_t min(const Image<pixel_t>& img)
	{
		pixel_t res = numeric_limits<pixel_t>::max();
		#pragma omp parallel if(img.pixelCount() > PARALLELIZATION_THRESHOLD)
		{
			pixel_t res_private = numeric_limits<pixel_t>::max();
			#pragma omp for nowait
			for (coord_t n = 0; n < img.pixelCount(); n++)
			{
				if (img(n) < res_private)
					res_private = img(n);

				// Showing progress info here would induce more processing than is done in the whole loop.
			}

			#pragma omp critical(minimum_reduction)
			{
				if (res_private < res)
					res = res_private;
			}
		}

		return res;

		// MSVC does not support min reduction
		//double res;
		//#pragma omp parallel for reduction(min:res)
		//for (coord_t n = 0; n < img.pixelCount(); n++)
		//{
		//	if (img(n) < res)
		//		res = img(n);
		//}
		//return res;

		//return projectAllDimensions<pixel_t, internals::minProjectionOp<pixel_t> >(img);
	}

	/**
	Calculates maximum of all pixels in the image.
	@param img Image to process.
	@return Maximum of all pixels.
	*/
	template<typename pixel_t> pixel_t max(const Image<pixel_t>& img)
	{
		pixel_t res = numeric_limits<pixel_t>::lowest();
		#pragma omp parallel if(img.pixelCount() > PARALLELIZATION_THRESHOLD)
		{
			pixel_t res_private = numeric_limits<pixel_t>::lowest();
			#pragma omp for nowait
			for (coord_t n = 0; n < img.pixelCount(); n++)
			{
				if (img(n) > res_private)
					res_private = img(n);

				// Showing progress info here would induce more processing than is done in the whole loop.
			}

			#pragma omp critical(maximum_reduction)
			{
				if (res_private > res)
					res = res_private;
			}
		}

		return res;

		// MSVC Does not support max reduction
		//double res;
		//#pragma omp parallel for reduction(max:res)
		//for (coord_t n = 0; n < img.pixelCount(); n++)
		//{
		//	if (img(n) > res)
		//		res = img(n);
		//}
		//return res;

		//return projectAllDimensions<pixel_t, internals::maxProjectionOp<pixel_t> >(img);
	}

	/**
	Calculates mean all pixels in the image.
	@param img Image to process.
	@return Mean of all pixels.
	*/
	template<typename pixel_t, typename out_t = double> out_t mean(const Image<pixel_t>& img)
	{
		return sum<pixel_t, out_t>(img) / (out_t)img.pixelCount();
	}

	/**
	Calculates mean of all pixels in the image, except those that have specified value.
	@param img Image to process.
	@param ignoreValue Value that should not be accounted for when calculating mean.
	@return Mean of all pixels.
	*/
	template<typename pixel_t, typename out_t = double> out_t maskedmean(const Image<pixel_t>& img, pixel_t ignoreValue)
	{
		out_t count;
		out_t s = maskedsum<pixel_t, out_t>(img, ignoreValue, count);
		return s / count;
	}

	/**
	Calculates mean and standard deviation of all pixels in the image.
	@param img Image to process.
	@return Vec2 containing mean at Vec2::x and standard deviation at Vec2::y.
	*/
	template<typename pixel_t, typename out_t = double> Vec2<out_t> meanAndStdDev(const Image<pixel_t>& img)
	{
		// Non-threaded version
		//out_t total = out_t();
		//out_t total2 = out_t();

		//for (coord_t n = 0; n < img.pixelCount(); n++)
		//{
		//	out_t val = (out_t)img(n);
		//	total += val;
		//	total2 += val * val;
		//}


		//out_t avg = total / img.pixelCount();
		//out_t std = sqrt((total2 - (total * total / img.pixelCount())) / (img.pixelCount() - 1));

		//return Vec2<out_t>(avg, std);

		out_t total = out_t();
		out_t total2 = out_t();

#pragma omp parallel if(img.pixelCount() > PARALLELIZATION_THRESHOLD)
		{
			out_t total_private = out_t();
			out_t total2_private = out_t();

#pragma omp for nowait
			for (coord_t n = 0; n < img.pixelCount(); n++)
			{
				out_t val = (out_t)img(n);
				total_private += val;
				total2_private += val * val;

				// Showing progress info here would induce more processing than is done in the whole loop.
			}

#pragma omp critical(sum_reduction)
			{
				total += total_private;
				total2 += total2_private;
			}
		}

		out_t avg = total / img.pixelCount();
		out_t std = sqrt((total2 - (total * total / img.pixelCount())) / (img.pixelCount() - 1));

		return Vec2<out_t>(avg, std);
	}

	/**
	Calculates standard deviation of all pixels in the image.
	@param img Image to process.
	@return Standard deviation of all pixels.
	*/
	template<typename pixel_t, typename out_t = double> out_t stdDev(const Image<pixel_t>& img)
	{
		Vec2<out_t> v = meanAndStdDev(img);
		return v.y;
	}

	/**
	Sum projection in some coordinate direction.
	@param img Input image.
	@param dimension Zero-based dimension over which the projection should be made.
	@param out Output image.
	@param showProgressInfo Set to true to show a progress bar.
	*/
	template<typename pixel_t, typename out_t> void sum(const Image<pixel_t>& img, size_t dimension, Image<out_t>& out, bool showProgressInfo = true)
	{
		projectDimension<pixel_t, out_t, internals::sumProjectionOp<pixel_t>>(img, dimension, out, 0, showProgressInfo);
	}

	/**
	Minimum projection in some coordinate direction.
	@param img Input image.
	@param dimension Zero-based dimension over which the projection should be made.
	@param out Output image.
	@param showProgressInfo Set to true to show a progress bar.
	*/
	template<typename pixel_t, typename out_t> void min(const Image<pixel_t>& img, size_t dimension, Image<out_t>& out, bool showProgressInfo = true)
	{
		projectDimension<pixel_t, out_t, internals::minProjectionOp<pixel_t> >(img, dimension, out, numeric_limits<double>::infinity(), showProgressInfo);
	}

	/**
	Maximum projection in some coordinate direction.
	@param img Input image.
	@param dimension Zero-based dimension over which the projection should be made.
	@param out Output image.
	@param showProgressInfo Set to true to show a progress bar.
	*/
	template<typename pixel_t, typename out_t> void max(const Image<pixel_t>& img, size_t dimension, Image<out_t>& out, bool showProgressInfo = true)
	{
		projectDimension<pixel_t, out_t, internals::maxProjectionOp<pixel_t> >(img, dimension, out, -numeric_limits<double>::infinity(), showProgressInfo);
	}

	/**
	Mean projection in some coordinate direction.
	@param img Input image.
	@param dimension Zero-based dimension over which the projection should be made.
	@param out Output image.
	@param showProgressInfo Set to true to show a progress bar.
	*/
	template<typename pixel_t, typename out_t> void mean(const Image<pixel_t>& img, size_t dimension, Image<out_t>& out, bool showProgressInfo = true)
	{
		sum(img, dimension, out, showProgressInfo);
		divide(out, (double)img.dimension(dimension));
	}






	/*
	Special projections
	*/


	template<typename pixel_t, void process(pixel_t, double&)> double projectAllDimensions(const Image<pixel_t>& img)
	{
		// Single-threaded version for use when no reduction is available.
		double result = 0;
		for (coord_t z = 0; z < img.depth(); z++)
		{
			for (coord_t y = 0; y < img.height(); y++)
			{
				for (coord_t x = 0; x < img.width(); x++)
				{
					process(img(x, y, z), result);
				}
			}
		}

		return result;
	}

	template<typename pixel_t, void process(pixel_t, pixel_t, double&)> double projectAllDimensions(const Image<pixel_t>& img1, const Image<pixel_t>& img2)
	{
		// Single-threaded version for use when no reduction is available.
		img1.checkSize(img2);
		double result = 0;
		for (coord_t z = 0; z < img1.depth(); z++)
		{
			for (coord_t y = 0; y < img1.height(); y++)
			{
				for (coord_t x = 0; x < img1.width(); x++)
				{
					process(img1(x, y, z), img2(x, y, z), result);
				}
			}
		}

		return result;
	}

	template<typename pixel_t> double twoNormSquared(const Image<pixel_t>& img)
	{
		return projectAllDimensions<pixel_t, internals::squareSumProjectionOp<pixel_t> >(img);
	}

	template<typename pixel_t> double crossSum(const Image<pixel_t>& img1, const Image<pixel_t>& img2)
	{
		return projectAllDimensions<pixel_t, internals::crossSumProjectionOp<pixel_t> >(img1, img2);
	}

	/**
	Calculates total absolute difference between two images, i.e. sum_r(abs(img1(r) - img2(r))), where the sum is taken over all pixels.
	*/
	template<typename pixel_t> double absDifference(const Image<pixel_t>& img1, const Image<pixel_t>& img2)
	{
		return projectAllDimensions<pixel_t, internals::absDifferenceProjectionOp<pixel_t> >(img1, img2);
	}


	namespace tests
	{
		void projections();
		void projections2();
	}
}
