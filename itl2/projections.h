#pragma once

#include "image.h"
#include "utilities.h"
#include "math/mathutils.h"
#include "pointprocess.h"

#include <set>
#include <unordered_set>
#include <limits>

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



	/**
	Project one dimension of the image.
	Use to create x, y and z projections that pick output value from second image.
	*/
	template<typename pixel_t, typename val_t, typename out_t, typename outval_t, void process(pixel_t, double&, val_t, val_t&)>
	void projectDimension(const Image<pixel_t>& img, const Image<val_t>& valImg, size_t dimension, Image<out_t>& out, Image<val_t>& outVal, double initialValue, val_t initialVal, bool showProgressInfo = true)
	{
		img.checkSize(valImg);

		size_t counter = 0;
		if (dimension == 2)
		{
			// Z project
			out.init(img.width(), img.height());
			outVal.init(img.width(), img.height());

#pragma omp parallel for if(img.pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
			for (coord_t y = 0; y < img.height(); y++)
			{
				for (coord_t x = 0; x < img.width(); x++)
				{
					double res = initialValue;
					val_t val = initialVal;
					for (coord_t z = 0; z < img.depth(); z++)
					{
						process(img(x, y, z), res, valImg(x, y, z), val);
					}
					out(x, y) = pixelRound<out_t>(res);
					outVal(x, y) = val;
				}

				showThreadProgress(counter, img.height(), showProgressInfo);
			}
		}
		else if (dimension == 1)
		{

			// Y project
			out.init(img.width(), img.depth());
			outVal.init(img.width(), img.depth());

#pragma omp parallel for if(img.pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
			for (coord_t z = 0; z < img.depth(); z++)
			{
				for (coord_t x = 0; x < img.width(); x++)
				{
					double res = initialValue;
					val_t val = initialVal;
					for (coord_t y = 0; y < img.height(); y++)
					{
						process(img(x, y, z), res, valImg(x, y, z), val);
					}
					out(x, z) = pixelRound<out_t>(res);
					outVal(x, z) = val;
				}

				showThreadProgress(counter, img.depth(), showProgressInfo);
			}
		}
		else if (dimension == 0)
		{
			// X project. Swap z and y in output to make the image more logical (in some sense...)
			out.init(img.depth(), img.height());
			outVal.init(img.depth(), img.height());

#pragma omp parallel for if(img.pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
			for (coord_t y = 0; y < img.height(); y++)
			{
				for (coord_t z = 0; z < img.depth(); z++)
				{
					double res = initialValue;
					val_t val = initialVal;
					for (coord_t x = 0; x < img.width(); x++)
					{
						process(img(x, y, z), res, valImg(x, y, z), val);
					}
					out(out.width() - z - 1, y) = pixelRound<out_t>(res);
					outVal(out.width() - z - 1, y) = val;
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

		template<typename pixel_t, typename val_t> void minProjectionValOp(pixel_t p, double& result, val_t val, val_t& outVal)
		{
			if (p < result)
			{
				result = (double)p;
				outVal = val;
			}
		}

		template<typename pixel_t, typename val_t> void maxProjectionValOp(pixel_t p, double& result, val_t val, val_t& outVal)
		{
			if (p > result)
			{
				result = (double)p;
				outVal = val;
			}
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
	Determines if all pixel values equal value v.
	@param img Image to compare
	@param v value to compare
	*/
	template<typename pixel_t> bool allEquals(const Image<pixel_t>& img, const pixel_t& v)
	{
		bool eq = true;
		#pragma omp parallel for if(img.pixelCount() > PARALLELIZATION_THRESHOLD)
		for (coord_t n = 0; n < img.pixelCount(); n++)
		{
			if (eq && img(n) != v)
				eq = false;
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
	Determines if pixel values in the two images are equal.
	@param a, b Images to compare
	*/
	template<typename pixel_t, typename = std::enable_if_t<std::is_floating_point_v<pixel_t> > > bool equals(const Image<pixel_t>& a, const Image<pixel_t>& b, pixel_t tolerance = NumberUtils<pixel_t>::tolerance())
	{
		a.checkSize(b);

		bool eq = true;
#pragma omp parallel for if(a.pixelCount() > PARALLELIZATION_THRESHOLD)
		for (coord_t n = 0; n < a.pixelCount(); n++)
		{
			if (eq && !NumberUtils<pixel_t>::equals(a(n), b(n), tolerance))
				eq = false;

			// Showing progress info here would induce more processing than is done in the whole loop.
		}

		return eq;
	}

	/**
	Determines if pixel values in the two images are not equal.
	@param a, b Images to compare
	*/
	template<typename pixel_t, typename = std::enable_if_t<std::is_floating_point_v<pixel_t> > > bool differs(const Image<pixel_t>& a, const Image<pixel_t>& b, pixel_t tolerance = NumberUtils<pixel_t>::tolerance())
	{
		return !equals(a, b, tolerance);
	}


	/**
	Finds all unique pixel values in the image and adds them to the set.
	*/
	template<typename pixel_t> void unique(const Image<pixel_t>& img, std::set<pixel_t>& values)
	{
		#pragma omp parallel if(img.pixelCount() > PARALLELIZATION_THRESHOLD)
		{
			std::unordered_set<pixel_t> privateValues;
			#pragma omp for nowait
			for (coord_t n = 0; n < img.pixelCount(); n++)
			{
				privateValues.emplace(img(n));
			}

			#pragma omp critical(unique_reduction)
			{
				for (pixel_t v : privateValues)
					values.emplace(v);
			}
		}
	}

	/**
	This is used to select suitable accumulator type for sum projections.
	*/
	template<class pixel_t> struct sum_intermediate_type {
		using type = typename std::conditional <
			std::is_floating_point_v<pixel_t>,
			double,			// floating point pixels -> double accumulator
			typename std::conditional<std::is_signed_v<pixel_t>,
				int64_t,	// signed integer pixels -> int64 accumulator
				uint64_t	// unsigned integer pixels -> uint64 accumulator
			>::type
		>::type;
	};

	/**
	Calculates sum of all pixels in the image.
	@param img Image to process.
	@return Sum of all pixels.
	*/
	template<typename pixel_t, typename out_t = double> out_t sum(const Image<pixel_t>& img)
	{
		using temp_t = typename sum_intermediate_type<pixel_t>::type;
		temp_t res = temp_t();
		#pragma omp parallel if(img.pixelCount() > PARALLELIZATION_THRESHOLD)
		{
			temp_t res_private = temp_t();
			#pragma omp for nowait
			for (coord_t n = 0; n < img.pixelCount(); n++)
			{
				res_private = NumberUtils<temp_t>::saturatingAdd(res_private, (temp_t)img(n));

				// Showing progress info here would induce more processing than is done in the whole loop.
			}

			#pragma omp critical(sum_reduction)
			{
				res = NumberUtils<temp_t>::saturatingAdd(res, res_private);
			}
		}

		return pixelRound<out_t>(res);

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
	Calculates sum of squares of all pixels in the image.
	@param img Image to process.
	@return Sum of all pixels.
	*/
	template<typename pixel_t, typename out_t = double> out_t squareSum(const Image<pixel_t>& img)
	{
		using temp_t = typename sum_intermediate_type<pixel_t>::type;
		temp_t res = temp_t();
		#pragma omp parallel if(img.pixelCount() > PARALLELIZATION_THRESHOLD)
		{
			temp_t res_private = temp_t();
			#pragma omp for nowait
			for (coord_t n = 0; n < img.pixelCount(); n++)
			{
				temp_t val = (temp_t)img(n);
				res_private = NumberUtils<temp_t>::saturatingAdd(res_private, NumberUtils<temp_t>::saturatingMultiply(val, val));

				// Showing progress info here would induce more processing than is done in the whole loop.
			}

			#pragma omp critical(sum2_reduction)
			{
				//res += res_private;
				res = NumberUtils<temp_t>::saturatingAdd(res, res_private);
			}
		}

		return pixelRound<out_t>(res);
	}

	/**
	Calculates sum of all pixels and sum of squares of all pixels in the image.
	*/
	template<typename pixel_t, typename out_t = double> Vec2<out_t> sumAndSquareSum(const Image<pixel_t>& img)
	{
		using temp_t = typename sum_intermediate_type<pixel_t>::type;

		temp_t total = temp_t();
		temp_t total2 = temp_t();

#pragma omp parallel if(img.pixelCount() > PARALLELIZATION_THRESHOLD)
		{
			temp_t total_private = temp_t();
			temp_t total2_private = temp_t();

#pragma omp for nowait
			for (coord_t n = 0; n < img.pixelCount(); n++)
			{
				temp_t val = (temp_t)img(n);
				total_private = NumberUtils<temp_t>::saturatingAdd(total_private, val);
				total2_private = NumberUtils<temp_t>::saturatingAdd(total2_private, NumberUtils<temp_t>::saturatingMultiply(val, val));

				// Showing progress info here would induce more processing than is done in the whole loop.
			}

#pragma omp critical(meanandstddev_reduction)
			{
				total = NumberUtils<temp_t>::saturatingAdd(total, total_private);
				total2 = NumberUtils<temp_t>::saturatingAdd(total2, total2_private);
			}
		}

		return Vec2<out_t>(pixelRound<out_t>(total), pixelRound<out_t>(total2));
	}


	/*
	Calculates sum of all pixels except those whose value is ignoreValue.
	Counts pixels that contributes to the sum and places the count to count argument.
	@param img Image to process.
	@param ignoreValue Value that should be ignored while calculating sum.
	@param count Count of pixels that do not have ignoreValue.
	@return Sum of all pixels except those whose value is ignoreValue.
	*/
	template<typename pixel_t, typename out_t = double> out_t maskedSum(const Image<pixel_t>& img, pixel_t ignoreValue, out_t& count)
	{
		using temp_t = typename sum_intermediate_type<pixel_t>::type;

		if (!NumberUtils<pixel_t>::isnan(ignoreValue))
		{
			temp_t res = temp_t();
			temp_t tempCount = temp_t();
#pragma omp parallel if(img.pixelCount() > PARALLELIZATION_THRESHOLD)
			{
				temp_t res_private = temp_t();
				temp_t count_private = temp_t();
#pragma omp for nowait
				for (coord_t n = 0; n < img.pixelCount(); n++)
				{
					pixel_t p = img(n);
					if (p != ignoreValue)
					{
						res_private = NumberUtils<temp_t>::saturatingAdd(res_private, (temp_t)p);
						count_private = NumberUtils<temp_t>::saturatingAdd(count_private, (temp_t)1);
					}

					// Showing progress info here would induce more processing than is done in the whole loop.
				}

#pragma omp critical(maskedsum_reduction)
				{
					res = NumberUtils<temp_t>::saturatingAdd(res, res_private);
					tempCount = NumberUtils<temp_t>::saturatingAdd(tempCount, count_private);
				}
			}

			count = pixelRound<out_t>(tempCount);
			return pixelRound<out_t>(res);
		}
		else
		{
			// ignoreValue is NaN, comparison logic must be different
			temp_t res = temp_t();
			temp_t tempCount = temp_t();
#pragma omp parallel if(img.pixelCount() > PARALLELIZATION_THRESHOLD)
			{
				temp_t res_private = temp_t();
				temp_t count_private = temp_t();
#pragma omp for nowait
				for (coord_t n = 0; n < img.pixelCount(); n++)
				{
					pixel_t p = img(n);
					if (!NumberUtils<pixel_t>::isnan(p))
					{
						res_private = NumberUtils<temp_t>::saturatingAdd(res_private, (temp_t)p);
						count_private = NumberUtils<temp_t>::saturatingAdd(count_private, (temp_t)1);
					}

					// Showing progress info here would induce more processing than is done in the whole loop.
				}

#pragma omp critical(maskedsum_reduction)
				{
					res = NumberUtils<temp_t>::saturatingAdd(res, res_private);
					tempCount = NumberUtils<temp_t>::saturatingAdd(tempCount, count_private);
				}
			}

			count = pixelRound<out_t>(tempCount);
			return pixelRound<out_t>(res);
		}
	}


	/*
	Calculates sum and sum of squares of all pixels except those whose value is ignoreValue.
	Counts pixels that contributes to the sum and places the count to count argument.
	@param img Image to process.
	@param ignoreValue Value that should be ignored while calculating sum.
	@param count Count of pixels that do not have ignoreValue.
	@return Sum and sum of squares of all pixels except those whose value is ignoreValue, in Vec2::x and Vec2::y, respectively.
	*/
	template<typename pixel_t, typename out_t = double> Vec2<out_t> maskedSumAndSquareSum(const Image<pixel_t>& img, pixel_t ignoreValue, out_t& count)
	{
		using temp_t = typename sum_intermediate_type<pixel_t>::type;

		if (!NumberUtils<pixel_t>::isnan(ignoreValue))
		{
			temp_t res = temp_t();
			temp_t res2 = temp_t();
			temp_t tempCount = temp_t();
#pragma omp parallel if(img.pixelCount() > PARALLELIZATION_THRESHOLD)
			{
				temp_t res_private = temp_t();
				temp_t res2_private = temp_t();
				temp_t count_private = temp_t();
#pragma omp for nowait
				for (coord_t n = 0; n < img.pixelCount(); n++)
				{
					pixel_t p = img(n);
					if (p != ignoreValue)
					{
						res_private = NumberUtils<temp_t>::saturatingAdd(res_private, (temp_t)p);
						res2_private = NumberUtils<temp_t>::saturatingAdd(res2_private, NumberUtils<temp_t>::saturatingMultiply((temp_t)p, (temp_t)p));
						count_private = NumberUtils<temp_t>::saturatingAdd(count_private, (temp_t)1);
					}

					// Showing progress info here would induce more processing than is done in the whole loop.
				}

#pragma omp critical(maskedsum2_reduction)
				{
					res = NumberUtils<temp_t>::saturatingAdd(res, res_private);
					res2 = NumberUtils<temp_t>::saturatingAdd(res2, res2_private);
					tempCount = NumberUtils<temp_t>::saturatingAdd(tempCount, count_private);
				}
			}

			count = pixelRound<out_t>(tempCount);
			return Vec2<out_t>(pixelRound<out_t>(res), pixelRound<out_t>(res2));
		}
		else
		{
			// ignoreValue is NaN, comparison logic must be different
			temp_t res = temp_t();
			temp_t res2 = temp_t();
			temp_t tempCount = temp_t();
#pragma omp parallel if(img.pixelCount() > PARALLELIZATION_THRESHOLD)
			{
				temp_t res_private = temp_t();
				temp_t res2_private = temp_t();
				temp_t count_private = temp_t();
#pragma omp for nowait
				for (coord_t n = 0; n < img.pixelCount(); n++)
				{
					pixel_t p = img(n);
					if (!NumberUtils<pixel_t>::isnan(p))
					{
						res_private = NumberUtils<temp_t>::saturatingAdd(res_private, (temp_t)p);
						res2_private = NumberUtils<temp_t>::saturatingAdd(res2_private, NumberUtils<temp_t>::saturatingMultiply((temp_t)p, (temp_t)p));
						count_private = NumberUtils<temp_t>::saturatingAdd(count_private, (temp_t)1);
					}

					// Showing progress info here would induce more processing than is done in the whole loop.
				}

#pragma omp critical(maskedsum2_reduction)
				{
					res = NumberUtils<temp_t>::saturatingAdd(res, res_private);
					res2 = NumberUtils<temp_t>::saturatingAdd(res2, res2_private);
					tempCount = NumberUtils<temp_t>::saturatingAdd(tempCount, count_private);
				}
			}

			count = pixelRound<out_t>(tempCount);
			return Vec2<out_t>(pixelRound<out_t>(res), pixelRound<out_t>(res2));
		}
	}

	/**
	Calculates minimum of all pixels in the image.
	@param img Image to process.
	@return Minimum of all pixels.
	*/
	template<typename pixel_t, typename result_t = pixel_t> result_t min(const Image<pixel_t>& img)
	{
		pixel_t res = std::numeric_limits<pixel_t>::max();
		#pragma omp parallel if(img.pixelCount() > PARALLELIZATION_THRESHOLD)
		{
			pixel_t res_private = std::numeric_limits<pixel_t>::max();
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

		return pixelRound<result_t>(res);

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
	template<typename pixel_t, typename result_t = pixel_t> result_t max(const Image<pixel_t>& img)
	{
		pixel_t res = std::numeric_limits<pixel_t>::lowest();
		#pragma omp parallel if(img.pixelCount() > PARALLELIZATION_THRESHOLD)
		{
			pixel_t res_private = std::numeric_limits<pixel_t>::lowest();
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

		return pixelRound<result_t>(res);

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
		return NumberUtils<out_t>::saturatingDivide(sum<pixel_t, out_t>(img), (out_t)img.pixelCount());
	}

	/**
	Calculates mean of all pixels in the image, except those that have specified value.
	@param img Image to process.
	@param ignoreValue Value that should not be accounted for when calculating mean.
	@param count At output, count of pixels whose value is not ignoreValue.
	@return Mean of all pixels.
	*/
	template<typename pixel_t, typename out_t = double> out_t maskedMean(const Image<pixel_t>& img, pixel_t ignoreValue, out_t& count)
	{
		out_t s = maskedSum<pixel_t, out_t>(img, ignoreValue, count);
		return NumberUtils<out_t>::saturatingDivide(s, count);
	}

	/**
	Calculates mean of all pixels in the image, except those that have specified value.
	@param img Image to process.
	@param ignoreValue Value that should not be accounted for when calculating mean.
	@return Mean of all pixels.
	*/
	template<typename pixel_t, typename out_t = double> out_t maskedMean(const Image<pixel_t>& img, pixel_t ignoreValue)
	{
		out_t dummy;
		return maskedMean<pixel_t, out_t>(img, ignoreValue, dummy);
	}

	namespace internals
	{
		/**
		Converts sum, square sum, and count to mean and standard deviation.
		*/
		template<typename out_t> Vec2<out_t> sumAndSquareSumToMeanAndStdDev(const Vec2<out_t>& stats, out_t count)
		{
			out_t total = stats[0];
			out_t total2 = stats[1];

			out_t avg = total / count;
			out_t std = sqrt((total2 - (total * total / count)) / (count - 1));

			return Vec2<out_t>(avg, std);
		}
	}

	/**
	Calculates mean and standard deviation of all pixels in the image.
	@param img Image to process.
	@return Vec2 containing mean at Vec2::x and standard deviation at Vec2::y.
	*/
	template<typename pixel_t, typename out_t = double> Vec2<out_t> meanAndStdDev(const Image<pixel_t>& img)
	{
		Vec2<out_t> stats = sumAndSquareSum<pixel_t, out_t>(img);
		return internals::sumAndSquareSumToMeanAndStdDev(stats, (out_t)img.pixelCount());
	}

	/**
	Calculates mean and standard deviation of all pixels in the image, except those whose value equals to ignoreValue.
	@param img Image to process.
	@return Vec2 containing mean at Vec2::x and standard deviation at Vec2::y.
	*/
	template<typename pixel_t, typename out_t = double> Vec2<out_t> maskedMeanAndStdDev(const Image<pixel_t>& img, pixel_t ignoreValue)
	{
		out_t count = 0;
		Vec2<out_t> stats = maskedSumAndSquareSum<pixel_t, out_t>(img, ignoreValue, count);
		return internals::sumAndSquareSumToMeanAndStdDev(stats, count);
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
		projectDimension<pixel_t, out_t, internals::minProjectionOp<pixel_t> >(img, dimension, out, std::numeric_limits<double>::infinity(), showProgressInfo);
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
		projectDimension<pixel_t, out_t, internals::maxProjectionOp<pixel_t> >(img, dimension, out, -std::numeric_limits<double>::infinity(), showProgressInfo);
	}


	/**
	Minimum projection in some coordinate direction.
	Picks output value from second image from the location of the minimum.
	@param img Input image.
	@param val Image where second output value is picked.
	@param dimension Zero-based dimension over which the projection should be made.
	@param out Output image.
	@param outVal Output where picked values will be placed.
	@param showProgressInfo Set to true to show a progress bar.
	*/
	template<typename pixel_t, typename val_t, typename out_t> void min(const Image<pixel_t>& img, const Image<val_t>& valImg, size_t dimension, Image<out_t>& out, Image<val_t>& outVal, bool showProgressInfo = true)
	{
		projectDimension<pixel_t, val_t, out_t, val_t, internals::minProjectionValOp<pixel_t, val_t> >(img, valImg, dimension, out, outVal, std::numeric_limits<double>::infinity(), val_t(), showProgressInfo);
	}

	/**
	Maximum projection in some coordinate direction.
	Picks output value from second image from the location of the maximum.
	@param img Input image.
	@param val Image where second output value is picked.
	@param dimension Zero-based dimension over which the projection should be made.
	@param out Output image.
	@param outVal Output where picked values will be placed.
	@param showProgressInfo Set to true to show a progress bar.
	*/
	template<typename pixel_t, typename val_t, typename out_t> void max(const Image<pixel_t>& img, const Image<val_t>& valImg, size_t dimension, Image<out_t>& out, Image<val_t>& outVal, bool showProgressInfo = true)
	{
		projectDimension<pixel_t, val_t, out_t, val_t, internals::maxProjectionValOp<pixel_t, val_t> >(img, valImg, dimension, out, outVal, -std::numeric_limits<double>::infinity(), val_t(), showProgressInfo);
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
