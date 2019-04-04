#pragma once

/*
The distance transform code is adapted from ITK from file
ITK/Modules/Filtering/DistanceMap/include/itkSignedMaurerDistanceMapImageFilter.hxx
The original file has been licensed under the Apache License, Version 2.0.

The algorithm used here is published in
Calvin R. Maurer, Vijay Raghavan, A Linear Time Algorithm for Computing Exact Euclidean Distance Transforms of Binary Images in Arbitrary Dimensions, IEEE TRANSACTIONS ON PATTERN ANALYSIS AND MACHINE INTELLIGENCE, VOL. 25, NO. 2, FEBRUARY 2003
*/

#include <vector>
#include <algorithm>
#include <omp.h>

#include "image.h"
#include "math/mathutils.h"
#include "pointprocess.h"
#include "utilities.h"


namespace itl2
{
	/*
	Helpers for distance map calculations.
	*/
	namespace internals
	{
		template<typename pixel_t> bool remove(pixel_t d1, pixel_t d2, pixel_t df, pixel_t x1, pixel_t x2, pixel_t xf)
		{
			pixel_t a = x2 - x1;
			pixel_t b = xf - x2;
			pixel_t c = xf - x1;

			pixel_t value = (c * ::abs(d2) - b * ::abs(d1) - a * ::abs(df) - a * b * c);

			return value > 0;
		}

		/*
		Helper for distance map calculation.
		Processes one row in dimension d whose start point is idx.
		Optimized version that does not store nearest object point for each dmap point.
		g and h are temporary images whose size must be output.dimension(d) x 1 x 1
		*/
		template<typename pixel_t> void voronoi(size_t d, math::Vec3c idx, Image<pixel_t>& output, Image<pixel_t>& g, Image<pixel_t>& h)
		{
			coord_t nd = output.dimension(d);

			//Image<float32_t> g(nd);
			//Image<float32_t> h(nd);

			int l = -1;

			for (coord_t i = 0; i < nd; i++)
			{
				idx[d] = i;

				pixel_t di = output(idx);

				pixel_t iw = static_cast<pixel_t>(i);

				if (di < numeric_limits<pixel_t>::max())
				{
					if (l < 1)
					{
						l++;
						g(l) = di;
						h(l) = iw;
					}
					else
					{
						while ((l >= 1) && remove(g(l - 1), g(l), di, h(l - 1), h(l), iw))
						{
							l--;
						}
						l++;
						g(l) = di;
						h(l) = iw;
					}
				}
			}

			if (l == -1)
			{
				return;
			}

			int ns = l;

			l = 0;

			for (coord_t i = 0; i < nd; i++)
			{
				pixel_t iw = static_cast<pixel_t>(i);

				pixel_t d1 = ::abs(g(l)) + (h(l) - iw) * (h(l) - iw);

				while (l < ns)
				{
					// be sure to compute d2 *only* if l < ns
					pixel_t d2 = ::abs(g(l + 1)) + (h(l + 1) - iw) * (h(l + 1) - iw);
					// then compare d1 and d2
					if (d1 <= d2)
					{
						break;
					}
					l++;
					d1 = d2;
				}
				idx[d] = i;

				output(idx) = d1;
			}

		}

		/*
		Helper for distance map calculation.
		Processes one row in dimension d whose start point is idx.
		Fills also nearestObjectPoint image by locations of nearest object point for each dmap point.
		g, P, and h are temporary images whose size must be output.dimension(d) x 1 x 1
		*/
		template<typename pixel_t> void voronoi(size_t d, math::Vec3c idx, Image<pixel_t>& output, Image<math::Vec3c>& nearestObjectPoint, Image<pixel_t>& g, Image<math::Vec3c>& P, Image<pixel_t>& h)
		{
			coord_t nd = output.dimension(d);

			//Image<float32_t> g(nd);
			//Image<Vec3c> P(nd);
			//Image<float32_t> h(nd);

			int l = -1;

			for (coord_t i = 0; i < nd; i++)
			{
				idx[d] = i;

				pixel_t di = output(idx);

				pixel_t iw = static_cast<pixel_t>(i);
				
				if (di < numeric_limits<pixel_t>::max())
				{
					if (l < 1)
					{
						l++;
						g(l) = di;
						h(l) = iw;
						P(l) = nearestObjectPoint(idx);
					}
					else
					{
						while ((l >= 1)	&& remove(g(l - 1), g(l), di, h(l - 1), h(l), iw))
						{
							l--;
						}
						l++;
						g(l) = di;
						h(l) = iw;
						P(l) = nearestObjectPoint(idx);
					}
				}
			}

			if (l == -1)
			{
				return;
			}

			int ns = l;

			l = 0;

			for (coord_t i = 0; i < nd; i++)
			{
				pixel_t iw = static_cast<pixel_t>(i);

				pixel_t d1 = ::abs(g(l)) + (h(l) - iw) * (h(l) - iw);
				math::Vec3c Pc = P(l);

				while (l < ns)
				{
					// be sure to compute d2 *only* if l < ns
					pixel_t d2 = ::abs(g(l + 1)) + (h(l + 1) - iw) * (h(l + 1) - iw);
					// then compare d1 and d2
					if (d1 <= d2)
					{
						break;
					}
					l++;
					d1 = d2;
					Pc = P(l);
				}
				idx[d] = i;

				output(idx) = d1;

				nearestObjectPoint(idx) = Pc;
			}

		}


		template<typename pixel_t>  void processDimension(Image<pixel_t>& output, size_t currentDimension, Image<math::Vec3c>* nearestObjectPoint, bool showProgressInfo = false)
		{
			// Determine count of pixels to process
			math::Vec3c reducedDimensions = output.dimensions();
			reducedDimensions[currentDimension] = 1;
			coord_t rowCount = reducedDimensions.x * reducedDimensions.y * reducedDimensions.z;

			size_t counter = 0;
			#pragma omp parallel if(!omp_in_parallel() && output.pixelCount() > PARALLELIZATION_THRESHOLD)
			{
				// Temporary buffer
				coord_t nd = output.dimension(currentDimension);
				Image<pixel_t> g(nd);
				Image<math::Vec3c> P;
				if (nearestObjectPoint)
					P.ensureSize(nd);
				Image<pixel_t> h(nd);

				// Determine start points of pixel rows and process each row
				#pragma omp for
				for (coord_t n = 0; n < rowCount; n++)
				{
					math::Vec3c start = indexToCoords(n, reducedDimensions);

					// Process the current row
					if (!nearestObjectPoint)
					{
						voronoi(currentDimension, start, output, g, h);
					}
					else
					{
						voronoi(currentDimension, start, output, *nearestObjectPoint, g, P, h);
					}
					

					showThreadProgress(counter, rowCount, showProgressInfo);
				}
			}
		}

	//	inline void processDimension(Image<float32_t>& output, size_t currentDimension, Image<math::Vec3c>* nearestObjectPoint)
	//	{
	//		// compute the number of rows first, so we can setup a progress reporter
	//		std::vector<coord_t> numberOfRows;
	//		for (size_t i = 0; i < output.dimensionality(); i++)
	//		{
	//			numberOfRows.push_back(1);
	//			for (unsigned int d = 0; d < output.dimensionality(); d++)
	//			{
	//				if (d != i)
	//				{
	//					numberOfRows[i] *= output.dimension(d);
	//				}
	//			}
	//		}

	//		// This variable provides the amount by which to divide the dimensionless index in order to get the index for each dimension.
	//		std::vector<size_t> k(output.dimensionality() - 1);
	//		size_t count = 0;
	//		k[count] = 1;
	//		for (size_t d = currentDimension + output.dimensionality() - 1; d > currentDimension + 1; d--)
	//		{
	//			k[count + 1] = k[count] * output.dimension(d % output.dimensionality());
	//			count++;
	//		}
	//		std::reverse(k.begin(), k.end());

	//		coord_t rowCount = numberOfRows[currentDimension];
	//		

	//		//size_t counter = 0;

	//		// Process all rows
	//		#pragma omp parallel if(rowCount > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
	//		{
	//			// Temporary buffers. P is not needed if nearest object (zero-set) point is not determined.
	//			coord_t nd = output.dimension(currentDimension);
	//			Image<float32_t> g(nd);
	//			Image<math::Vec3c> P;
	//			if (nearestObjectPoint)
	//				P.ensureSize(nd);
	//			Image<float32_t> h(nd);

	//			#pragma omp for
	//			for (coord_t n = 0; n < rowCount; n++)
	//			{
	//				// Determine first point of the row that should be processed
	//				size_t index = n;
	//				size_t count = 0;
	//				math::Vec3c idx;
	//				for (size_t d = currentDimension + 1; d < currentDimension + output.dimensionality(); d++)
	//				{
	//					idx[d % output.dimensionality()] = static_cast<coord_t>(static_cast<double>(index) / static_cast<double>(k[count]));

	//					index %= k[count];
	//					count++;
	//				}

	//				// Process the current row
	//				if (!nearestObjectPoint)
	//				{
	//					voronoi(currentDimension, idx, output, g, h);
	//				}
	//				else
	//				{
	//					voronoi(currentDimension, idx, output, *nearestObjectPoint, g, P, h);
	//				}

	//				//showThreadProgress(counter, rowCount);
	//			}
	//		}
	//	}
	}

	/**
	Calculates squared distance transform of image in-place.
	Template argument pixel_t is usually float32_t or uint32_t.
	@param img Input and output image. Pixels belonging to the background must be set to zero, and all other pixels must be set to numeric_limits<pixel_t>::max(). See also prepareDistanceTransform.
	*/
	template<typename pixel_t> void distanceTransform2(Image<pixel_t>& img, Image<math::Vec3c>* nearestObjectPoint = 0)
	{
		if (nearestObjectPoint)
		{
			nearestObjectPoint->ensureSize(img);

			// Init each nearest point to point to itself if the point is in zero-distance set,
			// and to point to nothing if the point is outside the zero-distance set.
			for (coord_t z = 0; z < img.depth(); z++)
			{
				for (coord_t y = 0; y < img.height(); y++)
				{
					for (coord_t x = 0; x < img.width(); x++)
					{
						if(img(x, y, z) == 0)
							(*nearestObjectPoint)(x, y, z) = math::Vec3c(x, y, z);
					}
				}
			}
		}

		for (size_t n = 0; n < img.dimensionality(); n++)
		{
			internals::processDimension(img, n, nearestObjectPoint);
		}

	}

	/**
	Calculates distance transform of image in-place.
	@param img Input and output image. Pixels belonging to the background must be set to zero, and all other pixels must be set to numeric_limits<float32_t>::max(). See also prepareDistanceTransform.
	*/
	inline void distanceTransform(Image<float32_t>& img, Image<math::Vec3c>* nearestObjectPoint = 0)
	{
		distanceTransform2(img, nearestObjectPoint);
		squareRoot(img);
	}

	/**
	Prepares output image for in-place distance transform based on input image of another data type.
	Input and output can be the same image.
	*/
	template<typename pixel_t, typename dmap_t> void prepareDistanceTransform(const Image<pixel_t>& input, Image<dmap_t>& output, pixel_t backgroundValue = 0)
	{
		output.ensureSize(input);

		for (coord_t n = 0; n < output.pixelCount(); n++)
		{
			if (input(n) == backgroundValue)
				output(n) = 0;
			else
				output(n) = numeric_limits<dmap_t>::max();
		}
	}

	/**
	Calculates squared distance transform of input and places it to output.
	Template argument out_t is usually int32_t.
	If nearestObjectPoint is nonzero, feature transformation is placed into that.
	@param input Input image containing the geometry. Pixels whose value equals backgroundValue are assumed to belong to background (zero distance set).
	@param output Will contain the distance transform.
	@param nearestObjectPoint If nonzero, will contain feature transform.
	@param backgroundValue Value of background pixels in input image.
	*/
	template<typename pixel_t, typename out_t> void distanceTransform2(const Image<pixel_t>& input, Image<out_t>& output, Image<math::Vec3c>* nearestObjectPoint = 0, pixel_t backgroundValue = 0)
	{
		prepareDistanceTransform(input, output, backgroundValue);
		distanceTransform2(output, nearestObjectPoint);
	}

	/**
	Calculates distance transform of input and places it to output.
	If nearestObjectPoint is nonzero, feature transformation is placed into that.
	@param input Input image containing the geometry. Pixels whose value equals backgroundValue are assumed to belong to background (zero distance set).
	@param output Will contain the distance transform.
	@param nearestObjectPoint If nonzero, will contain feature transform.
	@param backgroundValue Value of background pixels in input image.
	*/
	template<typename pixel_t> void distanceTransform(const Image<pixel_t>& input, Image<float32_t>& output, Image<math::Vec3c>* nearestObjectPoint = 0, pixel_t backgroundValue = 0)
	{
		prepareDistanceTransform(input, output, backgroundValue);
		distanceTransform(output, nearestObjectPoint);
	}

	namespace tests
	{
		void dmap1();
	}

}
