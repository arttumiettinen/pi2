#pragma once

#include <iostream>
#include <cmath>

#include "image.h"
#include "utilities.h"
#include "math/vec2.h"
#include "math/vec3.h"

using math::Vec3c;
using math::pixelRound;

namespace itl2
{

	/**
	Calculates unweighted histogram of input image.
	@param img Image whose histogram is calculated.
	@param histogram Image that will store the histogram. Bin count is pixel count in this image. Make sure that the pixel data type of the image can hold big enough values; floating point image is recommended.
	@param range Gray value range for the histogram.
	@param edgeSkip This many pixels at the image edge are not considered in the histogram.
	*/
	template<typename pixel_t, typename hist_t> void histogram(const Image<pixel_t>& img, Image<hist_t>& histogram, const math::Vec2d& range, coord_t edgeSkip = 0)
	{

		coord_t dim = histogram.pixelCount();
		size_t counter = 0;

#pragma omp parallel if(img.pixelCount() > PARALLELIZATION_THRESHOLD)
		{
			Image<hist_t> privateHist(histogram.dimensions());
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
							coord_t bin = (coord_t)floor(((pix - range.min) / (range.max - range.min)) * (double)dim);

							// Problem with above expression is that the terms inside floor() may give, e.g. 0.2899999998 for pixel
							// that should go to bin 290-300. The floor makes it end in bin 280-290.
							// Check that the pixel really belongs to the bin determined using above expression, and adjust if necessary.
							// TODO: There is probably some numerically stable algorithm that does not need this check.
							pixel_t binMin = pixelRound<pixel_t>(range.min + (double)bin / (double)dim * (range.max - range.min));
							pixel_t binMax = pixelRound<pixel_t>(range.min + (double)(bin + 1) / (double)dim * (range.max - range.min));
							if (pix < binMin)
								bin--;
							else if (pix >= binMax)
								bin++;

							if (bin < 0)
								bin = 0;
							else if (bin >= dim)
								bin = dim - 1;

							privateHist(bin)++;
						}
					}
				}

				showThreadProgress(counter, img.depth());
			}

#pragma omp critical(sum_reduction)
			{
				for (coord_t n = 0; n < histogram.pixelCount(); n++)
				{
					histogram(n) += privateHist(n);
				}
			}
		}


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

		//	showThreadProgress(counter, img.depth());
		//}
	}


	namespace tests
	{
		void histogram();
	}

	///*
	//* Calculates n-dimensional histogram.
	//* @param channels One data image for each histogram dimension.
	//* @param weight Image containing weight for each pixel. Pass 0 for uniform weight of 1 for each pixel.
	//* @param ranges Data range for each histogram dimension. Data outside the range will be added to the nearest edge bin.
	//* @param histogram The histogram will be placed to this image. Dimensionality of this image must be count of items in channels array. The dimensions
	//* of the image will be counts of bins in each dimension.
	//* @param edgeSkip Pixels that are less than this many pixels from any edge of the image are not added to the histogram.
	//*/
	//template<typename Tin, typename Thist> void histogram(const vector<Image<Tin>*>& channels, Image<Tin>* weight, const vector<Vec2d>& ranges, BasicImage<Thist>& hist, size_t edgeSkip = 0)
	//{
	//	cout << "histogram..." << endl;

	//	if(channels.size() != hist.dimensionality())
	//		throw ITLException("Channel count and histogram dimensionality do not match.");

	//	if(ranges.size() != channels.size())
	//		throw ITLException("Range count is not channel count.");

	//	if(channels.size() <= 0)
	//		throw ITLException("At least one channel must be specified.");

	//	for(size_t n = 0; n < channels.size(); n++)
	//	{
	//		if(channels[n]->dimensions() != channels[0]->dimensions())
	//			throw ITLException("Individual channel dimensions are not equal.");
	//		if(weight)
	//		{
	//			if(channels[n]->dimensions() != weight->dimensions())
	//				throw ITLException("Weight dimensions do not equal channel dimensions.");
	//		}
	//	}

	//	// Create cursor for each channel.
	//	vector<LinearCursor<Tin>* > channelCursors;
	//	channelCursors.reserve(channels.size());
	//	for(size_t n = 0; n < channels.size(); n++)
	//		channelCursors.push_back(channels[n]->createLinearCursor());

	//	// Create cursor for weight, if any.
	//	LinearCursor<Tin>* weightCursor = 0;
	//	if(weight)
	//		weightCursor = weight->createLinearCursor();

	//	// Output image is accessed randomly.

	//	// Build the histogram.
	//	vector<coord_t> binvect = utilities::createParams<coord_t>(0, channels.size());

	//	vector<coord_t> position(channels[0]->dimensionality());
	//	vector<size_t> tempSpace(channels[0]->dimensionality());
	//	size_t k = 0;
	//	size_t dispStep = channels[0]->pixelCount() / 100;
	//	for(; channelCursors[0]->hasNext(); )
	//	{
	//		size_t edgeDist = 0;
	//		if(edgeSkip > 0)
	//		{
	//			channelCursors[0]->getPosition(position, tempSpace);
	//			edgeDist = edgeDistance(position, channels[0]->dimensions());
	//		}

	//		if(edgeDist >= edgeSkip)
	//		{
	//			// Get coordinates
	//			for(size_t n = 0; n < channelCursors.size(); n++)
	//			{
	//				double coord = **channelCursors[n];
	//				Vec2d range = ranges[n];
	//				size_t dim = hist.getDimension(n);

	//				coord_t bin = (coord_t)floor(((coord - range.min) / (range.max - range.min)) * (double)(dim));
	//				if(bin < 0)
	//					bin = 0;
	//				else if(bin >= (coord_t)dim)
	//					bin = (coord_t)dim - 1;

	//				binvect[n] = bin;
	//			}

	//			// Get weight value
	//			Thist w = 1;
	//			if(weightCursor)
	//				w = (Thist)(**weightCursor);

	//			// Add weight to proper bin
	//			Thist c = hist.getPixel(binvect);
	//			hist.setPixel(binvect, c + w);
	//		}

	//		// Proceed cursors
	//		for(size_t n = 0; n < channelCursors.size(); n++)
	//			channelCursors[n]->proceed();
	//		if(weightCursor)
	//			weightCursor->proceed();

	//		k++;
	//		if(k > dispStep)
	//		{
	//			k = 0;
	//			cout << mathutils::round(channelCursors[0]->getProgress() * 100) << " %\r" << flush;
	//		}
	//	}


	//	// Free cursors.
	//	for(size_t n = 0; n < channels.size(); n++)
	//		channels[n]->closeCursor(channelCursors[n]);

	//	if(weight)
	//		weight->closeCursor(weightCursor);
	//	
	//}

	///*
	//* Calculate histogram of an image.
	//* @param img Source image.
	//* @param weight Weight for each pixel or 0 for uniform weight of 1.
	//* @param histogram 1-dimensional image that receives the histogram.
	//*/
	//template<typename Tin, typename Thist> void histogram(Image<Tin>& img, BasicImage<Thist>& hist, Image<Tin>* weight, Vec2d& range, size_t edgeSkip = 0)
	//{
	//	vector<Image<Tin>*> imgv;
	//	imgv.push_back(&img);

	//	vector<Vec2d> ranges;
	//	ranges.push_back(range);

	//	histogram<Tin, Thist>(imgv, weight, ranges, hist, edgeSkip);
	//}

	///*
	//* Calculate histogram of an image.
	//* @param img Source image.
	//* @param weight Weight for each pixel or 0 for uniform weight of 1.
	//* @param histogram 1-dimensional image that receives the histogram.
	//*/
	//template<typename Tin, typename Thist> void histogram(Image<Tin>& img, BasicImage<Thist>& hist, Image<Tin>* weight = 0, size_t edgeSkip = 0)
	//{
	//	Vec2d range = Vec2d(numeric_limits<Tin>::lowest(), numeric_limits<Tin>::max());
	//	histogram<Tin, Thist>(img, hist, weight, range, edgeSkip);
	//}

	//template<typename Tin, typename Thist> void histogram(Image<Tin>** channels, Image<Tin>* weight, Vec2d* ranges, BasicImage<Thist>& hist, size_t channelCount, size_t edgeSkip = 0)
	//{
	//	vector<Image<Tin>*> imgv;
	//	vector<Vec2d> rangesv;

	//	for(size_t n = 0; n < channelCount; n++)
	//	{
	//		imgv.push_back(channels[n]);
	//		rangesv.push_back(ranges[n]);
	//	}

	//	histogram<Tin, Thist>(imgv, weight, rangesv, hist, edgeSkip);
	//}
}
