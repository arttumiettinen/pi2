#pragma once

#include <vector>

#include "command.h"
#include "math/vec2.h"
#include "histogram.h"
#include "distributable.h"
#include "pilibutilities.h"
#include "io/raw.h"
#include "pointprocess.h"
#include "utilities.h"
#include "specialcommands.h"
#include "distributedtempimage.h"

namespace pilib
{
	namespace internals
	{
		template<typename hist_t> void saveHist(std::string fname, const Image<hist_t>& histogram, Distributor::BLOCK_INDEX_ARG_TYPE blockIndex)
		{
			if (fname.length() > 0)
			{
				TimingFlag flag(TimeClass::IO);

				if (blockIndex >= 0)
					fname = fname + "_" + itl2::toString(blockIndex);

				raw::writed(histogram, fname);
			}
		}

		template<typename hist_t, typename raw_t> void combineHistograms(size_t fileCount, const std::string& tempFilename, DistributedImage<hist_t>& outputHist, Image<raw_t>& rawHist)
		{
			rawHist.ensureSize(outputHist.dimensions());
			setValue(rawHist, 0);
			
			raw_t summer = 0;
			
			std::cout << "Combining results..." << std::endl;
			for (size_t n = 0; n < fileCount; n++)
			{
				std::string fname = tempFilename + "_" + itl2::toString(n);
				raw::internals::expandRawFilename(fname);

				std::cout << "Reading " << fname << std::endl;

				
				Image<raw_t> currHist;
				{
					TimingFlag flag(TimeClass::IO);
					raw::read(currHist, fname);
				}

				add(rawHist, currHist);

				{
					TimingFlag flag(TimeClass::IO);
					fs::remove(fname);
				}
				
			}
			
			// Convert data to output format
			Image<hist_t> temp;
			setValue(temp, rawHist);

			// Push the histogram to the distributed image.
			outputHist.setData(temp);
		}

		inline void fillBins(Image<float32_t>& binStarts, double hmin, double hmax, coord_t binCount)
		{
			binStarts.ensureSize(binCount);
			for (coord_t n = 0; n < binCount; n++)
			{
				binStarts(n) = (float32_t)(hmin + (double)n / (double)binCount* (hmax - hmin));
			}
		}
	}

	template<typename pixel_t, typename hist_t> class HistogramCommand : public Command, public Distributable
	{
	protected:
		friend class CommandList;

		HistogramCommand() : Command("hist", "Calculate histogram of an image.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "input image", "Image whose histogram will be calculated."),
				CommandArgument<Image<hist_t> >(ParameterDirection::Out, "histogram", "Image that will contain the histogram on output. The number of pixels in this image will be set to match bin count. Please use pixel data type that can hold large enough values, e.g. float32 or uint64."),
				CommandArgument<Image<float32_t> >(ParameterDirection::Out, "bin starts", "Image that will contain the histogram bin start coordinates when the function finishes."),
				CommandArgument<double>(ParameterDirection::In, "histogram minimum", "Minimum value to be included in the histogram. Values less than minimum are included in the first bin.", 0),
				CommandArgument<double>(ParameterDirection::In, "histogram maximum", "The end of the last bin. Values greater than maximum are included in the last bin.", (double)std::numeric_limits<pixel_t>::max()),
				CommandArgument<coord_t>(ParameterDirection::In, "bin count", "Count of bins. The output image will be resized to contain this many pixels.", 256),
				CommandArgument<std::string>(ParameterDirection::In, "output file", "Name of file where the histogram data is to be saved. Specify empty string to disable saving. This argument is used internally in distributed processing, but can be used to (conveniently?) save the histogram in .raw format.", ""),
				CommandArgument<Distributor::BLOCK_INDEX_ARG_TYPE>(ParameterDirection::In, Distributor::BLOCK_INDEX_ARG_NAME, "Index of image block that we are currently processing. This argument is used internally in distributed processing and should normally be set to negative value. If positive, this number is appended to output file name.", -1)
			})
		{
		}

	public:

		virtual void run(std::vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& in = *pop<Image<pixel_t>* >(args);
			Image<hist_t>& hist = *pop<Image<hist_t>* >(args);
			Image<float32_t>& bins = *pop<Image<float32_t>*>(args);
			double hmin = pop<double>(args);
			double hmax = pop<double>(args);
			coord_t bincount = pop<coord_t>(args);
			std::string fname = pop<std::string>(args);
			Distributor::BLOCK_INDEX_ARG_TYPE blockIndex = pop<Distributor::BLOCK_INDEX_ARG_TYPE>(args);

			hist.ensureSize(bincount);
			histogram(in, hist, Vec2d(hmin, hmax));

			internals::saveHist(fname, hist, blockIndex);

			internals::fillBins(bins, hmin, hmax, bincount);

			// This saves in .csv format
			//if (fname.length() > 0)
			//{
			//	ofstream out(fname);
			//	out << "Bin start [pixel]\tCount" << std::endl;
			//	for (coord_t n = 0; n < hist.pixelCount(); n++)
			//	{
			//		pixel_t binStart = pixelRound<pixel_t>(hmin + (double)n / (double)hist.pixelCount() * (hmax - hmin));
			//		out << binStart << "\t" << hist(n) << std::endl;
			//	}
			//}
		}

		using Distributable::runDistributed;

		virtual std::vector<std::string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& in = *pop<DistributedImage<pixel_t>* >(args);
			DistributedImage<hist_t>& hist = *pop<DistributedImage<hist_t>* >(args);
			DistributedImage<float32_t>& bins = *pop<DistributedImage<float32_t>*>(args);
			double hmin = pop<double>(args);
			double hmax = pop<double>(args);
			coord_t bincount = pop<coord_t>(args);
			std::string fname = pop<std::string>(args);
			Distributor::BLOCK_INDEX_ARG_TYPE blockIndex = pop<Distributor::BLOCK_INDEX_ARG_TYPE>(args);

			hist.ensureSize(bincount);
			bins.ensureSize(bincount);


			std::string tempFilename = createTempFilename("histogram");

			//std::vector<ParamVariant> args2 = { &in, &hist, hmin, hmax, bincount, tempFilename, Distributor::BLOCK_INDEX_ARG_TYPE() };
			//std::vector<std::string> output = distributor.distribute(this, args2);


			// Run the block histograms using uint64_t data type
			DistributedTempImage<uint64_t> dummy(distributor, "histogram_dummy", Vec3c(1, 1, 1), DistributedImageStorageType::Raw);
			std::vector<ParamVariant> args2 = { &in, &dummy.get(), &bins, hmin, hmax, bincount, tempFilename, Distributor::BLOCK_INDEX_ARG_TYPE() };

			std::vector<std::string> output = distributor.distribute(&CommandList::get<HistogramCommand<pixel_t, uint64_t> >(), args2);

			// Combine everything
			Image<uint64_t> temp(hist.dimensions());
			internals::combineHistograms(output.size(), tempFilename, hist, temp);

			internals::saveHist(fname, temp, blockIndex);

			// Create bins
			Image<float32_t> tempBins(bincount);
			internals::fillBins(tempBins, hmin, hmax, bincount);
			bins.setData(tempBins);

			return std::vector<std::string>();
		}

		virtual void getCorrespondingBlock(const std::vector<ParamVariant>& args, size_t argIndex, Vec3c& readStart, Vec3c& readSize, Vec3c& writeFilePos, Vec3c& writeImPos, Vec3c& writeSize) const override
		{
			if (argIndex == 1)
			{
				// Always load output image and the bins image, but do not write.
				DistributedImage<hist_t>& hist = *std::get<DistributedImage<hist_t>* >(args[argIndex]);
				readStart = Vec3c(0, 0, 0);
				readSize = hist.dimensions();
				writeFilePos = Vec3c(0, 0, 0);
				writeImPos = Vec3c(0, 0, 0);
				writeSize = Vec3c(0, 0, 0);
			}
			else if(argIndex == 2)
			{
				// Always load output image and the bins image, but do not write.
				DistributedImage<float32_t>& hist = *std::get<DistributedImage<float32_t>* >(args[argIndex]);
				readStart = Vec3c(0, 0, 0);
				readSize = hist.dimensions();
				writeFilePos = Vec3c(0, 0, 0);
				writeImPos = Vec3c(0, 0, 0);
				writeSize = Vec3c(0, 0, 0);
			}
		}

		virtual size_t getRefIndex(const std::vector<ParamVariant>& args) const override
		{
			// Input image is the reference image.
			return 0;
		}
	};


	template<typename pixel1_t, typename pixel2_t> class Histogram2Command : public Command, public Distributable
	{
	protected:
		friend class CommandList;

		Histogram2Command() : Command("hist2", "Calculate bivariate histogram of two images. In the bivariate histogram position (i, j) counts number of locations where the value of the first input image is i and the value of the second input image is j (assuming minimum = 0, maximum = data type max, and bin size = 1). If image sizes are different, only the region available in both images is included in the histogram. In this case, a warning message is shown.",
			{
				CommandArgument<Image<pixel1_t> >(ParameterDirection::In, "first input image", "First image whose histogram will be calculated."),
				CommandArgument<double>(ParameterDirection::In, "min 1", "Minimum value of first image to be included in the histogram.", 0),
				CommandArgument<double>(ParameterDirection::In, "max 1", "The end of the last bin. This is the smallest value above minimum that is not included in the histogram.", (double)std::numeric_limits<pixel1_t>::max()),
				CommandArgument<coord_t>(ParameterDirection::In, "bin count 1", "Count of bins. The first dimension of the output image will be resized to contain this many pixels.", 256),

				CommandArgument<Image<pixel2_t> >(ParameterDirection::In, "second input image", "Second image whose histogram will be calculated."),
				CommandArgument<double>(ParameterDirection::In, "min 2", "Minimum value of second image to be included in the histogram.", 0),
				CommandArgument<double>(ParameterDirection::In, "max 2", "The end of the last bin. This is the smallest value above minimum that is not included in the histogram.", (double)std::numeric_limits<pixel2_t>::max()),
				CommandArgument<coord_t>(ParameterDirection::In, "bin count 2", "Count of bins. The second dimension of the output image will be resized to contain this many pixels.", 256),

				CommandArgument<Image<float32_t> >(ParameterDirection::Out, "histogram", "Image that will contain the histogram on output."),
				CommandArgument<Image<float32_t> >(ParameterDirection::Out, "bin starts 1", "Image that will contain the histogram bin start coordinates in the first dimension when the function finishes."),
				CommandArgument<Image<float32_t> >(ParameterDirection::Out, "bin starts 2", "Image that will contain the histogram bin start coordinates in the second dimension when the function finishes."),

				CommandArgument<std::string>(ParameterDirection::In, "output file", "Name of file where the histogram data is to be saved. Specify empty string to disable saving. This argument is used internally in distributed processing, but can be used to (conveniently?) save the histogram in .raw format.", ""),
				CommandArgument<Distributor::BLOCK_INDEX_ARG_TYPE>(ParameterDirection::In, Distributor::BLOCK_INDEX_ARG_NAME, "Index of image block that we are currently processing. This argument is used internally in distributed processing and should normally be set to negative value. If positive, this number is appended to output file name.", -1)
			})
		{
		}

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			Image<pixel1_t>& in1 = *pop<Image<pixel1_t>* >(args);
			double hmin1 = pop<double>(args);
			double hmax1 = pop<double>(args);
			coord_t bincount1 = pop<coord_t>(args);

			Image<pixel2_t>& in2 = *pop<Image<pixel2_t>* >(args);
			double hmin2 = pop<double>(args);
			double hmax2 = pop<double>(args);
			coord_t bincount2 = pop<coord_t>(args);

			Image<float32_t>& hist = *pop<Image<float32_t>* >(args);
			Image<float32_t>& bins1 = *pop<Image<float32_t>* >(args);
			Image<float32_t>& bins2 = *pop<Image<float32_t>* >(args);

			std::string fname = pop<std::string>(args);
			Distributor::BLOCK_INDEX_ARG_TYPE blockIndex = pop<Distributor::BLOCK_INDEX_ARG_TYPE>(args);

			Image<uint64_t> rawHist;
			itl2::multiHistogram(rawHist, 0, ImageAndRange(in1, Vec2d(hmin1, hmax1), bincount1), ImageAndRange(in2, Vec2d(hmin2, hmax2), bincount2));

			internals::saveHist(fname, rawHist, blockIndex);

			internals::fillBins(bins1, hmin1, hmax1, bincount1);
			internals::fillBins(bins2, hmin2, hmax2, bincount2);

			setValue(hist, rawHist);
		}

		virtual std::vector<std::string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel1_t>& in1 = *pop<DistributedImage<pixel1_t>* >(args);
			double hmin1 = pop<double>(args);
			double hmax1 = pop<double>(args);
			coord_t bincount1 = pop<coord_t>(args);

			DistributedImage<pixel2_t>& in2 = *pop<DistributedImage<pixel2_t>* >(args);
			double hmin2 = pop<double>(args);
			double hmax2 = pop<double>(args);
			coord_t bincount2 = pop<coord_t>(args);

			DistributedImage<float32_t>& hist = *pop<DistributedImage<float32_t>* >(args);
			DistributedImage<float32_t>& bins1 = *pop<DistributedImage<float32_t>* >(args);
			DistributedImage<float32_t>& bins2 = *pop<DistributedImage<float32_t>* >(args);

			std::string fname = pop<std::string>(args);
			Distributor::BLOCK_INDEX_ARG_TYPE blockIndex = pop<Distributor::BLOCK_INDEX_ARG_TYPE>(args);


			in1.checkSize(in2.dimensions());

			hist.ensureSize(bincount1, bincount2);
			bins1.ensureSize(bincount1);
			bins2.ensureSize(bincount2);

			// This does not work reliably. It could work reliably if reference dimensions were
			// min(in1.dimensions(), in2.dimensions()), but now that is not supported as we have
			// reference images and not reference dimensions.
			//if (in1.dimensions() != in2.dimensions())
			//	std::cout << "Warning: Image sizes are not equal. The histogram will only include region [0, 0, 0] - " << (componentwiseMin(in1.dimensions(), in2.dimensions())) << "." << std::endl;

			//hist.ensureSize(bincount1, bincount2);


			//std::string tempFilename = createTempFilename("histogram2");

			//std::vector<ParamVariant> args2 = { &in1, hmin1, hmax1, bincount1, &in2, hmin2, hmax2, bincount2, &rawHist, tempFilename, Distributor::BLOCK_INDEX_ARG_TYPE() };
			//std::vector<std::string> output = distributor.distribute(this, args2);

			//hist.ensureSize(bincount1, bincount2); // ensure size again as passing zeroes to writeSize in getCorrespondingBlock will set the size of the image to zero.

			//Image<uint64_t> localHist;
			//internals::combineHistograms(output.size(), tempFilename, hist, localHist);

			//internals::saveHist(fname, localHist, blockIndex);

			
			std::string tempFilename = createTempFilename("histogram2");

			// Calculate block histograms
			std::vector<ParamVariant> args2 = { &in1, hmin1, hmax1, bincount1, &in2, hmin2, hmax2, bincount2, &hist, &bins1, &bins2, tempFilename, Distributor::BLOCK_INDEX_ARG_TYPE() };
			std::vector<std::string> output = distributor.distribute(&CommandList::get<Histogram2Command<pixel1_t, pixel2_t> >(), args2);

			// Combine everything
			// ensure size again as passing zeroes to writeSize in getCorrespondingBlock will set the size of the image to zero.
			hist.ensureSize(bincount1, bincount2);
			Image<uint64_t> temp(hist.dimensions());
			internals::combineHistograms(output.size(), tempFilename, hist, temp);

			internals::saveHist(fname, temp, blockIndex);

			// Create bins
			Image<float32_t> tempBins(bincount1);
			internals::fillBins(tempBins, hmin1, hmax1, bincount1);
			bins1.setData(tempBins);

			tempBins.ensureSize(bincount2);
			internals::fillBins(tempBins, hmin2, hmax2, bincount2);
			bins2.setData(tempBins);

			return std::vector<std::string>();
		}

		virtual void getCorrespondingBlock(const std::vector<ParamVariant>& args, size_t argIndex, Vec3c& readStart, Vec3c& readSize, Vec3c& writeFilePos, Vec3c& writeImPos, Vec3c& writeSize) const override
		{
			if (argIndex == 8)
			{
				// Always load output image, but do not write.
				DistributedImage<float32_t>& hist = *std::get<DistributedImage<float32_t>* >(args[argIndex]);
				readStart = Vec3c(0, 0, 0);
				readSize = hist.dimensions();
				writeFilePos = Vec3c(0, 0, 0);
				writeImPos = Vec3c(0, 0, 0);
				writeSize = Vec3c(0, 0, 0);
			}
			else if (argIndex == 9 || argIndex == 10)
			{
				// Always load bins image, but do not write.
				DistributedImage<float32_t>& img = *std::get<DistributedImage<float32_t>* >(args[argIndex]);
				readStart = Vec3c(0, 0, 0);
				readSize = img.dimensions();
				writeFilePos = Vec3c(0, 0, 0);
				writeImPos = Vec3c(0, 0, 0);
				writeSize = Vec3c(0, 0, 0);
			}
		}

		virtual size_t getRefIndex(const std::vector<ParamVariant>& args) const override
		{
			// Input image is the reference image.
			return 0;
		}
	};


	template<typename pixel_t> class WeightedHistogramCommand : public Command, public Distributable
	{
	protected:
		friend class CommandList;

		WeightedHistogramCommand() : Command("whist", "Calculate weighted histogram of an image.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "input image", "Image whose histogram will be calculated."),
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "weight", "Weight image. This image must contain weight value for each pixel in the input images"),
				CommandArgument<Image<float32_t> >(ParameterDirection::Out, "histogram", "Image that will contain the histogram on output. The number of pixels in this image will be set to match bin count."),
				CommandArgument<Image<float32_t> >(ParameterDirection::Out, "bin starts", "Image that will contain the histogram bin start coordinates when the function finishes."),
				CommandArgument<double>(ParameterDirection::In, "histogram minimum", "Minimum value to be included in the histogram. Values less than minimum are included in the first bin.", 0),
				CommandArgument<double>(ParameterDirection::In, "histogram maximum", "The end of the last bin. Values greater than maximum are included in the last bin.", (double)std::numeric_limits<pixel_t>::max()),
				CommandArgument<coord_t>(ParameterDirection::In, "bin count", "Count of bins. The output image will be resized to contain this many pixels.", 256),
				CommandArgument<std::string>(ParameterDirection::In, "output file", "Name of file where the histogram data is to be saved. Specify empty string to disable saving. This argument is used internally in distributed processing, but can be used to (conveniently?) save the histogram in .raw format.", ""),
				CommandArgument<Distributor::BLOCK_INDEX_ARG_TYPE>(ParameterDirection::In, Distributor::BLOCK_INDEX_ARG_NAME, "Index of image block that we are currently processing. This argument is used internally in distributed processing and should normally be set to negative value. If positive, this number is appended to output file name.", -1)
			})
		{
		}

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& in = *pop<Image<pixel_t>* >(args);
			Image<float32_t>& weight = *pop<Image<float32_t>* >(args);
			Image<float32_t>& hist = *pop<Image<float32_t>* >(args);
			Image<float32_t>& bins = *pop<Image<float32_t>*>(args);
			double hmin = pop<double>(args);
			double hmax = pop<double>(args);
			coord_t bincount = pop<coord_t>(args);
			std::string fname = pop<std::string>(args);
			Distributor::BLOCK_INDEX_ARG_TYPE blockIndex = pop<Distributor::BLOCK_INDEX_ARG_TYPE>(args);

			using sum_t = typename histogram_intermediate_type<float32_t, float32_t>::type;
			Image<sum_t> rawHist;

			weight.checkSize(in);
			rawHist.ensureSize(bincount);
			itl2::histogram(in, rawHist, Vec2d(hmin, hmax), 0, &weight);

			internals::saveHist(fname, rawHist, blockIndex);

			internals::fillBins(bins, hmin, hmax, bincount);

			setValue(hist, rawHist);
		}

		virtual std::vector<std::string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& in = *pop<DistributedImage<pixel_t>* >(args);
			DistributedImage<float32_t>& weight = *pop<DistributedImage<float32_t>*>(args);
			DistributedImage<float32_t>& hist = *pop<DistributedImage<float32_t>* >(args);
			DistributedImage<float32_t>& bins = *pop<DistributedImage<float32_t>*>(args);
			double hmin = pop<double>(args);
			double hmax = pop<double>(args);
			coord_t bincount = pop<coord_t>(args);
			std::string fname = pop<std::string>(args);
			Distributor::BLOCK_INDEX_ARG_TYPE blockIndex = pop<Distributor::BLOCK_INDEX_ARG_TYPE>(args);

			weight.checkSize(in.dimensions());
			hist.ensureSize(bincount);
			bins.ensureSize(bincount);


			std::string tempFilename = createTempFilename("histogram");



			// Run the block histograms using suitable data type
			//using sum_t = typename histogram_intermediate_type<float32_t, float32_t>::type;
			//DistributedTempImage<sum_t> dummy(distributor, "histogram_dummy", Vec3c(1, 1, 1));
			std::vector<ParamVariant> args2 = { &in, &weight, &hist, &bins, hmin, hmax, bincount, tempFilename, Distributor::BLOCK_INDEX_ARG_TYPE() };

			std::vector<std::string> output = distributor.distribute(&CommandList::get<WeightedHistogramCommand<pixel_t> >(), args2);

			// Combine everything
			hist.ensureSize(bincount);
			using sum_t = typename histogram_intermediate_type<float32_t, float32_t>::type;
			Image<sum_t> temp(hist.dimensions());
			internals::combineHistograms(output.size(), tempFilename, hist, temp);

			internals::saveHist(fname, temp, blockIndex);

			// Create bins
			Image<float32_t> tempBins(bincount);
			internals::fillBins(tempBins, hmin, hmax, bincount);
			bins.setData(tempBins);

			return std::vector<std::string>();
		}

		virtual void getCorrespondingBlock(const std::vector<ParamVariant>& args, size_t argIndex, Vec3c& readStart, Vec3c& readSize, Vec3c& writeFilePos, Vec3c& writeImPos, Vec3c& writeSize) const override
		{
			if (argIndex == 2)
			{
				// Always load output image and the bins image, but do not write.
				DistributedImage<float32_t>& hist = *std::get<DistributedImage<float32_t>* >(args[argIndex]);
				readStart = Vec3c(0, 0, 0);
				readSize = hist.dimensions();
				writeFilePos = Vec3c(0, 0, 0);
				writeImPos = Vec3c(0, 0, 0);
				writeSize = Vec3c(0, 0, 0);
			}
			else if (argIndex == 3)
			{
				// Always load output image and the bins image, but do not write.
				DistributedImage<float32_t>& hist = *std::get<DistributedImage<float32_t>* >(args[argIndex]);
				readStart = Vec3c(0, 0, 0);
				readSize = hist.dimensions();
				writeFilePos = Vec3c(0, 0, 0);
				writeImPos = Vec3c(0, 0, 0);
				writeSize = Vec3c(0, 0, 0);
			}
		}

		virtual size_t getRefIndex(const std::vector<ParamVariant>& args) const override
		{
			// Input image is the reference image.
			return 0;
		}

		virtual size_t getDistributionDirection2(const std::vector<ParamVariant>& args) const override
		{
			return 1;
		}
	};


	template<typename pixel1_t, typename pixel2_t> class WeightedHistogram2Command : public Command, public Distributable
	{
	protected:
		friend class CommandList;

		WeightedHistogram2Command() : Command("whist2", "Calculate weighted bivariate histogram of two images. In the bivariate histogram position $(i, j)$ counts total weight of locations where the value of the first input image is $i$ and the value of the second input image is $j$ (assuming minimum = 0, maximum = data type max, and bin size = 1). If image sizes are different, only the region available in both images is included in the histogram. In this case, a warning message is shown.",
			{
				CommandArgument<Image<pixel1_t> >(ParameterDirection::In, "first input image", "First image whose histogram will be calculated."),
				CommandArgument<double>(ParameterDirection::In, "min 1", "Minimum value of first image to be included in the histogram.", 0),
				CommandArgument<double>(ParameterDirection::In, "max 1", "The end of the last bin. This is the smallest value above minimum that is not included in the histogram.", (double)std::numeric_limits<pixel1_t>::max()),
				CommandArgument<coord_t>(ParameterDirection::In, "bin count 1", "Count of bins. The first dimension of the output image will be resized to contain this many pixels.", 256),

				CommandArgument<Image<pixel2_t> >(ParameterDirection::In, "second input image", "Second image whose histogram will be calculated."),
				CommandArgument<double>(ParameterDirection::In, "min 2", "Minimum value of second image to be included in the histogram.", 0),
				CommandArgument<double>(ParameterDirection::In, "max 2", "The end of the last bin. This is the smallest value above minimum that is not included in the histogram.", (double)std::numeric_limits<pixel2_t>::max()),
				CommandArgument<coord_t>(ParameterDirection::In, "bin count 2", "Count of bins. The second dimension of the output image will be resized to contain this many pixels.", 256),

				CommandArgument<Image<float32_t> >(ParameterDirection::In, "weight", "Weight image. This image must contain weight value for each pixel in the input images"),

				CommandArgument<Image<float32_t> >(ParameterDirection::Out, "histogram", "Image that will contain the histogram on output."),
				CommandArgument<Image<float32_t> >(ParameterDirection::Out, "bin starts 1", "Image that will contain the histogram bin start coordinates in the first dimension when the function finishes."),
				CommandArgument<Image<float32_t> >(ParameterDirection::Out, "bin starts 2", "Image that will contain the histogram bin start coordinates in the second dimension when the function finishes."),

				CommandArgument<std::string>(ParameterDirection::In, "output file", "Name of file where the histogram data is to be saved. Specify empty string to disable saving. This argument is used internally in distributed processing, but can be used to (conveniently?) save the histogram in .raw format.", ""),
				CommandArgument<Distributor::BLOCK_INDEX_ARG_TYPE>(ParameterDirection::In, Distributor::BLOCK_INDEX_ARG_NAME, "Index of image block that we are currently processing. This argument is used internally in distributed processing and should normally be set to negative value. If positive, this number is appended to output file name.", -1)
			})
		{
		}

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			Image<pixel1_t>& in1 = *pop<Image<pixel1_t>* >(args);
			double hmin1 = pop<double>(args);
			double hmax1 = pop<double>(args);
			coord_t bincount1 = pop<coord_t>(args);

			Image<pixel2_t>& in2 = *pop<Image<pixel2_t>* >(args);
			double hmin2 = pop<double>(args);
			double hmax2 = pop<double>(args);
			coord_t bincount2 = pop<coord_t>(args);

			Image<float32_t>& weight = *pop<Image<float32_t>* >(args);

			Image<float32_t>& hist = *pop<Image<float32_t>* >(args);
			Image<float32_t>& bins1 = *pop<Image<float32_t>* >(args);
			Image<float32_t>& bins2 = *pop<Image<float32_t>* >(args);

			std::string fname = pop<std::string>(args);
			Distributor::BLOCK_INDEX_ARG_TYPE blockIndex = pop<Distributor::BLOCK_INDEX_ARG_TYPE>(args);

			using sum_t = typename histogram_intermediate_type<float32_t, float32_t>::type;
			Image<sum_t> rawHist;
			itl2::multiHistogram(rawHist, 0, &weight, ImageAndRange(in1, Vec2d(hmin1, hmax1), bincount1), ImageAndRange(in2, Vec2d(hmin2, hmax2), bincount2));

			internals::saveHist(fname, rawHist, blockIndex);

			internals::fillBins(bins1, hmin1, hmax1, bincount1);
			internals::fillBins(bins2, hmin2, hmax2, bincount2);

			setValue(hist, rawHist);
		}

		virtual std::vector<std::string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel1_t>& in1 = *pop<DistributedImage<pixel1_t>* >(args);
			double hmin1 = pop<double>(args);
			double hmax1 = pop<double>(args);
			coord_t bincount1 = pop<coord_t>(args);

			DistributedImage<pixel2_t>& in2 = *pop<DistributedImage<pixel2_t>* >(args);
			double hmin2 = pop<double>(args);
			double hmax2 = pop<double>(args);
			coord_t bincount2 = pop<coord_t>(args);

			DistributedImage<float32_t>& weight = *pop<DistributedImage<float32_t>* >(args);

			DistributedImage<float32_t>& hist = *pop<DistributedImage<float32_t>* >(args);
			DistributedImage<float32_t>& bins1 = *pop<DistributedImage<float32_t>* >(args);
			DistributedImage<float32_t>& bins2 = *pop<DistributedImage<float32_t>* >(args);

			std::string fname = pop<std::string>(args);
			Distributor::BLOCK_INDEX_ARG_TYPE blockIndex = pop<Distributor::BLOCK_INDEX_ARG_TYPE>(args);


			in1.checkSize(in2.dimensions());
			in1.checkSize(weight.dimensions());
			hist.ensureSize(bincount1, bincount2);
			bins1.ensureSize(bincount1);
			bins2.ensureSize(bincount2);

			std::string tempFilename = createTempFilename("whistogram2");

			// Calculate block histograms
			std::vector<ParamVariant> args2 = { &in1, hmin1, hmax1, bincount1, &in2, hmin2, hmax2, bincount2, &weight, &hist, &bins1, &bins2, tempFilename, Distributor::BLOCK_INDEX_ARG_TYPE() };
			std::vector<std::string> output = distributor.distribute(&CommandList::get<WeightedHistogram2Command<pixel1_t, pixel2_t> >(), args2);

			// Combine everything
			// ensure size again as passing zeroes to writeSize in getCorrespondingBlock will set the size of the image to zero.
			hist.ensureSize(bincount1, bincount2);
			using sum_t = typename histogram_intermediate_type<float32_t, float32_t>::type;
			Image<sum_t> temp(hist.dimensions());
			internals::combineHistograms(output.size(), tempFilename, hist, temp);

			internals::saveHist(fname, temp, blockIndex);

			// Create bins
			Image<float32_t> tempBins(bincount1);
			internals::fillBins(tempBins, hmin1, hmax1, bincount1);
			bins1.setData(tempBins);

			tempBins.ensureSize(bincount2);
			internals::fillBins(tempBins, hmin2, hmax2, bincount2);
			bins2.setData(tempBins);

			return std::vector<std::string>();
		}

		virtual void getCorrespondingBlock(const std::vector<ParamVariant>& args, size_t argIndex, Vec3c& readStart, Vec3c& readSize, Vec3c& writeFilePos, Vec3c& writeImPos, Vec3c& writeSize) const override
		{
			if (argIndex == 9)
			{
				// Always load output image, but do not write.
				DistributedImage<float32_t>& hist = *std::get<DistributedImage<float32_t>* >(args[argIndex]);
				readStart = Vec3c(0, 0, 0);
				readSize = hist.dimensions();
				writeFilePos = Vec3c(0, 0, 0);
				writeImPos = Vec3c(0, 0, 0);
				writeSize = Vec3c(0, 0, 0);
			}
			else if (argIndex == 10 || argIndex == 11)
			{
				// Always load bins image, but do not write.
				DistributedImage<float32_t>& img = *std::get<DistributedImage<float32_t>* >(args[argIndex]);
				readStart = Vec3c(0, 0, 0);
				readSize = img.dimensions();
				writeFilePos = Vec3c(0, 0, 0);
				writeImPos = Vec3c(0, 0, 0);
				writeSize = Vec3c(0, 0, 0);
			}
		}

		virtual size_t getRefIndex(const std::vector<ParamVariant>& args) const override
		{
			// Input image is the reference image.
			return 0;
		}

		virtual size_t getDistributionDirection2(const std::vector<ParamVariant>& args) const override
		{
			return 1;
		}
	};

}
