#pragma once

#include "autothreshold.h"
#include "conversions.h"

#include "commandsbase.h"
#include "distributedtempimage.h"

#include "pointprocesscommands.h"
#include "histogramcommands.h"

#include <string>


namespace pilib
{

	inline std::string autoThresholdMethodsHelp()
	{
		return
			"**Otsu**\n"
			"\n"
			"Otsu's thresholding algorithm.\n"
			"\n"
			"The threshold value is calculated from the histogram of the image. "
			"The first and the second arguments are values corresponding to the minimum and maximum of the histogram. "
			"Out-of-range values will be placed to the first and the last bins, respectively. "
			"The third argument is the count of bins in the histogram.\n"
			"\n"
			"Described in N. Otsu, A Threshold Selection Method from Gray-Level Histograms, IEEE Transactions on Systems, Man, and Cybernetics 9(1), 1979.\n"
			"\n"
			"Original C++ code by Jordan Bevik, ported to ImageJ plugin by Gabriel Landini, and finally ported to pi2.\n"
			"\n"
			"\n"
			"**Huang**\n"
			"\n"
			"Huang's fuzzy thresholding method.\n"
			"\n"
			"Uses Shannon's entropy function (one can also use Yager's entropy function).\n"
			"\n"
			"The threshold value is calculated from the histogram of the image. "
			"The first and the second arguments are values corresponding to the minimum and maximum of the histogram. "
			"Out-of-range values will be placed to the first and the last bins, respectively. "
			"The third argument is the count of bins in the histogram.\n"
			"\n"
			"Described in "
			"L.-K. Huang and M.-J.J. Wang, Image Thresholding by Minimizing the Measures of Fuzziness, Pattern Recognition 28(1), 1995.\n"
			"\n"
			"Original code by M. Emre Celebi, ported to ImageJ plugin by G. Landini from E. Celebi's fourier_0.8 routines, then ported to pi2.\n"
			"\n"
			"\n"
			"**Intermodes**\n"
			"\n"
			"Assumes a bimodal histogram. The histogram needs is iteratively smoothed using a "
			"running average of size 3 until there are only two local maxima at j and k. "
			"The threshold t is (j+k)/2.\n"
			"\n"
			"Images with histograms having extremely unequal peaks or a broad and "
			"flat valleys are unsuitable for this method.\n"
			"\n"
			"The threshold value is calculated from the histogram of the image. "
			"The first and the second arguments are values corresponding to the minimum and maximum of the histogram. "
			"Out-of-range values will be placed to the first and the last bins, respectively. "
			"The third argument is the count of bins in the histogram.\n"
			"\n"
			"Described in  "
			"J. M. S. Prewitt and M. L. Mendelsohn, The analysis of cell images, Annals of the New York Academy of Sciences 128, 1966.\n"
			"\n"
			"Ported to ImageJ plugin by Gabriel Landini from Antti Niemisto's Matlab code (GPL), and then to pi2.\n"
			"\n"
			"\n"
			"\n"
			"**IsoData**\n"
			"\n"
			"Iterative procedure based on the isodata algorithm described in "
			"T.W. Ridler and S. Calvard, Picture thresholding using an iterative selection method, IEEE Transactions on System, Man and Cybernetics, SMC-8, 1978.\n"
			"\n"
			"The procedure divides the image into objects and background by taking an initial threshold, "
			"then the averages of the pixels at or below the threshold and pixels above are computed. "
			"The averages of those two values are computed, the threshold is incremented and the "
			"process is repeated until the threshold is larger than the composite average. That is, "
			"threshold = (average background + average objects) / 2\n"
			"\n"
			"The threshold value is calculated from the histogram of the image. "
			"The first and the second arguments are values corresponding to the minimum and maximum of the histogram. "
			"Out-of-range values will be placed to the first and the last bins, respectively. "
			"The third argument is the count of bins in the histogram.\n"
			"\n"
			"The code implementing this method is probably originally from NIH Image, then ported to ImageJ, and then to pi2.\n"
			"\n"
			"\n"
			"**Li**\n"
			"\n"
			"Li's Minimum Cross Entropy thresholding method.\n"
			"\n"
			"The threshold value is calculated from the histogram of the image. "
			"The first and the second arguments are values corresponding to the minimum and maximum of the histogram. "
			"Out-of-range values will be placed to the first and the last bins, respectively. "
			"The third argument is the count of bins in the histogram.\n"
			"\n"
			"This implementation is based on the iterative version of the algorithm, described in "
			"C.H. Li and P.K.S. Tam, An Iterative Algorithm for Minimum Cross Entropy Thresholding, Pattern Recognition Letters 18(8), 1998.\n"
			"\n"
			"Ported to ImageJ plugin by G.Landini from E Celebi's fourier_0.8 routines, and then to pi2.\n"
			"\n"
			"\n"
			"**MaxEntropy**\n"
			"\n"
			"Implements Kapur-Sahoo-Wong (Maximum Entropy) thresholding method described in "
			"J.N. Kapur, P.K. Sahoo, and A.K.C Wong, A New Method for Gray-Level Picture Thresholding Using the Entropy of the Histogram, Graphical Models and Image Processing 29(3), 1985.\n"
			"\n"
			"The threshold value is calculated from the histogram of the image. "
			"The first and the second arguments are values corresponding to the minimum and maximum of the histogram. "
			"Out-of-range values will be placed to the first and the last bins, respectively. "
			"The third argument is the count of bins in the histogram.\n"
			"\n"
			"Original code by M. Emre Celebi, ported to ImageJ plugin by G.Landini, then ported to pi2.\n"
			"\n"
			"\n"
			"\n"
			"**Mean**\n"
			"\n"
			"The threshold is shifted mean of the greyscale data, i.e. $t = \\mu - c$, where $t$ is the threshold, $\\mu$ is the mean of the image or the neighbourhood, and $c$ is the shift.\n"
			"\n"
			"First argument is the shift value $c$.\n"
			"\n"
			"Described in C. A. Glasbey, An analysis of histogram-based thresholding algorithms, CVGIP: Graphical Models and Image Processing 55, 1993.\n"
			"\n"
			"\n"
			"\n"
			"**MinError**\n"
			"\n"
			"Implements minimum error thresholding method described in  "
			"J. Kittler and J. Illingworth, Minimum error thresholding, Pattern Recognition 19, 1986.\n"
			"\n"
			"The threshold value is calculated from the histogram of the image. "
			"The first and the second arguments are values corresponding to the minimum and maximum of the histogram. "
			"Out-of-range values will be placed to the first and the last bins, respectively. "
			"The third argument is the count of bins in the histogram.\n"
			"\n"
			"Code is originally from Antti Niemisto's Matlab code (GPL), then ported to ImageJ by Gabriel Landini, "
			"and then to pi2.\n"
			"\n"
			"\n"
			"\n"
			"**Minimum**\n"
			"\n"
			"Assumes a bimodal histogram. The histogram needs is iteratively smoothed (using a "
			"running average of size 3) until there are only two local maxima. "
			"Threshold $t$ is such that $y_{t-1} > y_t$ and $y_t \\leq y_{t+1}$. "
			"Images with histograms having extremely unequal peaks or a broad and "
			"flat valleys are unsuitable for this method.\n"
			"\n"
			"The threshold value is calculated from the histogram of the image. "
			"The first and the second arguments are values corresponding to the minimum and maximum of the histogram. "
			"Out-of-range values will be placed to the first and the last bins, respectively. "
			"The third argument is the count of bins in the histogram.\n"
			"\n"
			"Described in J. M. S. Prewitt and M. L. Mendelsohn, The analysis of cell images, Annals of the New York Academy of Sciences 128, 1966.\n"
			"\n"
			"Original Matlab code by Antti Niemisto, ported to ImageJ by Gabriel Landini, "
			"then ported to pi2.\n"
			"\n"
			"\n"
			"**Moments**\n"
			"\n"
			"Attempts to preserve the moments of the original image in the thresholded result.\n"
			"\n"
			"The threshold value is calculated from the histogram of the image. "
			"The first and the second arguments are values corresponding to the minimum and maximum of the histogram. "
			"Out-of-range values will be placed to the first and the last bins, respectively. "
			"The third argument is the count of bins in the histogram.\n"
			"\n"
			"Described in  W. Tsai, Moment-preserving thresholding: a new approach, Computer Vision,	Graphics, and Image Processing 29, 1985.\n"
			"\n"
			"Original code by M. Emre Celebi ported to ImageJ by Gabriel Landini, and then to pi2.\n"
			"\n"
			"\n"
			"**Percentile**\n"
			"\n"
			"Assumes the fraction of foreground pixels to be given value.\n"
			"\n"
			"The threshold value is calculated from the histogram of the image. "
			"The first and the second arguments are values corresponding to the minimum and maximum of the histogram. "
			"Out-of-range values will be placed to the first and the last bins, respectively. "
			"The third argument is the count of bins in the histogram. "
			"The fourth argument is the fraction of foreground pixels.\n"
			"\n"
			"Described in W. Doyle, Operation useful for similarity-invariant pattern recognition, Journal of the Association for Computing Machinery 9, 1962.\n"
			"\n"
			"Original code by Antti Niemisto, ported to ImageJ by Gabriel Landini, then to pi2.\n"
			"\n"
			"\n"
			"\n"
			"**RenyiEntropy**\n"
			"\n"
			"Similar to the MaxEntropy method, but using Renyi's entropy instead.\n"
			"\n"
			"The threshold value is calculated from the histogram of the image. "
			"The first and the second arguments are values corresponding to the minimum and maximum of the histogram. "
			"Out-of-range values will be placed to the first and the last bins, respectively. "
			"The third argument is the count of bins in the histogram.\n"
			"\n"
			"Described in J.N. Kapur, P.K. Sahoo, and A.K.C Wong, A New Method for Gray-Level Picture Thresholding Using the Entropy of the Histogram, Graphical Models and Image Processing 29(3), 1985.\n"
			"\n"
			"Original code by M. Emre Celebi, ported to ImageJ plugin by G.Landini, then ported to pi2.\n"
			"\n"
			"\n"
			"**Shanbhag**\n"
			"\n"
			"Described in A.G. Shanhbag, Utilization of Information Measure as a Means of Image Thresholding, Graphical Models and Image Processing 56(5), 1994.\n"
			"\n"
			"The threshold value is calculated from the histogram of the image. "
			"The first and the second arguments are values corresponding to the minimum and maximum of the histogram. "
			"Out-of-range values will be placed to the first and the last bins, respectively. "
			"The third argument is the count of bins in the histogram.\n"
			"\n"
			"Original code by M. Emre Celebi, ported to ImageJ plugin by Gabriel Landini, then to pi2.\n"
			"\n"
			"\n"
			"**Triangle**\n"
			"\n"
			"The triangle algorithm assumes a peak near either end of the histogram, finds minimum near the other end, "
			"draws a line between the minimum and maximum, and sets threshold $t$ to a value for which the point $(t, y(t))$ is "
			"furthest away from the line (where $y$ is the histogram).\n"
			"\n"
			"The threshold value is calculated from the histogram of the image. "
			"The first and the second arguments are values corresponding to the minimum and maximum of the histogram. "
			"Out-of-range values will be placed to the first and the last bins, respectively. "
			"The third argument is the count of bins in the histogram.\n"
			"\n"
			"Described in G.W. Zack, W.E. Rogers, and S.A. Latt, Automatic Measurement of Sister Chromatid Exchange Frequency, Journal of Histochemistry and Cytochemistry 25(7), 1977.\n"
			"\n"
			"Original code by Johannes Schindelin, modified by Gabriel Landini, then ported to pi2.\n"
			"\n"
			"\n"
			"**Yen**\n"
			"\n"
			"Yen thresholding method described in "
			"J.C. Yen, F.K. Chang, and S. Chang, A New Criterion for Automatic Multilevel Thresholding, IEEE Transactions on Image Processing 4(3), 1995.\n"
			"\n"
			"The threshold value is calculated from the histogram of the image. "
			"The first and the second arguments are values corresponding to the minimum and maximum of the histogram. "
			"Out-of-range values will be placed to the first and the last bins, respectively. "
			"The third argument is the count of bins in the histogram.\n"
			"\n"
			"Original code by M. Emre Celebi, ported to ImageJ plugin by Gabriel Landini, then to pi2.\n"
			"\n"
			"\n"
			"**Median**\n"
			"\n"
			"The threshold is shifted median of the greyscale data, i.e. $t = M - c$, where $t$ is the threshold, $M$ is the median of the image, and $c$ is the shift.\n"
			"\n"
			"The median is calculated from image histogram, so its accuracy might degrade if too small bin count is used. "
			"The first and the second arguments are values corresponding to the minimum and maximum of the histogram. "
			"Out-of-range values will be placed to the first and the last bins, respectively. "
			"The third argument is the count of bins in the histogram. "
			"The fourth argument is the shift value $c$.\n"
			"\n"
			"\n"
			"**MidGrey**\n"
			"\n"
			"The threshold is $(m + M) / 2 - c$, where $m$ and $M$ are the minimum and the maximum value of the image or neighbourhood, respectively. "
			"The value $c$ is a user-specified shift.\n"
			"\n"
			"The first argument is the shift value $c$.\n"
			"\n"
			"\n"
			"**Niblack**\n"
			"\n"
			"Niblack's thresholding method.\n"
			"\n"
			"This method is mostly suited for local thresholding.\n"
			"\n"
			"The threshold is $\\mu + k * \\sigma - c$, where $\\mu$ and $\\sigma$ are the mean and the standard deviation of the image or the neighbourhood, respectively. "
			"The values $k$ and $c$ are user-specified scaling constant and shift.\n"
			"\n"
			"The first argument is the scaling constant $k$. "
			"The second argument is the shift value $c$.\n"
			"\n"
			"Described in W. Niblack, An introduction to Digital Image Processing, Prentice-Hall, 1986.\n"
			"\n"
			"\n"
			"**Phansalkar**\n"
			"\n"
			"Phansalkar's thresholding method.\n"
			"\n"
			"This method is mostly suited for local thresholding.\n"
			"\n"
			"Threshold is $\\mu * (1.0 + p * \\exp(-q * \\mu) + k * (\\sigma / r - 1))$ where $\\mu$ and $\\sigma$ are the mean and the standard deviation of the image or the neighbourhood, respectively. "
			"The values $k$, $r$, $p$, and $q$ are user-specified parameters of the method. "
			"The default values are "
			"$k = 0.25$, "
			"$r = 0.5$,  "
			"$p = 2.0$, and "
			"$q = 10.0$.\n"
			"\n"
			"The first four arguments are the parameters $k$, $r$, $p$, and $q$.\n"
			"\n"
			"Described in N. Phansalskar, S. More, and A. Sabale, et al., Adaptive local thresholding for detection of nuclei in diversity stained cytology images, International Conference on Communications and Signal Processing (ICCSP), 2011.\n"
			"\n"
			"\n"
			"**Sauvola**\n"
			"\n"
			"Sauvola's thresholding method.\n"
			"\n"
			"This method is mostly suited for local thresholding.\n"
			"\n"
			"The threshold is $\\mu * (1 + k * (\\sigma / r - 1))$, where $\\mu$ and $\\sigma$ are the mean and the standard deviation of the image or the neighbourhood, respectively. "
			"The values $k$ and $r$ are user-specified scaling constants.\n"
			"\n"
			"The first argument is the scaling constant $k$. "
			"The second argument is the scaling constant $r$.\n"
			"\n"
			"Described in Sauvola, J and Pietaksinen, M, Adaptive Document Image Binarization, Pattern Recognition 33(2), 2000.\n"
			"\n"
			"\n"
			"**Bernsen**\n"
			"\n"
			"Finds Bernsen's thresholding method.\n"
			"\n"
			"This method is mostly suited for local thresholding.\n"
			"\n"
			"The method uses a user-provided contrast threshold. If the local contrast (max - min) is above or "
			"equal to the contrast threshold, the threshold is set at the local midgrey value (the mean of the minimum "
			"and maximum grey values in the local window (or whole image in the case of global thresholding)). "
			"If the local contrast is below the contrast threshold the neighbourhood is considered to consist only of one "
			"class and the pixel is set to object or background depending on the value of the midgrey.\n"
			"\n"
			"The first argument is the local contrast threshold.\n"
			"\n"
			"Described in J. Bernsen, Dynamic Thresholding of Grey-Level Images, Proceedings of the 8th International Conference on Pattern Recognition, 1986.\n"
			"\n"
			"The original code is written by Gabriel Landini, and it has been ported to pi2.\n"
			"\n";
	}

	template<typename pixel_t> class AutoThresholdCommand : public OneImageInPlaceCommand<pixel_t>, public Distributable
	{
	protected:
		friend class CommandList;

		AutoThresholdCommand() : OneImageInPlaceCommand<pixel_t>("autothreshold",
			"Thresholds image by automatically determined threshold value. The supported thresholding methods are\n\n" +
			autoThresholdMethodsHelp(),
			{
				CommandArgument<string>(ParameterDirection::In, "method", "Thresholding method that will be applied.", "Otsu"),
				CommandArgument<double>(ParameterDirection::In, "argument 1", "Argument for the thresholding method. The purpose of this argument depends on the method, see the list above. Specify nan in order to use a method-specific default value.", std::numeric_limits<double>::quiet_NaN()),
				CommandArgument<double>(ParameterDirection::In, "argument 2", "Argument for the thresholding method. The purpose of this argument depends on the method, see the list above. Specify nan in order to use a method-specific default value.", std::numeric_limits<double>::quiet_NaN()),
				CommandArgument<double>(ParameterDirection::In, "argument 3", "Argument for the thresholding method. The purpose of this argument depends on the method, see the list above. Specify nan in order to use a method-specific default value.", std::numeric_limits<double>::quiet_NaN()),
				CommandArgument<double>(ParameterDirection::In, "argument 4", "Argument for the thresholding method. The purpose of this argument depends on the method, see the list above. Specify nan in order to use a method-specific default value.", std::numeric_limits<double>::quiet_NaN()),
			},
			"localthreshold, threshold")
		{

		}

	public:
		virtual void run(Image<pixel_t>& img, std::vector<ParamVariant>& args) const override
		{
			string methods = pop<string>(args);
			double arg0 = pop<double>(args);
			double arg1 = pop<double>(args);
			double arg2 = pop<double>(args);
			double arg3 = pop<double>(args);

			AutoThresholdMethod method = fromString<AutoThresholdMethod>(methods);

			autoThreshold(img, method, arg0, arg1, arg2, arg3);
		}

		virtual std::vector<std::string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& img = *pop<DistributedImage<pixel_t>*>(args);
			string methods = pop<string>(args);
			double arg0 = pop<double>(args);
			double arg1 = pop<double>(args);
			double arg2 = pop<double>(args);
			double arg3 = pop<double>(args);

			AutoThresholdMethod method = fromString<AutoThresholdMethod>(methods);

			double th;
			
			if (method == AutoThresholdMethod::Otsu
				|| method == AutoThresholdMethod::Huang
				|| method == AutoThresholdMethod::Intermodes
				|| method == AutoThresholdMethod::IsoData
				|| method == AutoThresholdMethod::Li
				|| method == AutoThresholdMethod::MaxEntropy
				|| method == AutoThresholdMethod::MinError
				|| method == AutoThresholdMethod::Minimum
				|| method == AutoThresholdMethod::Moments
				|| method == AutoThresholdMethod::RenyiEntropy
				|| method == AutoThresholdMethod::Shanbhag
				|| method == AutoThresholdMethod::Triangle
				|| method == AutoThresholdMethod::Yen
				|| method == AutoThresholdMethod::Median
				|| method == AutoThresholdMethod::Percentile)
			{
				double rangeMin = arg0;
				double rangeMax = arg1;
				double binCount = arg2;
				double arg = arg3;

				itl2::internals::autothreshold::setDefault(rangeMin, itl2::internals::ThresholdDefaults<pixel_t>::range().x);
				itl2::internals::autothreshold::setDefault(rangeMax, itl2::internals::ThresholdDefaults<pixel_t>::range().y);
				itl2::internals::autothreshold::setDefault(binCount, (double)itl2::internals::ThresholdDefaults<pixel_t>::binCount());

				size_t count = pixelRound<size_t>(binCount);

				DistributedTempImage<float32_t> dhist(distributor, "autothreshold_hist", (coord_t)count);
				DistributedTempImage<float32_t> dbins(distributor, "autothreshold_bins", (coord_t)count);
				CommandList::get<HistogramCommand<pixel_t, float32_t>>().runDistributed(distributor, { &img, &dhist.get(), &dbins.get(), rangeMin, rangeMax, (coord_t)count, std::string(""), Distributor::BLOCK_INDEX_ARG_TYPE() });

				Image<float32_t> hist(count);
				dhist.get().readTo(hist);

				// Uh oh, format conversion
				Image<double> histDouble(count);
				convert(hist, histDouble);

				double bin = itl2::internals::histogramAutoThreshold(histDouble, method, arg);
				th = bin / count * (rangeMax - rangeMin) + rangeMin;
			}
			else if (method == AutoThresholdMethod::Mean)
			{
				DistributedTempImage<float32_t> temp(distributor, "autothreshold_mean");
				CommandList::get<MeanAllPixelsCommand<pixel_t>>().runDistributed(distributor, { &img, &temp.get(), false });
				double mu = (double)temp.get().getValue();
				th = itl2::internals::autothreshold::mean(mu, arg0);
			}
			else if (method == AutoThresholdMethod::MidGrey)
			{
				DistributedTempImage<pixel_t> tempMin(distributor, "autothreshold_min");
				DistributedTempImage<pixel_t> tempMax(distributor, "autothreshold_max");
				CommandList::get<MinAllPixelsCommand<pixel_t>>().runDistributed(distributor, { &img, &tempMin.get(), false });
				CommandList::get<MaxAllPixelsCommand<pixel_t>>().runDistributed(distributor, { &img, &tempMax.get(), false });
				double m = (double)tempMin.get().getValue();
				double M = (double)tempMax.get().getValue();
				th = itl2::internals::autothreshold::midgrey(m, M, arg0);
			}
			else if (method == AutoThresholdMethod::Niblack)
			{
				DistributedTempImage<float32_t> tempMean(distributor, "autothreshold_mean");
				DistributedTempImage<float32_t> tempStd(distributor, "autothreshold_std");
				CommandList::get<MeanAllPixelsCommand<pixel_t>>().runDistributed(distributor, { &img, &tempMean.get(), false });
				CommandList::get<StdDevAllPixelsCommand<pixel_t>>().runDistributed(distributor, { &img, &tempStd.get(), false });
				double mean = (double)tempMean.get().getValue();
				double std = (double)tempStd.get().getValue();
				
				th = itl2::internals::autothreshold::niblack(mean, std, arg0, arg1);
			}
			else if (method == AutoThresholdMethod::Phansalkar)
			{
				DistributedTempImage<float32_t> tempMean(distributor, "autothreshold_mean");
				DistributedTempImage<float32_t> tempStd(distributor, "autothreshold_std");
				CommandList::get<MeanAllPixelsCommand<pixel_t>>().runDistributed(distributor, { &img, &tempMean.get(), false });
				CommandList::get<StdDevAllPixelsCommand<pixel_t>>().runDistributed(distributor, { &img, &tempStd.get(), false });
				double mean = (double)tempMean.get().getValue();
				double std = (double)tempStd.get().getValue();

				th = itl2::internals::autothreshold::phansalkar(mean, std, arg0, arg1, arg2, arg3);
			}
			else if (method == AutoThresholdMethod::Sauvola)
			{
				DistributedTempImage<float32_t> tempMean(distributor, "autothreshold_mean");
				DistributedTempImage<float32_t> tempStd(distributor, "autothreshold_std");
				CommandList::get<MeanAllPixelsCommand<pixel_t>>().runDistributed(distributor, { &img, &tempMean.get(), false });
				CommandList::get<StdDevAllPixelsCommand<pixel_t>>().runDistributed(distributor, { &img, &tempStd.get(), false });
				double mean = (double)tempMean.get().getValue();
				double std = (double)tempStd.get().getValue();

				th = itl2::internals::autothreshold::sauvola(mean, std, arg0, arg1, arg2);
			}
			else if (method == AutoThresholdMethod::Bernsen)
			{
				DistributedTempImage<pixel_t> tempMin(distributor, "autothreshold_min");
				DistributedTempImage<pixel_t> tempMax(distributor, "autothreshold_max");
				CommandList::get<MinAllPixelsCommand<pixel_t>>().runDistributed(distributor, { &img, &tempMin.get(), false });
				CommandList::get<MaxAllPixelsCommand<pixel_t>>().runDistributed(distributor, { &img, &tempMax.get(), false });
				double m = (double)tempMin.get().getValue();
				double M = (double)tempMax.get().getValue();
				th = itl2::internals::autothreshold::bernsen(m, M, arg0, (m + M) / 2.0);
			}
			else
			{
				throw ITLException(string("Unsupported auto-thresholding method: ") + itl2::toString(method));
			}
			
			CommandList::get<ThresholdConstantCommand<pixel_t> >().runDistributed(distributor, { &img, th });
			
			return std::vector<std::string>();
		}
	};



	template<typename pixel_t> class LocalThresholdCommand : public OverlapDistributable<TwoImageInputOutputCommand<pixel_t>>
	{
	protected:
		friend class CommandList;

		LocalThresholdCommand() : OverlapDistributable<TwoImageInputOutputCommand<pixel_t>>("localthreshold",
			"Local thresholding. Determines threshold value for a pixel based on its neighbourhood. The supported thresholding methods are\n\n" +
			autoThresholdMethodsHelp(),
			{
				CommandArgument<Vec3c>(ParameterDirection::In, "radius", "Radius of neighbourhood. Diameter will be $2r+1$.", Vec3c(1, 1, 1)),
				CommandArgument<string>(ParameterDirection::In, "method", "Thresholding method that will be applied.", "Otsu"),
				CommandArgument<double>(ParameterDirection::In, "argument 1", "Argument for the thresholding method. The purpose of this argument depends on the method, see the list above. Specify nan in order to use a method-specific default value.", std::numeric_limits<double>::quiet_NaN()),
				CommandArgument<double>(ParameterDirection::In, "argument 2", "Argument for the thresholding method. The purpose of this argument depends on the method, see the list above. Specify nan in order to use a method-specific default value.", std::numeric_limits<double>::quiet_NaN()),
				CommandArgument<double>(ParameterDirection::In, "argument 3", "Argument for the thresholding method. The purpose of this argument depends on the method, see the list above. Specify nan in order to use a method-specific default value.", std::numeric_limits<double>::quiet_NaN()),
				CommandArgument<double>(ParameterDirection::In, "argument 4", "Argument for the thresholding method. The purpose of this argument depends on the method, see the list above. Specify nan in order to use a method-specific default value.", std::numeric_limits<double>::quiet_NaN()),
				CommandArgument<BoundaryCondition>(ParameterDirection::In, "boundary condition", string("Type of boundary condition. ") + boundaryConditionHelp(), BoundaryCondition::Nearest)
			},
			"autothreshold, threshold")
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, std::vector<ParamVariant>& args) const override
		{
			Vec3c r = pop<Vec3c>(args);
			string methods = pop<string>(args);
			double arg0 = pop<double>(args);
			double arg1 = pop<double>(args);
			double arg2 = pop<double>(args);
			double arg3 = pop<double>(args);
			BoundaryCondition bc = pop<BoundaryCondition>(args);

			AutoThresholdMethod method = fromString<AutoThresholdMethod>(methods);

			localThreshold(in, out, r, method, arg0, arg1, arg2, arg3, bc);
		}

		virtual Vec3c calculateOverlap(const vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& in = *std::get<DistributedImage<pixel_t>*>(args[0]);
			DistributedImage<pixel_t>& out = *std::get<DistributedImage<pixel_t>*>(args[1]);
			Vec3c r = std::get<Vec3c>(args[2]);

			in.mustNotBe(out);
			out.ensureSize(in);
			
			return 2 * r + Vec3c(1, 1, 1) + Vec3c(2, 2, 2); // Extra 2 pixels margin for safety.
		}
	};
}