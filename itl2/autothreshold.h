#pragma once

#include "histogram.h"
#include "pointprocess.h"
#include "math/vec2.h"
#include "projections.h"
#include "median.h"
#include "neighbourhood.h"
#include "filters.h"

#include <algorithm>

namespace itl2
{
	enum class AutoThresholdMethod
	{
		/**
		Otsu's thresholding algorithm.

		The threshold value is calculated from the histogram of the image.
		The first and the second arguments are values corresponding to the minimum and maximum of the histogram.
		Out-of-range values will be placed to the first and the last bins, respectively.
		The third argument is the count of bins in the histogram.

		Described in N. Otsu, A Threshold Selection Method from Gray-Level Histograms, IEEE Transactions on Systems, Man, and Cybernetics 9(1), 1979.

		Original C++ code by Jordan Bevik, ported to ImageJ plugin by Gabriel Landini, and finally ported to pi2.
		*/
		Otsu,

		/**
		Huang's fuzzy thresholding method.

		Uses Shannon's entropy function (one can also use Yager's entropy function).

		The threshold value is calculated from the histogram of the image.
		The first and the second arguments are values corresponding to the minimum and maximum of the histogram.
		Out-of-range values will be placed to the first and the last bins, respectively.
		The third argument is the count of bins in the histogram.

		Described in
		L.-K. Huang and M.-J.J. Wang, Image Thresholding by Minimizing the Measures of Fuzziness, Pattern Recognition 28(1), 1995.

		Original code by M. Emre Celebi, ported to ImageJ plugin by G. Landini from E. Celebi's fourier_0.8 routines, then ported to pi2.
		*/
		Huang,

		/**
		Assumes a bimodal histogram. The histogram needs is iteratively smoothed using a
		running average of size 3 until there are only two local maxima at j and k.
		The threshold t is (j+k)/2.

		Images with histograms having extremely unequal peaks or a broad and
		flat valleys are unsuitable for this method.

		The threshold value is calculated from the histogram of the image.
		The first and the second arguments are values corresponding to the minimum and maximum of the histogram.
		Out-of-range values will be placed to the first and the last bins, respectively.
		The third argument is the count of bins in the histogram.

		Described in
		J. M. S. Prewitt and M. L. Mendelsohn, The analysis of cell images, Annals of the New York Academy of Sciences 128, 1966.

		Ported to ImageJ plugin by Gabriel Landini from Antti Niemisto's Matlab code (GPL), and then to pi2.
		*/
		Intermodes,

		/**
		Iterative procedure based on the isodata algorithm described in
		T.W. Ridler and S. Calvard, Picture thresholding using an iterative selection method, IEEE Transactions on System, Man and Cybernetics, SMC-8, 1978.

		The procedure divides the image into objects and background by taking an initial threshold,
		then the averages of the pixels at or below the threshold and pixels above are computed.
		The averages of those two values are computed, the threshold is incremented and the
		process is repeated until the threshold is larger than the composite average. That is,
		threshold = (average background + average objects) / 2

		The threshold value is calculated from the histogram of the image.
		The first and the second arguments are values corresponding to the minimum and maximum of the histogram.
		Out-of-range values will be placed to the first and the last bins, respectively.
		The third argument is the count of bins in the histogram.

		The code implementing this method is probably originally from NIH Image, then ported to ImageJ, and then to pi2.
		*/
		IsoData,

		/**
		Li's Minimum Cross Entropy thresholding method.

		The threshold value is calculated from the histogram of the image.
		The first and the second arguments are values corresponding to the minimum and maximum of the histogram.
		Out-of-range values will be placed to the first and the last bins, respectively.
		The third argument is the count of bins in the histogram.

		This implementation is based on the iterative version of the algorithm, described in
		C.H. Li and P.K.S. Tam, An Iterative Algorithm for Minimum Cross Entropy Thresholding, Pattern Recognition Letters 18(8), 1998.

		Ported to ImageJ plugin by G.Landini from E Celebi's fourier_0.8 routines, and then to pi2.
		*/
		Li,

		/**
		Implements Kapur-Sahoo-Wong (Maximum Entropy) thresholding method described in
		J.N. Kapur, P.K. Sahoo, and A.K.C Wong, A New Method for Gray-Level Picture Thresholding Using the Entropy of the Histogram, Graphical Models and Image Processing 29(3), 1985.

		The threshold value is calculated from the histogram of the image.
		The first and the second arguments are values corresponding to the minimum and maximum of the histogram.
		Out-of-range values will be placed to the first and the last bins, respectively.
		The third argument is the count of bins in the histogram.

		Original code by M. Emre Celebi, ported to ImageJ plugin by G.Landini, then ported to pi2.
		*/
		MaxEntropy,

		/**
		The threshold is shifted mean of the greyscale data, i.e. $t = \mu - c$, where $t$ is the threshold, $\mu$ is the mean of the image or the neighbourhood, and $c$ is the shift.

		First argument is the shift value $c$.

		Described in C. A. Glasbey, An analysis of histogram-based thresholding algorithms, CVGIP: Graphical Models and Image Processing 55, 1993.
		*/
		Mean,

		/**
		Implements minimum error thresholding method described in
		J. Kittler and J. Illingworth, Minimum error thresholding, Pattern Recognition 19, 1986.

		The threshold value is calculated from the histogram of the image.
		The first and the second arguments are values corresponding to the minimum and maximum of the histogram.
		Out-of-range values will be placed to the first and the last bins, respectively.
		The third argument is the count of bins in the histogram.

		Code is originally from Antti Niemisto's Matlab code (GPL), then ported to ImageJ by Gabriel Landini,
		and then to pi2.
		*/
		MinError,

		/**
		Assumes a bimodal histogram. The histogram needs is iteratively smoothed (using a
		running average of size 3) until there are only two local maxima.
		Threshold $t$ is such that $y_{t-1} > y_t$ and $y_t <= y_{t+1}$.
		Images with histograms having extremely unequal peaks or a broad and
		flat valleys are unsuitable for this method.

		The threshold value is calculated from the histogram of the image.
		The first and the second arguments are values corresponding to the minimum and maximum of the histogram.
		Out-of-range values will be placed to the first and the last bins, respectively.
		The third argument is the count of bins in the histogram.

		Described in J. M. S. Prewitt and M. L. Mendelsohn, The analysis of cell images, Annals of the New York Academy of Sciences 128, 1966.

		Original Matlab code by Antti Niemisto, ported to ImageJ by Gabriel Landini,
		then ported to pi2.
		*/
		Minimum,

		/**
		Attempts to preserve the moments of the original image in the thresholded result.

		The threshold value is calculated from the histogram of the image.
		The first and the second arguments are values corresponding to the minimum and maximum of the histogram.
		Out-of-range values will be placed to the first and the last bins, respectively.
		The third argument is the count of bins in the histogram.

		Described in  W. Tsai, Moment-preserving thresholding: a new approach, Computer Vision,	Graphics, and Image Processing 29, 1985.

		Original code by M. Emre Celebi ported to ImageJ by Gabriel Landini, and then to pi2.
		*/
		Moments,

		/**
		Assumes the fraction of foreground pixels to be given value.

		The threshold value is calculated from the histogram of the image.
		The first and the second arguments are values corresponding to the minimum and maximum of the histogram.
		Out-of-range values will be placed to the first and the last bins, respectively.
		The third argument is the count of bins in the histogram.
		The fourth argument is the fraction of foreground pixels.

		Described in W. Doyle, Operation useful for similarity-invariant pattern recognition, Journal of the Association for Computing Machinery 9, 1962.

		Original code by Antti Niemisto, ported to ImageJ by Gabriel Landini, then to pi2.
		*/
		Percentile,

		/**
		Similar to the MaxEntropy method, but using Renyi's entropy instead.

		The threshold value is calculated from the histogram of the image.
		The first and the second arguments are values corresponding to the minimum and maximum of the histogram.
		Out-of-range values will be placed to the first and the last bins, respectively.
		The third argument is the count of bins in the histogram.

		Described in J.N. Kapur, P.K. Sahoo, and A.K.C Wong, A New Method for Gray-Level Picture Thresholding Using the Entropy of the Histogram, Graphical Models and Image Processing 29(3), 1985.

		Original code by M. Emre Celebi, ported to ImageJ plugin by G.Landini, then ported to pi2.
		*/
		RenyiEntropy,

		/**
		Described in A.G. Shanhbag, Utilization of Information Measure as a Means of Image Thresholding, Graphical Models and Image Processing 56(5), 1994.

		The threshold value is calculated from the histogram of the image.
		The first and the second arguments are values corresponding to the minimum and maximum of the histogram.
		Out-of-range values will be placed to the first and the last bins, respectively.
		The third argument is the count of bins in the histogram.

		Original code by M. Emre Celebi, ported to ImageJ plugin by Gabriel Landini, then to pi2.
		*/
		Shanbhag,

		/**
		The triangle algorithm assumes a peak near either end of the histogram, finds minimum near the other end,
		draws a line between the minimum and maximum, and sets threshold $t$ to a value for which the point $(t, y(t))$ is
		furthest away from the line (where $y$ is the histogram).

		The threshold value is calculated from the histogram of the image.
		The first and the second arguments are values corresponding to the minimum and maximum of the histogram.
		Out-of-range values will be placed to the first and the last bins, respectively.
		The third argument is the count of bins in the histogram.

		Described in G.W. Zack, W.E. Rogers, and S.A. Latt, Automatic Measurement of Sister Chromatid Exchange Frequency, Journal of Histochemistry and Cytochemistry 25(7), 1977.

		Original code by Johannes Schindelin, modified by Gabriel Landini, then ported to pi2.
		*/
		Triangle,

		/**
		Yen thresholding method described in
		J.C. Yen, F.K. Chang, and S. Chang, A New Criterion for Automatic Multilevel Thresholding, IEEE Transactions on Image Processing 4(3), 1995.

		The threshold value is calculated from the histogram of the image.
		The first and the second arguments are values corresponding to the minimum and maximum of the histogram.
		Out-of-range values will be placed to the first and the last bins, respectively.
		The third argument is the count of bins in the histogram.

		Original code by M. Emre Celebi, ported to ImageJ plugin by Gabriel Landini, then to pi2.
		*/
		Yen,

		/**
		The threshold is shifted median of the greyscale data, i.e. $t = M - c$, where $t$ is the threshold, $M$ is the median of the image, and $c$ is the shift.

		The median is calculated from image histogram, so its accuracy might degrade if too small bin count is used.
		The first and the second arguments are values corresponding to the minimum and maximum of the histogram.
		Out-of-range values will be placed to the first and the last bins, respectively.
		The third argument is the count of bins in the histogram.
		The fourth argument is the shift value $c$.
		*/
		Median,

		/**
		The threshold is $(m + M) / 2 - c$, where $m$ and $M$ are the minimum and the maximum value of the image or neighbourhood, respectively.
		The value $c$ is a user-specified shift.

		The first argument is the shift value $c$.
		*/
		MidGrey,

		/**
		Niblack's thresholding method.

		This method is mostly suited for local thresholding.

		The threshold is $\mu + k * \sigma - c$, where $\mu$ and $\sigma$ are the mean and the standard deviation of the image or the neighbourhood, respectively.
		The values $k$ and $c$ are user-specified scaling constant and shift.

		The first argument is the scaling constant $k$.
		The second argument is the shift value $c$.

		Described in W. Niblack, An introduction to Digital Image Processing, Prentice-Hall, 1986.
		*/
		Niblack,

		/**
		Phansalkar's thresholding method.

		This method is mostly suited for local thresholding.

		Threshold is $\mu * (1.0 + p * \exp(-q * \mu) + k * (\sigma / r - 1))$ where $\mu$ and $\sigma$ are the mean and the standard deviation of the image or the neighbourhood, respectively.
		The values $k$, $r$, $p$, and $q$ are user-specified parameters of the method.
		The default values are
		$k = 0.25$,
		$r = 0.5$,
		$p = 2.0$, and
		$q = 10.0$.

		The first four arguments are the parameters $k$, $r$, $p$, and $q$.

		Described in N. Phansalskar, S. More, and A. Sabale, et al., Adaptive local thresholding for detection of nuclei in diversity stained cytology images, International Conference on Communications and Signal Processing (ICCSP), 2011.
		*/
		Phansalkar,

		/**
		Sauvola's thresholding method.

		This method is mostly suited for local thresholding.

		The threshold is $\mu * (1 + k * (\sigma / r - 1))$, where $\mu$ and $\sigma$ are the mean and the standard deviation of the image or the neighbourhood, respectively.
		The values $k$ and $r$ are user-specified scaling constants.

		The first argument is the scaling constant $k$.
		The second argument is the scaling constant $r$.

		Described in Sauvola, J and Pietaksinen, M, Adaptive Document Image Binarization, Pattern Recognition 33(2), 2000.
		*/
		Sauvola,

		/**
		Finds Bernsen's thresholding method.

		This method is mostly suited for local thresholding.

		The method uses a user-provided contrast threshold. If the local contrast (max - min) is above or
		equal to the contrast threshold, the threshold is set at the local midgrey value (the mean of the minimum
		and maximum grey values in the local window (or whole image in the case of global thresholding)).
		If the local contrast is below the contrast threshold the neighbourhood is considered to consist only of one
		class and the pixel is set to object or background depending on the value of the midgrey.

		The first argument is the local contrast threshold.

		Described in J. Bernsen, Dynamic Thresholding of Grey-Level Images, Proceedings of the 8th International Conference on Pattern Recognition, 1986.

		The original code is written by Gabriel Landini, and it has been ported to pi2.
		*/
		Bernsen
	};

	template<>
	inline AutoThresholdMethod fromString(const string& str0)
	{
		string str = str0;
		toLower(str);
		if (str == "otsu")
			return AutoThresholdMethod::Otsu;
		if (str == "huang")
			return AutoThresholdMethod::Huang;
		if (str == "intermodes")
			return AutoThresholdMethod::Intermodes;
		if (str == "isodata")
			return AutoThresholdMethod::IsoData;
		if (str == "li")
			return AutoThresholdMethod::Li;
		if (str == "maxentropy")
			return AutoThresholdMethod::MaxEntropy;
		if (str == "mean")
			return AutoThresholdMethod::Mean;
		if (str == "minerror")
			return AutoThresholdMethod::MinError;
		if (str == "minimum")
			return AutoThresholdMethod::Minimum;
		if (str == "moments")
			return AutoThresholdMethod::Moments;
		if (str == "percentile")
			return AutoThresholdMethod::Percentile;
		if (str == "renyi" || str == "renyientropy")
			return AutoThresholdMethod::RenyiEntropy;
		if (str == "shanbhag")
			return AutoThresholdMethod::Shanbhag;
		if (str == "triangle")
			return AutoThresholdMethod::Triangle;
		if (str == "yen")
			return AutoThresholdMethod::Yen;

		if (str == "median")
			return AutoThresholdMethod::Median;
		if (str == "midgrey")
			return AutoThresholdMethod::MidGrey;
		if (str == "niblack")
			return AutoThresholdMethod::Niblack;
		if (str == "phansalkar")
			return AutoThresholdMethod::Phansalkar;
		if (str == "sauvola")
			return AutoThresholdMethod::Sauvola;
		if (str == "bernsen")
			return AutoThresholdMethod::Bernsen;

		throw ITLException("Invalid auto-threshold method: " + str0);
	}

	inline string toString(const AutoThresholdMethod type)
	{
		switch (type)
		{
		case AutoThresholdMethod::Otsu: return "Otsu";
		case AutoThresholdMethod::Huang: return "Huang";
		case AutoThresholdMethod::Intermodes: return "Intermodes";
		case AutoThresholdMethod::IsoData: return "IsoData";
		case AutoThresholdMethod::Li: return "Li";
		case AutoThresholdMethod::MaxEntropy: return "MaxEntropy";
		case AutoThresholdMethod::Mean: return "Mean";
		case AutoThresholdMethod::MinError: return "MinError";
		case AutoThresholdMethod::Minimum: return "Minimum";
		case AutoThresholdMethod::Moments: return "Moments";
		case AutoThresholdMethod::Percentile: return "Percentile";
		case AutoThresholdMethod::RenyiEntropy: return "Renyi";
		case AutoThresholdMethod::Shanbhag: return "Shanbhag";
		case AutoThresholdMethod::Triangle: return "Triangle";
		case AutoThresholdMethod::Yen: return "Yen";
		case AutoThresholdMethod::Median: return "Median";
		case AutoThresholdMethod::MidGrey: return "MidGrey";
		case AutoThresholdMethod::Niblack: return "Niblack";
		case AutoThresholdMethod::Phansalkar: return "Phansalkar";
		case AutoThresholdMethod::Sauvola: return "Sauvola";
		case AutoThresholdMethod::Bernsen: return "Bernsen";
		default: throw ITLException("Invalid auto-threshold method.");
		}
	}


	namespace internals
	{
		namespace autothreshold
		{
			/**
			If val is nan, replaces it by def.
			*/
			inline void setDefault(double& val, double def)
			{
				if (std::isnan(val))
					val = def;
			}

			/**
			These functions determine the index of bin corresponding to threshold as given by each of the methods
			indicated in the function name.
			The code is ported from ImageJ (public domain) and Fiji (GPL).
			The input parameter is the histogram.
			*/
			double otsu(const Image<double>& histogram);
			double huang(const Image<double>& data);
			double intermodes(const Image<double>& data);
			double isoData(const Image<double>& data);
			double li(const Image<double>& data);
			double maxEntropy(const Image<double>& data);
			//double mean(const Image<double>& data, double c);
			double minErrorI(const Image<double>& data);
			double minimum(const Image<double>& data);
			double moments(const Image<double>& data);
			double percentile(const Image<double>& data, double foregroundFraction);
			double renyiEntropy(const Image<double>& data);
			double shanbhag(const Image<double>& data);
			double triangle(const Image<double>& data);
			double yen(const Image<double>& data);

			/**
			Finds mean threshold.
			Threshold is M - c.
			@param M mean value in the image/neighbourhood.
			@param c Shift value.
			*/
			inline double mean(double M, double c = 0.0)
			{
				// Defaults: arg0 = c = 0
				setDefault(c, 0);
				return M - c;
			}

			/**
			Finds median threshold.
			Threshold is M - c.
			@param M Median value in the image/neighbourhood.
			@param c Shift value.
			*/
			inline double median(double M, double c = 0.0)
			{
				// Defaults: arg = c = 0
				setDefault(c, 0);

				return M - c;
			}

			/**
			Finds mid-grey threshold.
			Threshold is (min + max) / 2 - c.
			@param min Minimum value in the image/neighbourhood.
			@param max Maximum value in the image/neighbourhood.
			@param c Shift value.
			*/
			inline double midgrey(double min, double max, double c = 0.0)
			{
				// Defaults: arg0 = c = 0
				setDefault(c, 0);

				return (min + max) / 2.0 - c;
			}

			/**
			Finds Niblack's threshold.

			Threshold is mean + k * standard_deviation - c.

			Described in W. Niblack, An introduction to Digital Image Processing, Prentice-Hall, 1986.

			@param M mean value in the image/neighbourhood.
			@param stddev Standard deviation value in the image/neighbourhood.
			@param c Shift value.
			*/
			inline double niblack(double mean, double stddev, double k = 1.0, double c = 0.0)
			{
				// Defaults: arg0 = k = 1.0, arg1 = c = 0.0
				setDefault(k, 1.0);
				setDefault(c, 0);

				return mean + k * stddev - c;
			}

			/**
			Finds Phansalkar's threshold.

			Threshold is mean * (1.0 + p * exp(-q * mean) + k * (stddev / r - 1.0)).

			Described in N. Phansalskar, S. More, and A. Sabale, et al., Adaptive local thresholding for detection of nuclei in diversity stained cytology images, International Conference on Communications and Signal Processing (ICCSP), 2011.

			@param M mean value in the image/neighbourhood.
			@param stddev Standard deviation value in the image/neighbourhood.
			*/
			inline double phansalkar(double mean, double stddev, double k = 0.25, double r = 0.5, double p = 2.0, double q = 10.0)
			{
				// Defaults: arg0 = k = 0.25, arg1 = r = 0.5, arg2 = p = 2.0, arg3 = q = 10.0
				setDefault(k, 0.25);
				setDefault(r, 0.5);
				setDefault(p, 2.0);
				setDefault(q, 10.0);

				return mean * (1.0 + p * exp(-q * mean) + k * (stddev / r - 1.0));
			}

			/**
			Finds Sauvola's threshold.

			Threshold is mean * ( 1 + k * ( standard_deviation / r - 1 ) ).

			Described in Sauvola, J and Pietaksinen, M, Adaptive Document Image Binarization, Pattern Recognition 33(2), 2000.

			@param M mean value in the image/neighbourhood.
			@param stddev Standard deviation value in the image/neighbourhood.
			@param c Shift value.
			*/
			inline double sauvola(double mean, double stddev, double k = 1.0, double r = 1.0, double c = 0.0)
			{
				// Defaults: arg0 = k = 1.0, arg1 = r = 1.0, arg2 = c = 0.0
				setDefault(k, 1.0);
				setDefault(r, 1.0);
				setDefault(c, 0.0);

				return mean * (1 + k * (stddev / r - 1));
			}

			/**
			Finds Bernsen's threshold value.

			The method uses a user-provided contrast threshold. If the local contrast (max-min) is above or
			equal to the contrast threshold, the threshold is set at the local midgrey value (the mean of the minimum
			and maximum grey values in the local window).
			If the local contrast is below the contrast threshold the neighbourhood is considered to consist only of one
			class and the pixel is set to object or background depending on the value of the midgrey.

			Described in J. Bernsen, Dynamic Thresholding of Grey-Level Images, Proceedings of the 8th International Conference on Pattern Recognition, 1986.

			@param min Minimum in the image/neighbourhood.
			@param max Maximum in the image/neighbourhood.
			@param contrastThreshold Threshold in local contrast that is used to differentiate neighbourhoods/images containing only single phase from those containing two phases.
			@param rangeCenter Center of the overall range of values in the whole image.
			*/
			inline double bernsen(double min, double max, double contrastThreshold, double rangeCenter)
			{
				// Defaults: arg0 = contrastThreshold = range/2
				setDefault(contrastThreshold, (min + max) / 2.0);

				double localContrast = abs(max - min);
				double midGray = (max + min) / 2.0;
				if (localContrast < contrastThreshold)
				{
					// There is only one phase in this region.
					if (midGray >= rangeCenter)
						// Mid-gray value is in the high part of the overall value range. Set threshold
						// to minimum value so that all pixels in the image/neighbourhood are set to 'object'.
						return min;
					else
						// Mid-gray value is in the low part of the overall value range. Set threshold
						// to maximum value so that all pixels in the image/neighbourhood are set to 'background'.
						return max;
				}
				else
				{
					// There are multiple phases in this region. Use midGray as the threshold.
					return midGray;
				}
			}
		}


		template<typename pixel_t, typename Enable = void> class ThresholdDefaults
		{
		public:
			static inline Vec2d range()
			{
				pixel_t::ThresholdDefaults_undefined;
				return Vec2d();
			}

			static inline size_t binCount()
			{
				pixel_t::ThresholdDefaults_undefined;
				return 0;
			}
		};


		template<class pixel_t> class ThresholdDefaults<pixel_t, std::enable_if_t<std::is_integral_v<pixel_t> > >
		{
		public:
			static inline Vec2d range()
			{
				return Vec2d((double)std::numeric_limits<pixel_t>::lowest(), (double)std::numeric_limits<pixel_t>::max());
			}

			static inline size_t binCount()
			{
				return std::min((coord_t)500, std::max((coord_t)0, (coord_t)std::round((double)std::numeric_limits<pixel_t>::max() - (double)std::numeric_limits<pixel_t>::lowest() + 1)));
			}
		};

		template<class pixel_t> class ThresholdDefaults<pixel_t, std::enable_if_t<std::is_floating_point_v<pixel_t> > >
		{
		public:
			static inline Vec2d range()
			{
				return Vec2d(0, 1);
			}

			static inline size_t binCount()
			{
				return 500;
			}
		};



		inline double histogramAutoThreshold(Image<double>& hist, AutoThresholdMethod method, double arg)
		{
			double bin;
			switch (method)
			{
			case AutoThresholdMethod::Otsu: bin = internals::autothreshold::otsu(hist); break;
			case AutoThresholdMethod::Huang: bin = internals::autothreshold::huang(hist); break;
			case AutoThresholdMethod::Intermodes: bin = internals::autothreshold::intermodes(hist); break;
			case AutoThresholdMethod::IsoData: bin = internals::autothreshold::isoData(hist); break;
			case AutoThresholdMethod::Li: bin = internals::autothreshold::li(hist); break;
			case AutoThresholdMethod::MaxEntropy: bin = internals::autothreshold::maxEntropy(hist); break;
			case AutoThresholdMethod::MinError: bin = internals::autothreshold::minErrorI(hist); break;
			case AutoThresholdMethod::Minimum: bin = internals::autothreshold::minimum(hist); break;
			case AutoThresholdMethod::Moments: bin = internals::autothreshold::moments(hist); break;
			case AutoThresholdMethod::RenyiEntropy: bin = internals::autothreshold::renyiEntropy(hist); break;
			case AutoThresholdMethod::Shanbhag: bin = internals::autothreshold::shanbhag(hist); break;
			case AutoThresholdMethod::Triangle: bin = internals::autothreshold::triangle(hist); break;
			case AutoThresholdMethod::Yen: bin = internals::autothreshold::yen(hist); break;
			case AutoThresholdMethod::Median:
			{
				// uh oh - format conversion for median calculation
				std::vector<double> bins;
				std::vector<double> counts;
				bins.reserve(hist.pixelCount());
				counts.reserve(hist.pixelCount());
				for (coord_t n = 0; n < hist.pixelCount(); n++)
				{
					bins.push_back((double)n);
					counts.push_back(hist(n));
				}

				double med = weightedMedian(bins, counts);
				bin = internals::autothreshold::median(med, arg);
			}
			break;
			case AutoThresholdMethod::Percentile: bin = internals::autothreshold::percentile(hist, arg); break;
			default: throw std::logic_error("Invalid auto threshold mode passed to helper function.");
			}

			return bin;
		}

		/**
		Calculates threshold value for given image/neighbourhood automatically.
		arg* parameters are documented in AutoThresholdMethod enum docs.
		*/
		template<typename pixel_t> double autoThresholdValue(const Image<pixel_t>& img,
			AutoThresholdMethod method,
			double arg0,
			double arg1,
			double arg2,
			double arg3)
		{
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

				internals::autothreshold::setDefault(rangeMin, internals::ThresholdDefaults<pixel_t>::range().x);
				internals::autothreshold::setDefault(rangeMax, internals::ThresholdDefaults<pixel_t>::range().y);
				internals::autothreshold::setDefault(binCount, (double)internals::ThresholdDefaults<pixel_t>::binCount());

				size_t count = pixelRound<size_t>(binCount);

				Image<double> hist(count);
				histogram<pixel_t, double, double>(img, hist, Vec2d(rangeMin, rangeMax), 0, nullptr);

				double bin = internals::histogramAutoThreshold(hist, method, arg);
				th = bin / count * (rangeMax - rangeMin) + rangeMin;
			}
			else if (method == AutoThresholdMethod::Mean)
			{
				double mu = mean(img);
				th = internals::autothreshold::mean(mu, arg0);
			}
			else if (method == AutoThresholdMethod::MidGrey)
			{
				double m = (double)min(img);
				double M = (double)max(img);
				th = internals::autothreshold::midgrey(m, M, arg0);
			}
			else if (method == AutoThresholdMethod::Niblack)
			{
				Vec2d stats = meanAndStdDev(img);
				th = internals::autothreshold::niblack(stats[0], stats[1], arg0, arg1);
			}
			else if (method == AutoThresholdMethod::Phansalkar)
			{
				Vec2d stats = meanAndStdDev(img);
				th = internals::autothreshold::phansalkar(stats[0], stats[1], arg0, arg1, arg2, arg3);
			}
			else if (method == AutoThresholdMethod::Sauvola)
			{
				Vec2d stats = meanAndStdDev(img);
				th = internals::autothreshold::sauvola(stats[0], stats[1], arg0, arg1, arg2);
			}
			else if (method == AutoThresholdMethod::Bernsen)
			{
				double m = (double)min(img);
				double M = (double)max(img);
				th = internals::autothreshold::bernsen(m, M, arg0, (m + M) / 2.0);
			}
			else
			{
				throw ITLException(string("Unsupported auto-thresholding method: ") + toString(method));
			}

			return th;
		}
	}

	/**
	Thresholds the given image with selected automatic thresholding method.
	@param img Image to be thresholded.
	@param method The thresholding method to use.
	@param arg0, arg1, arg2, arg3 Method-specific arguments. See AutoThresholdMethod for details.
	@return The threshold value that was used.
	*/
	template<typename pixel_t> pixel_t autoThreshold(Image<pixel_t>& img,
		AutoThresholdMethod method = AutoThresholdMethod::Otsu,
		double arg0 = std::numeric_limits<double>::quiet_NaN(),
		double arg1 = std::numeric_limits<double>::quiet_NaN(),
		double arg2 = std::numeric_limits<double>::quiet_NaN(),
		double arg3 = std::numeric_limits<double>::quiet_NaN())
	{
		double th = internals::autoThresholdValue(img, method, arg0, arg1, arg2, arg3);

		pixel_t thi = pixelRound<pixel_t>(th);

		threshold(img, thi);

		return thi;
	}


	namespace internals
	{
		struct LocalThresholdSettings
		{
			AutoThresholdMethod method;
			double arg0, arg1, arg2, arg3;
		};

		/**
		Processes single neighbourhood in local thresholding.
		NOTE: This function supports only rectangular neighbourhoods.
		*/
		template<typename pixel_t> typename NumberUtils<pixel_t>::FloatType localThresholdProcessNeighbourhood(const Image<pixel_t>& nb, const Image<pixel_t>& mask, const LocalThresholdSettings& settings)
		{
			double th = autoThresholdValue(nb, settings.method, settings.arg0, settings.arg1, settings.arg2, settings.arg3);
			if (intuitive::gt(nb(nb.dimensions() / 2), th))
				return (typename NumberUtils<pixel_t>::FloatType)1;
			else
				return (typename NumberUtils<pixel_t>::FloatType)0;
		}
	}

	/**
	Applies local thresholding to the image.
	@param img Input image.
	@param out Output image. This will store the binarized image.
	@param method Thesholding method.
	@param arg0, arg1, arg2, arg3 Arguments for the thresholding method. These are documented in docs of AutoThresholdMethod enumeration.
	@param bc Boundary condition.
	*/
	template<typename pixel_t> void localThreshold(const Image<pixel_t>& img,
		Image<pixel_t>& out,
		const Vec3c& radius = Vec3c(5, 5, 5),
		AutoThresholdMethod method = AutoThresholdMethod::Otsu,
		double arg0 = std::numeric_limits<double>::quiet_NaN(),
		double arg1 = std::numeric_limits<double>::quiet_NaN(),
		double arg2 = std::numeric_limits<double>::quiet_NaN(),
		double arg3 = std::numeric_limits<double>::quiet_NaN(),
		BoundaryCondition bc = BoundaryCondition::Nearest)
	{
		img.mustNotBe(out);
		out.ensureSize(img);

		internals::LocalThresholdSettings settings = {method, arg0, arg1, arg2, arg3};

		// TODO: This can be done faster (using separable filtering) for non-histogram-based auto threshold methods.

		filter<pixel_t, pixel_t, const internals::LocalThresholdSettings&, internals::localThresholdProcessNeighbourhood<pixel_t>>(img, out, radius, settings, NeighbourhoodType::Rectangular, bc);
	}

	namespace tests
	{
		void autothreshold();
		void localThreshold();
	}

}
