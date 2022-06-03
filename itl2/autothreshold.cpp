
#include "autothreshold.h"
#include "pointprocess.h"
#include "io/raw.h"

#include <iostream>
#include <cmath>
using namespace std;

namespace itl2
{
	namespace internals
	{

		namespace autothreshold
		{

			/*
			All the functions below are ported (and some modified) from ImageJ (public domain) and Fiji (GPL).
			*/

			double otsu(const Image<double>& data)
			{
				// Otsu's threshold algorithm
				// C++ code by Jordan Bevik <Jordan.Bevic@qtiworld.com>
				// ported to ImageJ plugin by G.Landini

				double S = 0; // The total intensity of the image
				double N = 0; // The total number of points
				for (int k = 0; k < data.pixelCount(); k++)
				{
					S += (double)k * data(k);   // Total histogram intensity
					N += data(k);				// Total number of data points
				}

				double Sk = 0; // The total intensity for all histogram points <=k
				double N1 = data(0); // N1 = # points with intensity <= k
				double BCV = 0; // The current Between Class Variance
				double BCVmax = 0; // maximum BCV
				double kStar = 0; // Optimal threshold

				// Look at each possible threshold value,
				// calculate the between-class variance, and decide if it's a max
				for (int k = 1; k < data.pixelCount() - 1; k++) // No need to check endpoints k = 0 or k = L-1
				{
					Sk += (double)k * data(k);
					N1 += data(k);

					// The float casting here is to avoid compiler warning about loss of precision and
					// will prevent overflow in the case of large saturated images
					double denom = (double)(N1) * (N - N1); // Maximum value of denom is (N^2)/4 =  approx. 3E10

					if (denom != 0)
					{
						// Float here is to avoid loss of precision when dividing
						double num = ((double)N1 / N) * S - Sk;  // Maximum value of num =  255*N = approx 8E7
						BCV = (num * num) / denom;
					}
					else
					{
						BCV = 0;
					}

					if (BCV >= BCVmax)
					{
						// Assign the best threshold found so far
						BCVmax = BCV;
						kStar = k;
					}
				}
				// kStar += 1;  // Use QTI convention that intensity -> 1 if intensity >= k
				// (the algorithm was developed for I-> 1 if I <= k.)
				return kStar;
			}



			double huang(const Image<double>& data)
			{
				// Implements Huang's fuzzy thresholding method 
				// Uses Shannon's entropy function (one can also use Yager's entropy function) 
				// Huang L.-K. and Wang M.-J.J. (1995) "Image Thresholding by Minimizing  
				// the Measures of Fuzziness" Pattern Recognition, 28(1): 41-51
				// M. Emre Celebi  06.15.2007
				// Ported to ImageJ plugin by G. Landini from E Celebi's fourier_0.8 routines

				/* Determine the first non-zero bin */
				coord_t first_bin = 0;
				for (coord_t ih = 0; ih < data.pixelCount(); ih++) {
					if (data(ih) != 0) {
						first_bin = ih;
						break;
					}
				}

				/* Determine the last non-zero bin */
				coord_t last_bin = data.pixelCount() - 1;
				for (coord_t ih = data.pixelCount() - 1; ih >= first_bin; ih--) {
					if (data(ih) != 0) {
						last_bin = ih;
						break;
					}
				}

				double term = 1.0 / (double)(last_bin - first_bin);
				Image<double> mu_0(data.pixelCount());
				double sum_pix = 0;
				double num_pix = 0;
				for (coord_t ih = first_bin; ih < data.pixelCount(); ih++) {
					sum_pix += ih * data(ih);
					num_pix += data(ih);
					/* NUM_PIX cannot be zero ! */
					mu_0(ih) = sum_pix / (double)num_pix;
				}

				Image<double> mu_1(data.pixelCount());
				sum_pix = 0;
				num_pix = 0;
				for (coord_t ih = last_bin; ih > 0; ih--) {
					sum_pix += ih * data(ih);
					num_pix += data(ih);
					/* NUM_PIX cannot be zero ! */
					mu_1(ih - 1) = sum_pix / (double)num_pix;
				}

				/* Determine the threshold that minimizes the fuzzy entropy */
				coord_t threshold = -1;
				double min_ent = numeric_limits<double>::max();
				for (coord_t it = 0; it < data.pixelCount(); it++) {
					double ent = 0.0;
					for (coord_t ih = 0; ih <= it; ih++) {
						/* Equation (4) in Ref. 1 */
						double mu_x = 1.0 / (1.0 + term * abs(ih - mu_0(it)));
						if (!((mu_x < 1e-06) || (mu_x > 0.999999))) {
							/* Equation (6) & (8) in Ref. 1 */
							ent += data(ih) * (-mu_x * log(mu_x) - (1.0 - mu_x) * log(1.0 - mu_x));
						}
					}

					for (coord_t ih = it + 1; ih < data.pixelCount(); ih++) {
						/* Equation (4) in Ref. 1 */
						double mu_x = 1.0 / (1.0 + term * abs(ih - mu_1(it)));
						if (!((mu_x < 1e-06) || (mu_x > 0.999999))) {
							/* Equation (6) & (8) in Ref. 1 */
							ent += data(ih) * (-mu_x * log(mu_x) - (1.0 - mu_x) * log(1.0 - mu_x));
						}
					}
					/* No need to divide by NUM_ROWS * NUM_COLS * LOG(2) ! */
					if (ent < min_ent) {
						min_ent = ent;
						threshold = it;
					}
				}
				return (double)threshold;
			}



			bool bimodalTest(const Image<double>& y) {
				coord_t len = y.pixelCount();
				bool b = false;
				int modes = 0;

				for (int k = 1; k < len - 1; k++) {
					if (y(k - 1) < y(k) && y(k + 1) < y(k)) {
						modes++;
						if (modes > 2)
							return false;
					}
				}
				if (modes == 2)
					b = true;
				return b;
			}

			void smoothUntilBimodal(Image<double>& hist)
			{
				int iter = 0;
				while (!bimodalTest(hist)) {
					//smooth with a 3 point running mean filter
					double previous = 0, current = 0, next = hist(0);
					for (int i = 0; i < hist.pixelCount() - 1; i++) {
						previous = current;
						current = next;
						next = hist(i + 1);
						hist(i) = (previous + current + next) / 3;
					}
					hist(hist.pixelCount() - 1) = (current + next) / 3;
					iter++;
					if (iter > 10000)
						throw ITLException("The histogram does not smooth into bimodal shape in 10000 iterations.");
				}
			}

			double intermodes(const Image<double>& data)
			{
				// J. M. S. Prewitt and M. L. Mendelsohn, "The analysis of cell images," in
				// Annals of the New York Academy of Sciences, vol. 128, pp. 1035-1053, 1966.
				// ported to ImageJ plugin by G.Landini from Antti Niemisto's Matlab code (GPL)
				// Original Matlab code Copyright (C) 2004 Antti Niemisto
				// See http://www.cs.tut.fi/~ant/histthresh/ for an excellent slide presentation
				// and the original Matlab code.
				//
				// Assumes a bimodal histogram. The histogram needs is smoothed (using a
				// running average of size 3, iteratively) until there are only two local maxima.
				// j and k
				// Threshold t is (j+k)/2.
				// Images with histograms having extremely unequal peaks or a broad and
				// flat valleys are unsuitable for this method.

				coord_t minbin = -1;
				coord_t maxbin = -1;
				for (coord_t i = 0; i < data.pixelCount(); i++)
					if (data(i) > 0)
						maxbin = i;

				for (coord_t i = data.pixelCount() - 1; i >= 0; i--)
					if (data(i) > 0)
						minbin = i;

				coord_t length = (maxbin - minbin) + 1;
				Image<double> hist(length);
				for (coord_t i = minbin; i <= maxbin; i++)
					hist(i - minbin) = data(i);

				smoothUntilBimodal(hist);

				// The threshold is the mean between the two peaks.
				coord_t tt = 0;
				for (coord_t i = 1; i < length - 1; i++) {
					if (hist(i - 1) < hist(i) && hist(i + 1) < hist(i)) {
						tt += i;
					}
				}

				double threshold = tt / 2.0;
				return threshold + minbin;

			}


			double minimum(const Image<double>& data) {
				// J. M. S. Prewitt and M. L. Mendelsohn, "The analysis of cell images," in
				// Annals of the New York Academy of Sciences, vol. 128, pp. 1035-1053, 1966.
				// ported to ImageJ plugin by G.Landini from Antti Niemisto's Matlab code (GPL)
				// Original Matlab code Copyright (C) 2004 Antti Niemisto
				// See http://www.cs.tut.fi/~ant/histthresh/ for an excellent slide presentation
				// and the original Matlab code.
				//
				// Assumes a bimodal histogram. The histogram needs is smoothed (using a
				// running average of size 3, iteratively) until there are only two local maxima.
				// Threshold t is such that yt-1 > yt <= yt+1.
				// Images with histograms having extremely unequal peaks or a broad and
				// flat valleys are unsuitable for this method.

				// This is the original preprocessing:
				//Image<double> iHisto(data.pixelCount());
				//for (int i = 0; i < data.pixelCount(); i++)
				//	iHisto(i) = data(i);

				// This is the same thing done in intermodes method:
				coord_t minbin = -1;
				coord_t maxbin = -1;
				for (coord_t i = 0; i < data.pixelCount(); i++)
					if (data(i) > 0)
						maxbin = i;

				for (coord_t i = data.pixelCount() - 1; i >= 0; i--)
					if (data(i) > 0)
						minbin = i;

				coord_t length = (maxbin - minbin) + 1;
				Image<double> iHisto(length);
				for (coord_t i = minbin; i <= maxbin; i++)
					iHisto(i - minbin) = data(i);

				smoothUntilBimodal(iHisto);

				//Image<double> tHisto(iHisto.pixelCount());
				//while (!bimodalTest(iHisto)) {
				//	//smooth with a 3 point running mean filter
				//	for (int i = 1; i < 255; i++)
				//		tHisto(i) = (iHisto(i - 1) + iHisto(i) + iHisto(i + 1)) / 3;
				//	tHisto(0) = (iHisto(0) + iHisto(1)) / 3; //0 outside
				//	tHisto(255) = (iHisto(254) + iHisto(255)) / 3; //0 outside
				//	System.arraycopy(tHisto, 0, iHisto, 0, iHisto.pixelCount());
				//	iter++;
				//	if (iter > 10000)
				//		throw ITLException("Minimum threshold not found after 10000 iterations.");
				//}

				// The threshold is the minimum between the two peaks.
				int threshold = -1;
				for (int i = 1; i < data.pixelCount() - 1; i++) {
					if (iHisto(i - 1) > iHisto(i) && iHisto(i + 1) >= iHisto(i)) {
						threshold = i;
						break;
					}
				}
				
				//return (double)threshold;
				return (double)threshold + minbin;
			}


			double isoData(const Image<double>& data)
			{
				// Also called intermeans
				// Iterative procedure based on the isodata algorithm [T.W. Ridler, S. Calvard, Picture 
				// thresholding using an iterative selection method, IEEE Trans. System, Man and 
				// Cybernetics, SMC-8 (1978) 630-632.] 
				// The procedure divides the image into objects and background by taking an initial threshold,
				// then the averages of the pixels at or below the threshold and pixels above are computed. 
				// The averages of those two values are computed, the threshold is incremented and the 
				// process is repeated until the threshold is larger than the composite average. That is,
				//  threshold = (average background + average objects)/2
				// The code in ImageJ that implements this function is the getAutoThreshold() method in the ImageProcessor class. 
				//
				// From: Tim Morris (dtm@ap.co.umist.ac.uk)
				// Subject: Re: Thresholding method?
				// posted to sci.image.processing on 1996/06/24
				// The algorithm implemented in NIH Image sets the threshold as that grey
				// value, G, for which the average of the averages of the grey values
				// below and above G is equal to G. It does this by initialising G to the
				// lowest sensible value and iterating:

				// L = the average grey value of pixels with intensities < G
				// H = the average grey value of pixels with intensities > G
				// is G = (L + H)/2?
				// yes => exit
				// no => increment G and repeat
				//

				coord_t g = 0;
				for (coord_t i = 1; i < data.pixelCount(); i++) {
					if (data(i) > 0) {
						g = i + 1;
						break;
					}
				}

				while (true) {
					double l = 0;
					double totl = 0;
					for (coord_t i = 0; i < g; i++) {
						totl += data(i);
						l += (data(i) * i);
					}

					double h = 0;
					double toth = 0;
					for (coord_t i = g + 1; i < data.pixelCount(); i++) {
						toth += data(i);
						h += data(i) * i;
					}

					if (totl > 0 && toth > 0) {
						l /= totl;
						h /= toth;
						if (g == (coord_t)round((l + h) / 2.0))
							break;
					}
					g++;
					if (g > data.pixelCount())
						return 0;
				}
				return (double)g;
			}





			double li(const Image<double>& data) {
				// Implements Li's Minimum Cross Entropy thresholding method
				// This implementation is based on the iterative version (Ref. 2) of the algorithm.
				// 1) Li C.H. and Lee C.K. (1993) "Minimum Cross Entropy Thresholding" 
				//    Pattern Recognition, 26(4): 617-625
				// 2) Li C.H. and Tam P.K.S. (1998) "An Iterative Algorithm for Minimum 
				//    Cross Entropy Thresholding"Pattern Recognition Letters, 18(8): 771-776
				// 3) Sezgin M. and Sankur B. (2004) "Survey over Image Thresholding 
				//    Techniques and Quantitative Performance Evaluation" Journal of 
				//    Electronic Imaging, 13(1): 146-165 
				//    http://citeseer.ist.psu.edu/sezgin04survey.html
				// Ported to ImageJ plugin by G.Landini from E Celebi's fourier_0.8 routines


				// threshold tolerance
				double tolerance = 0.5;
				double num_pixels = 0;
				for (int ih = 0; ih < data.pixelCount(); ih++)
					num_pixels += data(ih);

				/* Calculate the mean gray-level */
				double mean = 0.0; /* mean gray-level in the image */
				for (int ih = 0 + 1; ih < data.pixelCount(); ih++) //0 + 1?
					mean += (double)ih * data(ih);
				mean /= num_pixels;

				/* Initial estimate */
				double new_thresh = mean;
				double old_thresh;
				coord_t threshold;
				do
				{
					old_thresh = new_thresh;
					threshold = (coord_t)round(old_thresh);   /* range */

					/* Calculate the means of background and object pixels */
					/* Background */
					double sum_back = 0;
					double num_back = 0;
					for (coord_t ih = 0; ih <= threshold; ih++) {
						sum_back += ih * data(ih);
						num_back += data(ih);
					}
					/* mean of the background pixels at a given threshold */
					double mean_back = (num_back == 0 ? 0.0 : (sum_back / (double)num_back));

					/* Object */
					double sum_obj = 0;
					double num_obj = 0;
					for (coord_t ih = threshold + 1; ih < data.pixelCount(); ih++) {
						sum_obj += (double)ih * data(ih);
						num_obj += data(ih);
					}
					/* mean of the object pixels at a given threshold */
					double mean_obj = (num_obj == 0 ? 0.0 : (sum_obj / (double)num_obj));

					/* Calculate the new threshold: Equation (7) in Ref. 2 */
					//new_thresh = simple_round ( ( mean_back - mean_obj ) / ( log ( mean_back ) - log ( mean_obj ) ) );
					//simple_round ( double x ) {
					// return ( int ) ( IS_NEG ( x ) ? x - .5 : x + .5 );
					//}
					//
					//#define IS_NEG( x ) ( ( x ) < -DBL_EPSILON ) 
					//DBL_EPSILON = 2.220446049250313E-16
					double temp = (mean_back - mean_obj) / (log(mean_back) - log(mean_obj));

					if (temp < -2.220446049250313E-16)
						new_thresh = (int)(temp - 0.5);
					else
						new_thresh = (int)(temp + 0.5);
					/*  Stop the iterations when the difference between the
					new and old threshold values is less than the tolerance */
				} while (abs(new_thresh - old_thresh) > tolerance);

				return (double)threshold;
			}

			double maxEntropy(const Image<double>& data) {
				// Implements Kapur-Sahoo-Wong (Maximum Entropy) thresholding method
				// Kapur J.N., Sahoo P.K., and Wong A.K.C. (1985) "A New Method for
				// Gray-Level Picture Thresholding Using the Entropy of the Histogram"
				// Graphical Models and Image Processing, 29(3): 273-285
				// M. Emre Celebi
				// 06.15.2007
				// Ported to ImageJ plugin by G.Landini from E Celebi's fourier_0.8 routines

				Image<double> norm_histo(data.pixelCount()); /* normalized histogram */
				Image<double> P1(data.pixelCount()); /* cumulative normalized histogram */
				Image<double> P2(data.pixelCount());

				double total = 0;
				for (coord_t ih = 0; ih < data.pixelCount(); ih++)
					total += data(ih);

				for (coord_t ih = 0; ih < data.pixelCount(); ih++)
					norm_histo(ih) = data(ih) / total;

				P1(0) = norm_histo(0);
				P2(0) = 1.0 - P1(0);
				for (coord_t ih = 1; ih < data.pixelCount(); ih++) {
					P1(ih) = P1(ih - 1) + norm_histo(ih);
					P2(ih) = 1.0 - P1(ih);
				}

				/* Determine the first non-zero bin */
				coord_t first_bin = 0;
				for (coord_t ih = 0; ih < data.pixelCount(); ih++) {
					if (!(abs(P1(ih)) < 2.220446049250313E-16)) {
						first_bin = ih;
						break;
					}
				}

				/* Determine the last non-zero bin */
				coord_t last_bin = data.pixelCount() - 1;
				for (coord_t ih = data.pixelCount() - 1; ih >= first_bin; ih--) {
					if (!(abs(P2(ih)) < 2.220446049250313E-16)) {
						last_bin = ih;
						break;
					}
				}

				// Calculate the total entropy each gray-level
				// and find the threshold that maximizes it 
				double max_ent = numeric_limits<double>::lowest();
				coord_t threshold = 0;
				for (coord_t it = first_bin; it <= last_bin; it++) {
					/* Entropy of the background pixels */
					double ent_back = 0.0;
					for (coord_t ih = 0; ih <= it; ih++) {
						if (data(ih) != 0) {
							ent_back -= (norm_histo(ih) / P1(it)) * log(norm_histo(ih) / P1(it));
						}
					}

					/* Entropy of the object pixels */
					double ent_obj = 0.0;
					for (coord_t ih = it + 1; ih < data.pixelCount(); ih++) {
						if (data(ih) != 0) {
							ent_obj -= (norm_histo(ih) / P2(it)) * log(norm_histo(ih) / P2(it));
						}
					}

					/* Total entropy */
					double tot_ent = ent_back + ent_obj;

					if (max_ent < tot_ent) {
						max_ent = tot_ent;
						threshold = it;
					}
				}
				return (double)threshold;
			}

			//double mean(const Image<double>& data, double shift)
			//{
			//	// Defaults: arg0 = c = 0
			//	setDefault(shift, 0);

			//	// C. A. Glasbey, "An analysis of histogram-based thresholding algorithms,"
			//	// CVGIP: Graphical Models and Image Processing, vol. 55, pp. 532-537, 1993.
			//	//
			//	// The threshold is the mean of the greyscale data

			//	double tot = 0;
			//	double sum = 0;
			//	for (coord_t i = 0; i < data.pixelCount(); i++) {
			//		tot += data(i);
			//		sum += (double)i * data(i);
			//	}
			//	double threshold = sum / tot;

			//	return threshold - shift;
			//}




			double A(const Image<double>& y, coord_t j) {
				if (j >= y.pixelCount())
					j = y.pixelCount() - 1;
				double x = 0;
				for (coord_t i = 0; i <= j; i++)
					x += y(i);
				return x;
			}

			double B(const Image<double>& y, coord_t j) {
				if (j >= y.pixelCount())
					j = y.pixelCount() - 1;
				double x = 0;
				for (coord_t i = 0; i <= j; i++)
					x += i * y(i);
				return x;
			}

			double C(const Image<double>& y, coord_t j) {
				if (j >= y.pixelCount())
					j = y.pixelCount() - 1;
				double x = 0;
				for (coord_t i = 0; i <= j; i++)
					x += i * i * y(i);
				return x;
			}


			double minErrorI(const Image<double>& data) {
				// Kittler and J. Illingworth, "Minimum error thresholding," Pattern Recognition, vol. 19, pp. 41-47, 1986.
				// C. A. Glasbey, "An analysis of histogram-based thresholding algorithms," CVGIP: Graphical Models and Image Processing, vol. 55, pp. 532-537, 1993.
				// Ported to ImageJ plugin by G.Landini from Antti Niemisto's Matlab code (GPL)
				// Original Matlab code Copyright (C) 2004 Antti Niemisto
				// See http://www.cs.tut.fi/~ant/histthresh/ for an excellent slide presentation
				// and the original Matlab code.

				coord_t threshold = (coord_t)std::round(itl2::mean(data)); //Initial estimate for the threshold is found with the MEAN algorithm.
				coord_t Tprev = -2;

				while (threshold != Tprev) {
					//Calculate some statistics.
					double mu = B(data, threshold) / A(data, threshold);
					double nu = (B(data, data.pixelCount() - 1) - B(data, threshold)) / (A(data, data.pixelCount() - 1) - A(data, threshold));
					double p = A(data, threshold) / A(data, data.pixelCount() - 1);
					double q = (A(data, data.pixelCount() - 1) - A(data, threshold)) / A(data, data.pixelCount() - 1);
					double sigma2 = C(data, threshold) / A(data, threshold) - (mu*mu);
					double tau2 = (C(data, data.pixelCount() - 1) - C(data, threshold)) / (A(data, data.pixelCount() - 1) - A(data, threshold)) - (nu*nu);

					//The terms of the quadratic equation to be solved.
					double w0 = 1.0 / sigma2 - 1.0 / tau2;
					double w1 = mu / sigma2 - nu / tau2;
					double w2 = (mu*mu) / sigma2 - (nu*nu) / tau2 + std::log10((sigma2*(q*q)) / (tau2*(p*p)));

					//If the next threshold would be imaginary, return with the current one.
					double sqterm = (w1*w1) - w0 * w2;
					if (sqterm < 0)
						throw ITLException("MinError(I) algorithm is not converging.");

					//The updated threshold is the integer part of the solution of the quadratic equation.
					Tprev = threshold;
					double temp = (w1 + sqrt(sqterm)) / w0;

					if (std::isnan(temp))
						threshold = Tprev;
					else
						threshold = (int)floor(temp);
				}
				return (double)threshold;
			}





			double moments(const Image<double>& data)
			{
				//  W. Tsai, "Moment-preserving thresholding: a new approach," Computer Vision,
				// Graphics, and Image Processing, vol. 29, pp. 377-393, 1985.
				// Ported to ImageJ plugin by G.Landini from the the open source project FOURIER 0.8
				// by  M. Emre Celebi , Department of Computer Science,  Louisiana State University in Shreveport
				// Shreveport, LA 71115, USA
				//  http://sourceforge.net/projects/fourier-ipal
				//  http://www.lsus.edu/faculty/~ecelebi/fourier.htm

				Image<double> histo(data.pixelCount());

				double total = 0;
				for (coord_t i = 0; i < data.pixelCount(); i++)
					total += data(i);

				for (coord_t i = 0; i < data.pixelCount(); i++)
					histo(i) = data(i) / total; //normalised histogram

				/* Calculate the first, second, and third order moments */
				double m0 = 1.0; // Zeroth moment is 1 as we just normalized the histogram.
				double m1 = 0.0;
				double m2 = 0.0;
				double m3 = 0.0;
				for (coord_t i = 0; i < data.pixelCount(); i++) {
					double di = (double)i;
					m1 += di * histo(i);
					m2 += di * di * histo(i);
					m3 += di * di * di * histo(i);
				}

				/*
				First 4 moments of the gray-level image should match the first 4 moments
				of the target binary image. This leads to 4 equalities whose solutions
				are given in the Appendix of Ref. 1
				*/
				double cd = m0 * m2 - m1 * m1;
				double c0 = (-m2 * m2 + m1 * m3) / cd;
				double c1 = (m0 * -m3 + m2 * m1) / cd;
				double z0 = 0.5 * (-c1 - sqrt(c1 * c1 - 4.0 * c0));
				double z1 = 0.5 * (-c1 + sqrt(c1 * c1 - 4.0 * c0));
				double p0 = (z1 - m1) / (z1 - z0);  /* Fraction of the object pixels in the target binary image */

				// The threshold is the gray-level closest  
				// to the p0-tile of the normalized histogram 
				double sum = 0;
				coord_t threshold = 0;
				for (coord_t i = 0; i < data.pixelCount(); i++) {
					sum += histo(i);
					if (sum > p0) {
						threshold = i;
						break;
					}
				}
				return (double)threshold;
			}



			double partialSum(const Image<double>& y, coord_t j)
			{
				double x = 0;
				for (coord_t i = 0; i <= j; i++)
					x += y(i);
				return x;
			}

			double percentile(const Image<double>& data, double foregroundFraction)
			{
				// Defaults: arg = foreground fraction = 0.5
				setDefault(foregroundFraction, 0.5);

				// W. Doyle, "Operation useful for similarity-invariant pattern recognition,"
				// Journal of the Association for Computing Machinery, vol. 9,pp. 259-267, 1962.
				// ported to ImageJ plugin by G.Landini from Antti Niemisto's Matlab code (GPL)
				// Original Matlab code Copyright (C) 2004 Antti Niemisto
				// See http://www.cs.tut.fi/~ant/histthresh/ for an excellent slide presentation
				// and the original Matlab code.

				double ptile = foregroundFraction; // default fraction of foreground pixels

				Image<double> avec(data.pixelCount());
				for (coord_t i = 0; i < data.pixelCount(); i++)
					avec(i) = 0.0;

				double total = partialSum(data, data.pixelCount());
				double temp = 1.0;
				coord_t threshold = 0;
				for (coord_t i = 0; i < data.pixelCount(); i++) {
					avec(i) = abs((partialSum(data, i) / total) - ptile);
					if (avec(i) < temp) {
						temp = avec(i);
						threshold = i;
					}
				}

				return (double)threshold;
			}





			double renyiEntropy(const Image<double>& data) {
				// Kapur J.N., Sahoo P.K., and Wong A.K.C. (1985) "A New Method for
				// Gray-Level Picture Thresholding Using the Entropy of the Histogram"
				// Graphical Models and Image Processing, 29(3): 273-285
				// M. Emre Celebi
				// 06.15.2007
				// Ported to ImageJ plugin by G.Landini from E Celebi's fourier_0.8 routines


				Image<double> norm_histo(data.pixelCount()); /* normalized histogram */
				Image<double> P1(data.pixelCount()); /* cumulative normalized histogram */
				Image<double> P2(data.pixelCount());

				double total = 0;
				for (coord_t ih = 0; ih < data.pixelCount(); ih++)
					total += data(ih);

				for (coord_t ih = 0; ih < data.pixelCount(); ih++)
					norm_histo(ih) = data(ih) / total;

				P1(0) = norm_histo(0);
				P2(0) = 1.0 - P1(0);
				for (coord_t ih = 1; ih < data.pixelCount(); ih++) {
					P1(ih) = P1(ih - 1) + norm_histo(ih);
					P2(ih) = 1.0 - P1(ih);
				}

				/* Determine the first non-zero bin */
				coord_t first_bin = 0;
				for (coord_t ih = 0; ih < data.pixelCount(); ih++) {
					if (!(abs(P1(ih)) < 2.220446049250313E-16)) {
						first_bin = ih;
						break;
					}
				}

				/* Determine the last non-zero bin */
				coord_t last_bin = data.pixelCount() - 1;
				for (coord_t ih = data.pixelCount() - 1; ih >= first_bin; ih--) {
					if (!(abs(P2(ih)) < 2.220446049250313E-16)) {
						last_bin = ih;
						break;
					}
				}

				/* Maximum Entropy Thresholding - BEGIN */
				/* ALPHA = 1.0 */
				/* Calculate the total entropy each gray-level
				and find the threshold that maximizes it
				*/
				coord_t threshold = 0; // was MIN_INT in original code, but if an empty image is processed it gives an error later on.
				double max_ent = 0.0;

				for (coord_t it = first_bin; it <= last_bin; it++) {
					/* Entropy of the background pixels */
					double ent_back = 0.0;
					for (coord_t ih = 0; ih <= it; ih++) {
						if (data(ih) != 0) {
							ent_back -= (norm_histo(ih) / P1(it)) * log(norm_histo(ih) / P1(it));
						}
					}

					/* Entropy of the object pixels */
					double ent_obj = 0.0;
					for (coord_t ih = it + 1; ih < data.pixelCount(); ih++) {
						if (data(ih) != 0) {
							ent_obj -= (norm_histo(ih) / P2(it)) * log(norm_histo(ih) / P2(it));
						}
					}

					/* Total entropy */
					double tot_ent = ent_back + ent_obj;

					// IJ.log(""+max_ent+"  "+tot_ent);

					if (max_ent < tot_ent) {
						max_ent = tot_ent;
						threshold = it;
					}
				}
				coord_t t_star2 = threshold;

				/* Maximum Entropy Thresholding - END */
				threshold = 0; //was MIN_INT in original code, but if an empty image is processed it gives an error later on.
				max_ent = 0.0;
				double alpha = 0.5; // Alpha parameter of the method
				double term = 1.0 / (1.0 - alpha);
				for (coord_t it = first_bin; it <= last_bin; it++) {
					/* Entropy of the background pixels */
					double ent_back = 0.0;
					for (coord_t ih = 0; ih <= it; ih++)
						ent_back += sqrt(norm_histo(ih) / P1(it));

					/* Entropy of the object pixels */
					double ent_obj = 0.0;
					for (coord_t ih = it + 1; ih < data.pixelCount(); ih++)
						ent_obj += sqrt(norm_histo(ih) / P2(it));

					/* Total entropy */
					double tot_ent = term * ((ent_back * ent_obj) > 0.0 ? log(ent_back * ent_obj) : 0.0);

					if (tot_ent > max_ent) {
						max_ent = tot_ent;
						threshold = it;
					}
				}

				coord_t t_star1 = threshold;

				threshold = 0; //was MIN_INT in original code, but if an empty image is processed it gives an error later on.
				max_ent = 0.0;
				alpha = 2.0;
				term = 1.0 / (1.0 - alpha);
				for (coord_t it = first_bin; it <= last_bin; it++) {
					/* Entropy of the background pixels */
					double ent_back = 0.0;
					for (coord_t ih = 0; ih <= it; ih++)
						ent_back += (norm_histo(ih) * norm_histo(ih)) / (P1(it) * P1(it));

					/* Entropy of the object pixels */
					double ent_obj = 0.0;
					for (coord_t ih = it + 1; ih < data.pixelCount(); ih++)
						ent_obj += (norm_histo(ih) * norm_histo(ih)) / (P2(it) * P2(it));

					/* Total entropy */
					double tot_ent = term * ((ent_back * ent_obj) > 0.0 ? log(ent_back * ent_obj) : 0.0);

					if (tot_ent > max_ent) {
						max_ent = tot_ent;
						threshold = it;
					}
				}

				coord_t t_star3 = threshold;

				/* Sort t_star values */
				if (t_star2 < t_star1)
					swap(t_star1, t_star2);

				if (t_star3 < t_star2)
					swap(t_star2, t_star3);

				if (t_star2 < t_star1)
					swap(t_star2, t_star1);

				/* Adjust beta values */
				int beta1, beta2, beta3;
				if (abs(t_star1 - t_star2) <= 5) {
					if (abs(t_star2 - t_star3) <= 5) {
						beta1 = 1;
						beta2 = 2;
						beta3 = 1;
					}
					else {
						beta1 = 0;
						beta2 = 1;
						beta3 = 3;
					}
				}
				else {
					if (abs(t_star2 - t_star3) <= 5) {
						beta1 = 3;
						beta2 = 1;
						beta3 = 0;
					}
					else {
						beta1 = 1;
						beta2 = 2;
						beta3 = 1;
					}
				}

				/* Determine the optimal threshold value */
				double omega = P1(t_star3) - P1(t_star1);
				double opt_threshold = (t_star1 * (P1(t_star1) + 0.25 * omega * beta1) + 0.25 * t_star2 * omega * beta2 + t_star3 * (P2(t_star3) + 0.25 * omega * beta3));

				return opt_threshold;
			}


			double shanbhag(const Image<double>& data) {
				// Shanhbag A.G. (1994) "Utilization of Information Measure as a Means of
				//  Image Thresholding" Graphical Models and Image Processing, 56(5): 414-419
				// Ported to ImageJ plugin by G.Landini from E Celebi's fourier_0.8 routines

				Image<double> norm_histo(data.pixelCount()); /* normalized histogram */
				Image<double> P1(data.pixelCount()); /* cumulative normalized histogram */
				Image<double> P2(data.pixelCount());

				double total = 0;
				for (coord_t ih = 0; ih < data.pixelCount(); ih++)
					total += data(ih);

				for (coord_t ih = 0; ih < data.pixelCount(); ih++)
					norm_histo(ih) = data(ih) / total;

				P1(0) = norm_histo(0);
				P2(0) = 1.0 - P1(0);
				for (coord_t ih = 1; ih < data.pixelCount(); ih++) {
					P1(ih) = P1(ih - 1) + norm_histo(ih);
					P2(ih) = 1.0 - P1(ih);
				}

				/* Determine the first non-zero bin */
				coord_t first_bin = 0;
				for (coord_t ih = 0; ih < data.pixelCount(); ih++) {
					if (!(abs(P1(ih)) < 2.220446049250313E-16)) {
						first_bin = ih;
						break;
					}
				}

				/* Determine the last non-zero bin */
				coord_t last_bin = data.pixelCount() - 1;
				for (coord_t ih = data.pixelCount() - 1; ih >= first_bin; ih--) {
					if (!(abs(P2(ih)) < 2.220446049250313E-16)) {
						last_bin = ih;
						break;
					}
				}

				// Calculate the total entropy each gray-level
				// and find the threshold that maximizes it 
				coord_t threshold = -1;
				double min_ent = numeric_limits<double>::max();

				for (coord_t it = first_bin; it <= last_bin; it++) {
					/* Entropy of the background pixels */
					double ent_back = 0.0;
					double term = 0.5 / P1(it);
					for (coord_t ih = 1; ih <= it; ih++) { //0+1?
						ent_back -= norm_histo(ih) * log(1.0 - term * P1(ih - 1));
					}
					ent_back *= term;

					/* Entropy of the object pixels */
					double ent_obj = 0.0;
					term = 0.5 / P2(it);
					for (coord_t ih = it + 1; ih < data.pixelCount(); ih++) {
						ent_obj -= norm_histo(ih) * log(1.0 - term * P2(ih));
					}
					ent_obj *= term;

					/* Total entropy */
					double tot_ent = abs(ent_back - ent_obj);

					if (tot_ent < min_ent) {
						min_ent = tot_ent;
						threshold = it;
					}
				}
				return (double)threshold;
			}


			double triangle(const Image<double>& histogram) {
				//  Zack, G. W., Rogers, W. E. and Latt, S. A., 1977,
				//  Automatic Measurement of Sister Chromatid Exchange Frequency,
				// Journal of Histochemistry and Cytochemistry 25 (7), pp. 741-753
				//
				//  modified from Johannes Schindelin plugin
				// 

				Image<double> data(histogram.dimensions());
				setValue(data, histogram);

				// find min and max
				coord_t min = 0, max = 0, min2 = 0;
				for (coord_t i = 0; i < data.pixelCount(); i++) {
					if (data(i) > 0) {
						min = i;
						break;
					}
				}
				if (min > 0)
					min--; // line to the (p==0) point, not to data(min)

				// The Triangle algorithm cannot tell whether the data is skewed to one side or another.
				// This causes a problem as there are 2 possible thresholds between the max and the 2 extremes
				// of the histogram.
				// Here I propose to find out to which side of the max point the data is furthest, and use that as
				//  the other extreme.
				for (coord_t i = data.pixelCount() - 1; i > 0; i--) {
					if (data(i) > 0) {
						min2 = i;
						break;
					}
				}
				if (min2 < data.pixelCount() - 1)
					min2++; // line to the (p==0) point, not to data(min)

				double dmax = 0;
				for (coord_t i = 0; i < data.pixelCount(); i++) {
					if (data(i) > dmax) {
						max = i;
						dmax = data(i);
					}
				}

				// find which is the furthest side
				//IJ.log(""+min+" "+max+" "+min2);
				bool inverted = false;
				if ((max - min) < (min2 - max)) {
					// reverse the histogram
					inverted = true;
					coord_t left = 0;    // index of leftmost element
					coord_t right = data.pixelCount() - 1; // index of rightmost element
					while (left < right) {
						// exchange the left and right elements
						swap(data(left), data(right));
						// move the bounds toward the center
						left++;
						right--;
					}
					min = data.pixelCount() - 1 - min2;
					max = data.pixelCount() - 1 - max;
				}

				if (min == max) {
					return (double)min;
				}

				// describe line by nx * x + ny * y - d = 0
				// nx is just the max frequency as the other point has freq=0
				double nx = data(max);   //-min; // data(min); //  lowest value bmin = (p=0)% in the image
				double ny = (double)(min - max);
				double d = sqrt(nx * nx + ny * ny);
				nx /= d;
				ny /= d;
				d = nx * min + ny * data(min);

				// find split point
				coord_t split = min;
				double splitDistance = 0;
				for (coord_t i = min + 1; i <= max; i++) {
					double newDistance = nx * i + ny * data(i) - d;
					if (newDistance > splitDistance) {
						split = i;
						splitDistance = newDistance;
					}
				}
				split--;

				if (inverted) {
					// The histogram might be used for something else, so let's reverse it back
					coord_t left = 0;
					coord_t right = data.pixelCount() - 1;
					while (left < right) {
						swap(data(left), data(right));
						left++;
						right--;
					}
					return (double)(data.pixelCount() - 1 - split);
				}
				else
				{
					return (double)split;
				}
			}


			double yen(const Image<double>& data) {
				// Implements Yen  thresholding method
				// 1) Yen J.C., Chang F.J., and Chang S. (1995) "A New Criterion 
				//    for Automatic Multilevel Thresholding" IEEE Trans. on Image 
				//    Processing, 4(3): 370-378
				// 2) Sezgin M. and Sankur B. (2004) "Survey over Image Thresholding 
				//    Techniques and Quantitative Performance Evaluation" Journal of 
				//    Electronic Imaging, 13(1): 146-165
				//    http://citeseer.ist.psu.edu/sezgin04survey.html
				//
				// M. Emre Celebi
				// 06.15.2007
				// Ported to ImageJ plugin by G.Landini from E Celebi's fourier_0.8 routines
				Image<double> norm_histo(data.pixelCount()); /* normalized histogram */
				Image<double> P1(data.pixelCount()); /* cumulative normalized histogram */
				Image<double> P1_sq(data.pixelCount());
				Image<double> P2_sq(data.pixelCount());

				double total = 0;
				for (coord_t ih = 0; ih < data.pixelCount(); ih++)
					total += data(ih);

				for (coord_t ih = 0; ih < data.pixelCount(); ih++)
					norm_histo(ih) = data(ih) / total;

				P1(0) = norm_histo(0);
				for (coord_t ih = 1; ih < data.pixelCount(); ih++)
					P1(ih) = P1(ih - 1) + norm_histo(ih);

				P1_sq(0) = norm_histo(0) * norm_histo(0);
				for (coord_t ih = 1; ih < data.pixelCount(); ih++)
					P1_sq(ih) = P1_sq(ih - 1) + norm_histo(ih) * norm_histo(ih);

				P2_sq(data.pixelCount() - 1) = 0.0;
				for (coord_t ih = data.pixelCount() - 2; ih >= 0; ih--)
					P2_sq(ih) = P2_sq(ih + 1) + norm_histo(ih + 1) * norm_histo(ih + 1);

				/* Find the threshold that maximizes the criterion */
				coord_t threshold = 0;
				double max_crit = numeric_limits<double>::lowest();
				for (coord_t it = 0; it < data.pixelCount(); it++) {
					double crit = -1.0 * ((P1_sq(it) * P2_sq(it)) > 0.0 ? log(P1_sq(it) * P2_sq(it)) : 0.0) + 2 * ((P1(it) * (1.0 - P1(it))) > 0.0 ? log(P1(it) * (1.0 - P1(it))) : 0.0);
					if (crit > max_crit) {
						max_crit = crit;
						threshold = it;
					}
				}
				return (double)threshold;
			}
		}
	}







	namespace tests
	{
		void autothreshold()
		{
			Image<uint16_t> img(256, 256, 129);
			raw::read(img, "../test_input_data/t1-head_256x256x129.raw");

			double th = autoThreshold(img);

			cout << "Threshold value = " << th << endl;

			raw::writed(img, "./autothreshold/otsu");

		}

		void localThreshold()
		{
			Image<uint16_t> img(256, 256, 129);
			raw::read(img, "../test_input_data/t1-head_256x256x129.raw");

			Image<uint16_t> out;
			itl2::localThreshold(img, out);

			raw::writed(out, "./autothreshold/local_otsu");
		}
	}
}