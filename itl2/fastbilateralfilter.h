#pragma once

#include "image.h"
#include "boundarycondition.h"
#include "projections.h"
#include "filters.h"
#include "math/mathutils.h"

#include <random>
#include <chrono>

namespace itl2
{

	//namespace internals
	//{
	//	/**
	//	Calculates estimate of suitable approximation degree for bilateralFilterGaussPolynomial.
	//	*/
	//	inline coord_t calculateDegree(double rangeSigma, double epsilon, double range)
	//	{
	//		if (rangeSigma >= 70)
	//			return 10;

	//		double lambda = range / rangeSigma;
	//		lambda *= lambda;

	//		double p = 1 + log(lambda);
	//		double q = -lambda - log(epsilon);
	//		double t = q / (exp(1) * lambda);
	//		double W0 = t - t * t + 3 * t*t*t / 2 - 8 * t*t*t*t / 4;
	//		double N0 = q / W0;
	//		if (rangeSigma < 30)
	//		{
	//			for (coord_t k = 1; k <= 3; k++)
	//			{
	//				N0 = N0 - (N0 * log(N0) - p * N0 - q) / (log(N0) + 1 - p);
	//			}
	//		}

	//		return (coord_t)ceil(N0);
	//	}
	//}

	///**
	//This is implementation of fast bilateral filtering using Gaussian Polynomial approximation for range kernel.
	//See
	//Chaudhury - Fast and Provably Accurate Bilateral Filtering
	//The approximation works best for large range sigmas, and requires a lot of iterations for small range sigmas.
	//TODO: Currently uses double temporary images as float temporary images seemed to not be accurate enough.
	//@param spatialSigma Standard deviation of Gaussian kernel used for spatial smoothing.
	//@param rangeSigma Standard deviation of Gaussian kernel used to avoid smoothing edges of features. Order of magnitude must be similar to difference between gray levels of background and objects.
	//@param degree Approximation degree. Typically in the range of 10-200. Pass negative value to estimate degree automatically.
	//@param tc Mean value of the interesting regions of the image.
	//@param bc Boundary condition (BoundaryCondition::Zero or Nearest).
	//*/
	//template<typename pixel_t> void bilateralFilterGaussPolynomial(Image<pixel_t>& img, double spatialSigma, double rangeSigma, double tc, coord_t degree = -1, BoundaryCondition bc = BoundaryCondition::Nearest)
	//{
	//	//double tc = (double)mean(img);
	//	//double tc = (min(img) + max(img)) / 2;
	//	//double tc = 240;
	//	//double tc = max(img) / 2;
	//	//double tc = 250;

	//	if (degree <= 0)
	//	{
	//		degree = internals::calculateDegree(rangeSigma, 0.1, 2 * tc);
	//		std::cout << "Using approximation degree " << degree << std::endl;
	//	}
	//	
	//	Image<double> G(img.dimensions(), 1);
	//	Image<double> F(img.dimensions());
	//	Image<double> P(img.dimensions(), 0);
	//	Image<double> Q(img.dimensions(), 0);
	//	//Image<double> H(img.dimensions());
	//	Image<double> Fbar(img.dimensions());

	//	for (coord_t i = 0; i < F.pixelCount(); i++)
	//	{
	//		double hi = (double)img(i) - tc;
	//		F(i) = exp(-(hi * hi) / (2 * rangeSigma * rangeSigma));
	//		//H(i) = hi / rangeSigma;
	//		Fbar(i) = F(i);
	//	}

	//	gaussFilter(Fbar, spatialSigma, bc);

	//	for (coord_t n = 1; n <= degree; n++)
	//	{
	//		std::cout << "Processing terms of degree " << n << "/" << degree << std::endl;

	//		for (coord_t i = 0; i < F.pixelCount(); i++)
	//		{
	//			Q(i) = Q(i) + G(i) * Fbar(i);
	//			//F(i) = H(i) * F(i);
	//			F(i) = (((double)img(i) - tc) / rangeSigma) * F(i);
	//			Fbar(i) = F(i);
	//		}

	//		gaussFilter(Fbar, spatialSigma, bc);

	//		for (coord_t i = 0; i < F.pixelCount(); i++)
	//		{
	//			P(i) = P(i) + G(i) * Fbar(i);
	//			//G(i) = H(i) * G(i) / (double)n;
	//			G(i) = (((double)img(i) - tc) / rangeSigma) * G(i) / (double)n;
	//		}
	//	}

	//	for (coord_t i = 0; i < F.pixelCount(); i++)
	//	{
	//		img(i) = pixelRound<pixel_t>(rangeSigma * P(i) / Q(i) + tc);
	//	}
	//}




	/**
	Sampling bilateral filter.
	Does not process all pixels in a neighbourhood but only a random subset of them.
	The idea is from
	Banterle - A Low-Memory, Straightforward and Fast Bilateral Filter Through Subsampling in Spatial Domain
	but here we take the samples from normal distribution around current pixel center and thus there's no need to calculate
	the spatial weight at all.
	@param in Input image.
	@param out Output image.
	@param spatialSigma Standard deviation of Gaussian kernel used for spatial smoothing.
	@param rangeSigma Standard deviation of Gaussian kernel used to avoid smoothing edges of features. Order of magnitude must be similar to difference between gray levels of background and objects.
	@param sampleCount Count of samples processed in each neighbourhood. Specify zero to determine sample count automatically.
	@param bc Boundary condition (BoundaryCondition::Zero or Nearest).
	@param randSeed Seed for random number generation. Pass zero to choose a seed automatically.
	*/
	template<typename pixel_t, typename out_t> void bilateralFilterSampling(const Image<pixel_t>& in, Image<out_t>& out, float32_t spatialSigma, float32_t rangeSigma, size_t sampleCount = 0, size_t randSeed = 0)
	{
		out.mustNotBe(in);
		out.ensureSize(in);

		
		if (randSeed == 0)
		{
			randSeed = (unsigned int)std::chrono::system_clock::now().time_since_epoch().count();
		}

		std::mt19937 gen((unsigned int)randSeed);
		std::normal_distribution<float32_t> d(0, spatialSigma);


		if (sampleCount <= 0)
		{
			// Determine sample count automatically
			coord_t nbRadius = (coord_t)ceil(3 * spatialSigma);
			sampleCount = 4 * (2 * nbRadius + 1) * in.dimensionality();
			if (sampleCount < 5)
				sampleCount = 5;
				
			std::cout << "Using " << sampleCount << " samples per neighbourhood." << std::endl;
		}

		// Generate random patterns so that no random generation needs to be done in each thread
		std::vector<std::vector<Vec3c> > patterns;
		for (size_t i = 0; i < 1024; i++)
		{
			std::vector<Vec3c> pattern;
			for (size_t n = 0; n < sampleCount; n++)
			{
				pattern.push_back(round(Vec3f(d(gen), d(gen), d(gen))));
			}
			patterns.push_back(pattern);
		}

		// Process all pixels in the image
		size_t counter = 0;
		#pragma omp parallel for if(!omp_in_parallel() && in.pixelCount() >= PARALLELIZATION_THRESHOLD)
		for (coord_t cz = 0; cz < in.depth(); cz++)
		{
			for (coord_t cy = 0; cy < in.height(); cy++)
			{
				for (coord_t cx = 0; cx < in.width(); cx++)
				{
					float32_t sum = 0.0;
					float32_t wsum = 0.0;

					Vec3c center(cx, cy, cz);
					float32_t centerVal = (float32_t)in(center);

					// Select pattern randomly and process all pixels in it.
					coord_t patternIndex = randc(patterns.size());
					const std::vector<Vec3c>& pattern = patterns[patternIndex];
					for (const Vec3c& p : pattern)
					{
						Vec3c wp = center + p;

						if (in.isInImage(wp))
						{
							float32_t pix = (float32_t)in(wp);

							//float32_t r2 = (wp - center).normSquared();
							float32_t r2 = p.normSquared();
							float32_t c = pix - centerVal;

							// Calculate only range weight, spatial weight is automatically accounted for by the sampling pattern.
							float32_t w = (1 / (rangeSigma * sqrtf(2 * PIf)) * ::expf(-(c * c) / (2 * rangeSigma * rangeSigma)));
							// NOTE: Using constant weight results in Gaussian filtering. There we require a lot more points to make a good-looking result.
							//float32_t w = 1;

							sum += pix * w;
							wsum += w;
						}
					}

					out(cx, cy, cz) = pixelRound<out_t>(sum / wsum);
				}
			}

			showThreadProgress(counter, in.depth());
		}
	}



	namespace tests
	{
		//void fastBilateralGaussPolynomial();
		void fastBilateralSampling();
	}
}
