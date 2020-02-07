#pragma once

#include <random>
#include <functional>
#include <chrono>

#include "image.h"


namespace itl2
{
	/**
	Adds noise to pixels of the image. Noise values are drawn from the generator given as an argument.
	The () operator of the generator must return random noise value added to the image pixel.
	*/
	template<typename pixel_t, class F> void noise(Image<pixel_t>& img, F generator)
	{

		// Not parallelized as generator is usually not capable to parallel processing.
		for (coord_t n = 0; n < img.pixelCount(); n++)
		{
			img(n) = pixelRound<pixel_t>(img(n) + generator());

			// Showing progress info here would induce more processing than is done in the whole loop.
		}
	}

	/**
	Add Gaussian noise to the image.
	@param img Image where the noise is added to.
	@param mean Mean of the normal distribution where noise samples are drawn from.
	@param stddev Standard deviation of the distribution where the noise samples are drawn.
	@param seed Random seed. Set to zero to use time-based seed.
	*/
	template<typename pixel_t> void noise(Image<pixel_t>& img, double mean = 0, double stddev = 25, unsigned int seed = 0)
	{
		if (seed == 0)
		{
			seed = (unsigned int)std::chrono::system_clock::now().time_since_epoch().count();
		}

		std::mt19937 gen(seed);
		std::normal_distribution<double> dist(mean, stddev);
		auto dice = std::bind(dist, gen);
		noise(img, dice);
	}
}