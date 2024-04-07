#pragma once

#include <random>
#include <functional>
#include <chrono>

#include "image.h"


namespace itl2
{
	/**
	Adds additive noise to pixels of the image. Noise values are drawn from the generator given as an argument.
	@param generator The () operator of the generator must return random noise value added to the image pixel.
	@param region Region where the noise is to be added.
	*/
	template<typename pixel_t, class F> void noise(Image<pixel_t>& img, F generator, const AABoxc& region)
	{

		// Not parallelized as generator is usually not capable to parallel processing.
		forAllInBox(region, [&] (coord_t x, coord_t y, coord_t z)
			{
				img(x, y, z) = pixelRound<pixel_t>(img(x, y, z) + generator());
			});
	}

	/**
	Adds additive noise to pixels of the image. Noise values are drawn from the generator given as an argument.
	@param generator The () operator of the generator must return random noise value added to the image pixel.
	*/
	template<typename pixel_t, class F> void noise(Image<pixel_t>& img, F generator)
	{
		noise(img, generator, img.bounds());
	}


	/**
	Add additive Gaussian noise to the image.
	@param img Image where the noise is added to.
	@param mean Mean of the normal distribution where noise samples are drawn from.
	@param stddev Standard deviation of the distribution where the noise samples are drawn. If set to zero, 10 % of typical value range of the pixel data type is used.
	@param seed Random seed. Set to zero to use time-based seed.
	*/
	template<typename pixel_t> void noise(Image<pixel_t>& img, const AABoxc& region, double mean = 0, double stddev = 0, unsigned int seed = 0)
	{
		if (seed == 0)
		{
			seed = (unsigned int)std::chrono::system_clock::now().time_since_epoch().count();
		}

		if (stddev == 0)
		{
			stddev = 0.1 * NumberUtils<pixel_t>::scale();
		}

		std::mt19937 gen(seed);
		std::normal_distribution<double> dist(mean, stddev);
		auto dice = std::bind(dist, gen);
		noise(img, dice, region);
	}

	/**
	Add additive Gaussian noise to the image.
	@param img Image where the noise is added to.
	@param mean Mean of the normal distribution where noise samples are drawn from.
	@param stddev Standard deviation of the distribution where the noise samples are drawn. If set to zero, 10 % of typical value range of the pixel data type is used.
	@param seed Random seed. Set to zero to use time-based seed.
	*/
	template<typename pixel_t> void noise(Image<pixel_t>& img, double mean = 0, double stddev = 0, unsigned int seed = 0)
	{
		noise(img, img.bounds(), mean, stddev, seed);
	}
}