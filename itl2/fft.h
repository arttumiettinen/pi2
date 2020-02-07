#pragma once


#include <complex>
#include <omp.h>

#include "image.h"
#include "math/vec3.h"

namespace itl2
{
	/*
	Initializes FFTW if it has not been initialized.
	There is no need to call this function unless you make calls to FFTW library yourself.
	*/
	void initFFTW();

	/**
	Calculates Discrete Cosine Transform of the input image.
	Calculates 1D DCT if img is 1-dimensional, 2D DCT if img is 2-dimensional etc.
	The transformation is calculated in-place.
	*/
	void dct(Image<float32_t>& img);

	/**
	Calculates inverse Discrete Cosine Transform of the input image.
	Calculates 1D DCT if img is 1-dimensional, 2D DCT if img is 2-dimensional etc.
	The transformation is calculated in-place.
	*/
	void idct(Image<float32_t>& img);

	/**
	Calculates FFT of the input image and places it to the output image.
	Initializes output image to correct size.
	Calculates 1D FFT if img is 1-dimensional, 2D FFT if img is 2-dimensional etc.
	Input image is not modified.
	*/
	void fft(Image<float32_t>& img, Image<complex32_t>& out);

	/**
	Calculates inverse FFT of the input image and places it to the output image.
	Output image size must be set to the size of the original image where the FFT was calculated from.
	Calculates 1D FFT if img is 1-dimensional, 2D FFT if img is 2-dimensional etc.
	Input image data is used as a temporary buffer.
	*/
	void ifft(Image<complex32_t>& img, Image<float32_t>& out);

	/**
	Gaussian filtering.
	Performs filtering in Fourier domain.
	Produces the same results than ImageJ Gaussian Blur 3D, except at image edges.
	@param zeroEdges Set to true to avoid having remnants of left side of the image at the right side in the filtered image etc. Causes zero edges to the image.
	*/
	void gaussFilterFFT(Image<float32_t>& img, Vec3d sigma, bool zeroEdges = true);

	/**
	Gaussian filtering.
	Performs filtering in Fourier domain.
	Produces the same results than ImageJ Gaussian Blur 3D, except at image edges.
	@param zeroEdges Set to true to avoid having remnants of left side of the image at the right side in the filtered image etc. Causes zero edges to the image.
	*/
	inline void gaussFilterFFT(Image<float32_t>& img, double sigma, bool zeroEdges = true)
	{
		gaussFilterFFT(img, Vec3d(sigma, sigma, sigma), zeroEdges);
	}

	/**
	Bandpass filtering.
	Performs filtering in Fourier domain.
	@param min_size Passband minimum.
	@param max_size Passband maximum.
	@param zeroEdges Set to true to avoid having remnants of left side of the image at the right side in the filtered image etc. Causes zero edges to the image.
	*/
	void bandpassFilter(Image<float32_t>& img, double min_size, double max_size, bool zeroEdges = true);

	/**
	Calculates power spectrum of Fourier transform.
	@param img Fourier transform image.
	@param pwr Power spectrum will be stored in this image.
	*/
	void powerSpectrum(const Image<complex32_t>& img, Image<float32_t>& pwr);

	/**
	Calculates shift between the two images with phase correlation method.
	@param img1 Reference image. This image is used as temporary storage so it will be modified.
	@param img2 Shifted image. This image is not modified.
	@param maxShift Maximal shift that is to be recognized.
	@param goodness Estimate of goodness of fit between img1 and shifted img2.
	@return Shift between img1 and img2.
	*/
	Vec3d phaseCorrelation(Image<float32_t>& img1, Image<float32_t>& img2, const Vec3c& maxShift, double& goodness);

	namespace tests
	{
		void fourierTransformPair();
		void dctPair();
		void bandpass();
		void phaseCorrelation();
		void phaseCorrelation2();
		void modulo();
	}
}