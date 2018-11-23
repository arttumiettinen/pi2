#pragma once

#include <cmath>
#include "image.h"
#include "filters.h"
#include "math/matrix3x3.h"
#include "projections.h"

using math::Vec3;
using math::Vec3d;

namespace itl2
{
	/**
	Assumes that the pixels of the given three images form vectors (xi, yi, zi) and normalizes them
	by dividing it by its norm. If norm is zero, no division is done.
	*/
	template<typename pixel_t> void normalize(Image<pixel_t>& x, Image<pixel_t>& y, Image<pixel_t>& z)
	{
		x.checkSize(y);
		x.checkSize(z);

#pragma omp parallel for if(x.pixelCount() > PARALLELIZATION_THRESHOLD)
		for (coord_t n = 0; n < x.pixelCount(); n++)
		{
			Vec3<pixel_t> v(x(n), y(n), z(n));
			v.normalize();
			x(n) = v.x;
			y(n) = v.y;
			z(n) = v.z;
		}
	}

	/**
	Assumes that the pixels of the given three images form vectors (xi, yi, zi) and calculates norm of each pixel.
	*/
	template<typename pixel_t> void norm(const Image<pixel_t>& x, const Image<pixel_t>& y, const Image<pixel_t>& z, Image<pixel_t>& output)
	{
		x.checkSize(y);
		x.checkSize(z);
		output.ensureSize(x);

#pragma omp parallel for if(x.pixelCount() > PARALLELIZATION_THRESHOLD)
		for (coord_t n = 0; n < x.pixelCount(); n++)
		{
			Vec3<pixel_t> v(x(n), y(n), z(n));
			output(n) = v.norm();
		}
	}

	/**
	Calculates outScale-normalized derivative according to Lindeberg.
	See https://en.wikipedia.org/wiki/Scale_space#Scale_selection
	@param img Image whose derivative will be calculated.
	@param df Normalized derivative is placed to this image.
	@param sigma Spatial outScale (standard deviation of Gaussian).
	@param dim1, dim2 Dimensions where the derivative should be calculated. If dim2 is negative, calculates df/dx_i where f is img and i is dim1; otherwise, calculates d^2f/(dx_i dx_j), where j is dim2.
	@param gamma Scaling exponent. Set to zero to disable scaling.
	*/
	template<typename pixel_t, typename out_t> void normalizedDerivative(const Image<pixel_t>& img, Image<out_t>& df, double sigma, coord_t dim1, coord_t dim2, double gamma)
	{
		gaussDerivative(img, df, sigma, dim1, dim2, Nearest);
		if (gamma != 0)
		{
			double m = (double)dim1;
			double n = 0;
			if (dim2 >= 0)
				n = (double)dim2;

			double scale = std::pow(sigma, (m + n)  * gamma / 2.0);
			multiply(df, scale);
		}
	}

	/**
	Calculates gradient of img.
	@param gamma Scale-space scaling exponent. Set to zero to disable scaling.
	*/
	template<typename pixel_t, typename out_t> void gradient(const Image<pixel_t>& img, Image<out_t>& dfdx, Image<out_t>& dfdy, Image<out_t>& dfdz, double sigma, double gamma = 0)
	{
		cout << "df / dx..." << endl;
		normalizedDerivative(img, dfdx, sigma, 0, -1, gamma);
		cout << "df / dy..." << endl;
		normalizedDerivative(img, dfdy, sigma, 1, -1, gamma);
		cout << "df / dz..." << endl;
		normalizedDerivative(img, dfdz, sigma, 2, -1, gamma);
	}

	/**
	Calculates Hessian of img.
	@param gamma Scale-space scaling exponent. Set to zero to disable scaling.
	*/
	template<typename pixel_t, typename out_t> void hessian(const Image<pixel_t>& img, Image<out_t>& Fxx, Image<out_t>& Fyy, Image<out_t>& Fzz, Image<out_t>& Fxy, Image<out_t>& Fxz, Image<out_t>& Fyz, double sigma, double gamma = 0)
	{
		cout << "d^2 f / dx^2..." << endl;
		normalizedDerivative(img, Fxx, sigma, 0, 0, gamma);
		cout << "d^2 f / dy^2..." << endl;
		normalizedDerivative(img, Fyy, sigma, 1, 1, gamma);
		cout << "d^2 f / dz^2..." << endl;
		normalizedDerivative(img, Fzz, sigma, 2, 2, gamma);
		cout << "d^2 f / dxdy..." << endl;
		normalizedDerivative(img, Fxy, sigma, 0, 1, gamma);
		cout << "d^2 f / dxdz..." << endl;
		normalizedDerivative(img, Fxz, sigma, 0, 2, gamma);
		cout << "d^2 f / dydz..." << endl;
		normalizedDerivative(img, Fyz, sigma, 1, 2, gamma);
	}

	/**
	Calculates gradientMagnitude of image.
	Allocates 3 temporary images of the same size and pixel type than original.
	@param gradientMagnitude Output image, can be the same than the input image.
	*/
	template<typename pixel_t> void gradientMagnitude(const Image<pixel_t>& img, Image<pixel_t>& gradientMagnitude, double derivativeSigma, double gamma = 0)
	{
		Image<pixel_t> dfdx, dfdy, dfdz;
		gradient(img, dfdx, dfdy, dfdz, derivativeSigma, gamma);
		norm(dfdx, dfdy, dfdz, gradientMagnitude);
	}

	/**
	Calculates mean curvature of surfaces defined by img(x, y, z) == 0.
	Allocates 9 temporary images of the same size than the original.
	The curvature values are valid only near surfaces of objects, elsewhere the values tend to infinity.
	See https://en.wikipedia.org/wiki/Mean_curvature#Implicit_form_of_mean_curvature for formula that is used to determine the curvature.
	(The corresponding 2D formula is https://en.wikipedia.org/wiki/Implicit_curve#Curvature)
	@param curvature Output image, can equal to input image.
	*/
	template<typename pixel_t> void meanCurvature(const Image<pixel_t>& img, double sigma, Image<pixel_t>& curvature)
	{
		/*
		// 2D only
		Image<pixel_t> Fx, Fy, Fxx, Fyy, Fxy;
		gaussDerivative(img, Fx, sigma, 0, -1, Nearest);
		gaussDerivative(img, Fy, sigma, 1, -1, Nearest);
		gaussDerivative(img, Fxx, sigma, 0, 0, Nearest);
		gaussDerivative(img, Fyy, sigma, 1, 1, Nearest);
		gaussDerivative(img, Fxy, sigma, 0, 1, Nearest);

		curvature.ensureSize(img);

		#pragma omp parallel for if(Fx.pixelCount() > PARALLELIZATION_THRESHOLD)
		for (coord_t n = 0; n < Fx.pixelCount(); n++)
		{
			curvature(n) = math::pixelRound<pixel_t>((-Fy(n) * Fy(n) * Fxx(n) + 2 * Fx(n) * Fy(n) * Fxy(n) - Fx(n) * Fx(n) * Fyy(n)) / ::pow(Fx(n) * Fx(n) + Fy(n) * Fy(n), 3.0/2.0));
		}
		*/


		Image<pixel_t> Fx, Fy, Fz, Fxx, Fyy, Fzz, Fxy, Fxz, Fyz;
		gradient(img, Fx, Fy, Fz, sigma, 0.0);
		hessian(img, Fxx, Fyy, Fzz, Fxy, Fxz, Fyz, sigma, 0.0);
		
		curvature.ensureSize(img);

		size_t counter = 0;
#pragma omp parallel for if(Fx.pixelCount() > PARALLELIZATION_THRESHOLD)
		for (coord_t n = 0; n < Fx.pixelCount(); n++)
		{
			math::Matrix3x3<pixel_t> Hess(
				Fxx(n), Fxy(n), Fxz(n),
				Fxy(n), Fyy(n), Fyz(n),
				Fxz(n), Fyz(n), Fzz(n));
			math::Vec3<pixel_t> nablaF(Fx(n), Fy(n), Fz(n));
			pixel_t nablaFNorm = nablaF.norm();

			pixel_t c = (Hess.bilinear(nablaF, nablaF) - nablaFNorm * nablaFNorm * Hess.trace()) / (2 * nablaFNorm * nablaFNorm * nablaFNorm);
			curvature(n) = math::pixelRound<pixel_t>(c);

			showThreadProgress(counter, Fx.pixelCount());
		}
	}


	/**
	Calculate quantities from the structure tensor.
	Set outputs corresponding to desired quantities to pointers to images and set other pointers to zeros.
	Creates 6 temporary images of the same size than the original.
	All output images can equal to the input image.
	@param img Original image.
	@param pl1, pl2, pl3 Eigenvalues of the structure tensor.
	@param phi1, theta1, phi2, theta2, phi3, theta3 Directions of eigenvectors in spherical coordinates.
	@param pcylindricality Cylindricality value.
	@param pplanarity Planarity value.
	@param penergy Energy.
	@param gamma Scale-space scaling exponent. Set to zero to disable.
	*/
	template<typename pixel_t> void structureTensor(Image<pixel_t>& img,
		double sigmad, double sigmat,
		Image<pixel_t>* pl1 = 0, Image<pixel_t>* pl2 = 0, Image<pixel_t>* pl3 = 0,
		Image<pixel_t>* pphi1 = 0, Image<pixel_t>* ptheta1 = 0,
		Image<pixel_t>* pphi2 = 0, Image<pixel_t>* ptheta2 = 0,
		Image<pixel_t>* pphi3 = 0, Image<pixel_t>* ptheta3 = 0,
		Image<pixel_t>* pcylindricality = 0, Image<pixel_t>* pplanarity = 0, Image<pixel_t>* penergy = 0,
		double gamma = 0
	)
	{
		if (pl1)
			pl1->ensureSize(img);
		if (pl2)
			pl2->ensureSize(img);
		if (pl3)
			pl3->ensureSize(img);
		if (pphi1)
			pphi1->ensureSize(img);
		if (pphi2)
			pphi2->ensureSize(img);
		if (pphi3)
			pphi3->ensureSize(img);
		if (ptheta1)
			ptheta1->ensureSize(img);
		if (ptheta2)
			ptheta2->ensureSize(img);
		if (ptheta3)
			ptheta3->ensureSize(img);
		if (pcylindricality)
			pcylindricality->ensureSize(img);
		if (pplanarity)
			pplanarity->ensureSize(img);
		if (penergy)
			penergy->ensureSize(img);

		Image<pixel_t> dx2;
		Image<pixel_t> dy2;
		Image<pixel_t> dz2;
		Image<pixel_t> dxdy;
		Image<pixel_t> dxdz;
		Image<pixel_t> dydz;

		gradient(img, dx2, dy2, dz2, sigmad, gamma);


		cout << "Six multiplications..." << endl;
		setValue(dxdy, dx2);
		setValue(dxdz, dx2);
		setValue(dydz, dy2);

		multiply(dxdy, dy2);
		multiply(dxdz, dz2);
		multiply(dydz, dz2);
		multiply(dx2, dx2);
		multiply(dy2, dy2);
		multiply(dz2, dz2);

		cout << "Blur 1/6..." << endl;
		gaussFilter(dx2, sigmat, Nearest);
		cout << "Blur 2/6..." << endl;
		gaussFilter(dy2, sigmat, Nearest);
		cout << "Blur 3/6..." << endl;
		gaussFilter(dz2, sigmat, Nearest);
		cout << "Blur 4/6..." << endl;
		gaussFilter(dxdy, sigmat, Nearest);
		cout << "Blur 5/6..." << endl;
		gaussFilter(dxdz, sigmat, Nearest);
		cout << "Blur 6/6..." << endl;
		gaussFilter(dydz, sigmat, Nearest);

		cout << "Solving eigenvalues and outputs..." << endl;
		size_t counter = 0;
		#pragma omp parallel for if(dx2.pixelCount() > PARALLELIZATION_THRESHOLD)
		for (coord_t n = 0; n < dx2.pixelCount(); n++)
		{
			math::Matrix3x3d ST(
				dx2(n), dxdy(n), dxdz(n),
				dxdy(n), dy2(n), dydz(n),
				dxdz(n), dydz(n), dz2(n));

			double lambda1, lambda2, lambda3;
			math::Vec3d v1, v2, v3;

			ST.eig(v1, v2, v3, lambda1, lambda2, lambda3);

			double r;
			double phi1, theta1, phi2, theta2, phi3, theta3;
			math::toSpherical(v1, r, phi1, theta1);
			math::toSpherical(v2, r, phi2, theta2);
			math::toSpherical(v3, r, phi3, theta3);

			double planarity = (lambda1 - lambda2) / lambda1;
			double cylindricality = (lambda2 - lambda3) / lambda1;
			double energy = lambda1 + lambda2 + lambda3;

			// Assign outputs
			if (pl1)
				(*pl1)(n) = math::pixelRound<pixel_t>(lambda1);
			if (pl2)
				(*pl2)(n) = math::pixelRound<pixel_t>(lambda2);
			if (pl3)
				(*pl3)(n) = math::pixelRound<pixel_t>(lambda3);
			if (pphi1)
				(*pphi1)(n) = math::pixelRound<pixel_t>(phi1);
			if (pphi2)
				(*pphi2)(n) = math::pixelRound<pixel_t>(phi2);
			if (pphi3)
				(*pphi3)(n) = math::pixelRound<pixel_t>(phi3);
			if (ptheta1)
				(*ptheta1)(n) = math::pixelRound<pixel_t>(theta1);
			if (ptheta2)
				(*ptheta2)(n) = math::pixelRound<pixel_t>(theta2);
			if (ptheta3)
				(*ptheta3)(n) = math::pixelRound<pixel_t>(theta3);
			if (pcylindricality)
				(*pcylindricality)(n) = math::pixelRound<pixel_t>(cylindricality);
			if (pplanarity)
				(*pplanarity)(n) = math::pixelRound<pixel_t>(planarity);
			if (penergy)
				(*penergy)(n) = math::pixelRound<pixel_t>(energy);


			showThreadProgress(counter, dx2.pixelCount());
		}
	}

	/**
	Calculate filters that respond to lines and tubes.
	Creates 6 temporary images of the same size than the original.
	All output images can equal to the input image.
	@param plambda123 Pointer to image where Sato's lambda123 line filter results should be saved.
	@param pV Pointer to image where Frangi's Vo line filter results should be saved.
	@param gamma Scale-space scaling exponent. Set to zero to disable.
	@param outScale The output values are scaled by this value. Pass zero to outScale to numberutils<pixel_t>::outScale().
	*/
	template<typename pixel_t> void lineFilter(Image<pixel_t>& img,
		double sigma,
		Image<pixel_t>* plambda123 = 0, double gamma23 = 1, double gamma12 = 1, double alpha_sato = 0.5,
		Image<pixel_t>* pV = 0, double c = 0.25, double alpha = 0.5, double beta = 0.5,
		double gamma = 0,
		double outScale = 0)
	{
		if (plambda123)
			plambda123->ensureSize(img);
		if (pV)
			pV->ensureSize(img);

		if (outScale == 0)
			outScale = (double)math::NumberUtils<pixel_t>::scale();
		
		typedef typename math::NumberUtils<pixel_t>::FloatType real_t;
		
		Image<real_t> Fxx, Fyy, Fzz, Fxy, Fxz, Fyz;
		hessian(img, Fxx, Fyy, Fzz, Fxy, Fxz, Fyz, sigma, gamma);
		
		size_t counter = 0;
		#pragma omp parallel for if(Fxx.pixelCount() > PARALLELIZATION_THRESHOLD)
		for (coord_t n = 0; n < Fxx.pixelCount(); n++)
		{
			math::Matrix3x3d Hess(
				Fxx(n), Fxy(n), Fxz(n),
				Fxy(n), Fyy(n), Fyz(n),
				Fxz(n), Fyz(n), Fzz(n));

			double lambda1, lambda2, lambda3;
			math::Vec3d v1, v2, v3;

			Hess.eig(v1, v2, v3, lambda1, lambda2, lambda3);

			if (plambda123)
			{
				// Calculate lambda123 from
				// Sato - Three-dimensional multi-outScale line filter for segmentation and visualization of curvilinear structures in medical images

				double lambda123;
				if (lambda3 < lambda2 && lambda2 < lambda1 && lambda1 <= 0)
					lambda123 = std::abs(lambda3) * std::pow(lambda2 / lambda3, gamma23) * std::pow(1 + lambda1 / std::abs(lambda2), gamma12);
				else if (lambda3 < lambda2 && lambda2 < 0 && 0 < lambda1 && lambda1 < std::abs(lambda2) / alpha_sato)
					lambda123 = std::abs(lambda3) * std::pow(lambda2 / lambda3, gamma23) * std::pow(1 - alpha_sato * lambda1 / std::abs(lambda2), gamma12);
				else
					lambda123 = 0;

				if (std::isnan(lambda123))
					lambda123 = 0;

				(*plambda123)(n) = math::pixelRound<pixel_t>(lambda123 * outScale);
			}

			if (pV)
			{
				// Calculate Vo from
				// Frangi - Multiscale vessel enhancement filtering

				// Re-order eigenvalues according to Frangi's sorting order: |lambda1| < |lambda2| < |lambda3|
				double al1 = std::abs(lambda1);
				double al2 = std::abs(lambda2);
				double al3 = std::abs(lambda3);
				if (al1 > al3)
				{
					swap(al1, al3);
					swap(lambda1, lambda3);
					swap(v1, v3);
				}
				if (al1 > al2)
				{
					swap(al1, al2);
					swap(lambda1, lambda2);
					swap(v1, v2);
				}
				if (al2 > al3)
				{
					swap(al2, al3);
					swap(lambda2, lambda3);
					swap(v2, v3);
				}

				double Rb = al1 / sqrt(al2 * al3);
				double Ra = al2 / al3;
				double S = sqrt(al1 * al1 + al2 * al2 + al3 * al3);
				double Vo;
				if (lambda2 > 0 || lambda3 > 0)
					Vo = 0;
				else
					Vo = (1 - exp(-(Ra * Ra) / (2 * alpha * alpha))) * exp(-(Rb * Rb) / (2 * beta * beta)) * (1 - exp(-(S * S) / (2 * c * c)));

				if (std::isnan(Vo))
					Vo = 0;

				(*pV)(n) = math::pixelRound<pixel_t>(Vo * outScale);
			}
		}
	}

	namespace tests
	{
		void curvature();
		void structureTensor();
		void lineFilter();
	}

}
