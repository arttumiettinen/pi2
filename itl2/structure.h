#pragma once

#include <cmath>
#include "image.h"
#include "filters.h"
#include "math/matrix3x3.h"
#include "projections.h"
#include "interpolation.h"
#include "floodfill.h"
#include "iteration.h"

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
	Output image can be any of the input images.
	*/
	template<typename pixel_t, typename out_t> void norm(const Image<pixel_t>& x, const Image<pixel_t>& y, const Image<pixel_t>& z, Image<out_t>& output)
	{
		x.checkSize(y);
		x.checkSize(z);
		output.ensureSize(x);

#pragma omp parallel for if(x.pixelCount() > PARALLELIZATION_THRESHOLD)
		for (coord_t n = 0; n < x.pixelCount(); n++)
		{
			Vec3<pixel_t> v(x(n), y(n), z(n));
			output(n) = pixelRound<out_t>(v.norm());
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
		gaussDerivative(img, df, sigma, dim1, dim2, BoundaryCondition::Nearest);
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
		dfdx.ensureSize(img);	// Allocate memory here to fail fast if there is not enough memory.
		dfdy.ensureSize(img);
		dfdz.ensureSize(img);

		//if(showProgressInfo)
		//	std::cout << "df / dx..." << std::endl;
		normalizedDerivative(img, dfdx, sigma, 0, -1, gamma);

		//if (showProgressInfo)
		//	std::cout << "df / dy..." << std::endl;
		normalizedDerivative(img, dfdy, sigma, 1, -1, gamma);

		//if (showProgressInfo)
		//	std::cout << "df / dz..." << std::endl;
		normalizedDerivative(img, dfdz, sigma, 2, -1, gamma);
	}

	/**
	Calculates Hessian of img.
	@param gamma Scale-space scaling exponent. Set to zero to disable scaling.
	*/
	template<typename pixel_t, typename out_t> void hessian(const Image<pixel_t>& img, Image<out_t>& Fxx, Image<out_t>& Fyy, Image<out_t>& Fzz, Image<out_t>& Fxy, Image<out_t>& Fxz, Image<out_t>& Fyz, double sigma, double gamma = 0)
	{
		Fxx.ensureSize(img);	// Allocate memory here to fail fast if there is not enough memory.
		Fyy.ensureSize(img);
		Fzz.ensureSize(img);
		Fxy.ensureSize(img);
		Fxz.ensureSize(img);
		Fyz.ensureSize(img);

		std::cout << "d^2 f / dx^2..." << std::endl;
		normalizedDerivative(img, Fxx, sigma, 0, 0, gamma);
		std::cout << "d^2 f / dy^2..." << std::endl;
		normalizedDerivative(img, Fyy, sigma, 1, 1, gamma);
		std::cout << "d^2 f / dz^2..." << std::endl;
		normalizedDerivative(img, Fzz, sigma, 2, 2, gamma);
		std::cout << "d^2 f / dxdy..." << std::endl;
		normalizedDerivative(img, Fxy, sigma, 0, 1, gamma);
		std::cout << "d^2 f / dxdz..." << std::endl;
		normalizedDerivative(img, Fxz, sigma, 0, 2, gamma);
		std::cout << "d^2 f / dydz..." << std::endl;
		normalizedDerivative(img, Fyz, sigma, 1, 2, gamma);
	}

	/**
	Calculates gradientMagnitude of image.
	Allocates 3 temporary images of the same size and pixel type than original.
	@param gradientMagnitude Output image, can be the same than the input image.
	*/
	template<typename pixel_t, typename out_t> void gradientMagnitude(const Image<pixel_t>& img, Image<out_t>& gradientMagnitude, double derivativeSigma, double gamma = 0)
	{
		Image<typename NumberUtils<pixel_t>::FloatType> dfdx, dfdy, dfdz;
		gradient(img, dfdx, dfdy, dfdz, derivativeSigma, gamma);
		norm(dfdx, dfdy, dfdz, gradientMagnitude);
	}

	/**
	Calculates mean curvature of surfaces defined by img(x, y, z) == 0 multiplied by the norm of the gradient of img.
	Without multiplication by ||nabla f||, curvature of f is very hard to interpret as the values are only valid near the zero contour from
	image processing point of view. See also surfaceCurvature functionality, where the principal curvatyres of surfaces are being calculated.
	See https://en.wikipedia.org/wiki/Mean_curvature#Implicit_form_of_mean_curvature for formulas that can be used to determine the curvature.
	(The corresponding 2D formula is at https://en.wikipedia.org/wiki/Implicit_curve#Curvature)
	Note that this function might produce invalid results on surfaces of structures whose thickness is of the same order than sigma.
	TODO: Gaussian curvature can be calculated with similar means.
	@param curvature Output image, can equal to input image.
	*/
	template<typename pixel_t, typename out_t> void meanCurvature(const Image<pixel_t>& img, double sigma, Image<out_t>& curvature)
	{
		curvature.ensureSize(img);

		// This calculates curvature from gradient and Hessian
		//if (img.dimensionality() == 2)
		//{
		//	// 2D only
		//	Image<pixel_t> Fx, Fy, Fxx, Fyy, Fxy;
		//	gaussDerivative(img, Fx, sigma, 0, -1, BoundaryCondition::Nearest);
		//	gaussDerivative(img, Fy, sigma, 1, -1, BoundaryCondition::Nearest);
		//	gaussDerivative(img, Fxx, sigma, 0, 0, BoundaryCondition::Nearest);
		//	gaussDerivative(img, Fyy, sigma, 1, 1, BoundaryCondition::Nearest);
		//	gaussDerivative(img, Fxy, sigma, 0, 1, BoundaryCondition::Nearest);

		//	#pragma omp parallel for if(Fx.pixelCount() > PARALLELIZATION_THRESHOLD)
		//	for (coord_t n = 0; n < Fx.pixelCount(); n++)
		//	{
		//		pixel_t L = Vec2<pixel_t>(Fx(n), Fy(n)).norm();
		//		curvature(n) = pixelRound<pixel_t>(L * (-Fy(n) * Fy(n) * Fxx(n) + 2 * Fx(n) * Fy(n) * Fxy(n) - Fx(n) * Fx(n) * Fyy(n)) / ::pow((pixel_t)1e-8 + Fx(n) * Fx(n) + Fy(n) * Fy(n), 3.0 / 2.0));
		//	}

		//}
		//else
		//{
		//	// 2D or 3D, here used for the 3D case only.

		//	Image<pixel_t> Fx, Fy, Fz, Fxx, Fyy, Fzz, Fxy, Fxz, Fyz;
		//	Fx.ensureSize(img);	// Allocate memory here to fail fast if there is not enough memory.
		//	Fy.ensureSize(img);
		//	Fz.ensureSize(img);
		//	Fxx.ensureSize(img);
		//	Fyy.ensureSize(img);
		//	Fzz.ensureSize(img);
		//	Fxy.ensureSize(img);
		//	Fxz.ensureSize(img);
		//	Fyz.ensureSize(img);

		//	gradient(img, Fx, Fy, Fz, sigma, 0.0);
		//	hessian(img, Fxx, Fyy, Fzz, Fxy, Fxz, Fyz, sigma, 0.0);

		//	size_t counter = 0;
		//	#pragma omp parallel for if(Fx.pixelCount() > PARALLELIZATION_THRESHOLD)
		//	for (coord_t n = 0; n < Fx.pixelCount(); n++)
		//	{
		//		Matrix3x3<pixel_t> Hess(
		//			Fxx(n), Fxy(n), Fxz(n),
		//			Fxy(n), Fyy(n), Fyz(n),
		//			Fxz(n), Fyz(n), Fzz(n));
		//		Vec3<pixel_t> nablaF(Fx(n), Fy(n), Fz(n));
		//		pixel_t nablaFNorm = nablaF.norm();

		//		pixel_t c = nablaFNorm * (Hess.bilinear(nablaF, nablaF) - nablaFNorm * nablaFNorm * Hess.trace()) / ((pixel_t)1e-8 + 2 * nablaFNorm * nablaFNorm * nablaFNorm);
		//		curvature(n) = pixelRound<pixel_t>(c);

		//		showThreadProgress(counter, Fx.pixelCount());
		//	}
		//}


		// This calculates curvature as divergence of unit normal (times nabla f norm), [-0.5 nabla . (nabla f / ||nabla f||)] * ||nabla f||
		// This version requires less temporaries than the Hessian method above.
		// TODO: At least two temporaries (L and phix) can still be easily removed.

		Image<float32_t> phix, phiy, phiz, L;

		gradient(img, phix, phiy, phiz, sigma, 0);
		norm(phix, phiy, phiz, L);

		add(L, 1e-8); // Add small number to stabilize regions where L is zero.
		divide(phix, L);
		divide(phiy, L);
		divide(phiz, L);

		gaussDerivative(phix, sigma, 0, -1, BoundaryCondition::Nearest);
		gaussDerivative(phiy, sigma, 1, -1, BoundaryCondition::Nearest);
		gaussDerivative(phiz, sigma, 2, -1, BoundaryCondition::Nearest);

		setValue(curvature, phix);
		add(curvature, phiy);
		add(curvature, phiz);

		multiply(curvature, L);
		multiply(curvature, -0.5);
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
	template<typename pixel_t> void structureTensor(const Image<pixel_t>& img,
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

		dx2.ensureSize(img);
		dy2.ensureSize(img);
		dz2.ensureSize(img);
		dxdy.ensureSize(img);
		dxdz.ensureSize(img);
		dydz.ensureSize(img);

		gradient(img, dx2, dy2, dz2, sigmad, gamma);


		std::cout << "Six multiplications..." << std::endl;
		setValue(dxdy, dx2);
		setValue(dxdz, dx2);
		setValue(dydz, dy2);

		multiply(dxdy, dy2);
		multiply(dxdz, dz2);
		multiply(dydz, dz2);
		multiply(dx2, dx2);
		multiply(dy2, dy2);
		multiply(dz2, dz2);

		std::cout << "Blur 1/6..." << std::endl;
		gaussFilter(dx2, sigmat, BoundaryCondition::Nearest);
		std::cout << "Blur 2/6..." << std::endl;
		gaussFilter(dy2, sigmat, BoundaryCondition::Nearest);
		std::cout << "Blur 3/6..." << std::endl;
		gaussFilter(dz2, sigmat, BoundaryCondition::Nearest);
		std::cout << "Blur 4/6..." << std::endl;
		gaussFilter(dxdy, sigmat, BoundaryCondition::Nearest);
		std::cout << "Blur 5/6..." << std::endl;
		gaussFilter(dxdz, sigmat, BoundaryCondition::Nearest);
		std::cout << "Blur 6/6..." << std::endl;
		gaussFilter(dydz, sigmat, BoundaryCondition::Nearest);

		std::cout << "Solving eigenvalues and outputs..." << std::endl;
		#pragma omp parallel for if(dx2.pixelCount() > PARALLELIZATION_THRESHOLD)
		for (coord_t n = 0; n < dx2.pixelCount(); n++)
		{
		// Progress reporting has been removed as it causes a lot of overhead when there are many threads.
		//size_t counter = 0;
		//#pragma omp parallel for if(dx2.depth() > PARALLELIZATION_THRESHOLD)
		//for (coord_t z = 0; z < dx2.depth(); z++)
		//{
		//	coord_t n0 = dx2.getLinearIndex(0, 0, z);
		//	coord_t n1 = dx2.getLinearIndex(0, 0, z + 1);
		//	for (coord_t n = n0; n < n1; n++)
		//	{
				Matrix3x3d ST(
					dx2(n), dxdy(n), dxdz(n),
					dxdy(n), dy2(n), dydz(n),
					dxdz(n), dydz(n), dz2(n));

				double lambda1, lambda2, lambda3;
				Vec3d v1, v2, v3;

				ST.eigsym(v1, v2, v3, lambda1, lambda2, lambda3);

				double r;
				double phi1, theta1, phi2, theta2, phi3, theta3;
				cartesianToSpherical(v1, r, phi1, theta1);
				cartesianToSpherical(v2, r, phi2, theta2);
				cartesianToSpherical(v3, r, phi3, theta3);

				double energy = lambda1 + lambda2 + lambda3;
				double planarity = (lambda1 - lambda2) / lambda1;
				double cylindricality = ((lambda2 - lambda3) / lambda1) * energy;


				// Assign outputs
				if (pl1)
					(*pl1)(n) = pixelRound<pixel_t>(lambda1);
				if (pl2)
					(*pl2)(n) = pixelRound<pixel_t>(lambda2);
				if (pl3)
					(*pl3)(n) = pixelRound<pixel_t>(lambda3);
				if (pphi1)
					(*pphi1)(n) = pixelRound<pixel_t>(phi1);
				if (pphi2)
					(*pphi2)(n) = pixelRound<pixel_t>(phi2);
				if (pphi3)
					(*pphi3)(n) = pixelRound<pixel_t>(phi3);
				if (ptheta1)
					(*ptheta1)(n) = pixelRound<pixel_t>(theta1);
				if (ptheta2)
					(*ptheta2)(n) = pixelRound<pixel_t>(theta2);
				if (ptheta3)
					(*ptheta3)(n) = pixelRound<pixel_t>(theta3);
				if (pcylindricality)
					(*pcylindricality)(n) = pixelRound<pixel_t>(cylindricality);
				if (pplanarity)
					(*pplanarity)(n) = pixelRound<pixel_t>(planarity);
				if (penergy)
					(*penergy)(n) = pixelRound<pixel_t>(energy);
			//}
			//showThreadProgress(counter, dx2.depth());
		}
	}

	/**
	Calculates orientation of cylindrical or tubular structures.
	@param img At input, must contain the original image, possibly binary. At output, contains 'energy' as defined in the structure tensor documentation.
	@param phi The azimuthal angle of orientation direction will be stored in this image. The angle is given in radians and measured from positive $x$-axis towards positive $y$-axis and is given in range $[-\pi, \pi]$.
	@param theta The polar angle of orientation direction will be stored in this image. The angle is given in radians and measured from positive $z$-axis towards $xy$-plane. The values are in range $[0, \pi]$.
	@param derSigma Scale parameter. Set to the preferred scale of edges that define the cylinders. Derivatives required in the stucture tensor are calculated using convolutions with derivative of a Gaussian function, and this parameter defines the standard deviation of the Gaussian.
	@param smoothSigma The structure tensor is smoothed by convolution with a Gaussian. This parameter defines the standard deviation of the smoothing Gaussian.
	*/
	inline void cylinderOrientation(Image<float32_t>& img,
		Image<float32_t>& phi, Image<float32_t>& theta,
		double derSigma = 1.0, double smoothSigma = 1.0)
	{
		structureTensor<float32_t>(img, derSigma, smoothSigma, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, &phi, &theta, nullptr, nullptr, &img, 0.0);
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
			outScale = (double)NumberUtils<pixel_t>::scale();
		
		typedef typename NumberUtils<pixel_t>::FloatType real_t;
		
		Image<real_t> Fxx, Fyy, Fzz, Fxy, Fxz, Fyz;
		hessian(img, Fxx, Fyy, Fzz, Fxy, Fxz, Fyz, sigma, gamma);
		
		size_t counter = 0;
		#pragma omp parallel for if(Fxx.pixelCount() > PARALLELIZATION_THRESHOLD)
		for (coord_t n = 0; n < Fxx.pixelCount(); n++)
		{
			Matrix3x3d Hess(
				Fxx(n), Fxy(n), Fxz(n),
				Fxy(n), Fyy(n), Fyz(n),
				Fxz(n), Fyz(n), Fzz(n));

			double lambda1, lambda2, lambda3;
			Vec3d v1, v2, v3;

			Hess.eigsym(v1, v2, v3, lambda1, lambda2, lambda3);

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

				(*plambda123)(n) = pixelRound<pixel_t>(lambda123 * outScale);
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
					std::swap(al1, al3);
					std::swap(lambda1, lambda3);
					std::swap(v1, v3);
				}
				if (al1 > al2)
				{
					std::swap(al1, al2);
					std::swap(lambda1, lambda2);
					std::swap(v1, v2);
				}
				if (al2 > al3)
				{
					std::swap(al2, al3);
					std::swap(lambda2, lambda3);
					std::swap(v2, v3);
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

				(*pV)(n) = pixelRound<pixel_t>(Vo * outScale);
			}
		}
	}

	/**
	Non-maximum suppression (edge thinning / local maxima finding) step combined to dual thresholding step of Canny edge detection.
	Performing dual thresholding at the same time with non-maximum suppression removes need to store maximal gradient values, and thus memory usage is decreased a little bit.
	@param dx, dy, dz Partial derivatives of the image.
	@param out Empty image whose value will be set to value of gradient if gradient has local maxima at that point; and to zero otherwise.
	@param lowerThreshold, upperThreshold The minimum and maximum threshold values. If gradient magnitude is between lowerThreshold and upperThreshold, output will be assigned to value 1. If gradient magnitude is above upperThreshold, output will be assigned to value 2. Otherwise, output will be zeroed.
	*/
	template<typename real_t, typename out_t> void nonMaximumSuppression(const Image<real_t>& dx, const Image<real_t>& dy, const Image<real_t>& dz, Image<out_t>& out, real_t minThreshold, real_t maxThreshold)
	{
		// Here we need double precision.
		// If using single precision, the results differ between local and distributed processing,
		// as p0 is different in the same pixel in the two operation modes, and as a result the f1 and f2
		// values will become slightly different. Double precision seems to mitigate this at least somehow.
		using real2_t = double;

		dx.checkSize(dy);
		dz.checkSize(dy);
		out.ensureSize(dx);

		LinearInterpolator<real2_t, real_t> interp(BoundaryCondition::Nearest);
		ProgressIndicator progress(out.depth());
		#pragma omp parallel for if(out.pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		for (coord_t z = 0; z < out.depth(); z++)
		{
			for (coord_t y = 0; y < out.height(); y++)
			{
				for (coord_t x = 0; x < out.width(); x++)
				{
					// Calculate gradient value at the current pixel
					Vec3<real2_t> p0((real2_t)x, (real2_t)y, (real2_t)z);
					Vec3<real2_t> dir(dx(x, y, z), dy(x, y, z), dz(x, y, z));
					real2_t f0;
					dir.normalize(f0);
					
					// Calculate gradient value in gradient direction
					Vec3<real2_t> p1 = p0 + dir;
					real2_t dx1 = interp(dx, p1.x, p1.y, p1.z);
					real2_t dy1 = interp(dy, p1.x, p1.y, p1.z);
					real2_t dz1 = interp(dz, p1.x, p1.y, p1.z);
					real2_t f1 = (Vec3<real2_t>(dx1, dy1, dz1)).norm();

					// Calculate gradient value in reverse gradient direction
					Vec3<real2_t> p2 = p0 - dir;
					real2_t dx2 = interp(dx, p2.x, p2.y, p2.z);
					real2_t dy2 = interp(dy, p2.x, p2.y, p2.z);
					real2_t dz2 = interp(dz, p2.x, p2.y, p2.z);
					real2_t f2 = (Vec3<real2_t>(dx2, dy2, dz2)).norm();

					// Store to output values only if gradient at current location is maximal compared to the other two positions
					if(NumberUtils<real2_t>::greaterThan(f0, f1) && NumberUtils<real2_t>::greaterThan(f0, f2))
					{
						// Dual thresholding
						if(NumberUtils<real2_t>::lessThan(f0, minThreshold))
							out(x, y, z) = 0;
						else if(NumberUtils<real2_t>::lessThan(f0, maxThreshold))
							out(x, y, z) = 1;
						else
							out(x, y, z) = 2;
					}
					else
						out(x, y, z) = 0;
				}
			}

			progress.step();
		}
	}



	namespace internals
	{
		template<typename pixel_t> void cannyPart1(Image<pixel_t>& img, double derivativeSigma, double lowerThreshold, double upperThreshold)
		{
			// 1. Pre-smoothing (skipped in this version as we are calculating derivatives with Gaussian convolution anyway)
			//if(preSmoothingSigma > 0)
			//	gaussFilter(img, preSmoothingSigma, BoundaryCondition::Nearest);

			// 2. Gradient
			std::cout << "Partial derivatives..." << std::endl;
			Image<float32_t> dx, dy, dz;
			gradient<pixel_t, float32_t>(img, dx, dy, dz, derivativeSigma);


			// 3. Non-maximum suppression - find local maxima of gradient combined to
			// 4. Dual threshold - edges with gradient value above upper threshold are "surely" edges,
			// and edges with gradient value between lower and upper threshold are edges only if
			// they touch "sure" edge.
			std::cout << "Non-maximum suppression and dual thresholding..." << std::endl;
			nonMaximumSuppression<float32_t, pixel_t>(dx, dy, dz, img, pixelRound<float32_t>(lowerThreshold), pixelRound<float32_t>(upperThreshold));

			// 4. Dual threshold - edges with gradient value above upper threshold are "surely" edges,
			// and edges with gradient value between lower and upper threshold are edges only if
			// they touch "sure" edge.
			// NOTE: This step is combined to nonMaximumSuppression so that output image does not need to store gradient (real) values, but only labels.
			//vector<float32_t> th = { lowerThreshold, upperThreshold };
			//multiThreshold(img, th);
		}

		template<typename pixel_t> void cannyPart2(Image<pixel_t>& img)
		{
			// 5. Edge tracking - Convert all those edges to "sure" that touch a "sure" edge.
			std::cout << "Edge tracking..." << std::endl;
			size_t changed = grow<pixel_t>(img, pixelRound<pixel_t>(2), pixelRound<pixel_t>(1));
			std::cout << changed << " pixels changed." << std::endl;
		}
	}

	/**
	Canny edge detection.
	Does not perform the first pre-smoothing Gaussian filtering step. Do it separately if you want to perform it.
	Allocates 3 float32 temporary images.
	@param img Image that contains the structure. The contents are replaced by the edge detection result.
	@param derivativeSigma Scale parameter for derivative calculation. Set to the preferred scale of edges that should be detected.
	@param lowerThreshold, upperThreshold Edges that have gradient value above upperThreshold are always considered as edges. Edges that have gradient value between lowerThreshold and upperThreshold are considered edges only if they touch some edge with gradient value above upperThreshold.
	*/
	template<typename pixel_t> void canny(Image<pixel_t>& img, double derivativeSigma, double lowerThreshold, double upperThreshold)
	{
		// The Canny operation is separated into 3 parts so that distributed implementation is easier.
		// In the distributed implementation (using overlapped block distribution) part 1 must be run once
		// and part 2 until convergence, and part 3 (the thresholding) once.
		internals::cannyPart1(img, derivativeSigma, lowerThreshold, upperThreshold);
		internals::cannyPart2(img);
		
		// 6. Remove all other non-sure edges.
		threshold<pixel_t>(img, 1);
	}


	/**
	Convert a color value from (hue, saturation, value) representation to (red, green, blue).
	Hue must be given in radians, saturation and value must be in range [0, 1].
	Output (r, g, b) values are in range [0, 1].
	*/
	inline void hsv2rgb(double h, double s, double v, double& r, double& g, double& b)
	{
		// Convert hue to degrees.
		h *= 57.2957795;

		if (s <= 0.0)
		{
			r = v;
			g = v;
			b = v;
			return;
		}

		double hh = h;
		while (hh < 0)
			hh += 360;
		while (hh > 360)
			hh -= 360;
		hh /= 60.0;
		long i = (long)hh;
		double ff = hh - i;
		double p = v * (1.0 - s);
		double q = v * (1.0 - (s * ff));
		double t = v * (1.0 - (s * (1.0 - ff)));

		switch (i) {
		case 0:
			r = v;
			g = t;
			b = p;
			break;
		case 1:
			r = q;
			g = v;
			b = p;
			break;
		case 2:
			r = p;
			g = v;
			b = t;
			break;

		case 3:
			r = p;
			g = q;
			b = v;
			break;
		case 4:
			r = t;
			g = p;
			b = v;
			break;
		case 5:
		default:
			r = v;
			g = p;
			b = q;
			break;
		}
	}

	/**
	Calculates angle between orientation direction and given main direction for each pixel in the image.
	For pixel at x, calculates angle between (phi(x), theta(x)) and (phim, thetam), and assigns it to alpha(x).
	@param phi, theta Azimuthal (phi) and polar (theta) angles as described in the documentation of structureTensor function.
	@param alpha Output image. This image CAN be either phi or theta.
	@param phim, thetam Azimuthal and polar component of comparison angle.
	*/
	inline void orientationDifference(const Image<float32_t>& phi, const Image<float32_t>& theta, Image<float32_t>& alpha, double phim, double thetam)
	{
		phi.checkSize(theta);
		alpha.ensureSize(phi);
		if (&alpha != &phi)
			setValue(alpha, phi);
		forAll(alpha, theta, [=](float32_t p, float32_t t)
			{
				double alpha = (cos(p) * cos(phim) + sin(p) * sin(phim)) * sin(t) * sin(thetam) + cos(t) * cos(thetam);
				if (alpha > 1)
					alpha = 1;
				else if (alpha < -1)
					alpha = -1;
				alpha = acos(abs(alpha));
				return alpha;
			});
	}

	/**
	Color codes orientation data according to deviation from a given main orientation.
	In the output, hue describes angle from main orientation (phim, thetam).
	Saturation is always 1 and value is the pixel value in geometry image normalized
	such that maximum value of the geometry image is mapped to 1, and zero to 0.
	@param geometry Original geometry.
	@param phi, theta Azimuthal (phi) and polar (theta) angles as described in the documentation of structureTensor function.
	@param r, g, b Output color components.
	@param phim, thetam Azimuthal (phim) and polar (thetam) angles corresponding to the main orientation.
	*/
	template<typename pixel_t> void mainOrientationColoring(const Image<pixel_t>& geometry, const Image<float32_t>& phi, const Image<float32_t>& theta,
		Image<uint8_t>& r, Image<uint8_t>& g, Image<uint8_t>& b,
		double phim, double thetam)
	{
		geometry.checkSize(phi);
		geometry.checkSize(theta);

		r.ensureSize(geometry);
		g.ensureSize(geometry);
		b.ensureSize(geometry);

		double geomMax = (double)max(geometry);

		#pragma omp parallel for if(geometry.pixelCount() > PARALLELIZATION_THRESHOLD)
		for (coord_t n = 0; n < geometry.pixelCount(); n++)
		{
			double t = theta(n);
			double p = phi(n);

			// Calculate angle between main orientation and (phi, theta)
			double alpha = (cos(p) * cos(phim) + sin(p) * sin(phim)) * sin(t) * sin(thetam) + cos(t) * cos(thetam);
			if (alpha > 1)
				alpha = 1;
			else if (alpha < -1)
				alpha = -1;
			alpha = acos(abs(alpha));

			double hue = 2 * alpha;
			double saturation = 1;
			//double value = 1;
			double value = geometry(n) / geomMax;

			double rr, gg, bb;
			hsv2rgb(hue, saturation, value, rr, gg, bb);

			rr *= 255;
			gg *= 255;
			bb *= 255;

			r(n) = pixelRound<uint8_t>(rr);
			g(n) = pixelRound<uint8_t>(gg);
			b(n) = pixelRound<uint8_t>(bb);
		}
	}




	/**
	Color coding of orientation data used in Axelsson - Estimating 3D fibre orientation in volume images.
	This color coding is most suited to materials where most orientations are in the xy-plane, e.g. paper or cardboard.
	In the output, hue describes angle between the positive x-axis and the projection of the orientation vector to the xy-plane,
	i.e. azimuthal component of the orientation direction.
	Absolute value of the z-coordinate of the orientation direction is mapped to saturation, maximum being at the xy-plane.
	Value is mapped to the pixel value in the geometry image normalized
	such that maximum value of the geometry image is mapped to 1, and zero to 0.
	@param geometry Original geometry.
	@param phi, theta Azimuthal (phi) and polar (theta) angles as described in the documentation of structureTensor function.
	@param r, g, b Output color components.
	*/
	template<typename pixel_t> void axelssonColoring(const Image<pixel_t>& geometry, const Image<float32_t>& phi, const Image<float32_t>& theta,
		Image<uint8_t>& r, Image<uint8_t>& g, Image<uint8_t>& b)
	{
		geometry.checkSize(phi);
		geometry.checkSize(theta);

		r.ensureSize(geometry);
		g.ensureSize(geometry);
		b.ensureSize(geometry);

		double geomMax = (double)max(geometry);

#pragma omp parallel for if(geometry.pixelCount() > PARALLELIZATION_THRESHOLD)
		for (coord_t n = 0; n < geometry.pixelCount(); n++)
		{
			double t = theta(n);
			double p = phi(n);

			// Orientations in xy plane are in full color.
			// Increasing orientation in z decreases saturation.
			double hue = 2 * p;
			double z = cos(t);
			double saturation = 1 - abs(z);
			double value = geometry(n) / geomMax;


			double rr, gg, bb;
			hsv2rgb(hue, saturation, value, rr, gg, bb);

			rr *= 255;
			gg *= 255;
			bb *= 255;

			r(n) = pixelRound<uint8_t>(rr);
			g(n) = pixelRound<uint8_t>(gg);
			b(n) = pixelRound<uint8_t>(bb);
		}
	}




	namespace tests
	{
		void structureTensor();
		void lineFilter();
		void canny();
	}

}
