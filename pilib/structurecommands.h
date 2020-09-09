#pragma once

#include "commandsbase.h"
#include "overlapdistributable.h"

#include "structure.h"
#include "surfacecurvature.h"

namespace pilib
{
	inline std::string orientationSeeAlso()
	{
		return "cylindricality, cylinderorientation, plateorientation, mainorientationcolor, axelssoncolor, orientationdifference";
	}

	class CylindricalityCommand : public OverlapDistributable<OneImageInPlaceCommand<float32_t> >
	{
	protected:
		friend class CommandList;

		CylindricalityCommand() : OverlapDistributable<OneImageInPlaceCommand<float32_t> >("cylindricality", "Estimates likelihood of structures being cylinders, based on eigenvalues of the local structure tensor. The input image is replaced with the cylindricality values.",
			{
				CommandArgument<double>(ParameterDirection::In, "derivative sigma", "Scale parameter. Set to the preferred scale of edges that define the cylinders. Derivatives required in the stucture tensor are calculated using convolutions with derivative of a Gaussian function, and this parameter defines the standard deviation of the Gaussian.", 1.0),
				CommandArgument<double>(ParameterDirection::In, "smoothing sigma", "The structure tensor is smoothed by convolution with a Gaussian. This parameter defines the standard deviation of the smoothing Gaussian.", 1.0)
			},
			orientationSeeAlso())
		{
		}

	public:
		virtual Vec3c calculateOverlap(const std::vector<ParamVariant>& args) const override
		{
			double derSigma = std::get<double>(args[1]);
			double smoothSigma = std::get<double>(args[2]);

			coord_t margin = itl2::round(3 * (derSigma + smoothSigma)) + 4;

			return Vec3c(margin, margin, margin);
		}

		virtual double calculateExtraMemory(const std::vector<ParamVariant>& args) const override
		{
			return 6.0;
		}

		virtual void run(Image<float32_t>& in, std::vector<ParamVariant>& args) const override
		{
			double derSigma = pop<double>(args);
			double smoothSigma = pop<double>(args);

			structureTensor<float32_t>(in, derSigma, smoothSigma, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, &in, nullptr, nullptr, 0.0);
		}

		virtual JobType getJobType(const std::vector<ParamVariant>& args) const override
		{
			return JobType::Slow;
		}
	};


	// Structure tensor planarity does not work very well so the command is disabled for now.
	//class PlanarityCommand : public OneImageInPlaceCommand<float32_t>
	//{
	//protected:
	//	friend class CommandList;
	//	PlanarityCommand() : OneImageInPlaceCommand<float32_t>("planarity", "Estimates likelihood of structures being planar, based on eigenvalues of the local structure tensor.",
	//		{
	//			CommandArgument<double>(ParameterDirection::In, "derivative sigma", "Scale parameter for derivative calculation. Set to the preferred scale of edges that should be detected. Derivatives are calculated using convolutions with derivative of Gaussian function, and this parameter defines the standard deviation of the Gaussian.", 1.0),
	//			CommandArgument<double>(ParameterDirection::In, "smoothing sigma", "The structure tensor is smoothed by convolution with Gaussian. This parameter defines the standard deviation of the smoothing Gaussian.", 1.0)
	//		})
	//	{
	//	}
	//public:
	//	virtual void run(Image<float32_t>& in, std::vector<ParamVariant>& args) const override
	//	{
	//		double derSigma = pop<double>(args);
	//		double smoothSigma = pop<double>(args);

	//		structureTensor<float32_t>(in, derSigma, smoothSigma, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &in, 0, 0.0);
	//	}
	//};



	class CylinderOrientationCommand : public OverlapDistributable<Command>
	{
	protected:
		friend class CommandList;

		CylinderOrientationCommand() : OverlapDistributable<Command>("cylinderorientation", "Estimates local orientation of cylindrical structures using the structure tensor method.",
			{
				CommandArgument<Image<float32_t> >(ParameterDirection::InOut, "geometry", "At input, an image that contains the geometry for which local orientation should be determined. At output, the orientation 'energy' will be stored in this image. It equals the sum of the eigenvalues of the structure tensor, and can be used to distinguish regions without any interfaces (i.e. no orientation, low energy value) from regions with interfaces (i.e. orientation available, high energy value)."),
				CommandArgument<Image<float32_t> >(ParameterDirection::Out, "phi", R"(The azimuthal angle of orientation direction will be stored in this image. The angle is given in radians and measured from positive $x$-axis towards positive $y$-axis and is given in range $[-\pi, \pi]$.)"),
				CommandArgument<Image<float32_t> >(ParameterDirection::Out, "theta", R"(The polar angle of orientation direction will be stored in this image. The angle is given in radians and measured from positive $z$-axis towards $xy$-plane. The values are in range $[0, \pi]$.)"),
				CommandArgument<double>(ParameterDirection::In, "derivative sigma", "Scale parameter. Set to the preferred scale of edges that define the cylinders. Derivatives required in the stucture tensor are calculated using convolutions with derivative of a Gaussian function, and this parameter defines the standard deviation of the Gaussian.", 1.0),
				CommandArgument<double>(ParameterDirection::In, "smoothing sigma", "The structure tensor is smoothed by convolution with a Gaussian. This parameter defines the standard deviation of the smoothing Gaussian.", 1.0),
			},
			orientationSeeAlso())
		{
		}

	public:
		virtual Vec3c calculateOverlap(const std::vector<ParamVariant>& args) const override
		{
			DistributedImage<float32_t>& in = *std::get<DistributedImage<float32_t>*>(args[0]);
			DistributedImage<float32_t>& phi = *std::get<DistributedImage<float32_t>*>(args[1]);
			DistributedImage<float32_t>& theta = *std::get<DistributedImage<float32_t>*>(args[2]);
			double derSigma = std::get<double>(args[3]);
			double smoothSigma = std::get<double>(args[4]);

			coord_t margin = itl2::round(3 * (derSigma + smoothSigma)) + 4;

			phi.ensureSize(in);
			theta.ensureSize(in);

			return Vec3c(margin, margin, margin);
		}

		virtual double calculateExtraMemory(const std::vector<ParamVariant>& args) const override
		{
			return 6.0 / 3.0;
		}

		virtual void run(std::vector<ParamVariant>& args) const override
		{
			Image<float32_t>& in = *pop<Image<float32_t>*>(args);
			Image<float32_t>& phi = *pop<Image<float32_t>*>(args);
			Image<float32_t>& theta = *pop<Image<float32_t>*>(args);
			double derSigma = pop<double>(args);
			double smoothSigma = pop<double>(args);

			structureTensor<float32_t>(in, derSigma, smoothSigma, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, &phi, &theta, nullptr, nullptr, &in, 0.0);
		}

		virtual JobType getJobType(const std::vector<ParamVariant>& args) const override
		{
			return JobType::Slow;
		}
	};



	class PlateOrientationCommand : public OverlapDistributable<Command>
	{
	protected:
		friend class CommandList;

		PlateOrientationCommand() : OverlapDistributable<Command>("plateorientation", "Estimates local orientation of planar structures (normal of the plane) using the structure tensor method.",
			{
				CommandArgument<Image<float32_t> >(ParameterDirection::InOut, "geometry", "At input, an image that contains the geometry for which local orientation should be determined. At output, the orientation 'energy' will be stored in this image. It equals the sum of the eigenvalues of the structure tensor, and can be used to distinguish regions without any interfaces (i.e. no orientation, low energy value) from regions with interfaces (i.e. orientation available, high energy value)."),
				CommandArgument<Image<float32_t> >(ParameterDirection::Out, "phi", R"(The azimuthal angle of orientation direction will be stored in this image. The angle is given in radians and measured from positive $x$-axis towards positive $y$-axis and is given in range $[-\pi, \pi]$.)"),
				CommandArgument<Image<float32_t> >(ParameterDirection::Out, "theta", R"(The polar angle of orientation direction will be stored in this image. The angle is given in radians and measured from positive $z$-axis towards $xy$-plane. The values are in range $[0, \pi]$.)"),
				CommandArgument<double>(ParameterDirection::In, "derivative sigma", "Scale parameter. Set to the preferred scale of edges that define the cylinders. Derivatives required in the stucture tensor are calculated using convolutions with derivative of a Gaussian function, and this parameter defines the standard deviation of the Gaussian.", 1.0),
				CommandArgument<double>(ParameterDirection::In, "smoothing sigma", "The structure tensor is smoothed by convolution with a Gaussian. This parameter defines the standard deviation of the smoothing Gaussian.", 1.0),
			},
			orientationSeeAlso())
		{
		}

	public:
		virtual Vec3c calculateOverlap(const std::vector<ParamVariant>& args) const override
		{
			DistributedImage<float32_t>& in = *std::get<DistributedImage<float32_t>*>(args[0]);
			DistributedImage<float32_t>& phi = *std::get<DistributedImage<float32_t>*>(args[1]);
			DistributedImage<float32_t>& theta = *std::get<DistributedImage<float32_t>*>(args[2]);
			double derSigma = std::get<double>(args[3]);
			double smoothSigma = std::get<double>(args[4]);

			coord_t margin = itl2::round(3 * (derSigma + smoothSigma)) + 4;

			phi.ensureSize(in);
			theta.ensureSize(in);

			return Vec3c(margin, margin, margin);
		}

		virtual double calculateExtraMemory(const std::vector<ParamVariant>& args) const override
		{
			return 6.0 / 3.0;
		}

		virtual void run(std::vector<ParamVariant>& args) const override
		{
			Image<float32_t>& in = *pop<Image<float32_t>*>(args);
			Image<float32_t>& phi = *pop<Image<float32_t>*>(args);
			Image<float32_t>& theta = *pop<Image<float32_t>*>(args);
			double derSigma = pop<double>(args);
			double smoothSigma = pop<double>(args);

			structureTensor<float32_t>(in, derSigma, smoothSigma, nullptr, nullptr, nullptr, &phi, &theta, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, &in, 0.0);
		}

		virtual JobType getJobType(const std::vector<ParamVariant>& args) const override
		{
			return JobType::Slow;
		}
	};



	class OrientationDifferenceCommand : public Command, public Distributable
	{
	protected:
		friend class CommandList;

		OrientationDifferenceCommand() : Command("orientationdifference",
			R"(For each pixel $x$, calculates angle between $(\phi(x), \theta(x))$ and $(\phi_m, \theta_m)$ and assigns that to the output image.)"
			"This command does the same calculation than what `mainorientationcode` command does, but outputs the angular difference values instead of color map, and does not weight the output values by the original geometry in any way.",
			{
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "phi", R"(The azimuthal angle of the local orientation direction. The angle is given in radians and measured from positive $x$-axis towards positive $y$-axis and is given in range $[-\pi, \pi]$.)"),
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "theta", R"(The polar angle of the local orientation direction. The angle is given in radians and measured from positive $z$-axis towards $xy$-plane. The values are in range $[0, \pi]$.)"),
				CommandArgument<Image<float32_t> >(ParameterDirection::Out, "alpha", R"(The output image. The values are in range $[0, \pi/2]$.)"),
				CommandArgument<double>(ParameterDirection::In, "phim", "The azimuthal angle of the main orientation direction in radians."),
				CommandArgument<double>(ParameterDirection::In, "thetam", "The polar angle of the main orientation direction in radians."),
			},
			orientationSeeAlso())
		{
		}

	public:

		virtual void run(std::vector<ParamVariant>& args) const override
		{
			Image<float32_t>& phi = *pop<Image<float32_t>*>(args);
			Image<float32_t>& theta = *pop<Image<float32_t>*>(args);
			Image<float32_t>& out = *pop<Image<float32_t>*>(args);
			double phim = pop<double>(args);
			double thetam = pop<double>(args);

			orientationDifference(phi, theta, out, phim, thetam);
		}

		virtual std::vector<string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const override
		{
			DistributedImage<float32_t>& phi = *std::get<DistributedImage<float32_t>*>(args[0]);
			DistributedImage<float32_t>& theta = *std::get<DistributedImage<float32_t>*>(args[1]);
			DistributedImage<float32_t>& out = *std::get<DistributedImage<float32_t>*>(args[2]);

			phi.checkSize(theta.dimensions());
			out.ensureSize(phi);

			return distributor.distribute(this, args);
		}

		using Distributable::runDistributed;

		virtual size_t getDistributionDirection2(const std::vector<ParamVariant>& args) const override
		{
			return 1;
		}

		virtual JobType getJobType(const std::vector<ParamVariant>& args) const override
		{
			return JobType::Fast;
		}

		virtual bool canDelay(const std::vector<ParamVariant>& args) const override
		{
			return true;
		}
	};




	template<typename pixel_t> class MainOrientationColoringCommand : public Command, public Distributable
	{
	protected:
		friend class CommandList;

		MainOrientationColoringCommand() : Command("mainorientationcolor",
			"Color codes orientation data according to deviation from a given main orientation. "
			"In the output, hue describes angle between the local orientation (phi and theta arguments) and the main orientation (phim and thetam arguments). "
			"Saturation is always 1 and value is the pixel value in the geometry image normalized "
			"such that maximum value of the geometry image is mapped to 1, and zero to 0.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "geometry", "An image that contains the geometry for which local orientation has been determined."),
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "phi", R"(The azimuthal angle of the local orientation direction. The angle is given in radians and measured from positive $x$-axis towards positive $y$-axis and is given in range $[-\pi, \pi]$.)"),
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "theta", R"(The polar angle of the local orientation direction. The angle is given in radians and measured from positive $z$-axis towards $xy$-plane. The values are in range $[0, \pi]$.)"),
				CommandArgument<double>(ParameterDirection::In, "phim", "The azimuthal angle of the main orientation direction in radians."),
				CommandArgument<double>(ParameterDirection::In, "thetam", "The polar angle of the main orientation direction in radians."),
				CommandArgument<Image<uint8_t> >(ParameterDirection::Out, "r", "Red channel of the result will be stored into this image."),
				CommandArgument<Image<uint8_t> >(ParameterDirection::Out, "g", "Green channel of the result will be stored into this image."),
				CommandArgument<Image<uint8_t> >(ParameterDirection::Out, "b", "Blue channel of the result will be stored into this image."),
			},
			orientationSeeAlso())
		{
		}

	public:

		virtual void run(std::vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& in = *pop<Image<pixel_t>*>(args);
			Image<float32_t>& phi = *pop<Image<float32_t>*>(args);
			Image<float32_t>& theta = *pop<Image<float32_t>*>(args);
			double phim = pop<double>(args);
			double thetam = pop<double>(args);
			Image<uint8_t>& r = *pop<Image<uint8_t>*>(args);
			Image<uint8_t>& g = *pop<Image<uint8_t>*>(args);
			Image<uint8_t>& b = *pop<Image<uint8_t>*>(args);

			mainOrientationColoring(in, phi, theta, r, g, b, phim, thetam);
		}

		virtual std::vector<string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& in = *std::get<DistributedImage<pixel_t>*>(args[0]);
			DistributedImage<float32_t>& phi = *std::get<DistributedImage<float32_t>*>(args[1]);
			DistributedImage<float32_t>& theta = *std::get<DistributedImage<float32_t>*>(args[2]);

			DistributedImage<uint8_t>& r = *std::get<DistributedImage<uint8_t>*>(args[5]);
			DistributedImage<uint8_t>& g = *std::get<DistributedImage<uint8_t>*>(args[6]);
			DistributedImage<uint8_t>& b = *std::get<DistributedImage<uint8_t>*>(args[7]);

			in.checkSize(phi.dimensions());
			in.checkSize(theta.dimensions());
			r.ensureSize(in);
			g.ensureSize(in);
			b.ensureSize(in);

			return distributor.distribute(this, args);
		}

		using Distributable::runDistributed;

		virtual size_t getDistributionDirection2(const std::vector<ParamVariant>& args) const override
		{
			return 1;
		}

		virtual JobType getJobType(const std::vector<ParamVariant>& args) const override
		{
			return JobType::Fast;
		}

		virtual bool canDelay(const std::vector<ParamVariant>& args) const override
		{
			return true;
		}
	};




	template<typename pixel_t> class AxelssonColoringCommand : public Command, public Distributable
	{
	protected:
		friend class CommandList;

		AxelssonColoringCommand() : Command("axelssoncolor",
			"Color coding of orientation data used in Axelsson - Estimating 3D fibre orientation in volume images. "
			"This color coding is most suited to materials where most orientations are in the $xy$-plane, e.g. paper or cardboard. "
			"In the output, hue describes angle between the positive $x$-axis and the projection of the orientation vector to the $xy$-plane, "
			"i.e. the azimuthal component of the orientation direction. "
			"Absolute value of the $z$-coordinate of the orientation direction is mapped to saturation, maximum being at the $xy$-plane. "
			"Value is mapped to the pixel value in the geometry image normalized "
			"such that the maximum value of the geometry image is mapped to 1, and zero to 0.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "geometry", "An image that contains the geometry for which local orientation has been determined."),
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "phi", R"(The azimuthal angle of the local orientation direction. The angle is given in radians and measured from positive $x$-axis towards positive $y$-axis and is given in range $[-\pi, \pi]$.)"),
				CommandArgument<Image<float32_t> >(ParameterDirection::In, "theta", R"(The polar angle of the local orientation direction. The angle is given in radians and measured from positive $z$-axis towards $xy$-plane. The values are in range $[0, \pi]$.)"),
				CommandArgument<Image<uint8_t> >(ParameterDirection::Out, "r", "Red channel of the result will be stored into this image."),
				CommandArgument<Image<uint8_t> >(ParameterDirection::Out, "g", "Green channel of the result will be stored into this image."),
				CommandArgument<Image<uint8_t> >(ParameterDirection::Out, "b", "Blue channel of the result will be stored into this image."),
			},
			orientationSeeAlso())
		{
		}

	public:

		virtual void run(std::vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& in = *pop<Image<pixel_t>*>(args);
			Image<float32_t>& phi = *pop<Image<float32_t>*>(args);
			Image<float32_t>& theta = *pop<Image<float32_t>*>(args);
			Image<uint8_t>& r = *pop<Image<uint8_t>*>(args);
			Image<uint8_t>& g = *pop<Image<uint8_t>*>(args);
			Image<uint8_t>& b = *pop<Image<uint8_t>*>(args);

			axelssonColoring(in, phi, theta, r, g, b);
		}

		virtual std::vector<string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& in = *std::get<DistributedImage<pixel_t>*>(args[0]);
			DistributedImage<float32_t>& phi = *std::get<DistributedImage<float32_t>*>(args[1]);
			DistributedImage<float32_t>& theta = *std::get<DistributedImage<float32_t>*>(args[2]);
			DistributedImage<uint8_t>& r = *std::get<DistributedImage<uint8_t>*>(args[3]);
			DistributedImage<uint8_t>& g = *std::get<DistributedImage<uint8_t>*>(args[4]);
			DistributedImage<uint8_t>& b = *std::get<DistributedImage<uint8_t>*>(args[5]);

			in.checkSize(phi.dimensions());
			in.checkSize(theta.dimensions());
			r.ensureSize(in);
			g.ensureSize(in);
			b.ensureSize(in);

			return distributor.distribute(this, args);
		}

		using Distributable::runDistributed;

		virtual size_t getDistributionDirection2(const std::vector<ParamVariant>& args) const override
		{
			return 1;
		}

		virtual JobType getJobType(const std::vector<ParamVariant>& args) const override
		{
			return JobType::Fast;
		}

		virtual bool canDelay(const std::vector<ParamVariant>& args) const override
		{
			return true;
		}
	};

	



	template<typename pixel_t> class SurfaceCurvatureCommand : public OverlapDistributable<Command>
	{
	protected:
		friend class CommandList;

		SurfaceCurvatureCommand() : OverlapDistributable<Command>("curvature",
			"Calculates curvature of surfaces in the image. "
			"Uses quadratic surface fitting algorithms in Petitjean - A Survey of Methods for Recovering Quadrics in Triangle Meshes. "
			"Pointwise surface normal is determined using principal component analysis of the covariance matrix of surface points near the center point. "
			"The surface normal orientation is chosen so that it points toward background voxels. "
			"Curvature is determined by transforming surface points near center point to a coordinate system where the $z$-direction is "
			"parallel to the surface normal, and then fitting a surface "
			"$f(x, y) = a x^2 + b x y + c y^2 + d$ "
			"to the tranformed points. The curvature values and directions are calculated from the coefficients $a$, $b$ and $c$. "
			"Finally, directions are then transformed back to the original coordinates. ",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "geometry", "The image containing the geometry. Non-zero pixels are assumed to be foreground."),
				CommandArgument<double>(ParameterDirection::In, "radius", "Geodesic radius (radius on the surface) of neighbourhood around surface point considered when determining curvature. Typically e.g. 5 gives good results.", 5),
				CommandArgument<Image<float32_t> >(ParameterDirection::Out, "kappa1", "Largest principal curvature will be placed to this image. Set image size to (1, 1, 1) to skip calculation of this quantity."),
				CommandArgument<Image<float32_t> >(ParameterDirection::Out, "kappa2", "Smallest principal curvature will be placed to this image. Set image size to (1, 1, 1) to skip calculation of this quantity."),
				CommandArgument<BoundaryCondition>(ParameterDirection::In, "boundary condition", string("Type of boundary condition. ") + boundaryConditionHelp(), BoundaryCondition::Nearest),
				CommandArgument<double>(ParameterDirection::In, "non-surface value", "Value that is used to fill the non-surface points in the kappa1 and kappa2 images.", std::numeric_limits<double>::signaling_NaN())
			},
			"meancurvature")
		{
		}

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& geom = *pop<Image<pixel_t>* >(args);
			double radius = pop<double>(args);
			Image<float32_t>* kappa1 = pop<Image<float32_t>* >(args);
			Image<float32_t>* kappa2 = pop<Image<float32_t>* >(args);
			BoundaryCondition bc = pop<BoundaryCondition>(args);
			double nonSurfaceValue = pop<double>(args);

			if (kappa1->dimensions() != Vec3c(1, 1, 1))
				kappa1->ensureSize(geom);
			else
				kappa1 = nullptr;

			if (kappa2->dimensions() != Vec3c(1, 1, 1))
				kappa2->ensureSize(geom);
			else
				kappa2 = nullptr;

			itl2::surfaceCurvature<pixel_t, float32_t>(geom, (float32_t)radius, kappa1, kappa2, nullptr, nullptr, bc, pixelRound<float32_t>(nonSurfaceValue));
		}

		virtual Vec3c calculateOverlap(const std::vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& geom = *std::get<DistributedImage<pixel_t>* >(args[0]);
			double radius = std::get<double>(args[1]);
			DistributedImage<float32_t>* kappa1 = std::get<DistributedImage<float32_t>* >(args[2]);
			DistributedImage<float32_t>* kappa2 = std::get<DistributedImage<float32_t>* >(args[3]);
			BoundaryCondition bc = std::get<BoundaryCondition>(args[4]);
			double nonSurfaceValue = std::get<double>(args[5]);

			if (kappa1->dimensions() != Vec3c(1, 1, 1))
				kappa1->ensureSize(geom);

			if (kappa2->dimensions() != Vec3c(1, 1, 1))
				kappa2->ensureSize(geom);

			coord_t margin = itl2::ceil(radius) + 3;

			return Vec3c(margin, margin, margin);
		}

		virtual JobType getJobType(const std::vector<ParamVariant>& args) const override
		{
			return JobType::Slow;
		}

		virtual size_t getRefIndex(const std::vector<ParamVariant>& args) const
		{
			// Input image is the reference image
			return 0;
		}
	};


	template<typename pixel_t> class MeanCurvatureCommand : public OverlapDistributable<Command>
	{
	protected:
		friend class CommandList;

		MeanCurvatureCommand() : OverlapDistributable<Command>("meancurvature",
			R"(Calculates mean curvature of surfaces defined by $f(x, y, z) = 0$, multiplied by the norm of the gradient of the image $f$. )"
			R"(Without multiplication by $|| \nabla f ||$, the curvature of $f$ is very hard to interpret as the values are only valid near the zero contour of $f$ )"
			R"((at least from the image processing point-of-view). )"
			R"(This command calculates curvature as divergence of unit normal: $-0.5 \nabla \cdot (\nabla f / ||\nabla f||) ||\nabla f||$. )"
			R"(See also https://en.wikipedia.org/wiki/Mean_curvature#Implicit_form_of_mean_curvature for formulas that can be used to determine the mean curvature. )"
			R"((The corresponding 2D formula is found at https://en.wikipedia.org/wiki/Implicit_curve#Curvature.) )"
			R"(Note that this function might produce invalid results for surfaces of structures whose thickness is similar to the value of the sigma parameter. )"
			R"(See `curvature` command for a version that does not suffer from this deficiency and that can be used to calculate the principal curvatures of surfaces in binary images. )",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "geometry", "The image containing the geometry. Non-zero pixels are assumed to be foreground. The image does not need to be a binary image."),
				CommandArgument<Image<float32_t> >(ParameterDirection::Out, "mean curvature", "The mean curvature will be placed into this image."),
				CommandArgument<double>(ParameterDirection::In, "sigma", "Scale parameter for derivative calculation. Set to the preferred scale of edges that should be detected. Derivatives are calculated using convolutions with derivative of Gaussian function, and this parameter defines the standard deviation of the Gaussian.", 1.0),
			},
			"curvature, derivative")
		{
		}

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& geom = *pop<Image<pixel_t>* >(args);
			Image<float32_t>& curvature = *pop<Image<float32_t>* >(args);
			double sigma = pop<double>(args);

			itl2::meanCurvature<pixel_t>(geom, sigma, curvature);
		}

		virtual Vec3c calculateOverlap(const std::vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& geom = *std::get<DistributedImage<pixel_t>* >(args[0]);
			DistributedImage<float32_t>& curvature = *std::get<DistributedImage<float32_t>* >(args[1]);
			double sigma = std::get<double>(args[2]);

			coord_t margin = itl2::round(3 * sigma) + 4;

			curvature.ensureSize(geom);

			return Vec3c(margin, margin, margin);
		}

		virtual double calculateExtraMemory(const std::vector<ParamVariant>& args) const override
		{
			double current = sizeof(pixel_t) + sizeof(float32_t);
			double total = current + 4 * sizeof(float32_t);
			// current * (1 + x) = total => x = total / current - 1
			return std::max(0.0, total / current - 1.0);
		}

		virtual JobType getJobType(const std::vector<ParamVariant>& args) const override
		{
			return JobType::Slow;
		}
	};
}