#pragma once

#include "commandsbase.h"
#include "overlapdistributable.h"
#include "filters.h"
#include "structure.h"

using namespace itl2;

namespace pilib
{

	


	/**
	Base class for commands have input and output image, neighbourhood radius and neighbourhood type parameters.
	*/
	template<typename input_t, typename output_t = input_t> class NeighbourhoodFilterCommand : public OverlapDistributable<TwoImageInputOutputCommand<input_t, output_t> >
	{
	public:
		NeighbourhoodFilterCommand(const string& name, const string& help, vector<CommandArgumentBase> extraArgs = {}) :
			OverlapDistributable<TwoImageInputOutputCommand<input_t, output_t> >(name, help,
				concat(concat({ CommandArgument<Vec3c>(ParameterDirection::In, "radius", "Radius of neighbourhood. Diameter will be 2*r+1.", Vec3c(1, 1, 1)) },
					extraArgs),
					{ CommandArgument<NeighbourhoodType>(ParameterDirection::In, "neighbourhood type", "Type of neighbourhood, either Ellipsoidal or Rectangular.", NeighbourhoodType::Ellipsoidal),
					  CommandArgument<BoundaryCondition>(ParameterDirection::In, "boundary condition", "Type of boundary condition, either Zero or Nearest.", BoundaryCondition::Nearest)})
				)
		{
		}

		virtual void run(Image<input_t>& in, Image<output_t>& out, vector<ParamVariant>& args) const
		{
			Vec3c r = pop<Vec3c>(args);
			NeighbourhoodType nbtype = get<NeighbourhoodType>(args[args.size() - 2]);
			BoundaryCondition bc = get<BoundaryCondition>(args[args.size() - 1]);
			args.erase(args.end() - 1);
			args.erase(args.end() - 1);

			run(in, out, r, nbtype, bc, args);
		}

		virtual Vec3c calculateOverlap(vector<ParamVariant>& args) const
		{
			DistributedImage<input_t>& in = *get<DistributedImage<input_t>* >(args[0]);
			DistributedImage<output_t>& out = *get<DistributedImage<output_t>* >(args[1]);
			Vec3c r = get<Vec3c>(args[2]);
			Vec3c margin = r + Vec3c(1, 1, 1);

			out.ensureSize(in.dimensions());

			return margin;
		}

		virtual void run(Image<input_t>& in, Image<output_t>& out, const Vec3c& r, NeighbourhoodType nbtype, BoundaryCondition bc, vector<ParamVariant>& args) const = 0;
	};






	template<typename seed_t, typename mask_t> class MorphoRecCommand : public OverlapDistributable<TwoImageInputParamCommand<seed_t, mask_t> >
	{
	public:
		MorphoRecCommand() : OverlapDistributable<TwoImageInputParamCommand<seed_t, mask_t> >("morphorec", "Morphological reconstruction. Dilates input image (seed image) and after each dilation constraints changes to nonzero pixels of parameter image (mask image).")
		{
		}

		virtual void run(Image<seed_t>& in, Image<mask_t>& param, vector<ParamVariant>& args) const
		{
			morphoRec<seed_t, mask_t>(in, param);
		}

		virtual Vec3c calculateOverlap(vector<ParamVariant>& args) const
		{
			return Vec3c(2, 2, 2); // Actually, overlap 1 should be enough.
		}
	};



	template<typename pixel_t> class VarianceFilterCommand : public NeighbourhoodFilterCommand<pixel_t>
	{
	public:
		VarianceFilterCommand() : NeighbourhoodFilterCommand<pixel_t>("variancefilter", "Variance filter. Replaces pixel by variance of pixels in its neighbourhood. For accurate results, use on images with floating point data type. Has optimized implementation for rectangular neighbourhoods.")
		{
		}

		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, const Vec3c& r, NeighbourhoodType nbtype, BoundaryCondition bc, vector<ParamVariant>& args) const
		{
			varianceFilter<pixel_t, pixel_t>(in, out, r, nbtype, bc);
		}
	};

	template<typename pixel_t> class StddevFilterCommand : public NeighbourhoodFilterCommand<pixel_t>
	{
	public:
		StddevFilterCommand() : NeighbourhoodFilterCommand<pixel_t>("stddevfilter", "Standard deviation filter. Replaces pixel by standard deviation of pixels in its neighbourhood. For accurate results, use on images with floating point data type. Has optimized implementation for rectangular neighbourhoods.")
		{
		}

		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, const Vec3c& r, NeighbourhoodType nbtype, BoundaryCondition bc, vector<ParamVariant>& args) const
		{
			stddevFilter<pixel_t, pixel_t>(in, out, r, nbtype, bc);
		}
	};

	template<typename pixel_t> class VaWeFilterCommand : public NeighbourhoodFilterCommand<pixel_t>
	{
	public:
		VaWeFilterCommand() :
			NeighbourhoodFilterCommand<pixel_t>("vawefilter", "Variance Weighted mean filtering.",
				{
					CommandArgument<double>(ParameterDirection::In, "noise standard deviation", "Standard deviation of noise. For a rough order of magnitude estimate, measure standard deviation from a region that does not contain any features.")
				}
				)
		{
		}

		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, const Vec3c& r, NeighbourhoodType nbtype, BoundaryCondition bc, vector<ParamVariant>& args) const
		{
			double noisestd = pop<double>(args);
			vaweFilter<pixel_t, pixel_t>(in, out, r, noisestd, nbtype, bc);
		}

		virtual JobType getJobType() const
		{
			return JobType::Slow;
		}
	};








	template<typename pixel_t> class MinFilterCommand : public NeighbourhoodFilterCommand<pixel_t>
	{
	public:
		MinFilterCommand() : NeighbourhoodFilterCommand<pixel_t>("minfilter", "Minimum filter. Replaces pixel by minimum of pixels in its neighbourhood.",
			{
				CommandArgument<bool>(ParameterDirection::In, "allow opt", "Set to true to allow use of approximate decompositions of spherical structuring elements using periodic lines.", true)
			})
		{
		}

		virtual Vec3c calculateOverlap(vector<ParamVariant>& args) const
		{
			Vec3c r = get<Vec3c>(args[2]);
			NeighbourhoodType nbtype = get<NeighbourhoodType>(args[args.size() - 2]);
			
			if (nbtype == NeighbourhoodType::Ellipsoidal && r.x == r.y && r.x == r.z)
				internals::optimizeStructuringElementCached(r.x);

			return NeighbourhoodFilterCommand<pixel_t>::calculateOverlap(args);
		}

		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, const Vec3c& r, NeighbourhoodType nbtype, BoundaryCondition bc, vector<ParamVariant>& args) const
		{
			bool allowOpt = pop<bool>(args);
			minFilter<pixel_t, pixel_t>(in, out, r, nbtype, bc, allowOpt);
		}

		virtual JobType getJobType() const
		{
			return JobType::Slow;
		}
	};

	template<typename pixel_t> class MaxFilterCommand : public NeighbourhoodFilterCommand<pixel_t>
	{
	public:
		MaxFilterCommand() : NeighbourhoodFilterCommand<pixel_t>("maxfilter", "Maximum filter. Replaces pixel by maximum of pixels in its neighbourhood.",
			{
				CommandArgument<bool>(ParameterDirection::In, "allow opt", "Set to true to allow use of approximate decompositions of spherical structuring elements using periodic lines.", true)
			})
		{
		}

		virtual Vec3c calculateOverlap(vector<ParamVariant>& args) const
		{
			Vec3c r = get<Vec3c>(args[2]);
			NeighbourhoodType nbtype = get<NeighbourhoodType>(args[args.size() - 2]);

			if (nbtype == NeighbourhoodType::Ellipsoidal && r.x == r.y && r.x == r.z)
				internals::optimizeStructuringElementCached(r.x);

			return NeighbourhoodFilterCommand<pixel_t>::calculateOverlap(args);
		}

		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, const Vec3c& r, NeighbourhoodType nbtype, BoundaryCondition bc, vector<ParamVariant>& args) const
		{
			bool allowOpt = pop<bool>(args);
			maxFilter<pixel_t, pixel_t>(in, out, r, nbtype, bc, allowOpt);
		}

		virtual JobType getJobType() const
		{
			return JobType::Slow;
		}
	};

	//template<typename pixel_t> class OpeningFilter2ParamCommand : public NeighbourhoodFilterCommand<pixel_t>
	//{
	//public:
	//	OpeningFilter2ParamCommand() :
	//		NeighbourhoodFilterCommand<pixel_t>("openingfilter", "Opening filter. Widens gaps in bright objects and removes objects smaller than neighbourhood size. Has optimized implementation for rectangular neighbourhoods. NOTE: This function uses the first input image as temporary storage space.",
	//			{
	//				CommandArgument<bool>(ParameterDirection::In, "allow opt", "Set to true to allow use of approximate decompositions of spherical structuring elements using periodic lines.", true)
	//			})
	//	{
	//	}

	//	virtual Vec3c calculateOverlap(vector<ParamVariant>& args) const
	//	{
	//		DistributedImage<pixel_t>& in = *get<DistributedImage<pixel_t>* >(args[0]);
	//		DistributedImage<pixel_t>& out = *get<DistributedImage<pixel_t>* >(args[1]);
	//		Vec3c r = get<Vec3c>(args[2]);
	//		NeighbourhoodType nbtype = get<NeighbourhoodType>(args[args.size() - 2]);

	//		Vec3c margin = 2 * r + Vec3c(1, 1, 1);

	//		out.ensureSize(in.dimensions());

	//		if (nbtype == NeighbourhoodType::Ellipsoidal && r.x == r.y && r.x == r.z)
	//			internals::optimizeStructuringElementCached(r.x);

	//		return margin;
	//	}

	//	virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, const Vec3c& r, NeighbourhoodType nbtype, BoundaryCondition bc, vector<ParamVariant>& args) const
	//	{
	//		bool allowOpt = pop<bool>(args);
	//		openingFilter<pixel_t>(in, out, r, nbtype, bc, allowOpt);
	//		setValue<pixel_t, pixel_t>(out, in);
	//	}

	//	virtual JobType getJobType() const
	//	{
	//		return JobType::Slow;
	//	}
	//};

	//template<typename pixel_t> class ClosingFilter2ParamCommand : public NeighbourhoodFilterCommand<pixel_t>
	//{
	//public:
	//	ClosingFilter2ParamCommand() :
	//		NeighbourhoodFilterCommand<pixel_t>("closingfilter", "Closing filter. Closes gaps in bright objects. Has optimized implementation for rectangular neighbourhoods. NOTE: This function uses the first input image as temporary storage space.",
	//			{
	//				CommandArgument<bool>(ParameterDirection::In, "allow opt", "Set to true to allow use of approximate decompositions of spherical structuring elements using periodic lines.", true)
	//			})
	//	{
	//	}

	//	virtual Vec3c calculateOverlap(vector<ParamVariant>& args) const
	//	{
	//		DistributedImage<pixel_t>& in = *get<DistributedImage<pixel_t>* >(args[0]);
	//		DistributedImage<pixel_t>& out = *get<DistributedImage<pixel_t>* >(args[1]);
	//		Vec3c r = get<Vec3c>(args[2]);
	//		NeighbourhoodType nbtype = get<NeighbourhoodType>(args[args.size() - 2]);

	//		Vec3c margin = 2 * r + Vec3c(1, 1, 1);

	//		out.ensureSize(in.dimensions());

	//		if(nbtype == NeighbourhoodType::Ellipsoidal && r.x == r.y && r.x == r.z)
	//			internals::optimizeStructuringElementCached(r.x);

	//		return margin;
	//	}

	//	virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, const Vec3c& r, NeighbourhoodType nbtype, BoundaryCondition bc, vector<ParamVariant>& args) const
	//	{
	//		bool allowOpt = pop<bool>(args);
	//		closingFilter<pixel_t>(in, out, r, nbtype, bc, allowOpt);
	//		setValue<pixel_t, pixel_t>(out, in);
	//	}

	//	virtual JobType getJobType() const
	//	{
	//		return JobType::Slow;
	//	}
	//};

	template<typename pixel_t> class OpeningFilterCommand : public BasicOneImageNeighbourhoodCommand<pixel_t>
	{
	public:
		OpeningFilterCommand() : BasicOneImageNeighbourhoodCommand<pixel_t>("openingfilter", "Opening filter. Widens gaps in bright objects and removes objects smaller than neighbourhood size. Has optimized implementation for rectangular neighbourhoods. Creates one temporary image of same size than input.",
			{
				CommandArgument<bool>(ParameterDirection::In, "allow opt", "Set to true to allow use of approximate decompositions of spherical structuring elements using periodic lines.", true)
			})
		{
		}

		virtual void run(Image<pixel_t>& in, coord_t r, NeighbourhoodType nbtype, BoundaryCondition bc, vector<ParamVariant>& args) const
		{
			bool allowOpt = pop<bool>(args);
			Image<pixel_t> tmp;
			openingFilter<pixel_t>(in, tmp, r, nbtype, bc, allowOpt);
		}
	};

	template<typename pixel_t> class ClosingFilterCommand : public BasicOneImageNeighbourhoodCommand<pixel_t>
	{
	public:
		ClosingFilterCommand() : BasicOneImageNeighbourhoodCommand<pixel_t>("closingfilter", "Closing filter. Closes gaps in bright objects. Has optimized implementation for rectangular neighbourhoods. Creates one temporary image of same size than input.",
			{
				CommandArgument<bool>(ParameterDirection::In, "allow opt", "Set to true to allow use of approximate decompositions of spherical structuring elements using periodic lines.", true)
			})
		{
		}

		virtual void run(Image<pixel_t>& in, coord_t r, NeighbourhoodType nbtype, BoundaryCondition bc, vector<ParamVariant>& args) const
		{
			bool allowOpt = pop<bool>(args);
			Image<pixel_t> tmp;
			closingFilter<pixel_t>(in, tmp, r, nbtype, bc, allowOpt);
		}
	};





	template<typename pixel_t> class BilateralFilterCommand : public OverlapDistributable<TwoImageInputOutputCommand<pixel_t> >
	{
	public:
		BilateralFilterCommand() :
			OverlapDistributable<TwoImageInputOutputCommand<pixel_t> >("bilateralfilter", "Bilateral filtering.",
				{
					CommandArgument<double>(ParameterDirection::In, "spatial sigma", "Standard deviation of gaussian kernel used for spatial smoothing."),
					CommandArgument<double>(ParameterDirection::In, "radiometric sigma", "Standard deviation of gaussian kernel used to avoid smoothing edges of features. Order of magnitude must be similar to difference between gray levels of background and objects."),
					CommandArgument<BoundaryCondition>(ParameterDirection::In, "boundary condition", "Type of boundary condition, either Zero or Nearest.", BoundaryCondition::Nearest)
				}
				)
		{
		}

		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, vector<ParamVariant>& args) const
		{
			double noisestd = pop<double>(args);
			double radstd = pop<double>(args);
			BoundaryCondition bc = pop<BoundaryCondition>(args);

			bilateralFilter(in, out, noisestd, radstd, bc);
		}

		virtual Vec3c calculateOverlap(vector<ParamVariant>& args) const
		{
			DistributedImage<pixel_t>& in = *get<DistributedImage<pixel_t>* >(args[0]);
			DistributedImage<pixel_t>& out = *get<DistributedImage<pixel_t>* >(args[1]);
			double sigma = get<double>(args[2]);
			coord_t margin = math::round(3 * sigma + 4);

			out.ensureSize(in.dimensions());

			return Vec3c(margin, margin, margin);
		}

		virtual JobType getJobType() const
		{
			return JobType::Slow;
		}
	};

	template<typename pixel_t> class GaussianFilterCommand : public OverlapDistributable<TwoImageInputOutputCommand<pixel_t> >
	{
	public:
		GaussianFilterCommand() :
			OverlapDistributable<TwoImageInputOutputCommand<pixel_t> >("gaussfilter", "Gaussian blurring. If optimization flag is set to true, processes integer images with more than 8 bits of resolution with separable convolution and float32 images with FFT filtering. If optimization flag is set to false, processes all integer images with normal convolution and float32 images with separable convolution.",
				{
					CommandArgument<Vec3d>(ParameterDirection::In, "spatial sigma", "Standard deviation of gaussian kernel."),
					CommandArgument<BoundaryCondition>(ParameterDirection::In, "boundary condition", "Type of boundary condition, either Zero or Nearest.", BoundaryCondition::Nearest),
					CommandArgument<bool>(ParameterDirection::In, "allow opt", "Set to true to allow optimized processing of images with high dynamic range data types. Might cause small deviations from true filtration result, and edge artefacts.", true)
				}
				)
		{
		}

		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, vector<ParamVariant>& args) const
		{
			Vec3d std = pop<Vec3d>(args);
			BoundaryCondition bc = pop<BoundaryCondition>(args);
			bool opt = pop<bool>(args);

			gaussFilter(in, out, std, opt, bc);
		}

		virtual Vec3c calculateOverlap(vector<ParamVariant>& args) const
		{
			DistributedImage<pixel_t>& in = *get<DistributedImage<pixel_t>* >(args[0]);
			DistributedImage<pixel_t>& out = *get<DistributedImage<pixel_t>* >(args[1]);
			Vec3d sigma = get<Vec3d>(args[2]);
			Vec3c margin = round(3 * sigma) + Vec3c(4, 4, 4);

			out.ensureSize(in.dimensions());

			return margin;
		}
	};

	template<typename pixel_t, typename output_t> class DerivativeCommand : public OverlapDistributable<TwoImageInputOutputCommand<pixel_t, output_t> >
	{
	public:
		DerivativeCommand() :
			OverlapDistributable<TwoImageInputOutputCommand<pixel_t, output_t> >("derivative", "Calculates Gaussian partial derivative of image, either df / dx_i or d^2 f / (dx_i dx_j).",
				{
					CommandArgument<Vec3d>(ParameterDirection::In, "spatial sigma", "Standard deviation of gaussian kernel."),
					CommandArgument<coord_t>(ParameterDirection::In, "dimension 1", "Dimension where the first partial derivative should be calculated (index i in the example above)."),
					CommandArgument<coord_t>(ParameterDirection::In, "dimension 2", "Dimension where the second partial derivative should be calculated (index j in the example above). Pass negative value to calculate only the first derivative.", -1),
					CommandArgument<BoundaryCondition>(ParameterDirection::In, "boundary condition", "Type of boundary condition, either Zero or Nearest.", BoundaryCondition::Nearest)
				}
				)
		{
		}

		virtual void run(Image<pixel_t>& in, Image<output_t>& out, vector<ParamVariant>& args) const
		{
			Vec3d std = pop<Vec3d>(args);
			coord_t dim1 = pop<coord_t>(args);
			coord_t dim2 = pop<coord_t>(args);
			BoundaryCondition bc = pop<BoundaryCondition>(args);

			gaussDerivative(in, out, std, dim1, dim2, bc);
		}

		virtual Vec3c calculateOverlap(vector<ParamVariant>& args) const
		{
			DistributedImage<pixel_t>& in = *get<DistributedImage<pixel_t>* >(args[0]);
			DistributedImage<output_t>& out = *get<DistributedImage<output_t>* >(args[1]);
			Vec3d sigma = get<Vec3d>(args[2]);
			Vec3c margin = round(3 * sigma) + Vec3c(4, 4, 4);

			out.ensureSize(in.dimensions());

			return margin;
		}
	};


	template<typename pixel_t, typename out_t> class GradientMagnitudeCommand : public OverlapDistributable<TwoImageInputOutputCommand<pixel_t, out_t> >
	{
	public:
		GradientMagnitudeCommand() :
			OverlapDistributable<TwoImageInputOutputCommand<pixel_t, out_t> >("gradientmagnitude", "Calculates norm of Gaussian gradient of image. The output image can be the same than the input image.",
				{
					CommandArgument<double>(ParameterDirection::In, "spatial sigma", "Standard deviation of gaussian derivative kernel.", 1.0),
					CommandArgument<double>(ParameterDirection::In, "gamma", "Scale-space scaling exponent. Set to zero to disable scaling.", 0.0)
				}
				)
		{
		}

		virtual void run(Image<pixel_t>& img, Image<out_t>& out, vector<ParamVariant>& args) const
		{
			double std = pop<double>(args);
			double gamma = pop<double>(args);

			gradientMagnitude(img, out, std, gamma);
		}

		virtual Vec3c calculateOverlap(vector<ParamVariant>& args) const
		{
			DistributedImage<pixel_t>& img = *get<DistributedImage<pixel_t>* >(args[0]);
			DistributedImage<out_t>& out = *get<DistributedImage<out_t>* >(args[1]);
			double std = get<double>(args[2]);
			Vec3c margin = math::round(3 * std) * Vec3c(1, 1, 1) + Vec3c(4, 4, 4);

			out.ensureSize(img.dimensions());

			return margin;
		}
	};


	template<typename pixel_t> class HighpassFilterCommand : public OverlapDistributable<TwoImageInputOutputCommand<pixel_t> >
	{
	public:
		HighpassFilterCommand() :
			OverlapDistributable<TwoImageInputOutputCommand<pixel_t> >("highpassfilter", "Gaussian high-pass filtering. Subtracts a Gaussian filtered version of input from itself. If optimization flag is set to true, processes uint16 images with separable convolution and float32 images with FFT filtering. If optimization flag is set to false, processes uint16 images with normal convolution and float32 images with separable convolution.",
				{
					CommandArgument<Vec3d>(ParameterDirection::In, "spatial sigma", "Standard deviation of gaussian kernel."),
					CommandArgument<double>(ParameterDirection::In, "shift", "Constant added to pixel values. Use, e.g., to set the mean of the filtered image to a desired value. Useful especially when filtering unsigned images where non-shifted highpass filtering will lead to negative values that will be clipped to zero.", 0),
					CommandArgument<BoundaryCondition>(ParameterDirection::In, "boundary condition", "Type of boundary condition, either Zero or Nearest.", BoundaryCondition::Nearest),
					CommandArgument<bool>(ParameterDirection::In, "allow opt", "Set to true to allow optimized processing of images with high dynamic range data types. Might cause small deviations from true filtration result, and edge artefacts.", true)
				}
				)
		{
		}

		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, vector<ParamVariant>& args) const
		{
			Vec3d std = pop<Vec3d>(args);
			double shift = pop<double>(args);
			BoundaryCondition bc = pop<BoundaryCondition>(args);
			bool opt = pop<bool>(args);

			highpassFilter(in, out, std, pixelRound<pixel_t>(shift), opt, bc);
		}

		virtual Vec3c calculateOverlap(vector<ParamVariant>& args) const
		{
			DistributedImage<pixel_t>& in = *get<DistributedImage<pixel_t>* >(args[0]);
			DistributedImage<pixel_t>& out = *get<DistributedImage<pixel_t>* >(args[1]);
			Vec3d sigma = get<Vec3d>(args[2]);
			Vec3c margin = round(3 * sigma) + Vec3c(4, 4, 4);

			out.ensureSize(in.dimensions());

			return margin;
		}
	};



	template<typename pixel_t> class SatoLineFilterCommand : public OverlapDistributable<TwoImageInputOutputCommand<pixel_t> >
	{
	public:
		SatoLineFilterCommand() :
			OverlapDistributable<TwoImageInputOutputCommand<pixel_t> >("satofilter", "Calculates line-enhancing filter lambda123 according to Sato.",
				{
					CommandArgument<double>(ParameterDirection::In, "spatial sigma", "Standard deviation of Gaussian kernel, determines scale of structures that are probed."),
					CommandArgument<double>(ParameterDirection::In, "scale", "Output values are scaled by this number. Pass in zero to scale output value 1 to maximum of the data type.", 0),
					CommandArgument<double>(ParameterDirection::In, "gamma", "Scale-space scaling exponent. Set to zero to disable scaling.", 0),
					CommandArgument<double>(ParameterDirection::In, "gamma23", "gamma_23 parameter; controls the sharpness of the selectivity for the cross-section isotropy. Only non-negative values are valid.", 1),
					CommandArgument<double>(ParameterDirection::In, "gamma12", "gamma_12 parameter. Non-negative values are valid.", 1),
					CommandArgument<double>(ParameterDirection::In, "alpha", "Alpha parameter, must be between 0 and 1.", 0.5),
				}
				)
		{
		}

		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, vector<ParamVariant>& args) const
		{
			double std = pop<double>(args);
			double outScale = pop<double>(args);
			double gamma = pop<double>(args);
			double gamma23 = pop<double>(args);
			double gamma12 = pop<double>(args);
			double alpha = pop<double>(args);
			
			lineFilter<pixel_t>(in, std, &out, gamma23, gamma12, alpha, 0, 0, 0, 0, gamma, outScale);
		}

		virtual Vec3c calculateOverlap(vector<ParamVariant>& args) const
		{
			DistributedImage<pixel_t>& in = *get<DistributedImage<pixel_t>* >(args[0]);
			DistributedImage<pixel_t>& out = *get<DistributedImage<pixel_t>* >(args[1]);
			double sigma = get<double>(args[2]);
			coord_t margin = math::round(3 * sigma) + 4;

			out.ensureSize(in.dimensions());

			return Vec3c(margin, margin, margin);
		}

		virtual double calculateExtraMemory(vector<ParamVariant>& args) const
		{
			DistributedImage<pixel_t>& in = *get<DistributedImage<pixel_t>* >(args[0]);
			return (6.0 * in.pixelCount() * sizeof(typename NumberUtils<pixel_t>::FloatType) / sizeof(pixel_t)) / (2.0 * in.pixelCount());
		}

		virtual JobType getJobType() const
		{
			return JobType::Slow;
		}
	};

	template<typename pixel_t> class FrangiLineFilterCommand : public OverlapDistributable<TwoImageInputOutputCommand<pixel_t> >
	{
	public:
		FrangiLineFilterCommand() :
			OverlapDistributable<TwoImageInputOutputCommand<pixel_t> >("frangifilter", "Calculates line-enhancing filter Vo according to Frangi.",
				{
					CommandArgument<double>(ParameterDirection::In, "spatial sigma", "Standard deviation of Gaussian kernel, determines scale of structures that are probed."),
					CommandArgument<double>(ParameterDirection::In, "output scale", "Output values are scaled by this number. Pass in zero to scale output value 1 to maximum of the data type.", 0),
					CommandArgument<double>(ParameterDirection::In, "c", "Sensitivity of the filter to deviation from background noise. Typical value is quarter of the value of the maximum intensity of the lines.", 0.25),
					CommandArgument<double>(ParameterDirection::In, "gamma", "Scale-space scaling exponent. Set to zero to disable scaling.", 0),
					CommandArgument<double>(ParameterDirection::In, "alpha", "Sensitivity of the filter to deviation from plate-like structures.", 0.5),
					CommandArgument<double>(ParameterDirection::In, "beta", "Sensitivity of the filter to deviation from blob-like structures.", 0.5)
				}
				)
		{
		}

		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, vector<ParamVariant>& args) const
		{
			double std = pop<double>(args);
			double outScale = pop<double>(args);
			double c = pop<double>(args);
			double gamma = pop<double>(args);
			double alpha = pop<double>(args);
			double beta = pop<double>(args);

			lineFilter<pixel_t>(in, std, 0, 0, 0, 0, &out, c, alpha, beta, gamma, outScale);
		}

		virtual Vec3c calculateOverlap(vector<ParamVariant>& args) const
		{
			DistributedImage<pixel_t>& in = *get<DistributedImage<pixel_t>* >(args[0]);
			DistributedImage<pixel_t>& out = *get<DistributedImage<pixel_t>* >(args[1]);
			double sigma = get<double>(args[2]);
			coord_t margin = math::round(3 * sigma) + 4;

			out.ensureSize(in.dimensions());

			return Vec3c(margin, margin, margin);
		}

		virtual double calculateExtraMemory(vector<ParamVariant>& args) const
		{
			DistributedImage<pixel_t>& in = *get<DistributedImage<pixel_t>* >(args[0]);
			return (6.0 * in.pixelCount() * sizeof(typename NumberUtils<pixel_t>::FloatType) / sizeof(pixel_t)) / (2.0 * in.pixelCount());
		}

		virtual JobType getJobType() const
		{
			return JobType::Slow;
		}
	};





	class BandpassFilterCommand : public OneImageInPlaceCommand<float32_t>
	{
	public:
		BandpassFilterCommand() :
			OneImageInPlaceCommand<float32_t>("bandpassfilter", "Bandpass filtering. Removes gray value variations smaller or larger in spatial extent than specified.",
				{
					CommandArgument<double>(ParameterDirection::In, "minimum size", "Variations smaller than than this value are removed.", 3.0),
					CommandArgument<double>(ParameterDirection::In, "maximum size", "Variations larger than this value are removed.", 40.0)
				}
				)
		{
		}

		virtual void run(Image<float32_t>& img, vector<ParamVariant>& args) const
		{
			double minSize = pop<double>(args);
			double maxSize = pop<double>(args);

			bandpassFilter(img, minSize, maxSize);
		}
	};

	class FFTCommand : public TwoImageInputOutputCommand<float32_t, complex32_t>
	{
	public:
		FFTCommand() :
			TwoImageInputOutputCommand<float32_t, complex32_t>("fft", "Calculates Fourier transform."
				)
		{
		}

		virtual void run(Image<float32_t>& in, Image<complex32_t>& out, vector<ParamVariant>& args) const
		{
			fft(in, out);
		}
	};

	class InverseFFTCommand : public TwoImageInputOutputCommand<complex32_t, float32_t>
	{
	public:
		InverseFFTCommand() :
			TwoImageInputOutputCommand<complex32_t, float32_t>("ifft", "Calculates inverse Fourier transform."
				)
		{
		}

		virtual void run(Image<complex32_t>& in, Image<float32_t>& out, vector<ParamVariant>& args) const
		{
			ifft(in, out);
		}
	};
}
