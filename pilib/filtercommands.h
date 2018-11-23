#pragma once

#include "commandsbase.h"
#include "overlapdistributable.h"

#include "itl2.h"

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
				concat(concat({ CommandArgument<Vec3c>(In, "radius", "Radius of neighbourhood. Diameter will be 2*r+1.", Vec3c(1, 1, 1)) },
					extraArgs),
					{ CommandArgument<NeighbourhoodType>(In, "neighbourhood type", "Type of neighbourhood, either Ellipsoidal or Rectangular.", Ellipsoidal),
					  CommandArgument<BoundaryCondition>(In, "boundary condition", "Type of boundary condition, either Zero or Nearest.", Zero)})
				)
		{
		}

		virtual void run(Image<input_t>& in, Image<output_t>& out, vector<ParamVariant>& args) const
		{
			Vec3c r = pop<Vec3c>(args);
			NeighbourhoodType nbtype = pop<NeighbourhoodType>(args);
			BoundaryCondition bc = pop<BoundaryCondition>(args);

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







	template<typename pixel_t> class MinFilterCommand : public NeighbourhoodFilterCommand<pixel_t>
	{
	public:
		MinFilterCommand() : NeighbourhoodFilterCommand<pixel_t>("minfilter", "Minimum filter. Replaces pixel by minimum of pixels in its neighbourhood. Has optimized implementation for rectangular neighbourhoods.")
		{
		}

		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, const Vec3c& r, NeighbourhoodType nbtype, BoundaryCondition bc, vector<ParamVariant>& args) const
		{
			minFilter<pixel_t, pixel_t>(in, out, r, nbtype, bc);
		}
	};

	template<typename pixel_t> class MaxFilterCommand : public NeighbourhoodFilterCommand<pixel_t>
	{
	public:
		MaxFilterCommand() : NeighbourhoodFilterCommand<pixel_t>("maxfilter", "Maximum filter. Replaces pixel by maximum of pixels in its neighbourhood. Has optimized implementation for rectangular neighbourhoods.")
		{
		}

		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, const Vec3c& r, NeighbourhoodType nbtype, BoundaryCondition bc, vector<ParamVariant>& args) const
		{
			maxFilter<pixel_t, pixel_t>(in, out, r, nbtype, bc);
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
					CommandArgument<double>(In, "noise standard deviation", "Standard deviation of noise. For a rough order of magnitude estimate, measure standard deviation from a region that does not contain any features.")
				}
				)
		{
		}

		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, const Vec3c& r, NeighbourhoodType nbtype, BoundaryCondition bc, vector<ParamVariant>& args) const
		{
			double noisestd = pop<double>(args);
			vaweFilter<pixel_t, pixel_t>(in, out, r, noisestd, nbtype, bc);
		}
	};






	template<typename pixel_t> class OpeningFilter2ParamCommand : public NeighbourhoodFilterCommand<pixel_t>
	{
	public:
		OpeningFilter2ParamCommand() :
			NeighbourhoodFilterCommand<pixel_t>("openingfilter", "Opening filter. Widens gaps in bright objects and removes objects smaller than neighbourhood size. Has optimized implementation for rectangular neighbourhoods. NOTE: This function uses the first input image as temporary storage space.",
				{}
				)
		{
		}

		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, const Vec3c& r, NeighbourhoodType nbtype, BoundaryCondition bc, vector<ParamVariant>& args) const
		{
			openingFilter<pixel_t>(in, out, r, nbtype, bc);
			setValue<pixel_t, pixel_t>(out, in);
		}
	};

	template<typename pixel_t> class ClosingFilter2ParamCommand : public NeighbourhoodFilterCommand<pixel_t>
	{
	public:
		ClosingFilter2ParamCommand() :
			NeighbourhoodFilterCommand<pixel_t>("closingfilter", "Closing filter. Closes gaps in bright objects. Has optimized implementation for rectangular neighbourhoods. NOTE: This function uses the first input image as temporary storage space.",
				{}
				)
		{
		}

		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, const Vec3c& r, NeighbourhoodType nbtype, BoundaryCondition bc, vector<ParamVariant>& args) const
		{
			closingFilter<pixel_t>(in, out, r, nbtype, bc);
			setValue<pixel_t, pixel_t>(out, in);
		}
	};

	template<typename pixel_t> class OpeningFilterCommand : public BasicOneImageNeighbourhoodCommand<pixel_t>
	{
	public:
		OpeningFilterCommand() : BasicOneImageNeighbourhoodCommand<pixel_t>("openingfilter", "Opening filter. Widens gaps in bright objects and removes objects smaller than neighbourhood size. Has optimized implementation for rectangular neighbourhoods. Creates one temporary image of same size than input.")
		{
		}

		virtual void run(Image<pixel_t>& in, coord_t r, NeighbourhoodType nbtype, BoundaryCondition bc, vector<ParamVariant>& args) const
		{
			Image<pixel_t> tmp;
			openingFilter<pixel_t>(in, tmp, r, nbtype, bc);
		}
	};

	template<typename pixel_t> class ClosingFilterCommand : public BasicOneImageNeighbourhoodCommand<pixel_t>
	{
	public:
		ClosingFilterCommand() : BasicOneImageNeighbourhoodCommand<pixel_t>("closingfilter", "Closing filter. Closes gaps in bright objects. Has optimized implementation for rectangular neighbourhoods. Creates one temporary image of same size than input.")
		{
		}

		virtual void run(Image<pixel_t>& in, coord_t r, NeighbourhoodType nbtype, BoundaryCondition bc, vector<ParamVariant>& args) const
		{
			Image<pixel_t> tmp;
			closingFilter<pixel_t>(in, tmp, r, nbtype, bc);
		}
	};





	template<typename pixel_t> class BilateralFilterCommand : public OverlapDistributable<TwoImageInputOutputCommand<pixel_t> >
	{
	public:
		BilateralFilterCommand() :
			OverlapDistributable<TwoImageInputOutputCommand<pixel_t> >("bilateralfilter", "Bilateral filtering.",
				{
					CommandArgument<double>(In, "spatial sigma", "Standard deviation of gaussian kernel used for spatial smoothing."),
					CommandArgument<double>(In, "radiometric sigma", "Standard deviation of gaussian kernel used to avoid smoothing edges of features. Order of magnitude must be similar to difference between gray levels of background and objects."),
					CommandArgument<BoundaryCondition>(In, "boundary condition", "Type of boundary condition, either Zero or Nearest.", Zero)
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
	};

	template<typename pixel_t> class GaussianFilterCommand : public OverlapDistributable<TwoImageInputOutputCommand<pixel_t> >
	{
	public:
		GaussianFilterCommand() :
			OverlapDistributable<TwoImageInputOutputCommand<pixel_t> >("gaussfilter", "Gaussian blurring. If optimization flag is set to true, processes integer images with more than 8 bits of resolution with separable convolution and float32 images with FFT filtering. If optimization flag is set to false, processes all integer images with normal convolution and float32 images with separable convolution.",
				{
					CommandArgument<Vec3d>(In, "spatial sigma", "Standard deviation of gaussian kernel."),
					CommandArgument<BoundaryCondition>(In, "boundary condition", "Type of boundary condition, either Zero or Nearest.", Zero),
					CommandArgument<bool>(In, "allow opt", "Set to true to allow optimized processing of images with high dynamic range data types. Might cause small deviations from true filtration result, and edge artefacts.", "true")
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

	template<typename pixel_t> class DerivativeCommand : public OverlapDistributable<TwoImageInputOutputCommand<pixel_t> >
	{
	public:
		DerivativeCommand() :
			OverlapDistributable<TwoImageInputOutputCommand<pixel_t> >("derivative", "Calculates Gaussian partial derivative of image, either df / dx_i or d^2 f / (dx_i dx_j).",
				{
					CommandArgument<Vec3d>(In, "spatial sigma", "Standard deviation of gaussian kernel."),
					CommandArgument<coord_t>(In, "dimension 1", "Dimension where the first partial derivative should be calculated (index i in the example above)."),
					CommandArgument<coord_t>(In, "dimension 2", "Dimension where the second partial derivative should be calculated (index j in the example above). Pass negative value to calculate only the first derivative.", -1),
					CommandArgument<BoundaryCondition>(In, "boundary condition", "Type of boundary condition, either Zero or Nearest.", Zero)
				}
				)
		{
		}

		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, vector<ParamVariant>& args) const
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
			DistributedImage<pixel_t>& out = *get<DistributedImage<pixel_t>* >(args[1]);
			Vec3d sigma = get<Vec3d>(args[2]);
			Vec3c margin = round(3 * sigma) + Vec3c(4, 4, 4);

			out.ensureSize(in.dimensions());

			return margin;
		}
	};


	template<typename pixel_t> class HighpassFilterCommand : public OverlapDistributable<TwoImageInputOutputCommand<pixel_t> >
	{
	public:
		HighpassFilterCommand() :
			OverlapDistributable<TwoImageInputOutputCommand<pixel_t> >("highpassfilter", "Gaussian high-pass filtering. Subtracts a Gaussian filtered version of input from itself. If optimization flag is set to true, processes uint16 images with separable convolution and float32 images with FFT filtering. If optimization flag is set to false, processes uint16 images with normal convolution and float32 images with separable convolution.",
				{
					CommandArgument<Vec3d>(In, "spatial sigma", "Standard deviation of gaussian kernel."),
					CommandArgument<double>(In, "shift", "Constant added to pixel values. Use, e.g., to set the mean of the filtered image to a desired value. Useful especially when filtering unsigned images where non-shifted highpass filtering will lead to negative values that will be clipped to zero.", 0),
					CommandArgument<BoundaryCondition>(In, "boundary condition", "Type of boundary condition, either Zero or Nearest.", Zero),
					CommandArgument<bool>(In, "allow opt", "Set to true to allow optimized processing of images with high dynamic range data types. Might cause small deviations from true filtration result, and edge artefacts.", "true")
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
					CommandArgument<double>(In, "spatial sigma", "Standard deviation of Gaussian kernel, determines scale of structures that are probed."),
					CommandArgument<double>(In, "scale", "Output values are scaled by this number. Pass in zero to scale output value 1 to maximum of the data type.", 0),
					CommandArgument<double>(In, "gamma", "Scale-space scaling exponent. Set to zero to disable scaling.", 0),
					CommandArgument<double>(In, "gamma23", "gamma_23 parameter; controls the sharpness of the selectivity for the cross-section isotropy. Only non-negative values are valid.", 1),
					CommandArgument<double>(In, "gamma12", "gamma_12 parameter. Non-negative values are valid.", 1),
					CommandArgument<double>(In, "alpha", "Alpha parameter, must be between 0 and 1.", 0.5),
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
	};

	template<typename pixel_t> class FrangiLineFilterCommand : public OverlapDistributable<TwoImageInputOutputCommand<pixel_t> >
	{
	public:
		FrangiLineFilterCommand() :
			OverlapDistributable<TwoImageInputOutputCommand<pixel_t> >("frangifilter", "Calculates line-enhancing filter Vo according to Frangi.",
				{
					CommandArgument<double>(In, "spatial sigma", "Standard deviation of Gaussian kernel, determines scale of structures that are probed."),
					CommandArgument<double>(In, "output scale", "Output values are scaled by this number. Pass in zero to scale output value 1 to maximum of the data type.", 0),
					CommandArgument<double>(In, "c", "Sensitivity of the filter to deviation from background noise. Typical value is quarter of the value of the maximum intensity of the lines.", 0.25),
					CommandArgument<double>(In, "gamma", "Scale-space scaling exponent. Set to zero to disable scaling.", 0),
					CommandArgument<double>(In, "alpha", "Sensitivity of the filter to deviation from plate-like structures.", 0.5),
					CommandArgument<double>(In, "beta", "Sensitivity of the filter to deviation from blob-like structures.", 0.5)
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
	};





	class BandpassFilterCommand : public OneImageInPlaceCommand<float32_t>
	{
	public:
		BandpassFilterCommand() :
			OneImageInPlaceCommand<float32_t>("bandpassfilter", "Bandpass filtering. Removes gray value variations smaller or larger in spatial extent than specified.",
				{
					CommandArgument<double>(In, "minimum size", "Variations smaller than than this value are removed.", 3.0),
					CommandArgument<double>(In, "maximum size", "Variations larger than this value are removed.", 40.0)
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
