#pragma once

#include "commandsbase.h"
#include "overlapdistributable.h"
#include "filters.h"
#include "structure.h"
#include "fastbilateralfilter.h"
#include "pilibutilities.h"
#include "standardhelp.h"

using namespace itl2;

namespace pilib
{


	inline std::string filterSeeAlso()
	{
		return "gaussfilter, bilateralfilter, bilateralfilterapprox, vawefilter, openingfilter, closingfilter, minfilter, maxfilter, medianfilter, variancefilter, stddevfilter, bandpassfilter, highpassfilter";
	}


	/**
	Base class for commands have input and output image, neighbourhood radius and neighbourhood type parameters.
	*/
	template<typename input_t, typename output_t = input_t> class NeighbourhoodFilterCommand : public OverlapDistributable<TwoImageInputOutputCommand<input_t, output_t> >
	{
	protected:
		friend class CommandList;

		NeighbourhoodFilterCommand(const string& name, const string& help, vector<CommandArgumentBase> extraArgs = {}) :
			OverlapDistributable<TwoImageInputOutputCommand<input_t, output_t> >(name, help,
				concat(concat({ CommandArgument<Vec3c>(ParameterDirection::In, "radius", "Radius of neighbourhood. Diameter will be $2r+1$.", Vec3c(1, 1, 1)) },
					extraArgs),
					{ CommandArgument<NeighbourhoodType>(ParameterDirection::In, "neighbourhood type", string("Type of neighbourhood. ") + neighbourhoodTypeHelp(), NeighbourhoodType::Ellipsoidal),
					  CommandArgument<BoundaryCondition>(ParameterDirection::In, "boundary condition", string("Type of boundary condition. ") + boundaryConditionHelp(), BoundaryCondition::Nearest)
					}),
				filterSeeAlso()
				)
		{
		}

	public:
		virtual void run(Image<input_t>& in, Image<output_t>& out, vector<ParamVariant>& args) const override
		{
			Vec3c r = pop<Vec3c>(args);
			NeighbourhoodType nbtype = std::get<NeighbourhoodType>(args[args.size() - 2]);
			BoundaryCondition bc = std::get<BoundaryCondition>(args[args.size() - 1]);
			args.erase(args.end() - 1);
			args.erase(args.end() - 1);

			run(in, out, r, nbtype, bc, args);
		}

		virtual Vec3c calculateOverlap(const vector<ParamVariant>& args) const override
		{
			DistributedImage<input_t>& in = *std::get<DistributedImage<input_t>* >(args[0]);
			DistributedImage<output_t>& out = *std::get<DistributedImage<output_t>* >(args[1]);
			Vec3c r = std::get<Vec3c>(args[2]);
			Vec3c margin = r + Vec3c(1, 1, 1);

			out.ensureSize(in.dimensions());

			return margin;
		}

		virtual void run(Image<input_t>& in, Image<output_t>& out, const Vec3c& r, NeighbourhoodType nbtype, BoundaryCondition bc, vector<ParamVariant>& args) const = 0;
	};






	template<typename seed_t, typename mask_t> class MorphoRecCommand : public TwoImageInputParamCommand<seed_t, mask_t>, public Distributable
	{
	protected:
		friend class CommandList;

		MorphoRecCommand() : TwoImageInputParamCommand<seed_t, mask_t>("morphorec", "Morphological reconstruction. Dilates input image (seed image) and after each dilation constraints changes to nonzero pixels of parameter image (mask image). Continues until the image does not change anymore.")
		{
		}

	public:
		virtual void run(Image<seed_t>& in, Image<mask_t>& param, vector<ParamVariant>& args) const override
		{
			size_t changed = morphoRec<seed_t, mask_t>(in, param);
			std::cout << std::endl << changed << " pixels changed." << std::endl;
		}

		virtual Vec3c getMargin(const vector<ParamVariant>& args) const override
		{
			return Vec3c(2, 2, 2); // Actually, overlap 1 should be enough.
		}

		virtual JobType getJobType(const vector<ParamVariant>& args) const override
		{
			return JobType::Slow;
		}

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override
		{
			// TODO: This could be made with IterableDistributable
			coord_t changed;
			do
			{
				vector<string> output = distributor.distribute(this, args);

				changed = parseTotalCount(output, "pixels changed.");
				std::cout << std::endl << changed << " pixels changed." << std::endl;
			} while (changed > 0);

			return vector<string>();
		}
	};



	template<typename pixel_t> class MedianFilterCommand : public NeighbourhoodFilterCommand<pixel_t>
	{
	protected:
		friend class CommandList;

		MedianFilterCommand() : NeighbourhoodFilterCommand<pixel_t>("medianfilter", "Median filtering. Replaces pixel by median of pixels in its neighbourhood. Removes noise from the image while preserving sharp edges.")
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, const Vec3c& r, NeighbourhoodType nbtype, BoundaryCondition bc, vector<ParamVariant>& args) const override
		{
			medianFilter<pixel_t, pixel_t>(in, out, r, nbtype, bc);
		}
	};

	template<typename pixel_t> class VarianceFilterCommand : public NeighbourhoodFilterCommand<pixel_t>
	{
	protected:
		friend class CommandList;

		VarianceFilterCommand() : NeighbourhoodFilterCommand<pixel_t>("variancefilter", "Variance filter. Replaces pixel by variance of pixels in its neighbourhood. For accurate results, use on images with floating point data type. Has optimized implementation for rectangular neighbourhoods.")
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, const Vec3c& r, NeighbourhoodType nbtype, BoundaryCondition bc, vector<ParamVariant>& args) const override
		{
			varianceFilter<pixel_t, pixel_t>(in, out, r, nbtype, bc);
		}
	};

	template<typename pixel_t> class StddevFilterCommand : public NeighbourhoodFilterCommand<pixel_t>
	{
	protected:
		friend class CommandList;

		StddevFilterCommand() : NeighbourhoodFilterCommand<pixel_t>("stddevfilter", "Standard deviation filter. Replaces pixel by standard deviation of pixels in its neighbourhood. For accurate results, use on images with floating point data type. Has optimized implementation for rectangular neighbourhoods.")
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, const Vec3c& r, NeighbourhoodType nbtype, BoundaryCondition bc, vector<ParamVariant>& args) const override
		{
			stddevFilter<pixel_t, pixel_t>(in, out, r, nbtype, bc);
		}
	};

	template<typename pixel_t> class VaWeFilterCommand : public NeighbourhoodFilterCommand<pixel_t>
	{
	protected:
		friend class CommandList;

		VaWeFilterCommand() :
			NeighbourhoodFilterCommand<pixel_t>("vawefilter", "Variance Weighted mean filtering. Removes noise from the image while trying to preserve edges. Edge preservation is achieved by filtering less aggressively in regions where local variance of pixel values is high. The filter thus assumes that high local variance corresponds to details of interest, and low local variance corresponds to regions containing solely noise.",
				{
					CommandArgument<double>(ParameterDirection::In, "noise standard deviation", "Standard deviation of noise. For a rough order of magnitude estimate, measure standard deviation from a region that does not contain any features.")
				}
				)
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, const Vec3c& r, NeighbourhoodType nbtype, BoundaryCondition bc, vector<ParamVariant>& args) const override
		{
			double noisestd = pop<double>(args);
			vaweFilter<pixel_t, pixel_t>(in, out, r, noisestd, nbtype, bc);
		}

		virtual JobType getJobType(const vector<ParamVariant>& args) const override
		{
			return JobType::Slow;
		}
	};




	inline std::string periodicLinesHelp()
	{
		return "Set to true to allow use of approximate decompositions of spherical structuring elements using periodic lines. As a result of the approximation processing is much faster but the true shape of the structuring element is not sphere but a regular polyhedron. See van Herk - A fast algorithm for local minimum and maximum filters on rectangular and octagonal kernels and Jones - Periodic lines Definition, cascades, and application to granulometries.";
	}



	template<typename pixel_t> class MinFilterCommand : public NeighbourhoodFilterCommand<pixel_t>
	{
	protected:
		friend class CommandList;

		MinFilterCommand() : NeighbourhoodFilterCommand<pixel_t>("minfilter", "Minimum filter. Replaces pixel by minimum of pixels in its neighbourhood.",
			{
				CommandArgument<bool>(ParameterDirection::In, "allow optimization", periodicLinesHelp(), true)
			})
		{
		}

	public:
		virtual Vec3c calculateOverlap(const vector<ParamVariant>& args) const override
		{
			Vec3c r = std::get<Vec3c>(args[2]);
			NeighbourhoodType nbtype = std::get<NeighbourhoodType>(args[args.size() - 2]);
			
			if (nbtype == NeighbourhoodType::Ellipsoidal && r.x == r.y && r.x == r.z)
				internals::optimizeStructuringElementCached(r.x);

			return NeighbourhoodFilterCommand<pixel_t>::calculateOverlap(args);
		}

		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, const Vec3c& r, NeighbourhoodType nbtype, BoundaryCondition bc, vector<ParamVariant>& args) const override
		{
			bool allowOpt = pop<bool>(args);
			minFilter<pixel_t, pixel_t>(in, out, r, nbtype, bc, allowOpt);
		}

		virtual JobType getJobType(const vector<ParamVariant>& args) const override
		{
			Vec3c r = std::get<Vec3c>(args[2]);
			if (r.max() > 10)
				return JobType::Slow;
			else
				return JobType::Normal;
		}
	};

	template<typename pixel_t> class MaxFilterCommand : public NeighbourhoodFilterCommand<pixel_t>
	{
	protected:
		friend class CommandList;

		MaxFilterCommand() : NeighbourhoodFilterCommand<pixel_t>("maxfilter", "Maximum filter. Replaces pixel by maximum of pixels in its neighbourhood.",
			{
				CommandArgument<bool>(ParameterDirection::In, "allow optimization", periodicLinesHelp(), true)
			})
		{
		}

	public:
		virtual Vec3c calculateOverlap(const vector<ParamVariant>& args) const override
		{
			Vec3c r = std::get<Vec3c>(args[2]);
			NeighbourhoodType nbtype = std::get<NeighbourhoodType>(args[args.size() - 2]);

			if (nbtype == NeighbourhoodType::Ellipsoidal && r.x == r.y && r.x == r.z)
				internals::optimizeStructuringElementCached(r.x);

			return NeighbourhoodFilterCommand<pixel_t>::calculateOverlap(args);
		}

		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, const Vec3c& r, NeighbourhoodType nbtype, BoundaryCondition bc, vector<ParamVariant>& args) const override
		{
			bool allowOpt = pop<bool>(args);
			maxFilter<pixel_t, pixel_t>(in, out, r, nbtype, bc, allowOpt);
		}

		virtual JobType getJobType(const vector<ParamVariant>& args) const override
		{
			Vec3c r = std::get<Vec3c>(args[2]);
			if (r.max() > 10)
				return JobType::Slow;
			else
				return JobType::Normal;
		}
	};

	//template<typename pixel_t> class OpeningFilter2ParamCommand : public NeighbourhoodFilterCommand<pixel_t>
	//{
	//protected:
	//	friend class CommandList;

	//	OpeningFilter2ParamCommand() :
	//		NeighbourhoodFilterCommand<pixel_t>("openingfilter", "Opening filter. Widens gaps in bright objects and removes objects smaller than neighbourhood size. Has optimized implementation for rectangular neighbourhoods. NOTE: This function uses the first input image as temporary storage space.",
	//			{
	//				CommandArgument<bool>(ParameterDirection::In, "allow opt", "Set to true to allow use of approximate decompositions of spherical structuring elements using periodic lines.", true)
	//			})
	//	{
	//	}
	//public:
	//	virtual Vec3c calculateOverlap(const vector<ParamVariant>& args) const
	//	{
	//		DistributedImage<pixel_t>& in = *std::get<DistributedImage<pixel_t>* >(args[0]);
	//		DistributedImage<pixel_t>& out = *std::get<DistributedImage<pixel_t>* >(args[1]);
	//		Vec3c r = std::get<Vec3c>(args[2]);
	//		NeighbourhoodType nbtype = std::get<NeighbourhoodType>(args[args.size() - 2]);

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

	//	virtual JobType getJobType(const vector<ParamVariant>& args) const
	//	{
	//		return JobType::Slow;
	//	}
	//};

	//template<typename pixel_t> class ClosingFilter2ParamCommand : public NeighbourhoodFilterCommand<pixel_t>
	//{
	//protected:
	//	friend class CommandList;

	//	ClosingFilter2ParamCommand() :
	//		NeighbourhoodFilterCommand<pixel_t>("closingfilter", "Closing filter. Closes gaps in bright objects. Has optimized implementation for rectangular neighbourhoods. NOTE: This function uses the first input image as temporary storage space.",
	//			{
	//				CommandArgument<bool>(ParameterDirection::In, "allow opt", "Set to true to allow use of approximate decompositions of spherical structuring elements using periodic lines.", true)
	//			})
	//	{
	//	}

	//public:
	//	virtual Vec3c calculateOverlap(const vector<ParamVariant>& args) const
	//	{
	//		DistributedImage<pixel_t>& in = *std::get<DistributedImage<pixel_t>* >(args[0]);
	//		DistributedImage<pixel_t>& out = *std::get<DistributedImage<pixel_t>* >(args[1]);
	//		Vec3c r = std::get<Vec3c>(args[2]);
	//		NeighbourhoodType nbtype = std::get<NeighbourhoodType>(args[args.size() - 2]);

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

	//	virtual JobType getJobType(const vector<ParamVariant>& args) const
	//	{
	//		return JobType::Slow;
	//	}
	//};


	// This is not derived from BasicOneImageNeighbourhoodCommand as that supports only coord_t r, not Vec3c r.
	template<typename pixel_t, void operation(Image<pixel_t>&, Image<pixel_t>&, const Vec3c&, NeighbourhoodType, BoundaryCondition, bool)> class OpeningClosingFilterCommandBase : public OverlapDistributable<OneImageInPlaceCommand<pixel_t> >
	{
	protected:
		friend class CommandList;

		OpeningClosingFilterCommandBase(const string& name, const string& help) : OverlapDistributable<OneImageInPlaceCommand<pixel_t> >(name, help,
			{
				CommandArgument<Vec3c>(ParameterDirection::In, "radius", "Radius of neighbourhood. Diameter will be $2r+1$.", Vec3c(1, 1, 1)),
				CommandArgument<bool>(ParameterDirection::In, "allow optimization", periodicLinesHelp(), true),
				CommandArgument<NeighbourhoodType>(ParameterDirection::In, "neighbourhood type", string("Type of neighbourhood. ") + neighbourhoodTypeHelp(), NeighbourhoodType::Ellipsoidal),
				CommandArgument<BoundaryCondition>(ParameterDirection::In, "boundary condition", string("Type of boundary condition. ") + boundaryConditionHelp(), BoundaryCondition::Nearest)
			},
			filterSeeAlso())
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const override
		{
			Vec3c r = pop<Vec3c>(args);
			bool allowOpt = pop<bool>(args);
			NeighbourhoodType nbtype = pop<NeighbourhoodType>(args);
			BoundaryCondition bc = pop<BoundaryCondition>(args);
			
			Image<pixel_t> tmp;
			operation(in, tmp, r, nbtype, bc, allowOpt);
		}

		virtual Vec3c calculateOverlap(const vector<ParamVariant>& args) const override
		{
			Vec3c r = std::get<Vec3c>(args[1]);
			NeighbourhoodType nbtype = std::get<NeighbourhoodType>(args[3]);

			Vec3c margin = 2 * r + Vec3c(1, 1, 1);

			if (nbtype == NeighbourhoodType::Ellipsoidal && r.x == r.y && r.x == r.z)
				internals::optimizeStructuringElementCached(r.x);

			return margin;
		}

		virtual double calculateExtraMemory(const vector<ParamVariant>& args) const override
		{
			return 1;
		}

		virtual JobType getJobType(const vector<ParamVariant>& args) const override
		{
			Vec3c r = std::get<Vec3c>(args[2]);
			if (r.max() > 10)
				return JobType::Slow;
			else
				return JobType::Normal;
		}
	};

	template<typename pixel_t> class OpeningFilterCommand : public OpeningClosingFilterCommandBase<pixel_t, openingFilter<pixel_t> >
	{
	protected:
		friend class CommandList;

		OpeningFilterCommand() : OpeningClosingFilterCommandBase<pixel_t, openingFilter<pixel_t> >("openingfilter", "Opening filter. Widens gaps in bright objects and removes objects smaller than neighbourhood size. Has optimized implementation for rectangular neighbourhoods. Creates one temporary image of same size than input.")
		{
		}

		//virtual void run(Image<pixel_t>& in, coord_t r, NeighbourhoodType nbtype, BoundaryCondition bc, vector<ParamVariant>& args) const
		//{
		//	bool allowOpt = pop<bool>(args);
		//	Image<pixel_t> tmp;
		//	openingFilter<pixel_t>(in, tmp, r, nbtype, bc, allowOpt);
		//}
	};

	template<typename pixel_t> class ClosingFilterCommand : public OpeningClosingFilterCommandBase<pixel_t, closingFilter<pixel_t> >
	{
	protected:
		friend class CommandList;

		ClosingFilterCommand() : OpeningClosingFilterCommandBase<pixel_t, closingFilter<pixel_t> >("closingfilter", "Closing filter. Closes gaps in bright objects. Has optimized implementation for rectangular neighbourhoods. Creates one temporary image of same size than input.")
		{
		}

		//virtual void run(Image<pixel_t>& in, coord_t r, NeighbourhoodType nbtype, BoundaryCondition bc, vector<ParamVariant>& args) const
		//{
		//	bool allowOpt = pop<bool>(args);
		//	Image<pixel_t> tmp;
		//	closingFilter<pixel_t>(in, tmp, r, nbtype, bc, allowOpt);
		//}
	};






	template<typename pixel_t> class BilateralFilterCommand : public OverlapDistributable<TwoImageInputOutputCommand<pixel_t> >
	{
	protected:
		friend class CommandList;

		BilateralFilterCommand() :
			OverlapDistributable<TwoImageInputOutputCommand<pixel_t> >("bilateralfilter", "Bilateral filtering. Removes noise from the image while trying to preserve sharp edges. The filter is realized as a weighted local average, where weight value depends on both spatial and radiometric distance to the central pixel.",
				{
					CommandArgument<double>(ParameterDirection::In, "spatial sigma", "Standard deviation of Gaussian kernel used for spatial smoothing."),
					CommandArgument<double>(ParameterDirection::In, "radiometric sigma", "Standard deviation of Gaussian kernel used to avoid smoothing edges of features. Order of magnitude must be similar to difference between gray levels of background and objects."),
					CommandArgument<BoundaryCondition>(ParameterDirection::In, "boundary condition", string("Type of boundary condition. ") + boundaryConditionHelp(), BoundaryCondition::Nearest)
				},
				filterSeeAlso()
				)
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, vector<ParamVariant>& args) const override
		{
			double noisestd = pop<double>(args);
			double radstd = pop<double>(args);
			BoundaryCondition bc = pop<BoundaryCondition>(args);

			bilateralFilter(in, out, noisestd, radstd, bc);
		}

		virtual Vec3c calculateOverlap(const vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& in = *std::get<DistributedImage<pixel_t>* >(args[0]);
			DistributedImage<pixel_t>& out = *std::get<DistributedImage<pixel_t>* >(args[1]);
			double sigma = std::get<double>(args[2]);
			coord_t margin = itl2::round(3 * sigma + 4);

			out.ensureSize(in.dimensions());

			return Vec3c(margin, margin, margin);
		}

		virtual JobType getJobType(const vector<ParamVariant>& args) const override
		{
			return JobType::Slow;
		}
	};




	template<typename pixel_t> class ApproxBilateralFilterCommand : public OverlapDistributable<TwoImageInputOutputCommand<pixel_t> >
	{
	protected:
		friend class CommandList;

		ApproxBilateralFilterCommand() :
			OverlapDistributable<TwoImageInputOutputCommand<pixel_t> >("bilateralfilterapprox", "Approximate bilateral filtering. Removes noise from the image while trying to preserve sharp edges. The filter is realized as a weighted local average, where weight value depends on both spatial and radiometric distance to the central pixel. The difference to exact bilateral filter is that the approximate version does not process all pixels in a neighbourhood but only a random subset of them. As a result, the filtering operation is much faster but the output may contain some random noise.",
				{
					CommandArgument<double>(ParameterDirection::In, "spatial sigma", "Standard deviation of Gaussian kernel used for spatial smoothing."),
					CommandArgument<double>(ParameterDirection::In, "radiometric sigma", "Standard deviation of Gaussian kernel used to avoid smoothing edges of features. Order of magnitude must be similar to difference between gray levels of background and objects."),
					CommandArgument<size_t>(ParameterDirection::In, "sample count", "Count of samples processed in each neighbourhood. Specify zero to determine sample count automatically.", 0)
					// Random seed is disable as it is a bit hard to make it work correctly in the distributed case.
					//CommandArgument<size_t>(ParameterDirection::In, "random seed", "Seed for random number generation. Specify zero to determine seed automatically from current time.", 0),
				},
				filterSeeAlso()
				)
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, vector<ParamVariant>& args) const override
		{
			double noisestd = pop<double>(args);
			double radstd = pop<double>(args);
			size_t sampleCount = pop<size_t>(args);
			//size_t seed = pop<size_t>(args);

			bilateralFilterSampling(in, out, (float32_t)noisestd, (float32_t)radstd, sampleCount);
		}

		virtual Vec3c calculateOverlap(const vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& in = *std::get<DistributedImage<pixel_t>* >(args[0]);
			DistributedImage<pixel_t>& out = *std::get<DistributedImage<pixel_t>* >(args[1]);
			double sigma = std::get<double>(args[2]);
			coord_t margin = itl2::round(3 * sigma + 4);

			out.ensureSize(in.dimensions());

			return Vec3c(margin, margin, margin);
		}
		
		virtual JobType getJobType(const vector<ParamVariant>& args) const override
		{
			return JobType::Slow;
		}
	};



	inline std::string gaussianOptimizationHelp()
	{
		return "If optimization flag is set to true, processes integer images with more than 8 bits of resolution with separable convolution and floating point images with FFT filtering. If optimization flag is set to false, processes all integer images with normal convolution and floating point images with separable convolution.";
	}

	template<typename pixel_t> class GaussianFilterCommand : public OverlapDistributable<TwoImageInputOutputCommand<pixel_t> >
	{
	protected:
		friend class CommandList;

		GaussianFilterCommand() :
			OverlapDistributable<TwoImageInputOutputCommand<pixel_t> >("gaussfilter", "Gaussian blurring. Removes imaging noise but makes images look unsharp. " + gaussianOptimizationHelp(),
				{
					CommandArgument<Vec3d>(ParameterDirection::In, "spatial sigma", "Standard deviation of Gaussian kernel."),
					CommandArgument<BoundaryCondition>(ParameterDirection::In, "boundary condition", string("Type of boundary condition. ") + boundaryConditionHelp(), BoundaryCondition::Nearest),
					CommandArgument<bool>(ParameterDirection::In, "allow optimization", "Set to true to allow optimized processing of images with high dynamic range data types. Might cause small deviations from true filtering result and edge artefacts.", true)
				},
				filterSeeAlso()
				)
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, vector<ParamVariant>& args) const override
		{
			Vec3d std = pop<Vec3d>(args);
			BoundaryCondition bc = pop<BoundaryCondition>(args);
			bool opt = pop<bool>(args);

			gaussFilter(in, out, std, opt, bc);
		}

		virtual Vec3c calculateOverlap(const vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& in = *std::get<DistributedImage<pixel_t>* >(args[0]);
			DistributedImage<pixel_t>& out = *std::get<DistributedImage<pixel_t>* >(args[1]);
			Vec3d sigma = std::get<Vec3d>(args[2]);
			Vec3c margin = itl2::round(3 * sigma) + Vec3c(4, 4, 4);

			out.ensureSize(in.dimensions());

			return margin;
		}
	};

	template<typename pixel_t, typename output_t> class DerivativeCommand : public OverlapDistributable<TwoImageInputOutputCommand<pixel_t, output_t> >
	{
	protected:
		friend class CommandList;

		DerivativeCommand() :
			OverlapDistributable<TwoImageInputOutputCommand<pixel_t, output_t> >("derivative", R"(Calculates Gaussian partial derivative of image, either $\partial f / \partial x_i$ or $\partial^2 f / (\partial x_i \partial x_j)$.)",
				{
					CommandArgument<Vec3d>(ParameterDirection::In, "spatial sigma", "Standard deviation of Gaussian kernel."),
					CommandArgument<coord_t>(ParameterDirection::In, "dimension 1", "Dimension where the first partial derivative should be calculated (index $i$ in the example above)."),
					CommandArgument<coord_t>(ParameterDirection::In, "dimension 2", "Dimension where the second partial derivative should be calculated (index $j$ in the example above). Pass negative value to calculate only the first derivative.", -1),
					CommandArgument<BoundaryCondition>(ParameterDirection::In, "boundary condition", string("Type of boundary condition. ") + boundaryConditionHelp(), BoundaryCondition::Nearest)
				},
				"gaussfilter, gradient, gradientmagnitude"
				)
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, Image<output_t>& out, vector<ParamVariant>& args) const override
		{
			Vec3d std = pop<Vec3d>(args);
			coord_t dim1 = pop<coord_t>(args);
			coord_t dim2 = pop<coord_t>(args);
			BoundaryCondition bc = pop<BoundaryCondition>(args);

			gaussDerivative(in, out, std, dim1, dim2, bc);
		}

		virtual Vec3c calculateOverlap(const vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& in = *std::get<DistributedImage<pixel_t>* >(args[0]);
			DistributedImage<output_t>& out = *std::get<DistributedImage<output_t>* >(args[1]);
			Vec3d sigma = std::get<Vec3d>(args[2]);
			Vec3c margin = itl2::round(3 * sigma) + Vec3c(4, 4, 4);

			out.ensureSize(in.dimensions());

			return margin;
		}
	};


	template<typename pixel_t, typename output_t> class GradientCommand : public OverlapDistributable<Command>
	{
	protected:
		friend class CommandList;

		GradientCommand() :
			OverlapDistributable<Command>("gradient", R"(Calculates Gaussian gradient $(\partial f / \partial x, \partial f / \partial y, \partial f / \partial z)$ of an image $f$. Each of the derivatives is calculated by convolving the image with the corresponding derivative of the Gaussian function.)",
				{
					CommandArgument<Image<pixel_t> >(ParameterDirection::In, "f", "Image whose gradient is to be calculated."),
					CommandArgument<double>(ParameterDirection::In, "spatial sigma", "Standard deviation of Gaussian kernel."),
					CommandArgument<Image<output_t> >(ParameterDirection::Out, "dfdx", "Derivative in $x$-direction"),
					CommandArgument<Image<output_t> >(ParameterDirection::Out, "dfdy", "Derivative in $y$-direction"),
					CommandArgument<Image<output_t> >(ParameterDirection::Out, "dfdz", "Derivative in $z$-direction"),
					CommandArgument<double>(ParameterDirection::In, "gamma", "Scale-space scaling exponent according to Lindeberg. Set to zero to disable scaling.", 0.0),
				},
				"derivative, gradientmagnitude"
				)
		{
		}

	public:
		virtual void run(vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& img = *pop<Image<pixel_t>*>(args);
			double std = pop<double>(args);
			Image<output_t>& dfdx = *pop<Image<output_t>*>(args);
			Image<output_t>& dfdy = *pop<Image<output_t>*>(args);
			Image<output_t>& dfdz = *pop<Image<output_t>*>(args);
			double gamma = pop<double>(args);

			gradient(img, dfdx, dfdy, dfdz, std, gamma);
		}

		virtual Vec3c calculateOverlap(const vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& img = *std::get<DistributedImage<pixel_t>*>(args[0]);
			double sigma = std::get<double>(args[1]);
			DistributedImage<output_t>& dfdx = *std::get<DistributedImage<output_t>*>(args[2]);
			DistributedImage<output_t>& dfdy = *std::get<DistributedImage<output_t>*>(args[3]);
			DistributedImage<output_t>& dfdz = *std::get<DistributedImage<output_t>*>(args[4]);

			coord_t margin = itl2::round(3 * sigma) + 4;

			dfdx.ensureSize(img.dimensions());
			dfdy.ensureSize(img.dimensions());
			dfdz.ensureSize(img.dimensions());

			return Vec3c(margin, margin, margin);
		}
	};



	template<typename pixel_t, typename out_t> class GradientMagnitudeCommand : public OverlapDistributable<TwoImageInputOutputCommand<pixel_t, out_t> >
	{
	protected:
		friend class CommandList;

		GradientMagnitudeCommand() :
			OverlapDistributable<TwoImageInputOutputCommand<pixel_t, out_t> >("gradientmagnitude", "Calculates norm of Gaussian gradient of image. The output image can be the same than the input image.",
				{
					CommandArgument<double>(ParameterDirection::In, "spatial sigma", "Standard deviation of Gaussian derivative kernel.", 1.0),
					CommandArgument<double>(ParameterDirection::In, "gamma", "Scale-space scaling exponent. Set to zero to disable scaling.", 0.0)
				}
				)
		{
		}

	public:
		virtual void run(Image<pixel_t>& img, Image<out_t>& out, vector<ParamVariant>& args) const override
		{
			double std = pop<double>(args);
			double gamma = pop<double>(args);

			gradientMagnitude(img, out, std, gamma);
		}

		virtual Vec3c calculateOverlap(const vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& img = *std::get<DistributedImage<pixel_t>* >(args[0]);
			DistributedImage<out_t>& out = *std::get<DistributedImage<out_t>* >(args[1]);
			double std = std::get<double>(args[2]);
			Vec3c margin = itl2::round(3 * std) * Vec3c(1, 1, 1) + Vec3c(4, 4, 4);

			out.ensureSize(img.dimensions());

			return margin;
		}

		virtual double calculateExtraMemory(const vector<ParamVariant>& args) const override
		{
			return 3.0 * sizeof(typename NumberUtils<pixel_t>::FloatType) / (0.5 * (sizeof(pixel_t) + sizeof(out_t)));
		}
	};


	template<typename pixel_t> class HighpassFilterCommand : public OverlapDistributable<TwoImageInputOutputCommand<pixel_t> >
	{
	protected:
		friend class CommandList;

		HighpassFilterCommand() :
			OverlapDistributable<TwoImageInputOutputCommand<pixel_t> >("highpassfilter",
"Gaussian high-pass filtering. "
"Use to remove smooth, large-scale gray-scale variations from the image. "
"\n\n"
"Subtracts a Gaussian filtered version of input from itself. " + gaussianOptimizationHelp(),
				{
					CommandArgument<Vec3d>(ParameterDirection::In, "spatial sigma", "Standard deviation of Gaussian kernel."),
					CommandArgument<double>(ParameterDirection::In, "shift", "Constant to be added to the pixel values. Use to set the mean of the filtered image to a desired value. Useful especially when filtering unsigned images where non-shifted highpass filtering will lead to negative values that will be clipped to zero.", 0),
					CommandArgument<BoundaryCondition>(ParameterDirection::In, "boundary condition", string("Type of boundary condition. ") + boundaryConditionHelp(), BoundaryCondition::Nearest),
					CommandArgument<bool>(ParameterDirection::In, "allow optimization", "Set to true to allow optimized processing of images with high dynamic range data types. Might cause small deviations from true filtration result, and edge artefacts.", true)
				},
				filterSeeAlso()
				)
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, vector<ParamVariant>& args) const override
		{
			Vec3d std = pop<Vec3d>(args);
			double shift = pop<double>(args);
			BoundaryCondition bc = pop<BoundaryCondition>(args);
			bool opt = pop<bool>(args);

			highpassFilter(in, out, std, pixelRound<pixel_t>(shift), opt, bc);
		}

		virtual Vec3c calculateOverlap(const vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& in = *std::get<DistributedImage<pixel_t>* >(args[0]);
			DistributedImage<pixel_t>& out = *std::get<DistributedImage<pixel_t>* >(args[1]);
			Vec3d sigma = std::get<Vec3d>(args[2]);
			Vec3c margin = itl2::round(3 * sigma) + Vec3c(4, 4, 4);

			out.ensureSize(in.dimensions());

			return margin;
		}
	};



	template<typename pixel_t> class SatoLineFilterCommand : public OverlapDistributable<TwoImageInputOutputCommand<pixel_t> >
	{
	protected:
		friend class CommandList;

		SatoLineFilterCommand() :
			OverlapDistributable<TwoImageInputOutputCommand<pixel_t> >("satofilter", "Calculates line-enhancing filter lambda123 according to Sato - Three-dimensional multi-scale line filter for segmentation and visualization of curvilinear structures in medical images.",
				{
					CommandArgument<double>(ParameterDirection::In, "spatial sigma", "Standard deviation of Gaussian kernel, determines scale of structures that are probed."),
					CommandArgument<double>(ParameterDirection::In, "scale", "Output values are scaled by this number. Pass in zero to scale output value 1 to maximum of the data type for integer data types and to 1 for floating point data types.", 0),
					CommandArgument<double>(ParameterDirection::In, "gamma", "Scale-space scaling exponent. Set to zero to disable scaling.", 0),
					CommandArgument<double>(ParameterDirection::In, "gamma23", R"($\gamma_{23}$ parameter; controls the sharpness of the selectivity for the cross-section isotropy. Only non-negative values are valid.)", 1),
					CommandArgument<double>(ParameterDirection::In, "gamma12", R"($\gamma_{12}$ parameter. Non-negative values are valid.)", 1),
					CommandArgument<double>(ParameterDirection::In, "alpha", "Alpha parameter, must be between 0 and 1.", 0.5),
				},
				"frangifilter"
				)
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, vector<ParamVariant>& args) const override
		{
			double std = pop<double>(args);
			double outScale = pop<double>(args);
			double gamma = pop<double>(args);
			double gamma23 = pop<double>(args);
			double gamma12 = pop<double>(args);
			double alpha = pop<double>(args);
			
			lineFilter<pixel_t>(in, std, &out, gamma23, gamma12, alpha, 0, 0, 0, 0, gamma, outScale);
		}

		virtual Vec3c calculateOverlap(const vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& in = *std::get<DistributedImage<pixel_t>* >(args[0]);
			DistributedImage<pixel_t>& out = *std::get<DistributedImage<pixel_t>* >(args[1]);
			double sigma = std::get<double>(args[2]);
			coord_t margin = itl2::round(3 * sigma) + 4;

			out.ensureSize(in.dimensions());

			return Vec3c(margin, margin, margin);
		}

		virtual double calculateExtraMemory(const vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& in = *std::get<DistributedImage<pixel_t>* >(args[0]);
			return (6.0 * in.pixelCount() * sizeof(typename NumberUtils<pixel_t>::FloatType) / sizeof(pixel_t)) / (2.0 * in.pixelCount());
		}

		virtual JobType getJobType(const vector<ParamVariant>& args) const override
		{
			return JobType::Slow;
		}
	};

	template<typename pixel_t> class FrangiLineFilterCommand : public OverlapDistributable<TwoImageInputOutputCommand<pixel_t> >
	{
	protected:
		friend class CommandList;

		FrangiLineFilterCommand() :
			OverlapDistributable<TwoImageInputOutputCommand<pixel_t> >("frangifilter", "Calculates line-enhancing filter Vo according to Frangi - Multiscale vessel enhancement filtering.",
				{
					CommandArgument<double>(ParameterDirection::In, "spatial sigma", "Standard deviation of Gaussian kernel, determines scale of structures that are probed."),
					CommandArgument<double>(ParameterDirection::In, "output scale", "Output values are scaled by this number. Pass in zero to scale output value 1 to maximum of the data type for integer data types and to 1 for floating point data types.", 0),
					CommandArgument<double>(ParameterDirection::In, "c", "Sensitivity of the filter to deviation from background noise. Typical value is quarter of the value of the maximum intensity of the lines.", 0.25),
					CommandArgument<double>(ParameterDirection::In, "gamma", "Scale-space scaling exponent. Set to zero to disable scaling.", 0),
					CommandArgument<double>(ParameterDirection::In, "alpha", "Sensitivity of the filter to deviation from plate-like structures.", 0.5),
					CommandArgument<double>(ParameterDirection::In, "beta", "Sensitivity of the filter to deviation from blob-like structures.", 0.5)
				},
				"satofilter"
				)
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, vector<ParamVariant>& args) const override
		{
			double std = pop<double>(args);
			double outScale = pop<double>(args);
			double c = pop<double>(args);
			double gamma = pop<double>(args);
			double alpha = pop<double>(args);
			double beta = pop<double>(args);

			lineFilter<pixel_t>(in, std, 0, 0, 0, 0, &out, c, alpha, beta, gamma, outScale);
		}

		virtual Vec3c calculateOverlap(const vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& in = *std::get<DistributedImage<pixel_t>* >(args[0]);
			DistributedImage<pixel_t>& out = *std::get<DistributedImage<pixel_t>* >(args[1]);
			double sigma = std::get<double>(args[2]);
			coord_t margin = itl2::round(3 * sigma) + 4;

			out.ensureSize(in.dimensions());

			return Vec3c(margin, margin, margin);
		}

		virtual double calculateExtraMemory(const vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& in = *std::get<DistributedImage<pixel_t>* >(args[0]);
			return (6.0 * in.pixelCount() * sizeof(typename NumberUtils<pixel_t>::FloatType) / sizeof(pixel_t)) / (2.0 * in.pixelCount());
		}

		virtual JobType getJobType(const vector<ParamVariant>& args) const override
		{
			return JobType::Slow;
		}
	};



	inline std::string fftSeeAlso()
	{
		return "fft, ifft, bandpassfilter, highpassfilter";
	}


	class BandpassFilterCommand : public OneImageInPlaceCommand<float32_t>
	{
	protected:
		friend class CommandList;

		BandpassFilterCommand() :
			OneImageInPlaceCommand<float32_t>("bandpassfilter", "Bandpass filtering. Removes gray value variations smaller or larger in spatial extent than specified.",
				{
					CommandArgument<double>(ParameterDirection::In, "minimum size", "Variations smaller than than this value are removed.", 3.0),
					CommandArgument<double>(ParameterDirection::In, "maximum size", "Variations larger than this value are removed.", 40.0)
				},
				fftSeeAlso()
				)
		{
		}

	public:
		virtual void run(Image<float32_t>& img, vector<ParamVariant>& args) const override
		{
			double minSize = pop<double>(args);
			double maxSize = pop<double>(args);

			bandpassFilter(img, minSize, maxSize);
		}
	};

	class FFTCommand : public TwoImageInputOutputCommand<float32_t, complex32_t>
	{
	protected:
		friend class CommandList;

		FFTCommand() :
			TwoImageInputOutputCommand<float32_t, complex32_t>("fft", "Calculates Fourier transform.", {}, fftSeeAlso())
		{
		}

	public:
		virtual void run(Image<float32_t>& in, Image<complex32_t>& out, vector<ParamVariant>& args) const override
		{
			fft(in, out);
		}
	};

	class InverseFFTCommand : public TwoImageInputOutputCommand<complex32_t, float32_t>
	{
	protected:
		friend class CommandList;

		InverseFFTCommand() :
			TwoImageInputOutputCommand<complex32_t, float32_t>("ifft", "Calculates inverse Fourier transform. The size of the output image must be set to the size of the original data where the FFT was calculated from. The input image will be used as a temporary buffer and its contents are therefore destroyed.", {}, fftSeeAlso()
				)
		{
		}

	public:
		virtual void run(Image<complex32_t>& in, Image<float32_t>& out, vector<ParamVariant>& args) const override
		{
			ifft(in, out);
		}
	};
}
