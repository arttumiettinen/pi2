#pragma once

#include "command.h"
#include "commandsbase.h"
#include "commandlist.h"
#include "distributable.h"
#include "pointprocesscommands.h"

#include "dmap.h"
#include "sdmap.h"

namespace pilib
{

	inline std::string dmapSeeAlso()
	{
		return "dmap, dmap2, tmap, sdmap";
	}

	template<typename pixel_t, typename dmap_t> class PrepareDistanceMapCommand : public InputOutputPointProcess<pixel_t, dmap_t>
	{
	protected:
		friend class CommandList;

		PrepareDistanceMapCommand() : InputOutputPointProcess<pixel_t, dmap_t>("preparedmap", "Prepares image for distance map calculation. This command is used internally in distributed processing; consider using dmap or dmap2 commands instead of this one.",
			{
				CommandArgument<double>(ParameterDirection::In, "background value", "Pixels belonging to the background are marked with this value in the input image.", 0)
			})
		{
		}

	public:

		virtual bool isInternal() const override
		{
			return true;
		}

		using Distributable::runDistributed;

		virtual void run(Image<pixel_t>& in, Image<dmap_t>& out, vector<ParamVariant>& args) const override
		{
			double bgval = pop<double>(args);

			out.ensureSize(in);

			prepareDistanceTransform(in, out, pixelRound<pixel_t>(bgval));
		}

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& input = *std::get<DistributedImage<pixel_t>* >(args[0]);
			DistributedImage<dmap_t>& output = *std::get<DistributedImage<dmap_t>* >(args[1]);
			
			output.ensureSize(input);

			return distributor.distribute(this, args);
		}
	};


	template<typename dmap_t> class DistanceMapProcessDimensionCommand : public OneImageInPlaceCommand<dmap_t>, public Distributable
	{
	protected:
		friend class CommandList;

		DistanceMapProcessDimensionCommand() : OneImageInPlaceCommand<dmap_t>("processdmapdimension", "Performs processing related to one dimension in distance map calculation. The full distance map is calculated by first calling preparedmap command and then this command for each dimension of the image. This command is used internally in distributed processing; consider using dmap or dmap2 commands instead of this one.",
			{
				CommandArgument<size_t>(ParameterDirection::In, "dimension", "Dimension to process.")
			})
		{
		}

	public:
		virtual bool isInternal() const override
		{
			return true;
		}

		using Distributable::runDistributed;

		virtual void run(Image<dmap_t>& img, vector<ParamVariant>& args) const override
		{
			size_t dim = pop<size_t>(args);

			itl2::internals::processDimension(img, dim, nullptr, true);
		}

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override
		{
			return distributor.distribute(this, args);
		}

		virtual size_t getDistributionDirection1(const vector<ParamVariant>& args) const override
		{
			size_t dim = std::get<size_t>(args[1]);

			switch (dim)
			{
			case 0: return 2;
			case 1: return 2;
			case 2: return 1;
			default: throw ITLException("Unsupported dimension.");
			}
		}
	};





	
	template<typename pixel_t, typename dmap_t> class DistanceMap2Command : public Command, public Distributable
	{
	protected:
		friend class CommandList;

		DistanceMap2Command() : Command("dmap2", "Calculates squared distance map or squared distance transform of a binary image. In order to calculate (non-squared) distance map, use the `dmap` command or take a square root of the result using the `squareroot` command.  Uses algorithm from Maurer - A Linear Time Algorithm for Computing Exact Euclidean Distance Transforms of Binary Images in Arbitrary Dimensions.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "input image", "Input image where background is marked with background value given by the third argument."),
				CommandArgument<Image<dmap_t> >(ParameterDirection::Out, "output image", "Output image (squared distance map) where pixel value equals squared Euclidean distance to the nearest background pixel. Input and output can be the same image if in-place transformation is preferred and input data type is suitable. For exact result, use integer type. Please note that squared values require relatively large pixel data type, e.g. int32, depending on the magnitude of distances in the image. If floating point pixel type is used, the results might contain artifacts for very large images due to floating point inaccuracy. If that is encountered, consider calculating squared distance map using integer data type and then converting to floating point format, possibly followed by square root operation. That is, dmap2(img, img); convert(img, float32); squareroot(img); Notice that this sequence may require more memory than the standard operation of this command because of the explicit conversion step."),
				CommandArgument<double>(ParameterDirection::In, "background value", "Pixels belonging to the background are marked with this value in the input image.", 0)
			},
			dmapSeeAlso())
		{
		}

	public:
		using Distributable::runDistributed;

		void run(Image<pixel_t>& input, Image<dmap_t>& output, double bgval) const
		{
			distanceTransform2(input, output, nullptr, pixelRound<pixel_t>(bgval));
		}

		virtual void run(vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& input = *pop<Image<pixel_t>* >(args);
			Image<dmap_t>& output = *pop<Image<dmap_t>* >(args);
			double bgval = pop<double>(args);

			run(input, output, bgval);
		}

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& input = *pop<DistributedImage<pixel_t>* >(args);
			DistributedImage<dmap_t>& output = *pop<DistributedImage<dmap_t>* >(args);
			double bgval = pop<double>(args);

			output.ensureSize(input);

			// Prepare
			CommandList::get<PrepareDistanceMapCommand<pixel_t, dmap_t> >().runDistributed(distributor, { &input, &output, bgval });

			// Calculate in each dimension
			for (size_t n = 0; n < output.dimensionality(); n++)
			{
				CommandList::get<DistanceMapProcessDimensionCommand<dmap_t> >().runDistributed(distributor, { &output, n });
			}

			return vector<string>();
		}
	};

	template<typename pixel_t, typename dmap_t> class DistanceMapCommand : public Command, public Distributable
	{
	protected:
		friend class CommandList;

		DistanceMapCommand() : Command("dmap", "Calculates distance map or distance transform of a binary image. In order to calculate squared distance map, use `dmap2` command. Uses algorithm from Maurer - A Linear Time Algorithm for Computing Exact Euclidean Distance Transforms of Binary Images in Arbitrary Dimensions.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "input image", "Input image where background is marked with background value given by the third argument."),
				CommandArgument<Image<dmap_t> >(ParameterDirection::Out, "output image", "Output image (distance map) where pixel value equals Euclidean distance to the nearest background pixel. Input and output can be the same image if in-place transformation is preferred and input data type is suitable. Integer pixel data types result in approximate distance map only. If floating point pixel type is used, the results might contain artifacts for very large images due to floating point inaccuracy. If that is encountered, consider calculating squared distance map using integer data type and then converting to floating point format, followed by square root operation. That is, dmap2(img, img); convert(img, float32); squareroot(img); Notice that this sequence may require more memory than the standard operation of this command because of the explicit conversion step."),
				CommandArgument<double>(ParameterDirection::In, "background value", "Pixels belonging to the background are marked with this value in the input image.", 0)
			},
			dmapSeeAlso())
		{
		}

	public:
		using Distributable::runDistributed;

		virtual void run(vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& input = *pop<Image<pixel_t>* >(args);
			Image<dmap_t>& output = *pop<Image<dmap_t>* >(args);
			double bgval = pop<double>(args);

			distanceTransform(input, output, nullptr, pixelRound<pixel_t>(bgval));
		}

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override
		{
			DistributedImage<dmap_t>& output = *std::get<DistributedImage<dmap_t>* >(args[1]);

			// Squared distance map
			CommandList::get<DistanceMap2Command<pixel_t, dmap_t> >().runDistributed(distributor, args);

			// Square root
			CommandList::get<SquareRootCommand<float32_t> >().runDistributed(distributor, { &output });

			return vector<string>();
		}
	};



	template<typename pixel_t> class SeededDistanceMapCommand : public Command
	{
	protected:
		friend class CommandList;

		SeededDistanceMapCommand() : Command("sdmap", "Calculates seeded distance map of a binary image. For each pixel in the image, finds the shortest pixel-by-pixel path to nearest seed region such that the path does not pass through regions that are marked as obstacles. The output is the length of the shortest path, calculated as the number of pixel-to-pixel steps that must be taken to follow the path. This is in contrast to the `dmap` command where the distance to the seed region is direct Euclidean distance and no separate obstacles are supported. Please note that comparing the output of `dmap` and `sdmap` commands might lead to surprising results due to the different definition of distance.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "seeds", "Seed image that contains the set of pixels where the distance is zero. The set is marked with nonzero values, i.e. the distance map will propagate to pixels that have zero value in this image. This image is not modified."),
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "geometry", "Image containing the geometry in which to calculate the distances. Regions marked with zero pixel value are 'unpassable obstacles'. The distance transform will only proceed to pixels whose color in this image matches the color of the seed point (in this image). This image is not modified."),
				CommandArgument<Image<float32_t> >(ParameterDirection::Out, "output", "Output image that will contain the distance to the nearest seed region, not passing through zero regions in the geometry image."),
				CommandArgument<Connectivity>(ParameterDirection::In, "connectivity", string("Connectivity of the distance map. ") + connectivityHelp(), Connectivity::AllNeighbours),
			},
			dmapSeeAlso())
		{
		}

	public:

		virtual void run(vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& seeds = *pop<Image<pixel_t>* >(args);
			Image<pixel_t>& geometry = *pop<Image<pixel_t>* >(args);
			Image<float32_t>& output = *pop<Image<float32_t>* >(args);
			Connectivity connectivity = pop<Connectivity>(args);

			seededDistanceMap(seeds, geometry, output, connectivity);
		}
	};



	template<typename pixel_t> class PolySeededDistanceMapCommand : public Command
	{
	protected:
		friend class CommandList;

		PolySeededDistanceMapCommand() : Command("polysdmap", "Calculates seeded distance map of a binary image. For each pixel in the image, finds the shortest segmented polygonal path to nearest seed region such that the path does not pass through regions that are marked as obstacles, and ensuring that each segment that makes up the path is shorter than specified in 'max segment length' argument. This command is highly experimental.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "seeds", "Seed image that contains the set of pixels where the distance is zero. The set is marked with nonzero values, i.e. the distance map will propagate to pixels that have zero value in this image. This image is not modified."),
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "geometry", "Image containing the geometry in which to calculate the distances. Regions marked with zero pixel value are 'unpassable obstacles'. The distance transform will only proceed to pixels whose color in this image matches the color of the seed point (in this image). This image is not modified."),
				CommandArgument<Image<float32_t> >(ParameterDirection::Out, "output", "Output image that will contain the distance to the nearest seed region, not passing through zero regions in the geometry image."),
				CommandArgument<double>(ParameterDirection::In, "max segment length", "Maximum length of segments that make up the path to the nearest seed point."),
			},
			dmapSeeAlso())
		{
		}

	public:

		virtual void run(vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& seeds = *pop<Image<pixel_t>* >(args);
			Image<pixel_t>& geometry = *pop<Image<pixel_t>* >(args);
			Image<float32_t>& output = *pop<Image<float32_t>* >(args);
			float32_t l = (float32_t)pop<double>(args);

			seededDistanceMapPolyline(seeds, geometry, output, l);
		}
	};


}
