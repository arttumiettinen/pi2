#pragma once

#include "command.h"
#include "distributable.h"
#include "commandsbase.h"
#include "danielsson.h"
#include "thickmap.h"
#include "dmapcommands.h"
#include "overlapdistributable.h"
#include "io/imagedatatype.h"
#include "specialcommands.h"
#include "pisystem.h"
#include "projectioncommands.h"
#include "math/aabox.h"
#include "pilibutilities.h"
#include "timing.h"

#include <vector>
#include <type_traits>

namespace pilib
{
	template<typename pixel_t> class Danielsson2Command : public OverlapDistributable<TwoImageInputOutputCommand<pixel_t> >
	{
	protected:
		friend class CommandList;

		Danielsson2Command() : OverlapDistributable<TwoImageInputOutputCommand<pixel_t> >(
			"danielsson2",
			"Calculates an estimate of the set of centers of locally maximal spheres from squared Euclidean distance map. Uses Danielsson's algorithm. The output image can be interpreted as a medial axis or a distance ridge. Drawing a sphere on each nonzero pixel in the output image, with radius and color equal to square root of pixel value, results in a local thickness map of the structure. See e.g. Yaorong Ge and J. Michael Fitzpatrick - On the Generation of Skeletons from Discrete Euclidean Distance Maps, IEEE Transactions on Pattern Analysis and Machine Intelligence 18, 1996; and P.-E. Danielsson, Euclidean distance mapping, Computer Graphics and Image Processing 14(3), 1980.",
			{
			},
			"dmap2, dmap, tmap")
		{
		}

	public:
		using Distributable::runDistributed;

		void run(Image<pixel_t>& in, Image<pixel_t>& out) const
		{
			centersOfLocallyMaximalSpheres(in, out);
		}

		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, std::vector<ParamVariant>& args) const override
		{
			run(in, out);
		}

		virtual Vec3c calculateOverlap(const std::vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& input = *std::get<DistributedImage<pixel_t>* >(args[0]);
			DistributedImage<pixel_t>& output = *std::get<DistributedImage<pixel_t>* >(args[1]);
			output.ensureSize(input);
			return Vec3c(1, 1, 1);
		}
	};

	
	
	
	template<typename pixel_t> class DrawSpheres2ProcessDimensionCommand : public TwoImageInputOutputCommand<pixel_t>, public Distributable
	{
	protected:
		friend class CommandList;

		DrawSpheres2ProcessDimensionCommand() : TwoImageInputOutputCommand<pixel_t>("processtmapdimension", "Performs processing related to one dimension in thickness map sphere plotting phase. This command is used internally in distributed processing; consider using `tmap` command instead of this one. The input image must be squared distance map, and the output image is the squared local radius map.",
			{
				CommandArgument<size_t>(ParameterDirection::In, "dimension", "Dimension to process."),
				CommandArgument<string>(ParameterDirection::In, "ri prefix", "Prefix for temporary data files."),
				CommandArgument<string>(ParameterDirection::In, "full dmap2 file", "Prefix of the full dmap2 file."),
				CommandArgument<Vec3c>(ParameterDirection::In, "full dimensions", "Dimensions of the fulle dmap2 image."),
				CommandArgument<Distributor::BLOCK_ORIGIN_ARG_TYPE>(ParameterDirection::In, Distributor::BLOCK_ORIGIN_ARG_NAME, "Origin of current block in coordinates of the full image. This argument is used internally in distributed processing to load correct block of temporary data.", Distributor::BLOCK_ORIGIN_ARG_TYPE(0, 0, 0)),
				CommandArgument<Distributor::BLOCK_INDEX_ARG_TYPE>(ParameterDirection::In, Distributor::BLOCK_INDEX_ARG_NAME, "Index of image block that we are currently processing. This argument is used internally in distributed processing and should be set to the index of the block the system is currently processing.", 0),
				CommandArgument<double>(ParameterDirection::In, "mean r", "Mean value of non-zero points in the distance map.")
			})
		{
		}

		pixel_t maxR2(const Image<itl2::internals::RiStorageSet>& ri, const Image<pixel_t>& dmap2, const Image<pixel_t>& dmap2Full, const Vec3c& blockPos) const
		{
			pixel_t M = 0;

			#pragma omp parallel if(ri.pixelCount() > PARALLELIZATION_THRESHOLD)
			{
				pixel_t Mlocal = 0;

				#pragma omp for nowait
				for (coord_t z = 0; z < ri.depth(); z++)
				{
					for (coord_t y = 0; y < ri.height(); y++)
					{
						for (coord_t x = 0; x < ri.width(); x++)
						{
							Vec3c pos(x, y, z);
							itl2::internals::RiSuperSet<pixel_t> set;
							itl2::internals::toRiSet(ri(pos), set, pos, dmap2, Vec3sc(blockPos), dmap2Full);

							for (const auto& item : set)
							{
								if (item.R2 > Mlocal)
									Mlocal = item.R2;
							}
						}
					}
				}

				#pragma omp critical(maxR2)
				{
					if (Mlocal > M)
						M = Mlocal;
				}
			}

			return M;
		}

	public:
		virtual bool isInternal() const override
		{
			return true;
		}

		using Distributable::runDistributed;

		virtual void run(Image<pixel_t>& dmap2, Image<pixel_t>& tmap2, std::vector<ParamVariant>& args) const override
		{
			size_t dim = pop<size_t>(args);
			string riPrefix = pop<string>(args);
			string dmap2FullFile = pop<string>(args);
			Vec3c fullDimensions = pop<Vec3c>(args);
			Distributor::BLOCK_ORIGIN_ARG_TYPE blockOrigin = pop<Distributor::BLOCK_ORIGIN_ARG_TYPE>(args);
			Distributor::BLOCK_INDEX_ARG_TYPE blockIndex = pop<Distributor::BLOCK_INDEX_ARG_TYPE>(args);
			double meanr = pop<double>(args);

			// Map full image
			string prefix = itl2::getPrefix(dmap2FullFile);
			Image<pixel_t> dmap2Full(prefix, true, fullDimensions);


			// Initialize ri image
			Image<itl2::internals::RiStorageSet> ri(dmap2.dimensions());
			pixel_t M = 0;
			if (dim > 0)
			{
				// Read ri from previous dimension output
				{
					TimingFlag flag(TimeClass::IO);
					itl2::internals::readRiBlock(ri, riPrefix + "_dim" + itl2::toString(dim - 1), blockOrigin);
				}
				M = maxR2(ri, dmap2, dmap2Full, blockOrigin);
			}
			else
			{
				// Initialize ri from dmap2
				for (coord_t z = 0; z < ri.depth(); z++)
				{
					for (coord_t y = 0; y < ri.height(); y++)
					{
						for (coord_t x = 0; x < ri.width(); x++)
						{
							Vec3c p(x, y, z);
							pixel_t R2 = dmap2(p);
							if (R2 > 0)
							{
								ri(p).setValues(itl2::internals::RiStorageItem{ (int16_t)(x + blockOrigin.x), (int16_t)(y + blockOrigin.y) });
								if (R2 > M)
									M = R2;
							}
						}
					}
				}
			}

			// Build circle lookup
			if (intuitive::ge(M, std::numeric_limits<int32_t>::max()))
				throw ITLException("The squared distance map contains too large values. (This error is easily avoidable by changing buildCircleLookup functionality.)");
			itl2::internals::buildCircleLookup((int32_t)M);

			// Process
			// TODO: If reading from memory mapped file becomes a problem, we can make a version that saves more data to the temp files, but does not need dmap2Full.
			itl2::internals::processDimensionSuper<pixel_t>(ri, dim, dmap2, tmap2, AABox<coord_t>::fromMinMax(blockOrigin, dmap2.dimensions()), dmap2Full);


			// Write temporary file, if any
			if (dim < getDimensionality(fullDimensions) - 1)
			{
				TimingFlag flag(TimeClass::IO);
				itl2::internals::writeRiBlock(ri, riPrefix + "_dim" + itl2::toString(dim), (uint16_t)blockIndex, blockOrigin, fullDimensions);
			}
		}

		virtual std::vector<string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& in = *std::get<DistributedImage<pixel_t>* >(args[0]);
			DistributedImage<pixel_t>& out = *std::get<DistributedImage<pixel_t>* >(args[1]);
			out.ensureSize(in);
			return distributor.distribute(this, args);
		}

		
		virtual size_t getPreferredSubdivisions(const std::vector<ParamVariant>& args) const override
		{
			// Same number of subdivisions than smallest number of jobs in the disk-space saving algorithm.
			// This is done merely to make the benchmarks run in comparable time.
			Vec3c fullDimensions = std::get<Vec3c>(args[5]);
			size_t ddir = getDistributionDirection1(args);

			// Make sure that we don't divide too much for small images.
			size_t divs = 8;
			while (divs > 1 && fullDimensions[ddir] / divs < 5)
				divs--;

			return std::max((size_t)1, divs);
		}
		

		virtual size_t getDistributionDirection1(const std::vector<ParamVariant>& args) const override
		{
			size_t dim = std::get<size_t>(args[2]);

			switch (dim)
			{
			case 0: return 2;
			case 1: return 2;
			case 2: return 1;
			default: throw ITLException("Unsupported dimension.");
			}
		}

		virtual double calculateExtraMemory(const std::vector<ParamVariant>& args) const override
		{
			double meanr = std::get<double>(args.back());

			// Derivation for formula below:
			// (Size of inputs) * (1 + extraMemFactor) = (required memory)
			// ((input size) + (output size)) * (1 + extraMemFactor) = (required memory)
			// (2 * dmap2.pixelCount()) * (1 + extraMemFactor) = (4.5 + 0.2 * meanr) * dmap2.pixelCount()
			// 2 * (1 + extraMemFactor) = (4.5 + 0.2 * meanr)
			// 1 + extraMemFactor = (4.5 + 0.2 * meanr) / 2
			// extraMemFactor = (4.5 + 0.2 * meanr) / 2 - 1
			// finally, add safety factor of 2 to end in
			//extraMemFactor = (4.5 + 0.2 * meanr) - 1

			return (4.5 + 0.2 * meanr) - 1;
		}
	};
	
	
	
	


	template<typename pixel_t> class DrawSpheres2Command : public TwoImageInputOutputCommand<pixel_t>, public Distributable
	{
	protected:
		friend class CommandList;

			DrawSpheres2Command() : TwoImageInputOutputCommand<pixel_t>("drawspheres2", "Draws a sphere centered in each non-zero point of the input image. Radius of the sphere is given by the square root of the value of the pixel. Color is given by the value of the pixel. This command is used to generate squared local radius map from squared distance map or output of `danielsson2` command.",
				{
					CommandArgument<bool>(ParameterDirection::In, "save memory", "For non-distributed processing: Set to true to use slower algorithm that uses less memory than the faster one. If set to true, the input and output images can be the same. For distributed processing: Only false is supported at the moment.", false),
					CommandArgument<Vec3c>(ParameterDirection::In, "block size", "Block size to use in distributed processing. Set to zero to use a default value calculated based on available memory. If calculation with default value fails, use smaller block size and report it to the authors. This argument has no effect when running in non-distributed mode.", Vec3c(0, 0, 0))
				},
				"tmap")
		{
		}

		double calcNonZeroMeanRDistributed(Distributor& distributor, DistributedImage<pixel_t>& dmap2) const
		{
			// Create temporary image to hold the result
			DistributedTempImage<float32_t> resultImgPointer(distributor, "dmap2_mean", 1, DistributedImageStorageType::Raw);
			DistributedImage<float32_t>& resultImg = resultImgPointer.get();

			// Run masked mean and grab the result
			CommandList::get<MaskedMeanAllPixelsCommand<pixel_t> >().runDistributed(distributor, { &dmap2, &resultImg, 0.0, false, true });
			float32_t result = resultImg.getValue();

			return result;
		}

	public:
		using Distributable::runDistributed;

		void run(Image<pixel_t>& in, Image<pixel_t>& out, bool saveMemory) const
		{
			if (!saveMemory)
			{
				dimredsuper::thickmap2<pixel_t>(in, out, nullptr, nullptr);
			}
			else
			{
				setValue(out, in);
				optimized::thickmap2(out, nullptr);
			}
		}

		virtual void run(Image<pixel_t>& in, Image<pixel_t>& out, std::vector<ParamVariant>& args) const override
		{
			bool saveMemory = pop<bool>(args);
			Vec3c blockSize = pop<Vec3c>(args);

			run(in, out, saveMemory);
		}

		virtual std::vector<string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& dmap2 = *pop<DistributedImage<pixel_t>* >(args);
			DistributedImage<pixel_t>& result = *pop<DistributedImage<pixel_t>* >(args);
			bool saveMemory = pop<bool>(args);
			Vec3c blockSize = pop<Vec3c>(args);

			distributor.flush();

			//if (!dmap2.isRaw())
			if(dmap2.currentReadSourceType() != DistributedImageStorageType::Raw)
				throw ITLException("The squared distance map must be saved to a .raw file (the file must be memory mappable).");

			dmap2.mustNotBe(result);

			result.ensureSize(dmap2);

			if (!saveMemory)
			{

				// Fast algorithm, large disk space requirement.

				double meanr = calcNonZeroMeanRDistributed(distributor, dmap2);


				// DEBUG
				std::cout << "Mean r' = " << meanr << std::endl;


				string riPrefix = createTempFilename("ri");
				for (size_t dim = 0; dim < dmap2.dimensionality(); dim++)
				{
					std::cout << "Process dimension " << (dim + 1) << "..." << std::endl;

					// Process this dimension
					CommandList::get<DrawSpheres2ProcessDimensionCommand<pixel_t> >().runDistributed(distributor, { &dmap2, &result, dim, riPrefix, dmap2.currentReadSource(), dmap2.dimensions(), Distributor::BLOCK_ORIGIN_ARG_TYPE(), Distributor::BLOCK_INDEX_ARG_TYPE(), meanr });

					std::cout << "Remove unnecessary temporary files..." << std::endl;

					// Delete temporary files from previous round
					{
						TimingFlag flag(TimeClass::IO);
						auto items = itl2::buildFileList(riPrefix + "_dim" + itl2::toString(dim - 1) + "_*");
						for (const string& file : items)
							itl2::deleteFile(file);
					}
				}
			}
			else
			{
				// TODO: Maybe we want to implement the old Python script version here?
				throw ITLException("Memsave-version of this algorithm is currently not available in the distributed processing mode.");
			}

			return std::vector<string>();
		}
	};










	template<typename pixel_t, typename out_t> class FinalizeThicknessMapCommand : public InputOutputPointProcess<pixel_t, out_t>
	{
	protected:
		friend class CommandList;

		FinalizeThicknessMapCommand() : InputOutputPointProcess<pixel_t, out_t>("finalizetmap", "Converts squared local radius map to local thickness map.",
			{
			},
			"tmap")
		{
		}

	public:
		using Distributable::runDistributed;

		void run(Image<pixel_t>& in, Image<out_t>& out) const
		{
			finalizeThickmap(in, out);
		}

		virtual void run(Image<pixel_t>& in, Image<out_t>& out, std::vector<ParamVariant>& args) const override
		{
			run(in, out);
		}
	};



	template<typename pixel_t> class RoundDistanceRidge2Command : public InPlacePointProcess<pixel_t>
	{
	protected:
		friend class CommandList;

		RoundDistanceRidge2Command() : InPlacePointProcess<pixel_t>("rounddistanceridge2", "Rounds distance values to nearest integer. Inputs and outputs squared distance values.", {}, "tmap")
		{
		}

	public:

		void run(Image<pixel_t>& in) const
		{
			roundDistanceRidge2(in);
		}

		virtual void run(Image<pixel_t>& in, std::vector<ParamVariant>& args) const override
		{
			run(in);
		}
	};



	template<typename pixel_t, typename out_t> class ThicknessMapCommand : public Command, public Distributable
	{
	protected:
		friend class CommandList;

		ThicknessMapCommand() : Command("tmap",
"Calculates local thickness map (i.e. opening transform) of a binary image. "
"Input image is used as temporary storage space. "
"The pixel data type must be able to contain large enough values; usually uint8 is too small. Uint32 or uint64 are recommended. "
"\n\n"
"The algorithm selection depends on the value of the 'save memory' argument. "
"If 'save memory' is true, the algorithm introduced in Hildebrand - A New Method for the Model-Independent Assessment of Thickness in Three-Dimensional Images is used."
"If 'save memory' is false, the separable algorithm in Lovric - Separable distributed local thickness algorithm for efficient morphological characterization of terabyte-scale volume images is applied. "
"The separable algorithm is usually much faster than the Hildebrand algorithm, but requires much more RAM (normal mode) or temporary disk space (distributed mode). ",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "input image", "Input image where background is marked with background value given by the third argument."),
				CommandArgument<Image<out_t> >(ParameterDirection::Out, "output image", "Output image (thickness map) where pixel value equals diameter of object at that location, i.e. diameter of the largest sphere that fits in foreground points and contains the object."),
				CommandArgument<double>(ParameterDirection::In, "background value", "Pixels belonging to the background are marked with this value in the input image.", 0),
				CommandArgument<bool>(ParameterDirection::In, "round distance values", "Set to true to generate approximate distance map by rounding distance values to nearest integers. Calculation of rounded distance map is faster, but the result is only approximation both in distance values and shape of structures.", false),
				CommandArgument<bool>(ParameterDirection::In, "save memory", "For non-distributed processing: Set to true to use slower algorithm (Hildebrand) that uses less memory than the faster one (separable). If set to true, the input and output images can be the same. For distributed processing: Only false is supported at the moment.", false),
				CommandArgument<Vec3c>(ParameterDirection::In, "block size", "Block size to use in distributed processing in sphere plotting phase. Set to zero to use a default value calculated based on available memory. If calculation with default value fails, use smaller block size and report it to the authors. This argument has no effect when running in non-distributed mode.", Vec3c(0, 0, 0))
			},
			"dmap2, dmap")
		{
		}

	public:
		using Distributable::runDistributed;

		virtual void run(std::vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& input = *pop<Image<pixel_t>* >(args);
			Image<out_t>& output = *pop<Image<out_t>* >(args);
			double bgval = pop<double>(args);
			bool roundedTmap = pop<bool>(args);
			bool saveMemory = pop<bool>(args);
			Vec3c blockSize = pop<Vec3c>(args);

			{
				Image<pixel_t> temp(input.dimensions());
				CommandList::get<DistanceMap2Command<pixel_t, pixel_t> >().run(input, input, bgval);
				CommandList::get<Danielsson2Command<pixel_t> >().run(input, temp);
				if(roundedTmap)
					CommandList::get<RoundDistanceRidge2Command<pixel_t> >().run(temp);
				CommandList::get<DrawSpheres2Command<pixel_t> >().run(temp, input, saveMemory);
			} // temp image is deleted from memory here, and output image is created below if it does not exist.
			CommandList::get<FinalizeThicknessMapCommand<pixel_t, out_t> >().run(input, output);
		}

		virtual std::vector<string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& input = *pop<DistributedImage<pixel_t>* >(args);
			DistributedImage<out_t>& output = *pop<DistributedImage<out_t>* >(args);
			double bgval = pop<double>(args);
			bool roundedTmap = pop<bool>(args);
			bool saveMemory = pop<bool>(args);
			Vec3c blockSize = pop<Vec3c>(args);

			output.ensureSize(input);

			{
				// Distance map
				CommandList::get<DistanceMap2Command<pixel_t, pixel_t> >().runDistributed(distributor, { &input, &input, bgval });

				// Create temporary image with Raw (memory mappable) storage type
				DistributedTempImage<pixel_t> tempImg(distributor, "thickmap_temp", input.dimensions(), DistributedImageStorageType::Raw);
				DistributedImage<pixel_t>& temp = tempImg.get();

				// Danielsson
				CommandList::get<Danielsson2Command<pixel_t> >().runDistributed(distributor, { &input, &temp });

				// Rounding?
				if (roundedTmap)
					CommandList::get<RoundDistanceRidge2Command<pixel_t> >().runDistributed(distributor, { &temp });

				// Propagate
				CommandList::get<DrawSpheres2Command<pixel_t> >().runDistributed(distributor, { &temp, &input, saveMemory, blockSize });
			}
			// Finalize
			CommandList::get<FinalizeThicknessMapCommand<pixel_t, out_t> >().runDistributed(distributor, { &input, &output });

			return std::vector<string>();
		}
	};

}

