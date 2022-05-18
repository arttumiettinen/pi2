#pragma once

#include "floodfill.h"
#include "command.h"
#include "commandsbase.h"
#include "distributable.h"
#include "pilibutilities.h"
#include "io/vectorio.h"
#include "pointprocesscommands.h"

namespace pilib
{
	inline std::string floodfillSeeAlso()
	{
		return "grow, growlabels, floodfill, regionremoval, morphorec";
	}

	


	template<typename label_t, typename weight_t> class GrowPriorityCommand : public TwoImageInputParamCommand<label_t, weight_t>
	{
	protected:
		friend class CommandList;

		GrowPriorityCommand() : TwoImageInputParamCommand<label_t, weight_t>("grow",
			"Grows regions from seed points outwards. Seeds points are all nonzero pixels in the input image, pixel value defining region label. Each seed is grown towards surrounding zero pixels. Fill priority for each pixel is read from the corresponding pixel in the parameter image. Pixels for which priority is zero or negative are never filled. This process is equal to Meyer's watershed algorithm for given set of seeds, and watershed cuts are borders between filled regions in the output image.",
			{},
			floodfillSeeAlso())
		{
		}

	public:
		virtual void run(Image<label_t>& labels, Image<weight_t>& weights, std::vector<ParamVariant>& args) const override
		{
			grow(labels, weights);
		}
	};

	template<typename pixel_t> class GrowCommand : public OneImageInPlaceCommand<pixel_t>, public Distributable
	{
	protected:
		friend class CommandList;

		GrowCommand() : OneImageInPlaceCommand<pixel_t>("grow", "Grows regions with source color to regions with target color as much as possible.",
			{
				CommandArgument<double>(ParameterDirection::In, "source color", "Color that defines regions that are going to be grown."),
				CommandArgument<double>(ParameterDirection::In, "target color", "Color where the regions will grow."),
				CommandArgument<Connectivity>(ParameterDirection::In, "connectivity", string("Connectivity of the regions to grow. ") + connectivityHelp(), Connectivity::NearestNeighbours),
			},
			floodfillSeeAlso())
		{
		}

	public:
		using Distributable::runDistributed;

		virtual void run(Image<pixel_t>& in, std::vector<ParamVariant>& args) const override
		{
			double src = pop<double>(args);
			double target = pop<double>(args);
			Connectivity connectivity = pop<Connectivity>(args);

			size_t changed = grow(in, pixelRound<pixel_t>(src), pixelRound<pixel_t>(target), connectivity);
			std::cout << std::endl << changed << " pixels changed." << std::endl;
		}

		virtual Vec3c getMargin(const std::vector<ParamVariant>& args) const override
		{
			return Vec3c(3, 3, 3);
		}

		virtual JobType getJobType(const std::vector<ParamVariant>& args) const override
		{
			return JobType::Normal;
		}

		virtual double calculateExtraMemory(const std::vector<ParamVariant>& args) const override
		{
			// Allocate some extra memory for priority queue
			return 1.0;
		}

		virtual std::vector<std::string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const override
		{
			// TODO: This could be made with IterableDistributable
			coord_t changed;
			do
			{
				std::vector<std::string> output = distributor.distribute(this, args);

				changed = parseTotalCount(output, "pixels changed.");
				std::cout << std::endl << changed << " pixels changed." << std::endl;
			} while (changed > 0);

			return std::vector<std::string>();
		}
	};


	template<typename pixel_t> class DualThresholdCommand : public OneImageInPlaceCommand<pixel_t>, public Distributable
	{
	protected:
		friend class CommandList;

		DualThresholdCommand() : OneImageInPlaceCommand<pixel_t>("dualthreshold", "First sets all pixels with value over upper threshold to 1. Then sets all regions to 1 that have value over lower threshold and that are connected to region that has value over upper threshold.",
			{
				CommandArgument<double>(ParameterDirection::In, "lower threshold", "Regions that have value below lower threshold value are discarded. Regions that have value between lower and upper thresholds are included in the result only if they touch some region that has value above upper threshold."),
				CommandArgument<double>(ParameterDirection::In, "upper threshold", "Regions that have value above upper threshold value are always included in the result. Regions that have value between lower and upper thresholds are included in the result only if they touch some regoin that has value above upper threshold.")
			})
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, std::vector<ParamVariant>& args) const override
		{
			double loThreshold = pop<double>(args);
			double hiThreshold = pop<double>(args);

			dualThreshold(in, pixelRound<pixel_t>(loThreshold), pixelRound<pixel_t>(hiThreshold));
		}

		virtual std::vector<std::string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& img = *pop<DistributedImage<pixel_t>* >(args);
			double loThreshold = pop<double>(args);
			double hiThreshold = pop<double>(args);

			// Multi-threshold to two classes.
			auto& mt = CommandList::get<DoubleThresholdCommand<pixel_t> >();
			mt.runDistributed(distributor, { &img, loThreshold, hiThreshold });

			// Convert all those structures to "sure" that touch a "sure" structure.
			auto& grow = CommandList::get<GrowCommand<pixel_t> >();
			grow.runDistributed(distributor, { &img, 2.0, 1.0, Connectivity::NearestNeighbours });

			// Threshold so that only "sure" structures are left.
			auto& th = CommandList::get<ThresholdConstantCommand<pixel_t> >();
			th.runDistributed(distributor, { &img, 1.0 });

			return std::vector<std::string>();
		}
	};



	template<typename pixel_t> class GrowLabelsCommand : public OneImageInPlaceCommand<pixel_t>, public Distributable
	{
	protected:
		friend class CommandList;

		GrowLabelsCommand() : OneImageInPlaceCommand<pixel_t>("growlabels",
			"Grows all colored regions as much as possible into pixels that have a specific color. "
			"In practice, this command first finds all unique colors in the image, and uses each set of "
			"pixels having the same color as seed points for a flood fill that proceeds to pixels whose value is given in the 'allowed color' argument. "
			"\n\n"
			"This growing method is suited only for situations where separate parts of the original structure are labelled and "
			"the labels must be grown back to the original structure. **If there are multiple labels in "
			"a connected component, non-labeled pixels are assigned the smallest label in the non-distributed version "
			"and (mostly) random label among all the possibilities in the distributed version.** "
			"Therefore, **this function is suited only for images containing separate blobs or particles**, where each "
			"particle contains seed point(s) of only single value. "
			"\n\n"
			"An alternative to this command is `morphorec`. "
			"It works such that each pixel will get the label of the nearest labeled pixel.",
			{
				CommandArgument<double>(ParameterDirection::In, "allowed color", "Color where other colors will be grown into."),
				CommandArgument<double>(ParameterDirection::In, "background color", "Background color. Values of pixels having this color are not changed. Set to the same value than allowed color to fill to all pixels."),
				CommandArgument<Connectivity>(ParameterDirection::In, "connectivity", string("Connectivity of the regions to grow. ") + connectivityHelp(), Connectivity::NearestNeighbours),
			},
			floodfillSeeAlso())
		{
		}

	public:
		using Distributable::runDistributed;

		virtual void run(Image<pixel_t>& in, std::vector<ParamVariant>& args) const override
		{
			double allowed = pop<double>(args);
			double bg = pop<double>(args);
			Connectivity connectivity = pop<Connectivity>(args);

			size_t changed = growAll(in, pixelRound<pixel_t>(allowed), pixelRound<pixel_t>(bg), connectivity);
			std::cout << std::endl << changed << " pixels changed." << std::endl;
		}

		virtual Vec3c getMargin(const std::vector<ParamVariant>& args) const override
		{
			return Vec3c(3, 3, 3);
		}

		virtual double calculateExtraMemory(const std::vector<ParamVariant>& args) const override
		{
			// Allocate some extra memory for priority queue in filling
			return 2.0;
		}

		virtual std::vector<std::string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const override
		{
			// TODO: This could be made with IterableDistributable
			coord_t changed;
			do
			{
				std::vector<std::string> output = distributor.distribute(this, args);

				changed = parseTotalCount(output, "pixels changed.");
				std::cout << std::endl << changed << " pixels changed." << std::endl;
			} while (changed > 0);

			return std::vector<std::string>();
		}
	};


	//template<typename pixel_t> class FloodFillCommand;

	///**
	//Helper command for distributed flood fill.
	//NOTE: This command supports division only along z-axis!
	//*/
	//template<typename pixel_t> class FloodFillBlockCommand : public OneImageInPlaceCommand<pixel_t>, public Distributable
	//{
	//protected:
	//	friend class CommandList;
	//	friend class FloodFillCommand<pixel_t>;

	//	FloodFillBlockCommand() : OneImageInPlaceCommand<pixel_t>("floodfillblock", "Helper for distributed flood fill command. Performs flood fill starting from seed points defined in a file. Saves seed points outside of current block into target files.",
	//		{
	//			CommandArgument<string>(ParameterDirection::In, "seeds source filename prefix", "Filename prefix for seeds input files."),
	//			CommandArgument<string>(ParameterDirection::In, "seeds target filename prefix", "Filename prefix for seeds output files."),
	//			CommandArgument<Vec3c>(ParameterDirection::In, "start point", "Starting point for the fill. This is the initial start point for the entire distributed fill - not the start point for any block. Seeds for each block are given in seed input files."),
	//			CommandArgument<double>(ParameterDirection::In, "original color", "Original color that we are filling. (the color of the region where the fill is allowed to proceed)"),
	//			CommandArgument<double>(ParameterDirection::In, "fill color", "Fill color."),
	//			CommandArgument<Connectivity>(ParameterDirection::In, "connectivity", string("Connectivity of the region to fill. ") + connectivityHelp(), Connectivity::AllNeighbours),
	//			CommandArgument<Vec3c>(ParameterDirection::In, "original size", "Size of the original image."),
	//			CommandArgument<Distributor::BLOCK_ORIGIN_ARG_TYPE>(ParameterDirection::In, Distributor::BLOCK_ORIGIN_ARG_NAME, "Origin of current calculation block in coordinates of the full image. This argument is used internally in distributed processing. Set to zero in normal usage.", Distributor::BLOCK_ORIGIN_ARG_TYPE(0, 0, 0))
	//		})
	//	{
	//	}

	//private:
	//	static void readPoint(std::ifstream& in, Vec3sc& val, int32_t z)
	//	{
	//		in.read((char*)&val.x, sizeof(int32_t));
	//		in.read((char*)&val.y, sizeof(int32_t));
	//		val.z = z;
	//	}

	//	static void readInputFile(const string& infile, int32_t z, vector<Vec3sc>& output)
	//	{
	//		if (fs::exists(infile))
	//			readListFile(infile, output, [=](std::ifstream& in, Vec3sc& val) { readPoint(in, val, z); });
	//	}

	//	static void writePoint(std::ofstream& out, const Vec3sc& p)
	//	{
	//		out.write((char*)&p.x, sizeof(int32_t));
	//		out.write((char*)&p.y, sizeof(int32_t));
	//	}

	//	static void writeOutputFile(const string& outfile, const std::vector<Vec3sc>& list)
	//	{
	//		writeListFile(outfile, list, writePoint);
	//	}

	//public:
	//	virtual void run(Image<pixel_t>& in, std::vector<ParamVariant>& args) const override
	//	{
	//		string seedsSourcePrefix = pop<string>(args);
	//		string seedsTargetPrefix = pop<string>(args);
	//		Vec3c startPoint = pop<Vec3c>(args);
	//		pixel_t origColor = pixelRound<pixel_t>(pop<double>(args));
	//		pixel_t fillColor = pixelRound<pixel_t>(pop<double>(args));
	//		Connectivity connectivity = pop<Connectivity>(args);
	//		Vec3c origSize = pop<Vec3c>(args);
	//		Distributor::BLOCK_ORIGIN_ARG_TYPE blockOrigin = pop<Distributor::BLOCK_ORIGIN_ARG_TYPE>(args);

	//		int32_t blockStartZ = (int32_t)blockOrigin.z;
	//		int32_t blockEndZ = (int32_t)(blockStartZ + in.depth() - 1);

	//		// Read seed points (three files: block start z, block end z, start point z)
	//		vector<Vec3sc> seeds;
	//		readInputFile(seedsSourcePrefix + "_" + itl2::toString(blockStartZ), blockStartZ, seeds);
	//		readInputFile(seedsSourcePrefix + "_" + itl2::toString(startPoint.z), (int32_t)startPoint.z, seeds);
	//		readInputFile(seedsSourcePrefix + "_" + itl2::toString(blockEndZ), blockEndZ, seeds);

	//		// Convert seeds to block coordinates
	//		for (Vec3sc& seed : seeds)
	//			seed.z -= blockStartZ;

	//		// Iterate block edge pixel values and save to a list.
	//		std::vector<pixel_t> edgeValues;
	//		edgeValues.reserve(in.width() * in.height() * 2);
	//		coord_t zStep = std::max<coord_t>(1, in.depth() - 1);
	//		for (coord_t z = 0; z < in.depth(); z += zStep)
	//		{
	//			for (coord_t y = 0; y < in.height(); y++)
	//			{
	//				for (coord_t x = 0; x < in.width(); x++)
	//				{
	//					edgeValues.push_back(in(x, y, z));
	//				}
	//			}
	//		}

	//		// Flood fill
	//		floodfill(in, seeds, origColor, fillColor, fillColor, connectivity);



	//		// Iterate block edge values and compare to saved, make a list of points that changed.
	//		// Neighbours of the changed points that are not in the current block are new seeds for the neighbouring block.
	//		std::vector<Vec3sc> beginNewSeeds;
	//		std::vector<Vec3sc> endNewSeeds;
	//		size_t n = 0;
	//		for (coord_t z = 0; z < in.depth(); z += zStep)
	//		{
	//			int32_t dz = 0;
	//			std::vector<Vec3sc>* newSeeds;
	//			if (z <= 0)
	//			{
	//				dz = -1;
	//				newSeeds = &beginNewSeeds;
	//			}
	//			else if (z >= in.depth() - 1)
	//			{
	//				dz = 1;
	//				newSeeds = &endNewSeeds;
	//			}
	//			else
	//			{
	//				throw std::logic_error("Flood fill has been distributed along unsupported block shape (there is a face whose normal is not +-z).");
	//			}

	//			int32_t newZ = (int32_t)z + (int32_t)blockStartZ + dz;

	//			for (coord_t y = 0; y < in.height(); y++)
	//			{
	//				for (coord_t x = 0; x < in.width(); x++)
	//				{
	//					if (dz != 0 && newZ >= 0 && newZ < origSize.z)
	//					{
	//						if (edgeValues[n] != in(x, y, z))
	//						{
	//							Vec3sc p((int32_t)x, (int32_t)y, newZ);
	//							if (connectivity == Connectivity::NearestNeighbours)
	//							{
	//								newSeeds->push_back(p + Vec3sc(0, 0, 0));
	//							}
	//							else if (connectivity == Connectivity::AllNeighbours)
	//							{
	//								newSeeds->push_back(p + Vec3sc(-1, -1, 0));
	//								newSeeds->push_back(p + Vec3sc(0, -1, 0));
	//								newSeeds->push_back(p + Vec3sc(1, -1, 0));
	//								newSeeds->push_back(p + Vec3sc(-1, 0, 0));
	//								newSeeds->push_back(p + Vec3sc(0, 0, 0));
	//								newSeeds->push_back(p + Vec3sc(1, 0, 0));
	//								newSeeds->push_back(p + Vec3sc(-1, 1, 0));
	//								newSeeds->push_back(p + Vec3sc(0, 1, 0));
	//								newSeeds->push_back(p + Vec3sc(1, 1, 0));
	//							}
	//							else
	//							{
	//								throw std::logic_error("Unsupported connectivity value in flood fill.");
	//							}
	//						}
	//					}
	//					n++;
	//				}
	//			}
	//		}

	//		if (beginNewSeeds.size() > 0)
	//		{
	//			writeOutputFile(seedsTargetPrefix + "_" + itl2::toString(blockStartZ - 1), beginNewSeeds);
	//		}

	//		if (endNewSeeds.size() > 0)
	//		{
	//			writeOutputFile(seedsTargetPrefix + "_" + itl2::toString(blockEndZ + 1), endNewSeeds);
	//		}

	//	}

	//	bool isInternal() const override
	//	{
	//		return true;
	//	}

	//	virtual bool needsToRunBlock(const std::vector<ParamVariant>& args, const Vec3c& readStart, const Vec3c& readSize, const Vec3c& writeFilePos, const Vec3c& writeImPos, const Vec3c& writeSize, size_t blockIndex) const override
	//	{
	//		string seedsSourcePrefix = std::get<string>(args[1]);
	//		Vec3c startPoint = std::get<Vec3c>(args[3]);

	//		// We need to run a block only if input files corresponding to block start z, block end z, or start point z exist.

	//		int32_t blockStartZ = (int32_t)readStart.z;
	//		int32_t blockEndZ = (int32_t)(blockStartZ + readSize.z - 1);

	//		if (fs::exists(seedsSourcePrefix + "_" + itl2::toString(blockStartZ)) ||
	//			fs::exists(seedsSourcePrefix + "_" + itl2::toString(blockEndZ)))
	//			return true;

	//		if (blockStartZ <= startPoint.z && startPoint.z <= blockEndZ)
	//			return fs::exists(seedsSourcePrefix + "_" + itl2::toString(startPoint.z));

	//		return false;
	//	}

	//	using Distributable::runDistributed;

	//	virtual std::vector<string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const override
	//	{
	//		return distributor.distribute(this, args);
	//	}

	//	virtual double calculateExtraMemory(const std::vector<ParamVariant>& args) const override
	//	{
	//		// The amount of extra memory is impossible to know, so we make a bad estimate.
	//		// TODO: How to improve this?
	//		return 1.0;
	//	}

	//	virtual JobType getJobType(const std::vector<ParamVariant>& args) const override
	//	{
	//		return JobType::Normal;
	//	}

	//	virtual bool canDelay(const std::vector<ParamVariant>& args) const override
	//	{
	//		return false;
	//	}
	//};



	//template<typename pixel_t> class FloodFillCommand : public OneImageInPlaceCommand<pixel_t>, public Distributable
	//{
	//protected:
	//	friend class CommandList;

	//	FloodFillCommand() : OneImageInPlaceCommand<pixel_t>("floodfill", "Performs flood fill. Fills start point and all its neighbours and their neighbours etc. recursively as long as the color of the pixel to be filled equals color of the start point.",
	//		{
	//			CommandArgument<Vec3c>(ParameterDirection::In, "start point", "Starting point for the fill."),
	//			CommandArgument<double>(ParameterDirection::In, "fill value", "Fill color."),
	//			CommandArgument<Connectivity>(ParameterDirection::In, "connectivity", string("Connectivity of the region to fill. ") + connectivityHelp(), Connectivity::AllNeighbours),
	//		},
	//		floodfillSeeAlso())
	//	{
	//	}

	//public:
	//	virtual void run(Image<pixel_t>& in, std::vector<ParamVariant>& args) const override
	//	{
	//		Vec3c startPoint = pop<Vec3c>(args);
	//		pixel_t color = pixelRound<pixel_t>(pop<double>(args));
	//		Connectivity connectivity = pop<Connectivity>(args);

	//		floodfill(in, startPoint, color, color, connectivity);
	//	}

	//	virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override
	//	{
	//		// Algorithm:
	//		// make initial seeds file
	//		//	Instead of collecting everything into single file, output seeds such that seeds with constant z coordinate go
	//		//	to a single file (do not write z at all). Then gathering of seeds is not necessary, floodFillBlockCommand needs
	//		//	to be run only for blocks for which seed file exists, and each job writes to 2 seed files (separate ones).
	//		// distribute (run only blocks that contain seeds)
	//		//	each job:
	//		//		read seeds files in the block
	//		//		run flood fill
	//		//		write seed files if there are seeds
	//		// if there are any output seed files, repeat distribute

	//		DistributedImage<pixel_t>& img = *pop<DistributedImage<pixel_t>*>(args);
	//		Vec3c startPoint = pop<Vec3c>(args);
	//		pixel_t fillColor = pixelRound<pixel_t>(pop<double>(args));
	//		Connectivity connectivity = pop<Connectivity>(args);

	//		if (!img.isInImage(startPoint))
	//			return vector<string>();

	//		pixel_t origColor = img.getPixel(startPoint);

	//		if (origColor == fillColor)
	//			return vector<string>();

	//		string seedsSourceFilenamePrefix = createTempFilename("flood_fill_seeds1");
	//		string seedsTargetFilenamePrefix = createTempFilename("flood_fill_seeds2");

	//		// Create initial seed file
	//		string seedsSourceFilename = seedsSourceFilenamePrefix + itl2::toString(startPoint.z);
	//		vector<Vec3sc> seeds;
	//		seeds.push_back(Vec3sc(startPoint));
	//		FloodFillBlockCommand<pixel_t>::writeOutputFile(seedsSourceFilenamePrefix + "_" + itl2::toString(startPoint.z), seeds);

	//		size_t it = 0;
	//		while (true)
	//		{
	//			it++;
	//			std::cout << "Iteration " << it << std::endl;

	//			// Run distributed fill
	//			CommandList::get<FloodFillBlockCommand<pixel_t>>().runDistributed(distributor, { &img, seedsSourceFilenamePrefix, seedsTargetFilenamePrefix,
	//																			startPoint, (double)origColor, (double)fillColor, connectivity,
	//																			img.dimensions(), Distributor::BLOCK_ORIGIN_ARG_TYPE() });


	//			// Delete sources (seedsSourceFilenamePrefix*)
	//			auto items = itl2::buildFileList(seedsSourceFilenamePrefix + "_*");
	//			for (const string& file : items)
	//				itl2::deleteFile(file);

	//			// Check if there are any seedsTarget files (seedsTargetFilenamePrefix*)
	//			items = itl2::buildFileList(seedsTargetFilenamePrefix + "_*");
	//			if (items.size() <= 0)
	//				break; // No seeds target files means there are no more seeds to process.

	//			std::swap(seedsSourceFilenamePrefix, seedsTargetFilenamePrefix);
	//		}

	//		// At this points all seedsSourceFilenamePrefix* files should have been deleted and
	//		// no new seedsTargetFilenamePrefix* files should exist. ==> no clean-up necessary.

	//		return vector<string>();
	//	}

	//	virtual bool canDelay(const std::vector<ParamVariant>& args) const override
	//	{
	//		return false;
	//	}
	//};













	template<typename pixel_t> class FloodFillCommand;

	/**
	Helper command for distributed flood fill.
	*/
	template<typename pixel_t> class FloodFillBlockCommand : public OneImageInPlaceCommand<pixel_t>, public Distributable
	{
	protected:
		friend class CommandList;
		friend class FloodFillCommand<pixel_t>;

		FloodFillBlockCommand() : OneImageInPlaceCommand<pixel_t>("floodfillblock", "Helper for distributed flood fill command. Performs flood fill starting from seed points defined in a file. Saves seed points outside of current block into target files.",
			{
				CommandArgument<string>(ParameterDirection::In, "seeds source filename prefix", "Filename prefix for seeds input files."),
				CommandArgument<string>(ParameterDirection::In, "seeds target filename prefix", "Filename prefix for seeds output files."),
				CommandArgument<Vec3c>(ParameterDirection::In, "start point", "Starting point for the fill. This is the initial start point for the entire distributed fill - not the start point for any block. Seeds for each block are given in seed input files."),
				CommandArgument<double>(ParameterDirection::In, "original color", "Original color that we are filling. (the color of the region where the fill is allowed to proceed)"),
				CommandArgument<double>(ParameterDirection::In, "fill color", "Fill color."),
				CommandArgument<Connectivity>(ParameterDirection::In, "connectivity", string("Connectivity of the region to fill. ") + connectivityHelp(), Connectivity::AllNeighbours),
				CommandArgument<Vec3c>(ParameterDirection::In, "original size", "Size of the original image."),
				CommandArgument<Distributor::BLOCK_INDEX3_ARG_TYPE>(ParameterDirection::In, Distributor::BLOCK_INDEX3_ARG_NAME, "Index of current calculation block. This argument is used internally in distributed processing. Set to zero in normal usage.", Distributor::BLOCK_INDEX3_ARG_TYPE(0, 0, 0)),
				CommandArgument<Distributor::BLOCK_ORIGIN_ARG_TYPE>(ParameterDirection::In, Distributor::BLOCK_ORIGIN_ARG_NAME, "Origin of current calculation block in coordinates of the full image. This argument is used internally in distributed processing. Set to zero in normal usage.", Distributor::BLOCK_ORIGIN_ARG_TYPE(0, 0, 0))
			})
		{
		}

	private:
		static void readPoint(std::ifstream& in, Vec3sc& val)
		{
			in.read((char*)&val.x, sizeof(int32_t));
			in.read((char*)&val.y, sizeof(int32_t));
			in.read((char*)&val.z, sizeof(int32_t));
		}

		static void readInputFile(const string& infile, vector<Vec3sc>& output)
		{
			if (fs::exists(infile))
				readListFile(infile, output, readPoint);
		}

		static void writePoint(std::ofstream& out, const Vec3sc& p)
		{
			out.write((char*)&p.x, sizeof(int32_t));
			out.write((char*)&p.y, sizeof(int32_t));
			out.write((char*)&p.z, sizeof(int32_t));
		}

		static void writeOutputFile(const string& outfile, const std::vector<Vec3sc>& list)
		{
			writeListFile(outfile, list, writePoint);
		}

		static vector<string> getInputFiles(const string& prefix, const Vec3c& blockIndex)
		{
			return buildFileList(prefix + "_from_*_to_" + itl2::toString(blockIndex) + ".dat");
		}

		/**
		Read all seed input files that are routed to this block.
		*/
		static void readInputFiles(const string& prefix, const Vec3c& blockIndex, vector<Vec3sc>& seeds)
		{
			// Read all files whose name is in format [prefix]_from_[something]_to_[blockIndex].dat
			vector<string> files = getInputFiles(prefix, blockIndex);
			for (const string& file : files)
			{
				readInputFile(file, seeds);
			}
		}

	public:
		virtual void run(Image<pixel_t>& in, std::vector<ParamVariant>& args) const override
		{
			string seedsSourcePrefix = pop<string>(args);
			string seedsTargetPrefix = pop<string>(args);
			Vec3c startPoint = pop<Vec3c>(args);
			pixel_t origColor = pixelRound<pixel_t>(pop<double>(args));
			pixel_t fillColor = pixelRound<pixel_t>(pop<double>(args));
			Connectivity connectivity = pop<Connectivity>(args);
			Vec3c origSize = pop<Vec3c>(args);
			Distributor::BLOCK_INDEX3_ARG_TYPE blockIndex = pop<Distributor::BLOCK_INDEX3_ARG_TYPE >(args);
			Distributor::BLOCK_ORIGIN_ARG_TYPE blockOrigin = pop<Distributor::BLOCK_ORIGIN_ARG_TYPE>(args);

			//int32_t blockStartZ = (int32_t)blockOrigin.z;
			//int32_t blockEndZ = (int32_t)(blockStartZ + in.depth() - 1);
			Vec3c blockEnd = blockOrigin + in.dimensions() - Vec3c(1, 1, 1);

			// Read seed points (three files: block start z, block end z, start point z)
			vector<Vec3sc> seeds;
			readInputFiles(seedsSourcePrefix, blockIndex, seeds);

			// Add possible start point to the seeds list.
			if (startPoint.max() >= 0)
				seeds.push_back(Vec3sc(startPoint));

			// Convert seeds to block coordinates
			for (Vec3sc& seed : seeds)
				seed -= Vec3sc(blockOrigin);


			// Iterate block edge pixel values and save to a list.
			std::vector<pixel_t> edgeValues;
			edgeValues.reserve(in.width() * in.height() * 2 + in.width() * in.depth() * 2 + in.height() * in.depth() * 2);
			forEdges(in.bounds(), in.dimensionality(), [&](coord_t x, coord_t y, coord_t z)
				{
					edgeValues.push_back(in(x, y, z));
				});


			// Flood fill
			floodfill(in, seeds, origColor, fillColor, fillColor, connectivity);


			// Iterate block edge values and compare to saved, make a list of points that changed.
			// Neighbours of the changed points that are not in the current block are new seeds for the neighbouring blocks.
			Image<std::vector<Vec3sc>> newSeeds(3, 3, 3);
			AABoxsc fullBounds = AABoxsc::fromPosSize(Vec3sc(0, 0, 0), Vec3sc(origSize));
			size_t n = 0;
			forEdges(in.bounds(), in.dimensionality(), [&](coord_t x, coord_t y, coord_t z)
				{
					if (edgeValues[n] != in(x, y, z))
					{
						Vec3c p(x, y, z);
						forNeighbours(p, connectivity, [&](coord_t nx, coord_t ny, coord_t nz)
							{
								Vec3sc np((int32_t)nx, (int32_t)ny, (int32_t)nz);
								np += Vec3sc(blockOrigin);
								if (fullBounds.contains(np))
								{
									Vec3sc pos(1, 1, 1);
									if (nx >= in.width())
										pos.x++;
									else if (nx < 0)
										pos.x--;

									if (ny >= in.height())
										pos.y++;
									else if (ny < 0)
										pos.y--;

									if (nz >= in.depth())
										pos.z++;
									else if (nz < 0)
										pos.z--;

									if (pos != Vec3sc(1, 1, 1))
									{
										newSeeds(pos).push_back(np);
									}
								}
							});
					}
					n++;
				});

			// Write all 8 seed lists if they are not empty.
			forAllPixels(newSeeds, [&](coord_t x, coord_t y, coord_t z)
				{
					std::vector<Vec3sc> seedList = newSeeds(x, y, z);
					if (seedList.size() > 0)
					{
						Vec3c targetBlockIndex = blockIndex + Vec3c(x, y, z) - Vec3c(1, 1, 1); // -1's here as the current block corresponds to coordinates (1, 1, 1) in the newSeeds array.
						writeOutputFile(seedsTargetPrefix + "_from_" + itl2::toString(blockIndex) + "_to_" + itl2::toString(targetBlockIndex) + ".dat", seedList);
					}
				});

		}

		bool isInternal() const override
		{
			return true;
		}

		virtual bool needsToRunBlock(const std::vector<ParamVariant>& args, const Vec3c& readStart, const Vec3c& readSize, const Vec3c& writeFilePos, const Vec3c& writeImPos, const Vec3c& writeSize, size_t blockNumber, const Vec3c& blockIndex3) const override
		{
			DistributedImage<pixel_t>& in = *std::get<DistributedImage<pixel_t>*>(args[0]);
			string seedsSourcePrefix = std::get<string>(args[1]);
			Vec3c startPoint = std::get<Vec3c>(args[3]);

			//These values are not available here. Block origin equals readStart argument and block index 3 equals blockIndex3 argument.
			//Distributor::BLOCK_INDEX3_ARG_TYPE blockIndex = std::get<Distributor::BLOCK_INDEX3_ARG_TYPE >(args[8]);
			//Distributor::BLOCK_ORIGIN_ARG_TYPE blockOrigin = std::get<Distributor::BLOCK_ORIGIN_ARG_TYPE >(args[9]);

			// We need to run if we have start point(s) and if they are in this block
			AABoxc blockBounds = in.bounds();
			blockBounds.translate(readStart);
			if (startPoint.max() >= 0 && blockBounds.contains(startPoint))
				return true;

			// We need to run if there are input files available for this block.
			return getInputFiles(seedsSourcePrefix, blockIndex3).size() > 0;
		}

		using Distributable::runDistributed;

		virtual std::vector<string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const override
		{
			return distributor.distribute(this, args);
		}

		virtual double calculateExtraMemory(const std::vector<ParamVariant>& args) const override
		{
			// The amount of extra memory is impossible to know, so we make a bad estimate.
			// TODO: How to improve this?
			return 1.0;
		}

		virtual JobType getJobType(const std::vector<ParamVariant>& args) const override
		{
			return JobType::Normal;
		}

		virtual bool canDelay(const std::vector<ParamVariant>& args) const override
		{
			return false;
		}

		virtual size_t getDistributionDirection1(const std::vector<ParamVariant>& args) const override
		{
			// Enable distribution in z-direction.
			return 2;
		}

		virtual size_t getDistributionDirection2(const std::vector<ParamVariant>& args) const override
		{
			// Enable distribution in y-direction in addition to the default z-directional distribution.
			return 1;
		}

		virtual size_t getDistributionDirection3(const std::vector<ParamVariant>& args) const override
		{
			// Enable distribution in x-direction in addition to the default z-directional distribution.
			return 0;
		}
	};


	template<typename pixel_t> class FloodFillCommand : public OneImageInPlaceCommand<pixel_t>, public Distributable
	{
	protected:
		friend class CommandList;

		FloodFillCommand() : OneImageInPlaceCommand<pixel_t>("floodfill", "Performs flood fill. Fills start point and all its neighbours and their neighbours etc. recursively as long as the color of the pixel to be filled equals color of the start point.",
			{
				CommandArgument<Vec3c>(ParameterDirection::In, "start point", "Starting point for the fill."),
				CommandArgument<double>(ParameterDirection::In, "fill value", "Fill color."),
				CommandArgument<Connectivity>(ParameterDirection::In, "connectivity", string("Connectivity of the region to fill. ") + connectivityHelp(), Connectivity::AllNeighbours),
			},
			floodfillSeeAlso())
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, std::vector<ParamVariant>& args) const override
		{
			Vec3c startPoint = pop<Vec3c>(args);
			pixel_t color = pixelRound<pixel_t>(pop<double>(args));
			Connectivity connectivity = pop<Connectivity>(args);

			floodfill(in, startPoint, color, color, connectivity);
		}

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override
		{
			// Algorithm:
			// make initial seeds file
			//	Instead of collecting everything into single file, output seeds such that seeds with constant z coordinate go
			//	to a single file (do not write z at all). Then gathering of seeds is not necessary, floodFillBlockCommand needs
			//	to be run only for blocks for which seed file exists, and each job writes to 2 seed files (separate ones).
			// distribute (run only blocks that contain seeds)
			//	each job:
			//		read seeds files in the block
			//		run flood fill
			//		write seed files if there are seeds
			// if there are any output seed files, repeat distribute

			DistributedImage<pixel_t>& img = *pop<DistributedImage<pixel_t>*>(args);
			Vec3c startPoint = pop<Vec3c>(args);
			pixel_t fillColor = pixelRound<pixel_t>(pop<double>(args));
			Connectivity connectivity = pop<Connectivity>(args);

			if (!img.isInImage(startPoint))
				return vector<string>();

			pixel_t origColor = img.getPixel(startPoint);

			if (origColor == fillColor)
				return vector<string>();

			string seedsSourceFilenamePrefix = createTempFilename("flood_fill_seeds1");
			string seedsTargetFilenamePrefix = createTempFilename("flood_fill_seeds2");

			// Create initial seed file
			//string seedsSourceFilename = seedsSourceFilenamePrefix + itl2::toString(startPoint.z);
			//vector<Vec3sc> seeds;
			//seeds.push_back(Vec3sc(startPoint));
			//FloodFillBlockCommand<pixel_t>::writeOutputFile(seedsSourceFilenamePrefix + "_from_initiator", seeds);
			
			size_t it = 0;
			while (true)
			{
				it++;
				std::cout << "Iteration " << it << std::endl;

				// Run distributed fill
				CommandList::get<FloodFillBlockCommand<pixel_t>>().runDistributed(distributor, { &img, seedsSourceFilenamePrefix, seedsTargetFilenamePrefix,
																				startPoint, (double)origColor, (double)fillColor, connectivity,
																				img.dimensions(), Distributor::BLOCK_INDEX3_ARG_TYPE(), Distributor::BLOCK_ORIGIN_ARG_TYPE() });

				// Erase intial start point
				startPoint = Vec3c(-1, -1, -1);

				// Delete sources (seedsSourceFilenamePrefix*)
				auto items = itl2::buildFileList(seedsSourceFilenamePrefix + "_*");
				for (const string& file : items)
					itl2::deleteFile(file);

				// Check if there are any seedsTarget files (seedsTargetFilenamePrefix*)
				items = itl2::buildFileList(seedsTargetFilenamePrefix + "_*");
				if (items.size() <= 0)
					break; // No seeds target files means there are no more seeds to process.

				std::swap(seedsSourceFilenamePrefix, seedsTargetFilenamePrefix);
			}

			// At this points all seedsSourceFilenamePrefix* files should have been deleted and
			// no new seedsTargetFilenamePrefix* files should exist. ==> no clean-up necessary.

			return vector<string>();
		}

		virtual bool canDelay(const std::vector<ParamVariant>& args) const override
		{
			return false;
		}
	};
}