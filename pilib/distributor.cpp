
#include "distributor.h"
#include "command.h"
#include "distributedimage.h"
#include "utilities.h"
#include "pisystem.h"
#include "exeutils.h"
#include "math/vectoroperations.h"
#include "whereamicpp.h"

#include <tuple>
#include "filesystem.h"

using namespace std;


namespace pilib
{
	fs::path findPi2()
	{
#if defined(__linux__)
		string piName = "pi2";
#elif defined(_WIN32)
		string piName = "pi2.exe";
#else
	#error Please configure executable file name style for this platform.
#endif

		fs::path p1 = fs::current_path() / piName;
		fs::path p2 = getModulePath();
		fs::path p3 = fs::path(getExecutablePath());
		if(p2.has_filename())
			p2 = p2.replace_filename(piName);

		if (fs::exists(p1) && fs::is_regular_file(p1))
			return p1;
		if (fs::exists(p2) && fs::is_regular_file(p2))
			return p2;
		if (fs::exists(p3) && fs::is_regular_file(p3))
			return p3;

		string message = "Unable to find " + piName + " program. It was searched from " + p1.string();
		if (p2 != p1)
			message += " and " + p2.string();
		if (p3 != p2 && p3 != p1)
			message += " and " + p3.string();
		message += ".";

		throw ITLException(message);
	}

	Distributor::Distributor(PISystem* piSystem) : piSystem(piSystem)
	{
		piCommand = findPi2().string();
	}

	void Distributor::readSettings(INIReader& reader)
	{
		showSubmittedScripts = reader.get<bool>("show_submitted_scripts", false);
		allowDelaying = reader.get<bool>("allow_delaying", true);
	}


	void adjustBlockDimensions(coord_t start, coord_t size, coord_t margin, coord_t fullSize, coord_t& in_left, coord_t& in_width, coord_t& out_left, coord_t& out_width)
	{
		in_left = start - margin;
		out_left = margin;

		coord_t in_right = in_left + size + 2 * margin;
		coord_t out_right = out_left + size;

		//if (in_left < 0 && in_right > fullSize)
		//{
		//	throw ITLException("Unimplemented situation in distribution.");
		//}

		if (in_left < 0)
		{
			out_left += in_left;
			out_right += in_left;
			in_left = 0;
		}

		if (in_right > fullSize)
		{
			in_right = fullSize;

			if (start + (out_right - out_left) > fullSize)
				out_right = fullSize - start + out_left;

			//coord_t delta = in_right - fullSize;
			//in_right -= delta;
			//out_right -= delta;

			//in_right = fullSize;
			//out_right = fullSize - in_left;
		}

		in_width = in_right - in_left;
		out_width = out_right - out_left;
	}

	/**
	Divides input images into blocks for distributed processing.
	@return Vector of block lists. Each list contains corresponding block of each argument. Each item has five fields:
	position where image read operation starts, size of block to read, position in file where output data should be written,
	position in the image where output data region starts, and size of block to write.
	*/
	vector<vector<tuple<Vec3c, Vec3c, Vec3c, Vec3c, Vec3c> > > determineBlocks(size_t refIndex, const Vec3c& margin, const Vec3c& blockCounts, const Command* command, vector<ParamVariant>& args, const vector<ParamVariant>& argsForGetCorrespondingBlock)
	{
		const DistributedImageBase* img = getDistributedImage(args[refIndex]);
		Vec3c refDims = img->dimensions();

		// Block size in the reference image
		Vec3c refSize = Vec3c(itl2::ceil((double)refDims.x / (double)blockCounts.x), itl2::ceil((double)refDims.y / (double)blockCounts.y), itl2::ceil((double)refDims.z / (double)blockCounts.z));


		// allBlocks[parameter index][block index] = tuple<read start, read size, output position in file, output start position in image, output region size>
		vector<vector<tuple<Vec3c, Vec3c, Vec3c, Vec3c, Vec3c> > > allBlocks(args.size(), vector<tuple<Vec3c, Vec3c, Vec3c, Vec3c, Vec3c> >());

		for (coord_t iz = 0; iz < blockCounts.z; iz++)
		{
			for (coord_t iy = 0; iy < blockCounts.y; iy++)
			{
				for (coord_t ix = 0; ix < blockCounts.x; ix++)
				{
					// Select a block from reference image
					Vec3c refStart(ix * refSize.x, iy * refSize.y, iz * refSize.z);

					// ceiling in refSize calculation may cause bad blocks, skip those
					if (refStart.x >= refDims.x || refStart.y >= refDims.y || refStart.z >= refDims.z)
						break;

					// Add margins to the reference block
					Vec3c refStartMargins, refSizeMargins, refWriteImPos, refWriteImSize;
					adjustBlockDimensions(refStart.x, refSize.x, margin.x, refDims.x, refStartMargins.x, refSizeMargins.x, refWriteImPos.x, refWriteImSize.x);
					adjustBlockDimensions(refStart.y, refSize.y, margin.y, refDims.y, refStartMargins.y, refSizeMargins.y, refWriteImPos.y, refWriteImSize.y);
					adjustBlockDimensions(refStart.z, refSize.z, margin.z, refDims.z, refStartMargins.z, refSizeMargins.z, refWriteImPos.z, refWriteImSize.z);

					if ((refWriteImSize.x > refSize.x || refWriteImSize.y > refSize.y || refWriteImSize.z > refSize.z) ||
						(refStart.x + refWriteImSize.x > refDims.x || refStart.y + refWriteImSize.y > refDims.y || refStart.z + refWriteImSize.z > refDims.z))
						throw ITLException("Bug check 1 failed in distribution block size calculation.");

					if (refWriteImSize.x < 0 || refWriteImSize.y < 0 || refWriteImSize.z < 0)
						throw ITLException("Bug check 2 failed in distribution block size calculation.");
					
					//vector<tuple<Vec3c, Vec3c, Vec3c, Vec3c, Vec3c> > blocks;
					for (size_t n = 0; n < args.size(); n++)
					{
						// Calculate corresponding block in image n
						// for input image: file position where read operation should start
						Vec3c readStart = refStartMargins;
						// for input image: size of block to read
						Vec3c readSize = refSizeMargins;

						// for output image: file position where the data should be written
						Vec3c writeFilePos = refStart;
						// for output image: position of good data region in the output image
						Vec3c writeImPos = refWriteImPos;
						// for output image: size of block to write
						Vec3c writeSize = refWriteImSize;

						if (n != refIndex)
						{
							// Cast to reference so that exception is thrown on error.
							dynamic_cast<const Distributable&>(*command).getCorrespondingBlock(argsForGetCorrespondingBlock, n, readStart, readSize, writeFilePos, writeImPos, writeSize);
						}

						allBlocks[n].push_back(make_tuple(readStart, readSize, writeFilePos, writeImPos, writeSize));
					}
				}
			}
		}

		return allBlocks;
	}

	/**
	Calculates memory requirement given image blocks and command.
	allBlocks[parameter index][block index] = tuple<read start, read size, output position in file, output start position in image, output region size>
	*/
	size_t calcRequiredMemory(const vector<vector<tuple<Vec3c, Vec3c, Vec3c, Vec3c, Vec3c> > >& allBlocks, const Command* command)
	{
		// Calculate pixel size of each arg (zero for non-image args)
		vector<size_t> pixelSizes;
		for (size_t n = 0; n < command->args().size(); n++)
		{
			const CommandArgumentBase& def = command->args()[n];
			pixelSizes.push_back(pixelSize(def.dataType()));
		}

		// Calculate data size for all parameters in each block
		size_t blockCount = allBlocks[0].size();
		vector<size_t> blockMems(blockCount, 0);
		for (size_t parami = 0; parami < command->args().size(); parami++)
		{
			auto& blocks = allBlocks[parami];

			size_t req = 0;
			for (size_t n = 0; n < blocks.size(); n++)
			{
				const auto& tup = blocks[n];
				Vec3c readSize = get<1>(tup);
				blockMems[n] += readSize.x * readSize.y * readSize.z * pixelSizes[parami];
			}
		}

		// Return maximum data size
		return max(blockMems);
	}

	/**
	Calculates size of largest block in the given block list in pixels.
	*/
	size_t maxBlockSize(const vector<tuple<Vec3c, Vec3c, Vec3c, Vec3c, Vec3c> >& blocks)
	{
		size_t maxSize = 0;
		for (size_t n = 0; n < blocks.size(); n++)
		{
			const auto& tup = blocks[n];
			Vec3c readSize = get<1>(tup);
			size_t size = readSize.x * readSize.y * readSize.z;
			maxSize = std::max(maxSize, size);
		}
		return maxSize;
	}

	/**
	Gets all input images in the given list of delayed commands.
	*/
	set<DistributedImageBase*> getInputImages(const vector<Delayed>& list)
	{
		set<DistributedImageBase*> result;
		for (size_t n = 0; n < list.size(); n++)
		{
			auto items = list[n].getInputImages();
			result.insert(items.begin(), items.end());
		}

		return result;
	}

	/**
	Gets all output images in the given list of delayed commands.
	*/
	set<DistributedImageBase*> getOutputImages(const vector<Delayed>& list)
	{
		set<DistributedImageBase*> result;
		for (size_t n = 0; n < list.size(); n++)
		{
			auto items = list[n].getOutputImages();
			result.insert(items.begin(), items.end());
		}

		return result;
	}

	/**
	Creates set of all images in the arguments of all commands in the given list.
	*/
	set<DistributedImageBase*> getAllImages(const vector<Delayed>& list)
	{
		set<DistributedImageBase*> result = getInputImages(list);
		set<DistributedImageBase*> outputs = getOutputImages(list);
		result.insert(outputs.begin(), outputs.end());
		return result;
	}

	/**
	Tests if set 'all' contains all elements of set 'tofind'.
	*/
	bool containsAll(const set<DistributedImageBase*>& all, const vector<DistributedImageBase*>& tofind)
	{
		for (DistributedImageBase* item : tofind)
		{
			if (all.find(item) == all.end())
				return false;
		}

		return true;
	}

	bool Distributor::canDelay(const Delayed& delayed) const
	{
		if (!delayed.canDelay())
			return false;

		if (delayedCommands.empty())
			return true;

		if(!containsAll(getAllImages(delayedCommands), delayed.getInputImages()))
			return false;

		return true;
	}

	bool Distributor::tryDelay(Delayed& d)
	{
		if (!canDelay(d))
			return false;

		delayedCommands.push_back(d);

		Vec3c margin;
		set<DistributedImageBase*> inputImages;
		set<DistributedImageBase*> outputImages;
		JobType jobType;
		map<DistributedImageBase*, vector<tuple<Vec3c, Vec3c, Vec3c, Vec3c, Vec3c> > > blocksPerImage;
		size_t memoryReq;
		try
		{
			determineDistributionConfiguration(margin, inputImages, outputImages, jobType, blocksPerImage, memoryReq);
			return true;
		}
		catch (ITLException)
		{
			delayedCommands.pop_back();
			return false;
		}
	}

	vector<string> Distributor::distribute(const Command* command, vector<ParamVariant>& args, vector<ParamVariant>* argsForGetCorrespondingBlock)
	{
		if (args.size() != command->args().size())
			throw logic_error("Invalid number of arguments.");

		Delayed d(piSystem, command, args, argsForGetCorrespondingBlock);

		if (!allowDelaying)
		{
			// Delaying is not allowed so process the command right away.
			delayedCommands.push_back(d);
			runDelayedCommands();
			return lastOutput;
		}

		if (tryDelay(d))
		{
			// The command was succesfully delayed.
			return vector<string>();
		}
		else
		{
			// The command could not be delayed, so run commands in the queue and
			// try to delay again.

			runDelayedCommands();

			if (tryDelay(d))
			{
				// The command was succesfully delayed.
				return vector<string>();
			}
			else
			{
				// This command can't be delayed even with empty queue, so run it right away.
				delayedCommands.push_back(d);
				runDelayedCommands();
				return lastOutput;
			}
		}
		
	}

	/**
	Calculates maximum margin in all commands in the given list.
	*/
	Vec3c getMaxMargin(const vector<Delayed>& delayedCommands)
	{
		Vec3c margin(0, 0, 0);
		for (const Delayed& d : delayedCommands)
		{
			margin = max(margin, d.getMargin());
		}
		return margin;
	}

	/**
	Calculates suitable job type for the delayed commands in the list.
	*/
	JobType getCombinedJobType(const vector<Delayed>& delayedCommands)
	{
		if (delayedCommands.size() <= 0)
			return JobType::Normal;

		if (delayedCommands.size() == 1)
			return delayedCommands[0].getJobType();

		JobType slowest = JobType::Normal; // Begin from Normal type as we are combining multiple jobs.
		for (const Delayed& d : delayedCommands)
		{
			JobType curr = d.getJobType();
			if (curr == JobType::Slow)
				slowest = JobType::Slow;
			else if (curr == JobType::Normal && slowest == JobType::Fast)
				slowest = JobType::Normal;
		}

		return slowest;
	}

	/**
	Finds preferred amount of subdivisions.
	*/
	size_t getPreferredSubdivisions(const vector<Delayed>& delayedCommands)
	{
		if (delayedCommands.size() <= 0)
			throw logic_error("No commands.");

		vector<size_t> divs;
		for (size_t n = 0; n < delayedCommands.size(); n++)
		{
			divs.push_back(delayedCommands[n].getPreferredSubdivisions());
		}

		return max(divs);
	}

	/**
	Finds first distribution direction.
	*/
	size_t getDistributionDirection1(const vector<Delayed>& delayedCommands)
	{
		if (delayedCommands.size() <= 0)
			throw logic_error("No commands.");

		size_t result = delayedCommands[0].getDistributionDirection1();

		for (size_t n = 1; n < delayedCommands.size(); n++)
		{
			if (delayedCommands[n].getDistributionDirection1() != result)
				throw ITLException("Incompatible distribution directions.");
		}

		return result;
	}

	/**
	Finds second distribution direction.
	*/
	size_t getDistributionDirection2(const vector<Delayed>& delayedCommands)
	{
		if (delayedCommands.size() <= 0)
			throw logic_error("No commands.");

		size_t result = delayedCommands[0].getDistributionDirection2();

		if (result >= numeric_limits<size_t>::max())
			return numeric_limits<size_t>::max();

		for (size_t n = 1; n < delayedCommands.size(); n++)
		{
			size_t dir = delayedCommands[n].getDistributionDirection2();
			
			// Any command that does not allow distribution in second direction disables multi-direction distribution.
			if (dir >= numeric_limits<size_t>::max())
				return numeric_limits<size_t>::max();

			// If commands have different allowed second directions, we disable second direction altogether.
			if(dir != result)
				return numeric_limits<size_t>::max();
		}

		return result;
	}

	void Distributor::flush()
	{
		runDelayedCommands();
	}


	void Distributor::determineDistributionConfiguration(Vec3c& margin, set<DistributedImageBase*>& inputImages, set<DistributedImageBase*>& outputImages, JobType& jobType, map<DistributedImageBase*, vector<tuple<Vec3c, Vec3c, Vec3c, Vec3c, Vec3c> > >& blocksPerImage, size_t& memoryReq)
	{
		if (delayedCommands.size() <= 0)
			return;

		margin = getMaxMargin(delayedCommands);

		// Find all input and output images
		inputImages = getInputImages(delayedCommands);
		outputImages = getOutputImages(delayedCommands);

		// Get job type
		jobType = getCombinedJobType(delayedCommands);


		// Find distribution directions
		size_t distributionDirection1 = getDistributionDirection1(delayedCommands);
		size_t distributionDirection2 = getDistributionDirection2(delayedCommands);

		if (distributionDirection1 > 2)
			throw logic_error("Invalid distribution direction.");

		if (distributionDirection2 == distributionDirection1)
			throw logic_error("Both distribution directions can't be the same.");

		// If we need to distribute in two directions, check that all output images are (written to) .raw files.
		// Otherwise distribution is not allowed as sequences can't be written to in parallel.
		bool distributionDirection2Allowed = true;
		if (distributionDirection2 <= 2)
		{
			for (DistributedImageBase* img : outputImages)
			{
				if (!img->isOutputRaw())
				{
					distributionDirection2Allowed = false;
					break;
				}
			}
		}

        size_t preferredSubdivisions = getPreferredSubdivisions(delayedCommands);

		// Determine blocks that must be loaded, given amount of subdivisions in each direction.
		// blocksPerImage[image pointer][block index] = tuple<block definition>
		Vec3c subDivisions(1, 1, 1);
		subDivisions[distributionDirection1] = preferredSubdivisions;
		size_t lastMemoryReq = numeric_limits<size_t>::max();
		while (true)
		{
			//string debugText;
			//debugText += "delayedCommands.size() = " + itl2::toString(delayedCommands.size()) + "\n";

			// Calculate block sizes for each command.
			blocksPerImage.clear();
			vector<size_t> extraMemPerCommand;
			for (Delayed& d : delayedCommands)
			{
				// blocksPerParameter[parameter index][block index] = tuple<block definition>
				vector<vector<tuple<Vec3c, Vec3c, Vec3c, Vec3c, Vec3c> > > blocksPerParameter = determineBlocks(d.getRefIndex(), margin, subDivisions, d.getCommand(), d.getArgs(), d.getArgsForGetCorrespondingBlock());

				// Store block size for each image.
				// If some image gets multiple different block sizes, throw an error. The commands can't be executed together like that.
				const vector<CommandArgumentBase>& argDefs = d.getCommand()->args();
				for (size_t n = 0; n < argDefs.size(); n++)
				{
					if (isImage(argDefs[n].dataType()))
					{
						DistributedImageBase* img = getDistributedImage(d.getArgs()[n]);
						if (blocksPerImage.find(img) != blocksPerImage.end())
						{
							// There is already a blocks list for this image.
							// Make sure it is the same than what was generated now.
							if (blocksPerImage[img] != blocksPerParameter[n])
								throw ITLException("The delayed commands require different block sizes for the same images.");
						}
						else
						{
							blocksPerImage[img] = blocksPerParameter[n];
						}
					}
				}

				// Calculate total memory requirement
				//memoryReq += (size_t)ceil(calcRequiredMemory(blocksPerParameter, d.getCommand()) * (1 + d.calculateExtraMemory()));

				// Calculate extra memory requirement for current command

				//debugText += "blocksPerParameter = \n";
				//for(const vector<tuple<Vec3c, Vec3c, Vec3c, Vec3c, Vec3c> >& v : blocksPerParameter)
				//{
				//	for (const auto& tup : v)
				//	{
				//		debugText += itl2::toString(get<0>(tup)) + ", " + itl2::toString(get<1>(tup)) + ", " + itl2::toString(get<2>(tup)) + ", " + itl2::toString(get<3>(tup)) + ", " + itl2::toString(get<4>(tup)) + "\n";
				//	}
				//}
				//debugText += "calcRequiredMemory(blocksPerParameter, d.getCommand()) = " + itl2::toString(calcRequiredMemory(blocksPerParameter, d.getCommand())) + "\n";
				//debugText += "d.calculateExtraMemory() = " + itl2::toString(d.calculateExtraMemory()) + "\n";

				extraMemPerCommand.push_back(pixelRound<size_t>(std::ceil(calcRequiredMemory(blocksPerParameter, d.getCommand()) * d.calculateExtraMemory())));
			}


			// Extra memory required by commands
			memoryReq = max(extraMemPerCommand);

			// Add memory required by image blocks
			for (const auto& item : blocksPerImage)
			{
				memoryReq += maxBlockSize(item.second) * item.first->pixelSize();
			}


			if (memoryReq <= allowedMemory())
				break;

			// Check that distribution is sane
			//bool failed = false;

			// Condition 1: memory requirement must go down when we increase subdivisions.
			if (memoryReq > lastMemoryReq)
			{
				//cout << debugText << endl;

				throw ITLException(string("Unable to find suitable subdivision. Memory requirement does not decrease from ") + bytesToString((double)lastMemoryReq) + " when decreasing processing block size, and only " + bytesToString((double)allowedMemory()) + " is available for a single process. See also cluster configuration file max_memory setting.");
				//failed = true;
			}

			// Condition 2: Images must not be subdivided more than their dimensions allow.
			for (DistributedImageBase* img : inputImages)
			{
				// TODO: This condition is a bit "hacky" and will probably fail some time.
				if ((img->dimensions()[distributionDirection1] > 1 && subDivisions[distributionDirection1] >= img->dimensions()[distributionDirection1]) ||
					(distributionDirection2 <= 2 && img->dimensions()[distributionDirection2] > 1 && subDivisions[distributionDirection2] >= img->dimensions()[distributionDirection2])
					)
				{
					//cout << debugText << endl;

					//cout << "subdivisions in 1. dir = " << subDivisions[distributionDirection1] << endl;
					//cout << "subdivisions in 2. dir = " << subDivisions[distributionDirection2] << endl;

					throw ITLException(string("Unable to find suitable subdivision. The smallest possible blocks of the input and output images require ") + bytesToString((double)memoryReq) + " of memory, but only " + bytesToString((double)allowedMemory()) + " is available for a single process. See also cluster configuration file max_memory setting.");
					//failed = true;
					//break;
				}
			}

			//if(failed)
			//	throw ITLException(string("Unable to find suitable subdivision. The smallest possible blocks of the input and output images require ") + bytesToString((double)memoryReq) + " of memory, but only " + bytesToString((double)allowedMemory()) + " is available for a single process. See also cluster configuration file max_memory setting.");

			lastMemoryReq = memoryReq;

			// TODO: Try replacing the subDivisions[distributionDirection1] < N condition by a condition that compares total pixel count in margins to total size of image.
			if (distributionDirection2 > 2 || subDivisions[distributionDirection1] < 80)
			{
				// Distribution direction 2 is not available or there are only a few subdivisions in direction 1
				// Favor direction 1 as it is often z and distributing in that direction is faster and more compatible than other directions.
				subDivisions[distributionDirection1]++;
			}
			else
			{
				if (!distributionDirection2Allowed)
					throw ITLException("The input images are so large that they must be distributed in two coordinate directions. Distribution in two directions is currently not supported for non-.raw image files. Consider saving all the input images as .raw before calling this command.");

				// Subdivide in 2 directions
				if (subDivisions[distributionDirection2] < subDivisions[distributionDirection1])
					subDivisions[distributionDirection2]++;
				else
					subDivisions[distributionDirection1]++;
			}
		}
	}

	/**
	Convert argument to string.
	*/
	string argumentToString(const CommandArgumentBase& argument, const ParamVariant& value)
	{
		ArgumentDataType dt = argument.dataType();

		string s;
		switch (dt)
		{
		case ArgumentDataType::String: s = get<string>(value); return s;
		case ArgumentDataType::Double: return itl2::toString(get<double>(value));
		case ArgumentDataType::Int: return itl2::toString(get<coord_t>(value));
		case ArgumentDataType::Size: return itl2::toString(get<size_t>(value));
		case ArgumentDataType::NBType: return itl2::toString(get<NeighbourhoodType>(value));
		case ArgumentDataType::BoundaryCond: return itl2::toString(get<BoundaryCondition>(value));
		case ArgumentDataType::Connectiv: return itl2::toString(get<Connectivity>(value));
		case ArgumentDataType::InterpolationMode: return itl2::toString(get<InterpolationMode>(value));
		case ArgumentDataType::Bool: return itl2::toString(get<bool>(value));
		case ArgumentDataType::Vect3d: return itl2::toString(get<Vec3d>(value));
		case ArgumentDataType::Vect3c: return itl2::toString(get<Vec3c>(value));
		case ArgumentDataType::ImageUInt8:
		case ArgumentDataType::ImageUInt16:
		case ArgumentDataType::ImageUInt32:
		case ArgumentDataType::ImageUInt64:
		case ArgumentDataType::ImageInt8:
		case ArgumentDataType::ImageInt16:
		case ArgumentDataType::ImageInt32:
		case ArgumentDataType::ImageInt64:
		case ArgumentDataType::ImageFloat32:
		case ArgumentDataType::ImageComplex32:
			return getDistributedImage(value)->uniqueName();
		default: throw ITLException("Data type not configured.");
		}
	}

	void Distributor::runDelayedCommands()
	{
		if (delayedCommands.size() <= 0)
			return;

		// Find suitable block size and other parameters.
		Vec3c margin;
		set<DistributedImageBase*> inputImages;
		set<DistributedImageBase*> outputImages;
		JobType jobType;
		map<DistributedImageBase*, vector<tuple<Vec3c, Vec3c, Vec3c, Vec3c, Vec3c> > > blocksPerImage;
		size_t memoryReq;
		determineDistributionConfiguration(margin, inputImages, outputImages, jobType, blocksPerImage, memoryReq);


		// If overlap is nonzero, InOut images must be saved to different file from which they are loaded.
		// Changing write targets should not cause bad state of image objects even if writeComplete() is not called.
		if (margin != Vec3c(0, 0, 0))
		{
			for (DistributedImageBase* img : inputImages)
			{
				if (outputImages.find(img) != outputImages.end())
				{
					// This is InOut image
					if (img->currentReadSource() == img->currentWriteTarget())
					{
						img->newWriteTarget();
					}
				}
			}
		}


		size_t jobCount = blocksPerImage.begin()->second.size();
		cout << "Submitting " << jobCount << " jobs, each estimated to require at most " << bytesToString((double)memoryReq) << " of RAM..." << endl;
		for (size_t i = 0; i < jobCount; i++)
		{
			// Build job script:
			// readblock(Block of input image 1)
			// readblock(Block of input image 2)
			// ...(for all input and input/output images)
			//
			// command1
			// command2
			// command3
			// ...(for all commands)
			//
			// writeblock(Block of output image 1)
			// writeblock(Block of output image 2)
			// ...(for all output images)
			//
			// print("Everything done")

			stringstream script;

			// Init so that we always print something (required at least in SLURM distributor)
			script << "echo(true, false);" << endl;

			// Image read commands
			for(DistributedImageBase* img : inputImages)
			{
				Vec3c readStart = get<0>(blocksPerImage[img][i]);
				Vec3c readSize = get<1>(blocksPerImage[img][i]);
				script << img->emitReadBlock(readStart, readSize, true);
			}
			// Output image creation commands
			for (DistributedImageBase* img : outputImages)
			{
				if (inputImages.find(img) == inputImages.end())
				{
					Vec3c readStart = get<0>(blocksPerImage[img][i]);
					Vec3c readSize = get<1>(blocksPerImage[img][i]);
					script << img->emitReadBlock(readStart, readSize, false);
				}
			}

			// Processing commands
			for (size_t cmdi = 0; cmdi < delayedCommands.size(); cmdi++)
			{
				const Command* command = delayedCommands[cmdi].getCommand();
				vector<ParamVariant>& args = delayedCommands[cmdi].getArgs();
				size_t refIndex = delayedCommands[cmdi].getRefIndex();

				script << command->name() << "(";
				for (size_t n = 0; n < args.size(); n++)
				{
					// Value of argument whose type is Vec3c and name is "block origin" is replaced by the origin of current calculation block.
					// This functionality is needed at least in skeleton tracing command.
					const CommandArgumentBase& argDef = command->args()[n];
					ParamVariant argVal = args[n];
					if (argDef.dataType() == parameterType<BLOCK_ORIGIN_ARG_TYPE>() && argDef.name() == BLOCK_ORIGIN_ARG_NAME)
					{
						DistributedImageBase* refImage = getDistributedImage(args[refIndex]);
						Vec3c readStart = get<0>(blocksPerImage[refImage][i]);

						argVal = readStart;
					}
					else if (argDef.dataType() == parameterType<BLOCK_INDEX_ARG_TYPE>() && argDef.name() == BLOCK_INDEX_ARG_NAME)
					{
						argVal = (coord_t)i;
					}

					script << "\"" << argumentToString(argDef, argVal) << "\"";
					if (n < args.size() - 1)
						script << ", ";
				}
				script << ");" << endl;
			}

			// Image write commands
			for (DistributedImageBase* img : outputImages)
			{
				// Only write if the image is still visible from the main PI system object
				if (piSystem->isDistributedImage(img))
				{
					Vec3c writeFilePos = get<2>(blocksPerImage[img][i]);
					Vec3c writeImPos = get<3>(blocksPerImage[img][i]);
					Vec3c writeSize = get<4>(blocksPerImage[img][i]);

					// Only write if writing is requested by the command.
					if(writeSize.min() > 0)
						script << img->emitWriteBlock(writeFilePos, writeImPos, writeSize);
				}
			}

			if (showSubmittedScripts)
			{
				cout << "Submitting pi2 script:" << endl;
				cout << script.str() << endl;
			}

			submitJob(script.str(), jobType);
		}

		// Run jobs first and set writeComplete() only after the jobs have finished to make sure that
		// the output image exists before writeComplete() is called.
		// Additionally, this order ensures that if jobs fail, the input images still point to the correct files.
		try
		{
			cout << "Waiting for jobs to finish..." << endl;
			lastOutput = waitForJobs();

			for (DistributedImageBase* img : outputImages)
			{
				img->writeComplete();
			}
		}
		catch (...)
		{
			delayedCommands.clear();
			throw;
		}

		// This may deallocate images that are not in PISystem anymore.
		delayedCommands.clear();
	}

}
