
#include "distributor.h"
#include "command.h"
#include "distributedimage.h"
#include "utilities.h"
#include "pisystem.h"
#include "exeutils.h"
#include "math/vectoroperations.h"
#include "whereamicpp.h"
#if defined(__linux__) || defined(__APPLE__)
#include <sys/types.h>
#include <dirent.h>
#endif
#include <tuple>
#include "filesystem.h"
#include "timing.h"

using namespace std;


namespace pilib
{

	void flushCache(const string& filename)
	{
#if defined(__linux__) || defined(__APPLE__)
		// This may flush NFS cache on the files of the folder where the log file lives.
		int fd = open(filename.c_str(), O_RDONLY);
		fsync(fd);
		close(fd);

		fs::path dir(filename);
		dir = dir.parent_path();
		DIR* dr = opendir(dir.string().c_str());
		closedir(dr);
#endif
	}

	void showProgressBar(const string& bar, size_t& barLength)
	{
		// Remove old bar
		for (size_t n = 0; n < barLength; n++)
			cout << ' ';
		if (barLength > 0)
			cout << '\r';

		// Store length of new bar
		barLength = bar.length();

		// Show new bar
		cout << bar;
		if (barLength > 0)
			cout << '\r';

		cout << flush;
	}

	char getProgressChar(int progress)
	{
		if (progress == JOB_FAILED)
			return '!';
		if (progress == JOB_WAITING)
			return 'W';


		// 20	40	60	80	100
		// _	.	o	O	*

		if (progress < 20)
			return '_';
		if (progress < 40)
			return '.';
		if (progress < 60)
			return 'o';
		if (progress < 80)
			return 'O';

		return '*';
	}

	string createProgressBar(const vector<int>& progress)
	{
		int readyCount = 0;
		int waitingCount = 0;

		for (size_t n = 0; n < progress.size(); n++)
		{
			if (progress[n] == JOB_WAITING)
				waitingCount++;
			else if (progress[n] >= 100)
				readyCount++;
		}

		int runningCount = (int)progress.size() - waitingCount - readyCount;

		stringstream s;
		s << "waiting: " << waitingCount << ", ready: " << readyCount << ", running: " << runningCount << "; ";

		for (size_t n = 0; n < progress.size(); n++)
		{
			if (progress[n] != JOB_WAITING && progress[n] < 100)
				s << getProgressChar(progress[n]);
		}

		return s.str();
	}



	fs::path findPi2()
	{
#if defined(__linux__) || defined(__APPLE__)
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

	Distributor::Distributor(PISystem* piSystem) :
		piSystem(piSystem),
		nn5ChunkSize(nn5::DEFAULT_CHUNK_SIZE),
		showSubmittedScripts(false),
		allowDelaying(false)
	{
		configDir = findPi2().remove_filename();
		piCommand = "\"" + findPi2().string() + "\"";
	}


	void Distributor::readSettings(INIReader& reader)
	{
		piCommand = reader.get<string>("pi2_command", piCommand);
		showSubmittedScripts = reader.get<bool>("show_submitted_scripts", false);
		allowDelaying = reader.get<bool>("allow_delaying", true);
		chunkSize(reader.get<Vec3c>("chunk_size", nn5::DEFAULT_CHUNK_SIZE));
		maxSubmittedJobCount = reader.get<size_t>("max_parallel_submit_count", 0);
		promoteThreshold = reader.get<size_t>("promote_threshold", 3);
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


	vector<vector<DistributionChunk> > determineBlocksFromBlockSize(size_t refIndex, const Vec3c& margin, const Vec3c& blockSize, const Command* command, vector<ParamVariant>& args, const vector<ParamVariant>& argsForGetCorrespondingBlock)
	{
		// allBlocks[parameter index][block index] = tuple<read start, read size, output position in file, output start position in image, output region size>
		vector<vector<DistributionChunk> > allBlocks(args.size(), vector<DistributionChunk>());

		const DistributedImageBase* img = getDistributedImage(args[refIndex]);
		Vec3c refDims = img->dimensions();

		nn5::internals::forAllChunks(refDims, blockSize, false, [&](const Vec3c& chunkIndex, const Vec3c& refStart)
			{

				// Add margins to the reference block
				Vec3c refStartMargins, refSizeMargins, refWriteImPos, refWriteImSize;
				adjustBlockDimensions(refStart.x, blockSize.x, margin.x, refDims.x, refStartMargins.x, refSizeMargins.x, refWriteImPos.x, refWriteImSize.x);
				adjustBlockDimensions(refStart.y, blockSize.y, margin.y, refDims.y, refStartMargins.y, refSizeMargins.y, refWriteImPos.y, refWriteImSize.y);
				adjustBlockDimensions(refStart.z, blockSize.z, margin.z, refDims.z, refStartMargins.z, refSizeMargins.z, refWriteImPos.z, refWriteImSize.z);

				if ((refWriteImSize.x > blockSize.x || refWriteImSize.y > blockSize.y || refWriteImSize.z > blockSize.z) ||
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

					allBlocks[n].push_back(DistributionChunk(readStart, readSize, writeFilePos, writeImPos, writeSize, chunkIndex));
				}
			});

		return allBlocks;
	}

	Vec3c blockSizeFromBlockCounts(size_t refIndex, const Vec3c& blockCounts, const vector<ParamVariant>& args)
	{
		const DistributedImageBase* img = getDistributedImage(args[refIndex]);
		Vec3c refDims = img->dimensions();

		// Block size in the reference image
		Vec3c refSize = Vec3c(itl2::ceil((double)refDims.x / (double)blockCounts.x), itl2::ceil((double)refDims.y / (double)blockCounts.y), itl2::ceil((double)refDims.z / (double)blockCounts.z));
		
		return refSize;
	}

	/**
	Divides input images into blocks for distributed processing.
	@return Vector of block lists. Each list contains corresponding block of each argument. Each item has five fields:
	position where image read operation starts, size of block to read, position in file where output data should be written,
	position in the image where output data region starts, and size of block to write.
	*/
	vector<vector<DistributionChunk> > determineBlocks(size_t refIndex, const Vec3c& margin, const Vec3c& blockCounts, const Command* command, vector<ParamVariant>& args, const vector<ParamVariant>& argsForGetCorrespondingBlock)
	{
		Vec3c refSize = blockSizeFromBlockCounts(refIndex, blockCounts, args);
		return determineBlocksFromBlockSize(refIndex, margin, refSize, command, args, argsForGetCorrespondingBlock);
	}

	/**
	Calculates memory requirement given image blocks and command.
	allBlocks[parameter index][block index] = tuple<read start, read size, output position in file, output start position in image, output region size>
	*/
	size_t calcRequiredMemory(const vector<vector<DistributionChunk> >& allBlocks, const Command* command)
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
				Vec3c readSize = tup.readSize;
				blockMems[n] += readSize.x * readSize.y * readSize.z * pixelSizes[parami];
			}
		}

		// Return maximum data size
		return max(blockMems);
	}

	/**
	Calculates size of largest block in the given block list in pixels.
	*/
	size_t maxBlockSize(const vector<DistributionChunk>& blocks)
	{
		size_t maxSize = 0;
		for (size_t n = 0; n < blocks.size(); n++)
		{
			const auto& tup = blocks[n];
			Vec3c readSize = tup.readSize;
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
		map<DistributedImageBase*, vector<DistributionChunk> > blocksPerImage;
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

	/**
	Finds third distribution direction.
	*/
	size_t getDistributionDirection3(const vector<Delayed>& delayedCommands)
	{
		if (delayedCommands.size() <= 0)
			throw logic_error("No commands.");

		size_t result = delayedCommands[0].getDistributionDirection3();

		if (result >= numeric_limits<size_t>::max())
			return numeric_limits<size_t>::max();

		for (size_t n = 1; n < delayedCommands.size(); n++)
		{
			size_t dir = delayedCommands[n].getDistributionDirection3();

			// Any command that does not allow distribution in second direction disables multi-direction distribution.
			if (dir >= numeric_limits<size_t>::max())
				return numeric_limits<size_t>::max();

			// If commands have different allowed second directions, we disable second direction altogether.
			if (dir != result)
				return numeric_limits<size_t>::max();
		}

		return result;
	}

	void Distributor::flush()
	{
		runDelayedCommands();
	}

	/**
	Checks if the output images can be written to in arbitrary blocks.
	*/
	bool canWriteArbitaryBlocks(set<DistributedImageBase*>& outputImages)
	{
		for (DistributedImageBase* img : outputImages)
		{
			// NOTE: Here we list the files that support distribution in all directions. This way it is not possible
			// to induce a bug here by adding new DistributedImageStorageTypes.
			if (img->currentWriteTargetType() != DistributedImageStorageType::NN5 &&
				img->currentWriteTargetType() != DistributedImageStorageType::Raw)
			{
				return false;
			}
		}

		return true;
	}


	void Distributor::determineDistributionConfiguration(Vec3c& margin, set<DistributedImageBase*>& inputImages, set<DistributedImageBase*>& outputImages, JobType& jobType, map<DistributedImageBase*, vector<DistributionChunk> >& blocksPerImage, size_t& memoryReq)
	{
		if (delayedCommands.size() <= 0)
			return;

		margin = getMaxMargin(delayedCommands);

		// Find all input and output images
		inputImages = getInputImages(delayedCommands);
		outputImages = getOutputImages(delayedCommands);

		// Get job type
		jobType = getCombinedJobType(delayedCommands);

		// Find NN5 chunk size and use that as preferred block size.
		Vec3c preferredBlockSizeMultiple(1, 1, 1);
		for (DistributedImageBase* p : outputImages)
		{
			if (p->currentWriteTargetType() == DistributedImageStorageType::NN5)
			{
				preferredBlockSizeMultiple = p->getChunkSize();
			}
		}


		// Find distribution directions
		size_t distributionDirection1 = getDistributionDirection1(delayedCommands);
		size_t distributionDirection2 = getDistributionDirection2(delayedCommands);
		size_t distributionDirection3 = getDistributionDirection3(delayedCommands);

		if (distributionDirection1 > 2)
			throw logic_error("Invalid distribution direction.");

		if (distributionDirection2 == distributionDirection1)
			throw logic_error("Distribution directions 1 and 2 can't be the same.");

		if (distributionDirection3 <= 2 &&
			(distributionDirection3 == distributionDirection1 || distributionDirection3 == distributionDirection2))
			throw logic_error("Distribution direction 3 can't be the same than 1 or 2.");

		// If we need to distribute in two directions, check that all output images support it.
		// Otherwise distribution is not allowed as e.g. sequences can't be written to in parallel.
		bool distributionDirection2Allowed = distributionDirection2 <= 2;
		bool distributionDirection3Allowed = distributionDirection2Allowed && distributionDirection3 <= 2;
		if (distributionDirection2Allowed || distributionDirection3Allowed)
		{
			if (!canWriteArbitaryBlocks(outputImages))
			{
				distributionDirection2Allowed = false;
				distributionDirection3Allowed = false;
			}
		}

        size_t preferredSubdivisions = getPreferredSubdivisions(delayedCommands);

		// Determine blocks that must be loaded, given amount of subdivisions in each direction.
		// blocksPerImage[image pointer][block index] = tuple<block definition>
		Vec3c subDivisions(1, 1, 1);
		subDivisions[distributionDirection1] = preferredSubdivisions;
		size_t lastMemoryReq = numeric_limits<size_t>::max();
		Vec3<bool> preferredBlockSizeTested(false, false, false);
		while (true)
		{
			// Calculate block sizes for each command.
			blocksPerImage.clear();
			vector<size_t> extraMemPerCommand;
			Vec3c refSize;
			for (Delayed& d : delayedCommands)
			{
				// blocksPerParameter[parameter index][block index] = tuple<block definition>
				//vector<vector<tuple<Vec3c, Vec3c, Vec3c, Vec3c, Vec3c> > > blocksPerParameter = determineBlocks(d.getRefIndex(), margin, subDivisions, d.getCommand(), d.getArgs(), d.getArgsForGetCorrespondingBlock());

				// Find block size, taking into account subdivisions and preferred block size multiple.
				refSize = blockSizeFromBlockCounts(d.getRefIndex(), subDivisions, d.getArgs());

				// Round block size down to the nearest preferred block size multiple.
				Vec3c temp = preferredBlockSizeMultiple;
				for (size_t n = 0; n < temp.size(); n++)
				{
					if (subDivisions[n] != 1)
					{
						if (temp[n] < refSize[n])
						{
							// preferredBlockSizeMultiple is smaller than refSize, so we can try to increase block size
							// in multiples of preferredBlockSizeMultiple until we go above refSize.
							while (temp[n] + preferredBlockSizeMultiple[n] <= refSize[n])
								temp[n] += preferredBlockSizeMultiple[n];
						}
						else
						{
							// refSize is smaller than preferredBlockSizeMultiple
							
							// Find fraction of preferredBlockSizeMultiple that is smaller than the block size
							coord_t preferredFraction = preferredBlockSizeMultiple[n];
							while (preferredFraction > 1 && preferredFraction > refSize[n])
								preferredFraction /= 2;

							// Find how many fractions we can take until we go above estimated refSize.
							// Note that this way we might test e.g. size preferredBlockSizeMultiple[n]
							// 2 times but that should not hurt.
							coord_t val = 0;
							while (val < refSize[n])
								val += preferredFraction;
							temp[n] = val;
							
							//if (!preferredBlockSizeTested[n])
							//{
							//	preferredBlockSizeTested[n] = true;
							//	temp[n] = preferredBlockSizeMultiple[n];
							//}
							////  Reduce size until
							////// it is smaller than refSize.
							////while (temp[n] > refSize[n])
							////	temp[n] /= 2;

							////// If we fail, go back to original refSize.
							////if (temp[n] < 1)
							////	temp[n] = refSize[n];
							//temp[n] = refSize[n];
						}
					}
					else
					{
						// 1 subdivision equals the full size.
						temp[n] = refSize[n];
					}
				}
				
				refSize = temp;

				// blocksPerParameter[parameter index][block index] = <block definition>
				vector<vector<DistributionChunk> > blocksPerParameter = determineBlocksFromBlockSize(d.getRefIndex(), margin, refSize, d.getCommand(), d.getArgs(), d.getArgsForGetCorrespondingBlock());


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
			{
				cout << "Block size = " << refSize << ", preferred multiple = " << preferredBlockSizeMultiple << endl;
				break;
			}

			// Condition 1: memory requirement must go down when we increase subdivisions.
			if (memoryReq > lastMemoryReq)
			{
				throw ITLException(string("Unable to find suitable subdivision. Memory requirement does not decrease from ") + bytesToString((double)lastMemoryReq) + " when decreasing processing block size, and only " + bytesToString((double)allowedMemory()) + " is available for a single process. See also cluster configuration file max_memory setting.");
			}

			// Condition 2: Images must not be subdivided more than their dimensions allow.
			for (DistributedImageBase* img : inputImages)
			{
				// TODO: This condition is a bit "hacky" and will probably fail some time.
				if ((img->dimensions()[distributionDirection1] > 1 && subDivisions[distributionDirection1] >= img->dimensions()[distributionDirection1]) ||
					(distributionDirection2Allowed && img->dimensions()[distributionDirection2] > 1 && subDivisions[distributionDirection2] >= img->dimensions()[distributionDirection2]) ||
					(distributionDirection3Allowed && img->dimensions()[distributionDirection3] > 1 && subDivisions[distributionDirection3] >= img->dimensions()[distributionDirection3])
					)
				{
					throw ITLException(string("Unable to find suitable subdivision. The smallest possible blocks of the input and output images require ") + bytesToString((double)memoryReq) + " of memory, but only " + bytesToString((double)allowedMemory()) + " is available for a single process. See also cluster configuration file max_memory setting.");
				}
			}

			lastMemoryReq = memoryReq;

			// The logic below:
			// If there is preferredBlockSizeMultiple:
			// 	Divide in distributionDirection1 until refSize[distributionDirection1] <= preferredBlockSizeMultiple[distributionDirection1]
			// 	Then divide in distributionDirection2 if it is available, until refSize[distributionDirection2] <= preferredBlockSizeMultiple[distributionDirection2]
			//	Then divide in distributionDirection3 if it is available, until refSize[distributionDirection3] <= preferredBlockSizeMultiple[distributionDirection3]
			//	Then divide the direction that has the least subdivisions.
			// else:
			//	Divide in distributionDirection1 and 2 and 3 like before
			// => Potentially less unsafe blocks than with just the "else" part.

			Vec3c oldSubDivisions = subDivisions;
			if (preferredBlockSizeMultiple != Vec3c(1, 1, 1))
			{
				if (refSize[distributionDirection1] > preferredBlockSizeMultiple[distributionDirection1])
				{
					subDivisions[distributionDirection1]++;
				}
				else
				{
					if (distributionDirection2Allowed)
					{
						if (refSize[distributionDirection2] > preferredBlockSizeMultiple[distributionDirection2])
						{
							subDivisions[distributionDirection2]++;
						}
						else
						{
							if (distributionDirection3Allowed)
							{
								if (refSize[distributionDirection3] > preferredBlockSizeMultiple[distributionDirection3])
								{
									subDivisions[distributionDirection3]++;
								}
							}
						}
					}
				}
			}
			
			if (subDivisions == oldSubDivisions)
			{
				// No preferred block size or we could not continue subdivisions accounting for preferred block size.
				// Prefer distributing in first distribution direction.
				
				// TODO: Try replacing the subDivisions[distributionDirection1] < N condition by a condition that compares total pixel count in margins to total size of image.
				if (!distributionDirection2Allowed || subDivisions[distributionDirection1] < 80)
				{
					// Distribution direction 2 is not available or there are only a few subdivisions in direction 1
					// Favor direction 1 as it is often z and distributing in that direction is faster and more compatible than other directions.
					subDivisions[distributionDirection1]++;
				}
				else
				{
					if (!distributionDirection2Allowed)
						throw ITLException("The input images are so large that they must be distributed in two coordinate directions. Distribution in two directions is currently not supported for current input image type.");

					// Subdivide in 2 directions
					if (subDivisions[distributionDirection2] < subDivisions[distributionDirection1])
						subDivisions[distributionDirection2]++;
					else
						subDivisions[distributionDirection1]++;
					// TODO: Subdivide in the third direction, too.
				}
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

	/**
	Combines jobs such that their number is approximately halved.
	*/
	vector<tuple<string, JobType>> combineSmallJobs(const vector<tuple<string, JobType>>& jobsToSubmit)
	{
		vector<tuple<string, JobType>> newJobs;
		for (size_t n = 0; n < jobsToSubmit.size(); n += 2)
		{
			if (n + 1 < jobsToSubmit.size())
			{
				const auto& a = jobsToSubmit[n];
				const auto& b = jobsToSubmit[n + 1];
				string newScript = get<0>(a) + "\n\n\n\n" + get<0>(b);
				JobType at = get<1>(a);
				JobType bt = get<1>(b);
				// Never combine multiple slow jobs, but do combine faster jobs to slow ones.
				if (at != JobType::Slow || bt != JobType::Slow)
				{
					// New job type will be the slower of the two job types.
					JobType newType;
					if (at == JobType::Slow || bt == JobType::Slow)
						newType = JobType::Slow;
					else if (at == JobType::Normal || bt == JobType::Normal)
						newType = JobType::Normal;
					else
						newType = JobType::Fast;
					newJobs.push_back(make_tuple(newScript, newType));
				}
			}
			else
			{
				// Only one job left, so move that to the new job list.
				newJobs.push_back(jobsToSubmit[n]);
			}
		}

		return newJobs;
	}

	JobType promote(JobType t)
	{
		if (t == JobType::Fast)
			return JobType::Normal;
		if (t == JobType::Normal)
			return JobType::Slow;
		return JobType::Slow;
	}

	/**
	If multiple smaller jobs are combined into one, their output is also combined into one element in the job output array.
	This output extract output for each individual job from the combined output.
	*/
	vector<string> separateCombinedJobOutput(const vector<string>& outputs, size_t originalJobCount, const string& jobStartLine)
	{
		vector<string> newOutput(originalJobCount, "");
		for (string output : outputs)
		{
			vector<string> lines = itl2::split(output, true);

			// Remove all lines that print jobStartLine
			vector<string> lines2;
			for (size_t n = 0; n < lines.size(); n++)
			{
				if (!startsWith(lines[n], string("print(\"") + jobStartLine))
				{
					// Normal line
					lines2.push_back(lines[n]);
				}
				else
				{
					// This is a print statement that starts a new job. Remove preceding 'clear("")' if any
					// Erase clear("") before the print, if any
					if (lines2[lines2.size() - 1] == "clear(\"\")")
						lines2.erase(lines2.begin() + lines2.size() - 1);
				}
			}

			lines = lines2;

			// Each "------ start of job X" line starts a new block of output
			size_t currentLine = 0;
			while (currentLine < lines.size())
			{
				// Check that we start with jobStartLine
				if (!startsWith(lines[currentLine], jobStartLine))
				{
					// Skip over extra output, e.g. from user commands before starting pi2.
					currentLine++;
					//throw ITLException("Invalid combined job output. No job start marker at the first line of a new block.");
				}
				else
				{
					// Extract job index
					string indexStr = lines[currentLine].substr(jobStartLine.length());
					itl2::trim(indexStr);
					coord_t index = itl2::fromString<coord_t>(indexStr);

					currentLine++;

					// Extract lines until the next jobStartLine
					string content;
					while (currentLine < lines.size() && !startsWith(lines[currentLine], jobStartLine))
					{
						content += lines[currentLine] + "\n";
						currentLine++;
					}

					trim(content);

					if (index < 0 || (size_t)index >= originalJobCount)
						throw ITLException(string("Invalid job index ") + indexStr);

					if (newOutput[index] != "")
						throw ITLException(string("Multiple outputs for job ") + itl2::toString(index));

					newOutput[index] = content;
				}
			}
			
		}

		return newOutput;
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
		map<DistributedImageBase*, vector<DistributionChunk> > blocksPerImage;
		size_t memoryReq;
		determineDistributionConfiguration(margin, inputImages, outputImages, jobType, blocksPerImage, memoryReq);


		// If overlap is nonzero, InOut images must be saved to different file from which they are loaded.
		// Changing write targets should not cause bad state of image objects even if writeComplete() is not called.
		//if (margin != Vec3c(0, 0, 0))
		//{
		for (DistributedImageBase* img : inputImages)
		{
			if (outputImages.find(img) != outputImages.end())
			{
				// This is InOut image

				// TODO: NN5 should be able to support concurrent reading and writing. Is there some bug somewhere?
				if (//img->currentWriteTargetType() == DistributedImageStorageType::NN5 ||
					margin != Vec3c(0, 0, 0))
				{
					if (img->currentReadSource() == img->currentWriteTarget())
					{
						img->newWriteTarget();
					}
				}
			}
		}
		//}


		
		
		// No skipping jobs if there are InOut images for which
		// currentReadSource() and currentWriteTarget() data types are different, as in that case
		// we need to copy data from read source to write target with file type conversion.
		bool jobSkippingAllowed = true;
		for (DistributedImageBase* img : inputImages)
		{
			if (outputImages.find(img) != outputImages.end())
			{
				// This is InOut image; test input and output data types.

				if (!img->isSavedToDisk())
				{
					cout << "Job skipping is not allowed as there are in-place processed images that are not saved to the disk yet." << endl;
					jobSkippingAllowed = false;
					break;
				}

				if (img->currentReadSource() != img->currentWriteTarget())
				{
					cout << "Job skipping is not allowed as there are in-place processed images that need to be copied from the input file to the output file." << endl;
					jobSkippingAllowed = false;
					break;
				}
			}
		}

		size_t jobCount = blocksPerImage.begin()->second.size();
		cout << "Submitting " << jobCount << " jobs, each estimated to require at most " << bytesToString((double)memoryReq) << " of RAM per job..." << endl;
		vector<size_t> skippedJobs;
		vector<tuple<string, JobType>> jobsToSubmit;
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

			// Init so that we always print something (required at least in the SLURM distributor)
			script << "echo(true, false);" << endl;

			// Image read commands
			for(DistributedImageBase* img : inputImages)
			{
				Vec3c readStart = blocksPerImage[img][i].readStart;
				Vec3c readSize = blocksPerImage[img][i].readSize;
				script << img->emitReadBlock(readStart, readSize, true);
			}

			// Output image creation commands
			for (DistributedImageBase* img : outputImages)
			{
				if (inputImages.find(img) == inputImages.end())
				{
					Vec3c readStart = blocksPerImage[img][i].readStart;
					Vec3c readSize = blocksPerImage[img][i].readSize;
					script << img->emitReadBlock(readStart, readSize, false);
				}
			}

			// Processing commands
			bool hasCommandsToRun = false;
			for (size_t cmdi = 0; cmdi < delayedCommands.size(); cmdi++)
			{
				const Command* command = delayedCommands[cmdi].getCommand();
				vector<ParamVariant>& args = delayedCommands[cmdi].getArgs();
				size_t refIndex = delayedCommands[cmdi].getRefIndex();

				DistributedImageBase* refImage = getDistributedImage(args[refIndex]);
				Vec3c readStart = blocksPerImage[refImage][i].readStart;
				Vec3c readSize = blocksPerImage[refImage][i].readSize;
				Vec3c writeFilePos = blocksPerImage[refImage][i].writeFilePos;
				Vec3c writeImPos = blocksPerImage[refImage][i].writeImPos;
				Vec3c writeSize = blocksPerImage[refImage][i].writeSize;
				Vec3c chunkIndex3 = blocksPerImage[refImage][i].chunkIndex3;
				if (delayedCommands[cmdi].needsToRun(readStart, readSize, writeFilePos, writeImPos, writeSize, i, chunkIndex3))
				{
					hasCommandsToRun = true;
					script << command->name() << "(";
					for (size_t n = 0; n < args.size(); n++)
					{
						// Value of argument whose type is Vec3c and name is "block origin" is replaced by the origin of current calculation block.
						// This functionality is needed at least in skeleton tracing command.
						const CommandArgumentBase& argDef = command->args()[n];
						ParamVariant argVal = args[n];
						if (argDef.dataType() == parameterType<BLOCK_ORIGIN_ARG_TYPE>() && argDef.name() == BLOCK_ORIGIN_ARG_NAME)
						{
							argVal = readStart;
						}
						else if (argDef.dataType() == parameterType<BLOCK_INDEX_ARG_TYPE>() && argDef.name() == BLOCK_INDEX_ARG_NAME)
						{
							argVal = (coord_t)i;
						}
						else if (argDef.dataType() == parameterType<BLOCK_INDEX3_ARG_TYPE>() && argDef.name() == BLOCK_INDEX3_ARG_NAME)
						{
							argVal = chunkIndex3;
						}

						script << "\"" << argumentToString(argDef, argVal) << "\"";
						if (n < args.size() - 1)
							script << ", ";
					}
					script << ");" << endl;
				}
			}

			// Image write commands
			for (DistributedImageBase* img : outputImages)
			{
				// Only write if the image is still visible from the main PI system object
				if (piSystem->isDistributedImage(img))
				{
					Vec3c writeFilePos = blocksPerImage[img][i].writeFilePos;
					Vec3c writeImPos = blocksPerImage[img][i].writeImPos;
					Vec3c writeSize = blocksPerImage[img][i].writeSize;

					// Only write if writing is requested by the command.
					if(writeSize.min() > 0)
						script << img->emitWriteBlock(writeFilePos, writeImPos, writeSize);
				}
			}

			if (hasCommandsToRun || !jobSkippingAllowed)
			{
				jobsToSubmit.push_back(make_tuple(script.str(), jobType));
			}
			else
			{
				skippedJobs.push_back(i);
			}
		}

		if (skippedJobs.size() > 0)
		{
			if (skippedJobs.size() == 1)
			{
				cout << "Job " << skippedJobs[0] << " is skipped as it was determined that the job will do nothing." << endl;
			}
			else
			{
				cout << "Jobs ";
				for (size_t n = 0; n < skippedJobs.size(); n++)
				{
					if (n > 0)
						cout << ", ";
					cout << skippedJobs[n];
				}
				cout << " were skipped as it was determined that the jobs will do nothing." << endl;
			}
		}

		Timer timer;
		timer.start();

		// Prepare images for concurrent writing
		map<DistributedImageBase*, vector<nn5::NN5Process>> nn5processes;
		for (size_t i = 0; i < jobCount; i++)
		{
			if (std::find(skippedJobs.begin(), skippedJobs.end(), i) == skippedJobs.end())
			{
				// The job i is not skipped
				for (const auto& item : blocksPerImage)
				{
					DistributedImageBase* img = item.first;
					if (outputImages.find(img) != outputImages.end())
					{
						const auto& valueList = item.second;
						Vec3c readStart = valueList[i].readStart;
						Vec3c readSize = valueList[i].readSize;
						Vec3c writeFilePos = valueList[i].writeFilePos;
						Vec3c writeImPos = valueList[i].writeImPos;
						Vec3c writeSize = valueList[i].writeSize;
						if (img->currentReadSource() == img->currentWriteTarget() &&
							inputImages.find(img) != inputImages.end())
						{
							// We are reading and writing the same file.
							// AND the image is used as input and output.
							// => The NN5Process must contain both read and write region.
							nn5processes[img].push_back(nn5::NN5Process{ AABoxc::fromPosSize(readStart, readSize), AABoxc::fromPosSize(writeFilePos, writeSize) });
						}
						else
						{
							// We are reading and writing different file.
							// => The NN5Process must contain only the write region as no reads are made from the same file.
							nn5processes[img].push_back(nn5::NN5Process{ AABoxc::fromPosSize(Vec3c(-1, -1, -1), Vec3c(0, 0, 0)), AABoxc::fromPosSize(writeFilePos, writeSize) });
						}
					}
				}
			}
		}

		// Check for overlapping writes. Those must not happen.
		for (DistributedImageBase* img : outputImages)
		{
			const auto& list = nn5processes[img];
			for (size_t n = 0; n < list.size(); n++)
			{
				const auto& p1 = list[n];
				for (size_t m = n + 1; m < list.size(); m++)
				{
					const auto& p2 = list[m];
					if (p1.writeBlock.overlapsExclusive(p2.writeBlock))
						throw ITLException(string("Multiple jobs would write to the same block in image ") + img->varName());
				}
			}
		}

		// Call startConcurrentWrite for all output images.
		for (DistributedImageBase* img : outputImages)
		{
			const auto& list = nn5processes[img];
			size_t unsafeCount = img->startConcurrentWrite(list);
			if (unsafeCount > 0)
				cout << "Image " << img->varName() << " has " << unsafeCount << " unsafe chunks." << endl;
		}

		Timing::Add(TimeClass::WritePreparation, timer.lap());

		// Combine small jobs
		const string jobStartLine = "------ start of job";
		size_t combinationRounds = 0;
		size_t originalJobCount = jobsToSubmit.size();
		if (maxSubmittedJobCount > 0)
		{
			
			if (jobsToSubmit.size() > maxSubmittedJobCount)
			{
				// Add clear command and job start marker to all the job scripts.
				for (size_t n = 0; n < jobsToSubmit.size(); n++)
				{
					get<0>(jobsToSubmit[n]) = "clear();\nprint('" + jobStartLine + " " + itl2::toString(n) + "'); \n" + get<0>(jobsToSubmit[n]);
				}

				// Combine jobs until the desired job count is reached.
				while (jobsToSubmit.size() > maxSubmittedJobCount)
				{
					size_t oldCount = jobsToSubmit.size();
					jobsToSubmit = combineSmallJobs(jobsToSubmit);

					// Do not continue if no jobs could be combined.
					if (jobsToSubmit.size() >= oldCount)
						break;

					combinationRounds++;
				}
			}
		}


		double tasksPerJob = (double)originalJobCount / (double)jobsToSubmit.size();

		if (combinationRounds > 0)
			cout << "Small jobs were combined into " << jobsToSubmit.size() << " larger jobs (" << std::fixed << std::setprecision(1) << tasksPerJob << " small jobs per combined job)." << endl;

		// Submit jobs
		for (auto& tup : jobsToSubmit)
		{
			string& script = get<0>(tup);
			JobType type = get<1>(tup);

			if (showSubmittedScripts)
			{
				cout << "Submitting pi2 script:" << endl;
				cout << script << endl;
			}

			// Promote jobs to slower queues if they are combined a lot.
			if (tasksPerJob >= promoteThreshold)
				type = promote(type);

			submitJob(script, type);
		}

		// Run jobs first and set writeComplete() only after the jobs have finished to make sure that
		// the output image exists before writeComplete() is called.
		// Additionally, this order ensures that if jobs fail, the input images still point to the correct files.
		try
		{
			cout << "Waiting for jobs to finish..." << endl;
			lastOutput = waitForJobs();

			Timing::Add(TimeClass::JobsInclQueuing, timer.lap());

			// Submit endConcurrentWrite jobs
			// Limit the number of jobs to reduce queuing latency.
			cout << "Submitting write finalization jobs..." << endl;
			size_t finJobCount = 0;
			for (DistributedImageBase* img : outputImages)
			{
				vector<Vec3c> chunks = img->getChunksThatNeedEndConcurrentWrite();
				
				size_t maxFinJobCount = maxSubmittedJobCount;
				if (maxFinJobCount <= 0)
					maxFinJobCount = 8;	// Some sane default value here.

				vector<string> scripts(maxFinJobCount, "");
				size_t n = 0;
				for (Vec3c& chunk : chunks)
				{
					string script = img->emitEndConcurrentWrite(chunk);
					if (script != "")
					{
						scripts[n] += script + "\n";
						n++;
						if (n >= scripts.size())
							n = 0;
					}
				}
				for (const string& script : scripts)
				{
					if (script != "")
					{
						submitJob(script, JobType::Normal);
						finJobCount++;
					}
				}
			}

			if (finJobCount > 0)
			{
				cout << "Waiting for jobs to finish..." << endl;
				// NOTE: Do not assign last output here, as that would override the output of the real commands
				// by the output of these new housekeeping jobs.
				waitForJobs();
			}
			else
			{
				cout << "No write finalization jobs were necessary." << endl;
			}
			
			for (DistributedImageBase* img : outputImages)
			{
				img->writeComplete();
			}

			Timing::Add(TimeClass::WriteFinalizationInclQueuing, timer.lap());


			if (combinationRounds > 0)
			{
				// Multiple jobs were combined into one, so we need to separate job output
				// so that each job has its own output in the lastOutput array.
				lastOutput = separateCombinedJobOutput(lastOutput, originalJobCount, jobStartLine);
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
