
#include "distributor.h"
#include "command.h"
#include "distributedimage.h"
#include "utilities.h"
#include "pisystem.h"
#include "exeutils.h"

#include <tuple>
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;

using namespace std;

namespace pilib
{
	const string Distributor::BLOCK_INDEX_ARG_NAME = "block index";
	const string Distributor::BLOCK_ORIGIN_ARG_NAME = "block origin";

	fs::path findPi2()
	{
#if defined(__linux__)
		string piName = "pi2";
#elif defined(_WIN32)
		string piName = "pi2.exe";
#else
	#error Please configure file name style for this platform.
#endif

		fs::path p1 = fs::current_path() / piName;
		fs::path p2 = fs::path(getExecutablePath());
		if(p2.has_filename())
			p2 = p2.replace_filename(piName);

		if (fs::exists(p1))
			return p1;
		if (fs::exists(p2))
			return p2;

		throw ITLException("Unable to find " + piName + " program. It was searched from " + p1.string() + " and " + p2.string());
	}

	Distributor::Distributor(PISystem* piSystem) : piSystem(piSystem)
	{
		piCommand = findPi2().string();
	}

	/**
	Determine block size for each parameter image.
	The images should be divided to subDivisions[i] blocks in dimension i.
	*/
	vector<Vec3c> calcBlockSizes(const Command* command, const vector<ParamVariant>& args, const Vec3c& subDivisions)
	{
		vector<Vec3c> blocks;
		for (size_t n = 0; n < args.size(); n++)
		{
			const CommandArgumentBase& def = command->args()[n];
			size_t ps = pixelSize(def.dataType());
			if (ps > 0)
			{
				const DistributedImageBase* img = args[n].dimgval;
				Vec3c dims = img->dimensions();
				blocks.push_back(Vec3c((coord_t)ceil((double)dims.x / (double)subDivisions.x), (coord_t)ceil((double)dims.y / (double)subDivisions.y), (coord_t)ceil((double)dims.z / (double)subDivisions.z)));
			}
		}
		return blocks;
	}

	//Vec3c calcBlockSize(const Vec3c& imageDimensions, const Vec3c& subDivisions)
	//{
	//	return Vec3c((coord_t)ceil((double)imageDimensions.x / (double)subDivisions.x), (coord_t)ceil((double)imageDimensions.y / (double)subDivisions.y), (coord_t)ceil((double)imageDimensions.z / (double)subDivisions.z));
	//}

	Vec3c getFirstImageDimensions(const Command* command, const vector<ParamVariant>& args)
	{
		for (size_t n = 0; n < args.size(); n++)
		{
			const CommandArgumentBase& def = command->args()[n];
			size_t ps = pixelSize(def.dataType());
			if (ps > 0)
			{
				const DistributedImageBase* img = args[n].dimgval;
				Vec3c dims = img->dimensions();
				return dims;
			}
		}

		throw ITLException("Command parameters contain no images.");
	}

	/**
	Calculates memory required by all images in the parameter list.
	@param subDivisions Count of subdivisions in each coordinate direction.
	*/
	size_t calcRequiredMemory(const Command* command, const vector<ParamVariant>& args, const Vec3c& subDivisions, double extraMemFactor)
	{
		vector<Vec3c> blocks = calcBlockSizes(command, args, subDivisions);

		// Calculate pixel sizes
		vector<size_t> pixelSizes;
		for (size_t n = 0; n < args.size(); n++)
		{
			const CommandArgumentBase& def = command->args()[n];
			size_t ps = pixelSize(def.dataType());
			if (ps > 0)
				pixelSizes.push_back(ps);
		}

		size_t totalSize = 0;
		for (size_t n = 0; n < args.size(); n++)
			totalSize += (size_t)blocks[n].x * (size_t)blocks[n].y * (size_t)blocks[n].z * pixelSizes[n];


		//// Find images from parameter list and calculate memory required to load them
		//size_t totalSize = 0;
		//for (size_t n = 0; n < args.size(); n++)
		//{
		//	const CommandArgumentBase& def = command->args()[n];
		//	size_t ps = pixelSize(def.dataType());
		//	if (ps > 0)
		//	{
		//		const DistributedImageBase* img = args[n].dimgval;
		//		Vec3c dims = img->dimensions();

		//		if (totalSize <= 0)
		//			firstImDimensions = dims;

		//		Vec3c blockDims = calcBlockSize(dims, subDivisions);
		//		coord_t size = blockDims.x * blockDims.y * blockDims.z * ps;
		//		totalSize += (size_t)size;
		//	}
		//}

		//if (totalSize <= 0)
		//	throw ITLException("Command parameters contain no images.");

		return (size_t)ceil(totalSize * (1 + extraMemFactor));
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
	Finds image that is used as reference (first output image or first input image if there are no outputs) in the task distribution process.
	Returns the index of the image in the argument list.
	*/
	size_t findReferenceImage(const Command* command)
	{
		coord_t firstImageIndex = -1;
		for (size_t n = 0; n < command->args().size(); n++)
		{
			const CommandArgumentBase& def = command->args()[n];
			if (isImage(def.dataType()))
			{
				if (firstImageIndex < 0)
					firstImageIndex = n;

				if (def.direction() == Out || def.direction() == InOut)
					return n;
			}
		}

		if(firstImageIndex < 0)
			throw ITLException("Command parameters contain no images.");

		return (size_t)firstImageIndex;
	}

	/**
	Divides input images into blocks for distributed processing.
	@return Vector of block lists. Each list contains corresponding block of each argument. Each item has five fields:
	position where image read operation starts, size of block to read, position in file where output data should be written,
	position in the image where output data region starts, and size of block to write.
	*/
	vector<vector<tuple<Vec3c, Vec3c, Vec3c, Vec3c, Vec3c> > > determineBlocks(size_t refIndex, const Vec3c& refDims, const Vec3c& margin, const Vec3c& blockCounts, const Command* command, vector<ParamVariant>& args, vector<ParamVariant>& argsForGetCorrespondingBlock)
	{
		// Block size in the reference image
		Vec3c refSize = Vec3c((coord_t)ceil((double)refDims.x / (double)blockCounts.x), (coord_t)ceil((double)refDims.y / (double)blockCounts.y), (coord_t)ceil((double)refDims.z / (double)blockCounts.z));

		// first index = block index,
		// second index = parameter index
		vector<vector<tuple<Vec3c, Vec3c, Vec3c, Vec3c, Vec3c> > > allBlocks;
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
					
					vector<tuple<Vec3c, Vec3c, Vec3c, Vec3c, Vec3c> > blocks;
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
							dynamic_cast<const Distributable*>(command)->getCorrespondingBlock(argsForGetCorrespondingBlock, n, readStart, readSize, writeFilePos, writeImPos, writeSize);
						}

						blocks.push_back(make_tuple(readStart, readSize, writeFilePos, writeImPos, writeSize));
					}

					allBlocks.push_back(blocks);
				}
			}
		}

		return allBlocks;
	}

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
		vector<size_t> blockMems;
		for (const auto& blocks : allBlocks)
		{
			size_t req = 0;
			for(size_t n = 0; n < blocks.size(); n++)
			{
				const auto& tup = blocks[n];
				Vec3c readSize = get<1>(tup);
				req += readSize.x * readSize.y * readSize.z * pixelSizes[n];
			}
			blockMems.push_back(req);
		}

		// Return maximum data size
		return max(blockMems);
	}

	vector<string> Distributor::distribute(const Command* command, vector<ParamVariant>& args, size_t distributionDirection, const Vec3c& margin, const string* outputFile, size_t refIndex, vector<ParamVariant>* argsForGetCorrespondingBlock)
	{
		if (distributionDirection > 2)
			throw ITLException("Invalid distribution direction.");

		if (args.size() != command->args().size())
			throw logic_error("Invalid number of arguments.");

		if (!argsForGetCorrespondingBlock)
			argsForGetCorrespondingBlock = &args;

		// Uh oh, it is not nice to cast... but I have not figured out how to use e.g. virtual inheritance and avoid writing messy constructor code
		// in each and every Command class constructor.
		const Distributable* distributable = dynamic_cast<const Distributable*>(command);
		if (!distributable)
			throw logic_error("Unable to distribute non-distributable command.");

		double extraMemFactor = distributable->calculateExtraMemory(args);

		// Find first output image or first input if there are no output images.
		if(refIndex > command->args().size())
			refIndex = findReferenceImage(command);

		// Determine blocks that must be loaded, given amount of subdivisions in each direction.
		const DistributedImageBase* img = args[refIndex].dimgval;
		Vec3c refDims = img->dimensions();
		Vec3c subDivisions(1, 1, 1);
		size_t memoryReq;
		vector<vector<tuple<Vec3c, Vec3c, Vec3c, Vec3c, Vec3c> > > blocks;
		while(true)
		{
			blocks = determineBlocks(refIndex, refDims, margin, subDivisions, command, args, *argsForGetCorrespondingBlock);

			memoryReq = calcRequiredMemory(blocks, command);

			if (memoryReq <= allowedMemory())
				break;

			if (subDivisions[distributionDirection] >= img->dimensions()[distributionDirection])
				throw ITLException(string("Unable to find suitable subdivision. The smallest possible blocks of the input and output images require ") + bytesToString((double)memoryReq) + " of memory, but only " + bytesToString((double)allowedMemory()) + " is available for a single process. See also cluster configuration file max_memory setting.");

			subDivisions[distributionDirection]++;
		}

		// If overlap is nonzero, InOut images must be saved to different file from which they are loaded.
		// Changing write targets should not cause bad state of image objects even if writeComplete() is not called.
		if (margin != Vec3c(0, 0, 0))
		{
			for (size_t n = 0; n < args.size(); n++)
			{
				const CommandArgumentBase& def = command->args()[n];
				if (isImage(def.dataType()) && def.direction() == InOut)
				{
					DistributedImageBase& img = *args[n].dimgval;
					if (img.currentReadSource() == img.currentWriteTarget())
					{
						img.newWriteTarget();
					}
				}
			}
		}

		//size_t jobCount = 0;
		//for (coord_t sz = 0; sz < firstDims.z; sz += blockSizes[0].z)
		//{
		//	for (coord_t sy = 0; sy < firstDims.y; sy += blockSizes[0].y)
		//	{
		//		for (coord_t sx = 0; sx < firstDims.x; sx += blockSizes[0].x)
		//		{
		//			

		//			

		//			// Build job script:
		//			// readblock(Block of input image 1)
		//			// readblock(Block of input image 2)
		//			// ...(for all input images)
		//			//
		//			// command
		//			//
		//			// writeblock(Block of output image 1)
		//			// writeblock(Block of output image 2)
		//			// ...(for all output images)
		//			//
		//			stringstream script;
		//			
		//			// Init so that we always print something (required at least in SLURM distributor)
		//			script << "echo(true, false);" << endl;

		//			// Image read commands
		//			for (size_t n = 0; n < args.size(); n++)
		//			{
		//				const CommandArgumentBase& def = command->args()[n];

		//				// emitReadBlock is called for both input and output (and inout) images so that output images are created
		//				// if they don't exist yet. Note that if an image doesn't exist and the corresponding file does not exist,
		//				// emitReadBlock becomes newimage command.
		//				// If image type is Out, it does not have to be read (even if the storage file exists).
		//				if (isImage(def.dataType()))
		//				{
		//					bool dataNeeded = def.direction() == In || def.direction() == InOut;
		//					script << args[n].dimgval->emitReadBlock(inPos, inSize, dataNeeded);
		//				}
		//			}

		//			// Processing command
		//			script << command->name() << "(";
		//			for (size_t n = 0; n < args.size(); n++)
		//			{
		//				// Value of argument whose type is Vec3c and name is "block origin" is replaced by the origin of current calculation block.
		//				// This functionality is needed at least in skeleton tracing command.
		//				const CommandArgumentBase& argDef = command->args()[n];
		//				ParamVariant argVal = args[n];
		//				if (argDef.dataType() == parameterType<BLOCK_ORIGIN_ARG_TYPE>() && argDef.name() == BLOCK_ORIGIN_ARG_NAME)
		//				{
		//					argVal.vix = inPos.x;
		//					argVal.viy = inPos.y;
		//					argVal.viz = inPos.z;
		//				}
		//				else if (argDef.dataType() == parameterType<BLOCK_INDEX_ARG_TYPE>() && argDef.name() == BLOCK_INDEX_ARG_NAME)
		//				{
		//					argVal.ival = jobCount;
		//				}

		//				script << piSystem->argumentToString(argDef, argVal, true);
		//				if (n < args.size() - 1)
		//					script << ", ";
		//			}
		//			script << ");" << endl;

		//			// Image write commands
		//			for (size_t n = 0; n < args.size(); n++)
		//			{
		//				const CommandArgumentBase& def = command->args()[n];
		//				// Only output images are saved.
		//				if (isImage(def.dataType()) && (def.direction() == Out || def.direction() == InOut))
		//				{
		//					script << args[n].dimgval->emitWriteBlock(Vec3c(sx, sy, sz), outPos, outSize, outputFile);
		//				}
		//			}

		//			//cout << "Submitting job script:" << endl;
		//			//cout << script.str() << endl;

		//			submitJob(script.str());
		//			jobCount++;
		//		}
		//	}
		//}

		cout << "Submitting " << blocks.size() << " jobs..." << endl;
		for(size_t i = 0; i < blocks.size(); i++)
		{
			// Build job script:
			// readblock(Block of input image 1)
			// readblock(Block of input image 2)
			// ...(for all input images)
			//
			// command
			//
			// writeblock(Block of output image 1)
			// writeblock(Block of output image 2)
			// ...(for all output images)
			//

			stringstream script;

			// Init so that we always print something (required at least in SLURM distributor)
			script << "echo(true, false);" << endl;

			// Image read commands
			for (size_t n = 0; n < args.size(); n++)
			{
				const CommandArgumentBase& def = command->args()[n];

				// emitReadBlock is called for both input and output (and inout) images so that output images are created
				// if they don't exist yet. Note that if an image doesn't exist and the corresponding file does not exist,
				// emitReadBlock becomes newimage command.
				// If image type is Out, it does not have to be read (even if the storage file exists).
				if (isImage(def.dataType()))
				{
					bool dataNeeded = def.direction() == In || def.direction() == InOut;

					Vec3c readStart = get<0>(blocks[i][n]);
					Vec3c readSize = get<1>(blocks[i][n]);

					script << args[n].dimgval->emitReadBlock(readStart, readSize, dataNeeded);
				}
			}

			// Processing command
			script << command->name() << "(";
			for (size_t n = 0; n < args.size(); n++)
			{
				// Value of argument whose type is Vec3c and name is "block origin" is replaced by the origin of current calculation block.
				// This functionality is needed at least in skeleton tracing command.
				const CommandArgumentBase& argDef = command->args()[n];
				ParamVariant argVal = args[n];
				if (argDef.dataType() == parameterType<BLOCK_ORIGIN_ARG_TYPE>() && argDef.name() == BLOCK_ORIGIN_ARG_NAME)
				{
					Vec3c readStart = get<0>(blocks[i][refIndex]);

					argVal.vix = readStart.x;
					argVal.viy = readStart.y;
					argVal.viz = readStart.z;
				}
				else if (argDef.dataType() == parameterType<BLOCK_INDEX_ARG_TYPE>() && argDef.name() == BLOCK_INDEX_ARG_NAME)
				{
					argVal.ival = i;
				}

				script << piSystem->argumentToString(argDef, argVal, true);
				if (n < args.size() - 1)
					script << ", ";
			}
			script << ");" << endl;

			// Image write commands
			for (size_t n = 0; n < args.size(); n++)
			{
				const CommandArgumentBase& def = command->args()[n];
				// Only output images are saved.
				if (isImage(def.dataType()) && (def.direction() == Out || def.direction() == InOut))
				{
					Vec3c writeFilePos = get<2>(blocks[i][n]);
					Vec3c writeImPos = get<3>(blocks[i][n]);
					Vec3c writeSize = get<4>(blocks[i][n]);
					script << args[n].dimgval->emitWriteBlock(writeFilePos, writeImPos, writeSize, outputFile);
				}
			}

			//cout << "Submitting job script:" << endl;
			//cout << script.str() << endl;

			submitJob(script.str());
		}

        // Run jobs first and set writeComplete() only after the jobs have finished to make sure that
        // the output image exists before writeComplete() is called.
        // Additionally, this order ensures that if jobs fail, the input images still point to the correct files.
		//cout << "Waiting for " << jobCount << " jobs to finish..." << endl;
		cout << "Waiting for jobs to finish..." << endl;
        vector<string> results = waitForJobs();
		
		for (size_t n = 0; n < args.size(); n++)
		{
			const CommandArgumentBase& def = command->args()[n];
			if (isImage(def.dataType()) && (def.direction() == Out || def.direction() == InOut))
			{
				args[n].dimgval->writeComplete();
			}
		}
		
		return results;
	}
}
