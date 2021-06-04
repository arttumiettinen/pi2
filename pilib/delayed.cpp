
#include "delayed.h"
#include "pisystem.h"

#include "command.h"
#include "distributable.h"

using namespace std;

namespace pilib
{
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

				if (def.direction() == ParameterDirection::Out || def.direction() == ParameterDirection::InOut)
					return n;
			}
		}

		if (firstImageIndex < 0)
			throw ITLException("Command parameters contain no images.");

		return (size_t)firstImageIndex;
	}

	Delayed::Delayed(PISystem* system, const Command* command, std::vector<ParamVariant>& args, std::vector<ParamVariant>* argsForGetCorrespondingBlock) :
		command(command),
		args(args),
		distributable(dynamic_cast<const Distributable*>(command))
	{
		if (!distributable)
			throw logic_error("Unable to use non-distributable command.");

		if (argsForGetCorrespondingBlock == 0)
			this->argsForGetCorrespondingBlock = args;
		else
			this->argsForGetCorrespondingBlock.insert(this->argsForGetCorrespondingBlock.begin(), argsForGetCorrespondingBlock->begin(), argsForGetCorrespondingBlock->end());

		// Get shared_ptr for each image argument so that they are not destroyed while waiting for distributed processing to start.
		for (size_t n = 0; n < args.size(); n++)
		{
			const CommandArgumentBase& def = command->args()[n];
			if (isImage(def.dataType()))
			{
				DistributedImageBase* img = getDistributedImage(args[n]);
				if (def.direction() == ParameterDirection::In || def.direction() == ParameterDirection::InOut)
				{
					inputImagePointers.push_back(system->getDistributedImagePointer(img));
					inputImages.push_back(inputImagePointers.back().get());
				}

				if (def.direction() == ParameterDirection::Out || def.direction() == ParameterDirection::InOut)
				{
					outputImagePointers.push_back(system->getDistributedImagePointer(img));
					outputImages.push_back(outputImagePointers.back().get());
				}
			}
		}

		// Find first output image or first input if there are no output images.
		refIndex = distributable->getRefIndex(args);
		if (refIndex > args.size())
			refIndex = findReferenceImage(command);
	}

	bool Delayed::canDelay() const
	{
		return distributable->canDelay(args);
	}

	Vec3c Delayed::getMargin() const
	{
		return distributable->getMargin(args);
	}

	size_t Delayed::getPreferredSubdivisions() const
	{
		return distributable->getPreferredSubdivisions(args);
	}

	size_t Delayed::getDistributionDirection1() const
	{
		return distributable->getDistributionDirection1(args);
	}

	size_t Delayed::getDistributionDirection2() const
	{
		return distributable->getDistributionDirection2(args);
	}

	size_t Delayed::getRefIndex() const
	{
		return refIndex;
	}

	double Delayed::calculateExtraMemory() const
	{
		return distributable->calculateExtraMemory(args);
	}

	JobType Delayed::getJobType() const
	{
		return distributable->getJobType(args);
	}

	bool Delayed::needsToRun(const Vec3c& readStart, const Vec3c& readSize, const Vec3c& writeFilePos, const Vec3c& writeImPos, const Vec3c& writeSize, size_t blockIndex) const
	{
		return distributable->needsToRunBlock(args, readStart, readSize, writeFilePos, writeImPos, writeSize, blockIndex);
	}
}