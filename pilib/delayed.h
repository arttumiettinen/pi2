#pragma once

#include <vector>
#include <memory>

#include "argumentdatatype.h"
#include "distributedimage.h"
#include "jobtype.h"

namespace pilib
{
	class Command;
	class Distributable;
	class PISystem;

	/**
	Stores command and its arguments for delayed execution.
	Keeps shared_ptrs to all the image arguments.
	*/
	class Delayed
	{
	private:
		/**
		Data needed to run a command later.
		*/
		const Command* command;
		const Distributable* distributable;

		size_t refIndex;

		std::vector<ParamVariant> args;
		std::vector<ParamVariant> argsForGetCorrespondingBlock;

		std::vector<DistributedImageBase*> inputImages;
		std::vector<DistributedImageBase*> outputImages;
		std::vector<std::shared_ptr<DistributedImageBase> > inputImagePointers;
		std::vector<std::shared_ptr<DistributedImageBase> > outputImagePointers;

	public:

		Delayed(PISystem* system, const Command* command, std::vector<ParamVariant>& args, std::vector<ParamVariant>* argsForGetCorrespondingBlock);

		const Command* getCommand() const
		{
			return command;
		}

		std::vector<ParamVariant>& getArgs()
		{
			return args;
		}

		const std::vector<ParamVariant>& getArgs() const
		{
			return args;
		}

		const std::vector<ParamVariant>& getArgsForGetCorrespondingBlock() const
		{
			return argsForGetCorrespondingBlock;
		}

		const std::vector<DistributedImageBase*>& getInputImages() const
		{
			return inputImages;
		}

		const std::vector<DistributedImageBase*> getOutputImages() const
		{
			return outputImages;
		}

		bool canDelay() const;

		Vec3c getMargin() const;

		size_t getPreferredSubdivisions() const;

		size_t getDistributionDirection1() const;

		size_t getDistributionDirection2() const;

		size_t getRefIndex() const;

		double calculateExtraMemory() const;

		JobType getJobType() const;
	};

}