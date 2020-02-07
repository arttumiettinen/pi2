#pragma once

#include "distributor.h"

namespace pilib
{
	/**
	Runs tasks sequentially on the local computer.
	*/
	class LocalDistributor : public Distributor
	{
	private:
		size_t allowedMem;

		/**
		Stores output of each subprocess that has been run since last call to waitForJobs.
		*/
		std::vector<std::string> outputs;

	public:
		LocalDistributor(PISystem* system);

		virtual void submitJob(const std::string& piCode, JobType jobType) override;

		virtual std::vector<std::string> waitForJobs() override;

		virtual size_t allowedMemory() const override
		{
			return allowedMem;
		}

		virtual void allowedMemory(size_t maxMem) override;
	};
}
