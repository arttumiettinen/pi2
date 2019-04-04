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
		vector<string> outputs;

		virtual void submitJob(const string& piCode, JobType jobType);

		virtual vector<string> waitForJobs();

		virtual size_t allowedMemory() const
		{
			return allowedMem;
		}

	public:
		LocalDistributor(PISystem* system);
	};
}
