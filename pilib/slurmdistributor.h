#pragma once

#include "distributor.h"
#include "argumentdatatype.h"

namespace pilib
{
	/**
	Distributes tasks using SLURM cluster manager system.
	*/
	class SLURMDistributor : public Distributor
	{
	private:
		size_t allowedMem;

		virtual void submitJob(const string& piCode);

		virtual vector<string> waitForJobs();

		virtual size_t allowedMemory() const
		{
			return allowedMem;
		}

		/**
		Stores job ids of all submitted jobs since last call to waitForJobs.
		*/
		vector<size_t> submittedJobs;

		/**
		Extra arguments for sbatch and sinfo.
		*/
		string extraArgs;

		/**
		Commands run on each node before pi.
		*/
		string jobInitCommands;

		/**
		Calculate progress of all jobs in submittedJobs array.
		*/
		vector<int> getJobProgress() const;

		/**
		Checks if the given job has finished.
		*/
		bool isJobDone(size_t jobIndex) const;

		/**
		Gets log of given job.
		*/
		string getLog(size_t jobIndex, bool flush = false) const;

		/**
		Gets job progress from log file.
		*/
		int getJobProgressFromLog(size_t jobIndex) const;

		/**
		Get character that is used to report progress of a task.
		*/
		char getProgressChar(int progress) const;

	public:
		SLURMDistributor(PISystem* system);
	};
}
