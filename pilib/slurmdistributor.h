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

		virtual void submitJob(const string& piCode, JobType jobType);

		virtual vector<string> waitForJobs();

		virtual size_t allowedMemory() const
		{
			return allowedMem;
		}

		void cancelAll();

		string getErrorMessage(size_t jobIndex) const;

		/**
		Stores (job slurm id, queue type, submission count) of all submitted jobs since last call to waitForJobs.
		*/
		vector<tuple<size_t, JobType, size_t> > submittedJobs;

		/**
		Extra arguments for sbatch and sinfo, for fast jobs
		*/
		string extraArgsFastJobs;

		/**
		Extra arguments for sbatch and sinfo, for normal jobs
		*/
		string extraArgsNormalJobs;

		/**
		Extra arguments for sbatch and sinfo, for slow jobs
		*/
		string extraArgsSlowJobs;

		/**
		Returns suitable sbatch arguments given type of job.
		*/
		string extraArgs(JobType jobType) const
		{
			switch (jobType)
			{
			case JobType::Fast: return extraArgsFastJobs;
			case JobType::Normal: return extraArgsNormalJobs;
			case JobType::Slow: return extraArgsSlowJobs;
			default: throw logic_error("Invalid JobType value.");
			}
		}

		/**
		Maximum number of submissions per job.
		This limits the count of re-submissions if jobs fail.
		*/
		size_t maxSubmissions;

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
		Gets SLURM error log of given job.
		*/
		string getSlurmErrorLog(size_t jobIndex, bool flush = false) const;

		/**
		Gets job progress from log file.
		*/
		int getJobProgressFromLog(size_t jobIndex) const;

		/**
		Get character that is used to report progress of a task.
		*/
		char getProgressChar(int progress) const;

		/**
		Create a progress bar string.
		*/
		string createProgressBar(const vector<int>& progress);

		/**
		Submits job with given index again.
		*/
		void resubmit(size_t jobIndex);

	public:
		SLURMDistributor(PISystem* system);
	};
}
