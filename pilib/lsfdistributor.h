#pragma once

#include "distributor.h"
#include "argumentdatatype.h"

namespace pilib
{
	/**
	Distributes tasks using LSF cluster manager system.
	*/
	class LSFDistributor : public Distributor
	{
	private:
		size_t allowedMem;

		std::string getErrorMessage(size_t jobIndex) const;

		/**
		Stores (job LSF id, queue type, submission count) of all submitted jobs since last call to waitForJobs.
		*/
		std::vector<std::tuple<size_t, JobType, size_t> > submittedJobs;

		/**
		Extra arguments for sbatch and sinfo, for fast jobs
		*/
		std::string extraArgsFastJobs;

		/**
		Extra arguments for sbatch and sinfo, for normal jobs
		*/
		std::string extraArgsNormalJobs;

		/**
		Extra arguments for sbatch and sinfo, for slow jobs
		*/
		std::string extraArgsSlowJobs;

		/**
		Commands used to run sbatch, squeue, scancel and sinfo.
		*/
		std::string bsubCommand, bjobsCommand, bkillCommand;

		/**
		Returns suitable sbatch arguments given type of job.
		*/
		std::string extraArgs(JobType jobType) const
		{
			switch (jobType)
			{
			case JobType::Fast: return extraArgsFastJobs;
			case JobType::Normal: return extraArgsNormalJobs;
			case JobType::Slow: return extraArgsSlowJobs;
			default: throw std::logic_error("Invalid JobType value.");
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
		std::string jobInitCommands;

		/**
		Identifies this running instance from others so that multiple LSF distributor instances can run from the same working folder.
		*/
		std::string myName;

		/**
		Cancels job with given LSF id.
		*/
		void cancelJob(size_t id) const;

		/**
		Cancel all jobs submitted by this object.
		*/
		void cancelAll() const;

		/**
		Calculate progress of all jobs in submittedJobs array.
		*/
		std::vector<int> getJobProgress() const;

		/**
		Checks if the given job has finished.
		*/
		bool isJobDone(size_t jobIndex) const;

		/**
		Gets log of given job.
		*/
		std::string getLog(size_t jobIndex, bool flush = false) const;

		/**
		Gets LSF error log of given job.
		*/
		std::string getErrorLog(size_t jobIndex, bool flush = false) const;

		/**
		Gets job progress from log file.
		*/
		int getJobProgressFromLog(size_t jobIndex) const;

		/**
		Submits job with given index again.
		*/
		void resubmit(size_t jobIndex);

		/**
		Creates unique name for a job.
		*/
		std::string makeJobName(size_t jobIndex) const;

		/**
		Creates input file name.
		*/
		std::string makeInputName(size_t jobIndex) const;

		/**
		Creates output log file name.
		*/
		std::string makeOutputName(size_t jobIndex) const;

		/**
		Creates error log file name.
		*/
		std::string makeErrorName(size_t jobIndex) const;

	public:
		LSFDistributor(PISystem* system);

		virtual void submitJob(const std::string& piCode, JobType jobType) override;

		virtual std::vector<std::string> waitForJobs() override;

		virtual size_t allowedMemory() const override
		{
			return allowedMem;
		}

		virtual void allowedMemory(size_t maxMem) override;
	};
}

