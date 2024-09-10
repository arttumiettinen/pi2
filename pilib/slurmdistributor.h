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

		std::string getErrorMessage(size_t jobIndex) const;

		/**
		Stores (job slurm id, queue type, submission count) of all submitted jobs since last call to waitForJobs.
		*/
		std::vector<std::tuple<size_t, JobType, size_t> > submittedJobs;

		/**
		Extra arguments for sbatch and sinfo, for fast jobs
		*/
		std::string extraArgsFastJobsSBatch, extraArgsFastJobsSInfo;

		/**
		Extra arguments for sbatch and sinfo, for normal jobs
		*/
		std::string extraArgsNormalJobsSBatch, extraArgsNormalJobsSInfo;

		/**
		Extra arguments for sbatch and sinfo, for slow jobs
		*/
		std::string extraArgsSlowJobsSBatch, extraArgsSlowJobsSInfo;

		/**
		Commands used to run sbatch, squeue, scancel and sinfo.
		*/
		std::string sbatchCommand, squeueCommand, scancelCommand, sinfoCommand;

		/**
		Interval for job progress polls. Give in milliseconds.
		*/
		int progressPollInterval;

		/**
		Returns suitable sbatch arguments given type of job.
		*/
		std::string extraArgsSBatch(JobType jobType) const
		{
			switch (jobType)
			{
			case JobType::Fast: return extraArgsFastJobsSBatch;
			case JobType::Normal: return extraArgsNormalJobsSBatch;
			case JobType::Slow: return extraArgsSlowJobsSBatch;
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
		Cancels job with given SLURM id.
		*/
		void cancelJob(size_t slurmId) const;

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
		Gets SLURM error log of given job.
		*/
		std::string getSlurmErrorLog(size_t jobIndex, bool flush = false) const;

		/**
		Gets job progress from log file.
		*/
		int getJobProgressFromLog(size_t jobIndex) const;

		/**
		Retrieves job queueing time [s] and execution time [s] from SLURM.
		*/
		std::tuple<double, double> getJobTimes(size_t jobIndex) const;

		/**
		Submits job with given index again.
		*/
		void resubmit(size_t jobIndex);

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

		/**
		Creates sbatch file name.
		*/
		std::string makeSbatchName(size_t jobIndex) const;

	public:
		SLURMDistributor(PISystem* system);

		virtual void submitJob(const std::string& piCode, JobType jobType) override;

		virtual std::vector<std::string> waitForJobs() override;

		virtual size_t allowedMemory() const override
		{
			return allowedMem;
		}

		virtual void allowedMemory(size_t maxMem) override;
	};
}
