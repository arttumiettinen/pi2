#pragma once

#include "argumentdatatype.h"
#include "jobtype.h"
#include "delayed.h"
#include "io/inireader.h"
#include <string>
#include <vector>
#include <set>

namespace pilib
{
	class Command;
	class PISystem;

	void flushCache(const std::string& filename);
	
	constexpr int JOB_FAILED = -2;
	constexpr int JOB_WAITING = -1;

	/**
	Gets character describing progress of single task.
	*/
	char getProgressChar(int progress);

	/**
	Creates a string containing progress bar for multiple tasks.
	*/
	std::string createProgressBar(const std::vector<int>& progress);

	/**
	Shows progress bar in console.
	*/
	void showProgressBar(const std::string& bar, size_t& barLength);

	/**
	Base class for objects that are used to distribute commands to multiple processes.
	*/
	class Distributor
	{
	private:

		/**
		Indicates if submitted job scripts (pi code) should be printed to screen.
		*/
		bool showSubmittedScripts = false;

		/**
		Indicates if command execution can be delayed in order to combine execution of multiple (compatible) commands for reduced I/O and scratch disk space need.
		*/
		bool allowDelaying = false;


		/**
		Pointer to the PI system object.
		*/
		PISystem* piSystem;

		/**
		Path and filename of pi2 executable.
		*/
		std::string piCommand;
		
		/**
		List of commands whose execution has been delayed.
		*/
		std::vector<Delayed> delayedCommands;

		/**
		Output from last execution of delayed commands.
		*/
		std::vector<string> lastOutput;

		/**
		Chunk size for NN5 datasets used in distributed mode.
		*/
		Vec3c nn5ChunkSize;

		/**
		Determines suitable block size etc. for running commands in delayedCommands list.
		Throws exception if the commands cannot be run together.
		*/
		void determineDistributionConfiguration(Vec3c& margin, std::set<DistributedImageBase*>& inputImages, std::set<DistributedImageBase*>& outputImages, JobType& jobType, std::map<DistributedImageBase*, std::vector<std::tuple<Vec3c, Vec3c, Vec3c, Vec3c, Vec3c> > >& blocksPerImage, size_t& memoryReq);

		/**
		Runs commands that have been accumulated to the delayed command list.
		*/
		void runDelayedCommands();

		/**
		Test if the given command can be added to the delayed commands queue.
		Does not account for other commands.
		*/
		bool canDelay(const Delayed& delayed) const;

		/**
		Try to add the given command to the delayed commands queue.
		If the resulting queue cannot be run, returns false.
		*/
		bool tryDelay(Delayed& d);

	protected:

		Distributor(PISystem* piSystem);

		/**
		Returns shell command (including absolute path) that is used to run pi2 program.
		*/
		std::string getPiCommand() const
		{
			return piCommand;
		}

		/**
		Reads general distributor settings from ini file.
		Call this from derived class constructor when settings file has been found.
		*/
		void readSettings(INIReader& reader);


	public:

		virtual ~Distributor()
		{
			// Nothing to do here.
		}

		/**
		Block origin command parameter type and name.
		*/
		typedef Vec3c BLOCK_ORIGIN_ARG_TYPE;
		inline static const std::string BLOCK_ORIGIN_ARG_NAME = "block origin";

		/**
		Block index command parameter type and name.
		*/
		typedef coord_t BLOCK_INDEX_ARG_TYPE;
		inline static const std::string BLOCK_INDEX_ARG_NAME = "block index";

		/**
		Enables or disables delaying.
		*/
		void delaying(bool enable)
		{
			if (allowDelaying)
				flush();

			allowDelaying = enable;
		}

		/**
		Enables or disables printing of command scripts to console.
		*/
		void showScripts(bool enable)
		{
			showSubmittedScripts = enable;
		}

		/**
		Sets chunk size for NN5 datasets.
		*/
		void chunkSize(Vec3c chunkSize)
		{
			if (chunkSize.min() <= 0)
				throw ITLException("Chunk size cannot contain negative or zero elements.");
			nn5ChunkSize = chunkSize;
		}

		Vec3c getChunkSize() const
		{
			return nn5ChunkSize;
		}

		/**
		Process all delayed commands (if any).
		*/
		void flush();

		/**
		Gets pointer to PISystem object.
		*/
		PISystem* getSystem()
		{
			return piSystem;
		}

		/**
		Run the given command in the distributed framework.
		@param commands Commands to run on each block of the input image, and their arguments.
		@param args Command arguments. (Images must be of type DistributedImage)
		@return Output written by each subprocess.
		*/
		std::vector<std::string> distribute(const Command* command, std::vector<ParamVariant>& args, std::vector<ParamVariant>* argsForGetCorrespondingBlock = 0);



		/**
		Submits a job with the given pi2 code.
		This may be used for complex commands for which distribute(...) function is not suitable.
		See e.g. thickness map calculation.
		May block until the job is finished (the case with local processing).
		*/
		virtual void submitJob(const std::string& piCode, JobType jobType) = 0;

		/**
		Waits until all jobs have completed.
		Throws exception if any of the jobs fails or job output does not end in line "Everything done.".
		@return Output written by each job.
		*/
		virtual std::vector<string> waitForJobs() = 0;

		/**
		Returns the amount of memory in bytes that single process is allowed to use.
		*/
		virtual size_t allowedMemory() const = 0;

		/**
		Set the amount of allowed memory in bytes for single process.
		Set to 0 to determine the amount automatically.
		*/
		virtual void allowedMemory(size_t maxMem) = 0;
	};

}
