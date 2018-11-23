#pragma once

#include "argumentdatatype.h"
#include <string>
#include <vector>

using std::string;
using std::vector;


namespace pilib
{
	class Command;
	class PISystem;

	/**
	Base class for objects that are used to distribute commands to multiple processes.
	*/
	class Distributor
	{
	private:

		/**
		Pointer to the PI system object.
		*/
		PISystem* piSystem;

		/**
		Path and filename of pi2 executable.
		*/
		string piCommand;

		/**
		Submits a job with the given pi2 code.
		May block until the job is finished (the case with local processing).
		*/
		virtual void submitJob(const string& piCode) = 0;

		/**
		Waits until all jobs have completed.
		Throws exception if any of the jobs fails.
		@return Output written by each job.
		*/
		virtual vector<string> waitForJobs() = 0;

		/**
		Returns the amount of memory in bytes that single process is allowed to use.
		*/
		virtual size_t allowedMemory() const = 0;

	protected:

		Distributor(PISystem* piSystem);

		/**
		Returns shell command (including absolute path) that is used to run pi2 program.
		*/
		string getPiCommand() const
		{
			return piCommand;
		}


	public:

		typedef math::Vec3c BLOCK_ORIGIN_ARG_TYPE;
		static const string BLOCK_ORIGIN_ARG_NAME;

		typedef coord_t BLOCK_INDEX_ARG_TYPE;
		static const string BLOCK_INDEX_ARG_NAME;

		/**
		Run the given command in the distributed framework.
		@param commands Commands to run on each block of the input image, and their arguments.
		@param args Command arguments. (Images must be of type DistributedImage)
		@param distributionDirection Set to 0 to distribute in x-direction, 1 to distribute in y-direction, or 2 to distribute in z-direction (preferred).
		@param margin Distributed blocks must overlap this many pixels.
		@param outputFile File where the output of the command must be saved. This is used to override default behaviour of saving to temporary location in write* commands.
		@param refIndex Index of argument that is the reference image for determination of distribution block sizes. Set to maximum possible value to use the first output image or the first input image if there are no outputs as reference.
		@return Output written by each subprocess.
		*/
		vector<string> distribute(const Command* command, vector<ParamVariant>& args, size_t distributionDirection, const Vec3c& margin, const string* outputFile = 0, size_t refIndex = numeric_limits<size_t>::max(), vector<ParamVariant>* argsForGetCorrespondingBlock = 0);
	};

}
