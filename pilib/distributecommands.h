#pragma once

#include <string>
#include "command.h"
#include "distributable.h"

namespace pilib
{
	inline std::string submitSeeAlso()
	{
		return "submitjob, waitforjobs";
	}

	class SubmitJobCommand : public Command, public Distributable
	{
	protected:
		friend class CommandList;

		SubmitJobCommand() : Command("submitjob", "Submits a pi2 script job to the active cluster system, or runs it locally if distributed processing is not active. Does not wait for job completion. Typically you would not use this method of job submission, but just call `distribute` function and let the system deal with submitting jobs.",
			{
				CommandArgument<string>(ParameterDirection::In, "script", "Pi2 script to run."),
				CommandArgument<string>(ParameterDirection::In, "job type", "Job type, either fast, normal, or slow.", "normal"),
			},
			submitSeeAlso())
		{
		}

	public:
		virtual void runInternal(PISystem* system, std::vector<ParamVariant>& args) const override;

		virtual void run(std::vector<ParamVariant>& args) const override
		{
			// Not used - runInternal does all the work.
		}

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override;
	};


	class WaitForJobsCommand : public Command, public Distributable
	{
	protected:
		friend class CommandList;

		WaitForJobsCommand() : Command("waitforjobs", "Waits until all running cluster jobs are complete. Use this command to wait for jobs submitted using the submitjob command.",
			{
			},
			submitSeeAlso())
		{
		}

	public:
		virtual void run(vector<ParamVariant>& args) const override
		{

		}

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override
		{
			distributor.waitForJobs();
			return vector<string>();
		}
	};
}