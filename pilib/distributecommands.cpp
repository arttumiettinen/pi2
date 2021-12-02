
#include "distributecommands.h"
#include "commandmacros.h"
#include "pisystem.h"

namespace pilib
{
	void addDistributeCommands()
	{
		CommandList::add<SubmitJobCommand>();
		CommandList::add<WaitForJobsCommand>();
	}


	void SubmitJobCommand::runInternal(PISystem* system, std::vector<ParamVariant>& args) const
	{
		string script = pop<string>(args);
		JobType jobType = itl2::fromString<JobType>(pop<string>(args));

		system->run(script);
	}

	vector<string> SubmitJobCommand::runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
	{
		string script = pop<string>(args);
		JobType jobType = itl2::fromString<JobType>(pop<string>(args));
		distributor.flush();
		distributor.submitJob(script, jobType);
		return vector<string>();
	}

}