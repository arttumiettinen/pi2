#pragma once

#include "command.h"
#include "distributable.h"

namespace pilib
{
	/**
	Identifies command to be trivially distributable, i.e. it runs the same way regardless of distributed processing state.
	Used for example for info and help commands.
	*/
	class TrivialDistributable : public Distributable, virtual public Command
	{
	public:

		TrivialDistributable()
		{

		}

		virtual void runDistributedInternal(PISystem* system, Distributor& distributor, vector<ParamVariant>& args) const
		{
			runInternal(system, args);
		}

		virtual void runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
		{
		}
	};
}
