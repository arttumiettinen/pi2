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
	protected:
		friend class CommandList;

		TrivialDistributable()
		{

		}

	public:
		using Distributable::runDistributed;

		virtual std::vector<std::string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const override
		{
			runInternal(distributor.getSystem(), args);
			return std::vector<std::string>();
		}
	};
}
