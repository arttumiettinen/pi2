#pragma once

#include "command.h"
#include "trivialdistributable.h"

namespace pilib
{
	/**
	Command that shows program version and information.
	This command is in its own small file as it changes in each git commit (as it shows version info
	derived from git commit info).
	*/
	class InfoCommand : virtual public Command, public TrivialDistributable
	{
	protected:
		friend class CommandList;

		InfoCommand() : Command("info", "Displays information about the computer and the PI system.",
			{},
			"help, license")
		{
		}

	public:
		virtual void run(vector<ParamVariant>& args) const override;
	};

}