#pragma once

#include "command.h"
#include "commandsbase.h"
#include "logger.h"


namespace pilib
{

	inline std::string loggingSeeAlso()
	{
		return "-, -";
	}


	inline std::string loggingMethodsHelp()
	{
		return
			"**logging**\n"
			"\n"
			"\n"
			"\n";
	}


	template<typename pixel_t> class setLoggingCommand : public OneImageCommand<pixel_t>
	{
	protected:
		friend class CommandList;

		setLoggingCommand() : OneImageCommand<pixel_t>("setLogging",
			"Sets logging true/false.",
			{
				CommandArgument<bool>(ParameterDirection::In, "loggingValue", "Sets logging to this value true/false.")

			})
		{
		}

	public:
		virtual void run(const Image<pixel_t>& img, std::vector<ParamVariant>& args) const override
		{
			bool loggingValue = pop<bool>(args);
			setLogging(loggingValue);
		}
	};


}