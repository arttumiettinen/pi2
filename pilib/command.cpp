
#include "command.h"
#include "distributable.h"
#include "trivialdistributable.h"

#include <sstream>

using namespace std;

namespace pilib
{
	void Command::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		run(args);
	}

	string Command::toString() const
	{
		stringstream msg;
		msg << name() << "(";
		const vector<CommandArgumentBase>& candArgs = args();
		for (size_t m = 0; m < candArgs.size(); m++)
		{
			msg << candArgs[m].toString();

			if (m < candArgs.size() - 1)
				msg << ", ";
		}
		msg << ")";
		return msg.str();
	}

	string Command::toSimpleString() const
	{
		stringstream msg;
		msg << name() << "(";
		const vector<CommandArgumentBase>& candArgs = args();
		for (size_t m = 0; m < candArgs.size(); m++)
		{
			msg << candArgs[m].toSimpleString();

			if (m < candArgs.size() - 1)
				msg << ", ";
		}
		msg << ")";
		return msg.str();
	}

	void printTitle(stringstream& msg, const string& str, int level)
	{
		char uc = '-';
		if(level <= 1)
			uc = '=';
			
		msg << str << endl;
		for (size_t n = 0; n < str.length(); n++)
			msg << uc;
		msg << endl << endl;
	}

	bool Command::canDistribute() const
	{
		return dynamic_cast<const Distributable*>(this) != 0;
	}

	string Command::helpString() const
	{
		stringstream msg;

		// Command name and help
		printTitle(msg, toSimpleString(), 1);
		msg << help << endl;
		msg << endl;

		// Argument descriptions
		if (arguments.size() > 0)
		{
			printTitle(msg, "Arguments:", 2);
			for (size_t n = 0; n < arguments.size(); n++)
			{
				msg << arguments[n].helpString();
			}
		}

		// Distributed processing information
		printTitle(msg, "Distributed processing:", 2);
		if (dynamic_cast<const TrivialDistributable*>(this))
		{
			msg << "This command can be used in a distributed processing environment, but it does not participate in distributed processing." << endl;
		}
		else if (dynamic_cast<const Distributable*>(this))
		{
			msg << "This command can be used in a distributed processing environment. Use 'distribute' command to change processing mode from local to distributed." << endl;
		}
		else
		{
			msg << "This command cannot be used in a distributed processing environment." << endl;
		}
		msg << endl;


		return msg.str();
	}
}
