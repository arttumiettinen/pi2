
#include "command.h"
#include "distributable.h"
#include "trivialdistributable.h"
#include "pilibutilities.h"
#include "stringutils.h"

#include <sstream>
#include <regex>

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

	bool Command::canDistribute() const
	{
		return dynamic_cast<const Distributable*>(this) != 0;
	}

	string removeTags(const string& help)
	{
		string result = help;

		// `commandname` we leave as-is

		// $inline math$ -> remove $-signs
		result = std::regex_replace(result, std::regex("\\$(.*?)\\$"), "$1");

		// **text in bold** -> remove **
		result = std::regex_replace(result, std::regex("\\*\\*(.*?)\\*\\*"), "$1");
		
		return result;
	}

	string replaceTags(const string& help)
	{
		string result = help;

		// `commandname`
		result = std::regex_replace(result, std::regex("`(.*?)`"), ":ref:`$1`");

		// $inline math$
		result = std::regex_replace(result, std::regex("\\$(.*?)\\$"), ":math:`$1`");

		// No need to do anything for **text in bold**
		
		return result;
	}

	string processTags(const string& str, HelpFormat format)
	{
		switch (format)
		{
		case HelpFormat::Text: return removeTags(str);
		case HelpFormat::Rst: return replaceTags(str);
		default: throw logic_error("Unsupported help format.");
		}
	}


	string Command::helpString(HelpFormat format, const vector<string>* argumentDataTypeOverride, bool syntaxAsTitle) const
	{
		stringstream msg;

		// Command name
		switch (format)
		{
			case HelpFormat::Text:
				printTitle(msg, toSimpleString(), 1);
				break;
			case HelpFormat::Rst:
				if (!syntaxAsTitle)
				{
					msg << endl;
					msg << "**Syntax:** :code:`" << toSimpleString() << "`" << endl;
					msg << endl;
				}
				else
				{
					printTitle(msg, string(":code:`") + toSimpleString() + "`", 1);
				}
				break;
			default:
				throw logic_error("Unsupported help format.");
		}

		// Help text
		msg << processTags(help, format) << endl;

		// Distribution text
		if (dynamic_cast<const TrivialDistributable*>(this))
		{
			msg << endl << "This command can be used in a distributed processing mode, but it does not participate in distributed processing." << endl;
		}
		else if (dynamic_cast<const Distributable*>(this))
		{
			msg << endl << processTags("This command can be used in a distributed processing mode. Use 'distribute' command to change processing mode from local to distributed.", format) << endl;
		}
		else
		{
			msg << endl << "This command cannot be used in a distributed processing mode. If you need it, please contact the authors." << endl;
		}

		msg << endl;

		string empty = "";

		// Argument descriptions
		if (arguments.size() > 0)
		{
			printTitle(msg, "Arguments", 2);
			for (size_t n = 0; n < arguments.size(); n++)
			{
				if (argumentDataTypeOverride)
				{
					if(argumentDataTypeOverride->size() >= arguments.size())
						msg << processTags(arguments[n].helpString(format, &(*argumentDataTypeOverride)[n]), format);
					else
						msg << processTags(arguments[n].helpString(format, &empty), format);
				}
				else
				{
					msg << processTags(arguments[n].helpString(format), format);
				}
			}
		}

		if (seeAlso != "")
		{
			printTitle(msg, "See also", 2);
			
			vector<string> parts;
			switch (format)
			{
			case HelpFormat::Text:
				msg << seeAlso << endl;
				break;
			case HelpFormat::Rst:

				parts = itl2::split(seeAlso, false, ',', true);

				for(size_t n = 0; n < parts.size(); n++)
				{
					if (n > 0)
						msg << ", ";
					msg << ":ref:`" << parts[n] << "`";
				}
				msg << endl;
				
				break;
			default:
				throw logic_error("Unsupported help format.");
			}
		}


		return msg.str();
	}
}
