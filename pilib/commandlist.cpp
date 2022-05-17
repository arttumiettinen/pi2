
#include "commandlist.h"
#include "pilibutilities.h"

#include <tuple>

using namespace itl2;
using namespace std;

namespace pilib
{
	/**
	Creates commands and adds them to the CommandList.
	These are declared in .cpp files corresponding to command header files.
	*/
	void addFilterCommands();
	void addSpecialCommands();
	void addIOCommands();
	void addOtherCommands();
	void addPointProcessCommands();
	void addThinAndSkeletonCommands();
	void addSpecialSystemCommands();
	void addTransformCommands();
	void addParticlesCommands();
	void addGenerationCommands();
	void addTomoCommands();
	void addProjectionCommands();
	void addConvertCommands();
	void addHistogramCommands();
	void addDmapCommands();
	void addTmapCommands();
	void addMaximaCommands();
	void addCarpetCommands();
	void addStructureCommands();
	void addAutoThresholdCommands();
	void addInfoCommands();
	void addPathCommands();
	void addMetadataCommands();
	void addEvalCommands();
	void addDistributeCommands();
	void addInpaintCommands();
	void addFloodfillCommands();

	vector<unique_ptr<Command> > CommandList::commands;

	CommandList::constructor CommandList::cons;

	CommandList::constructor::constructor()
	{
		addFilterCommands();
		addSpecialCommands();
		addIOCommands();
		addOtherCommands();
		addFloodfillCommands();
		addPointProcessCommands();
		addThinAndSkeletonCommands();
		addSpecialSystemCommands();
		addTransformCommands();
		addParticlesCommands();
		addGenerationCommands();
		addTomoCommands();
		addProjectionCommands();
		addConvertCommands();
		addHistogramCommands();
		addDmapCommands();
		addTmapCommands();
		addMaximaCommands();
		addCarpetCommands();
		addStructureCommands();
		addAutoThresholdCommands();
		addInfoCommands();
		addPathCommands();
		addMetadataCommands();
		addEvalCommands();
		addDistributeCommands();
		addInpaintCommands();
	}


	string getCombinedHelp(HelpFormat format, const vector<Command*> &cmds, size_t start, size_t end)
	{
		// For each group, for each parameter, get data type, convert to string
		// Get help for first command in each group, override command data type strings.

		if (end - start <= 0)
			return "";

		// Show default values with data types only if mode is Rst and there are multiple commands and different default values
		vector<bool> showDefaultsWithDataTypes;
		for (size_t n = start; n < end; n++)
		{
			while (showDefaultsWithDataTypes.size() < cmds[n]->args().size())
				showDefaultsWithDataTypes.push_back(format == HelpFormat::Text);
		}

		if (format == HelpFormat::Rst)
		{
			vector<string> defaults;
			for (size_t n = start; n < end; n++)
			{
				while (defaults.size() < cmds[n]->args().size())
					defaults.push_back("");

				for (size_t m = 0; m < cmds[n]->args().size(); m++)
				{
					bool defAllowed = cmds[n]->args()[m].defaultAllowed();
					if (defAllowed)
					{
						string def = cmds[n]->args()[m].defaultValue();
						if (n <= start)
						{
							defaults[m] = def;
						}
						else if (def != defaults[m])
						{
							defaults[m] = def;
							showDefaultsWithDataTypes[m] = true;
						}
					}
				}
			}
		}

		// For each parameter, for each command, get data type and convert to string
		vector<string> parameterDataTypes;
		for (size_t n = start; n < end; n++)
		{
			while (parameterDataTypes.size() < cmds[n]->args().size())
				parameterDataTypes.push_back("");

			for (size_t m = 0; m < cmds[n]->args().size(); m++)
			{
				string dt;

				if(showDefaultsWithDataTypes[m])
					dt = cmds[n]->args()[m].dataTypeHelpString();
				else
					dt = cmds[n]->args()[m].dataTypeHelpStringNoDefault();

				if (!contains(parameterDataTypes[m], string(", ") + dt))
					parameterDataTypes[m] += ", " + dt;
			}
		}

		for (size_t n = 0; n < parameterDataTypes.size(); n++)
			parameterDataTypes[n] = parameterDataTypes[n].substr(2);

		return cmds[start]->helpString(format, &parameterDataTypes, format == HelpFormat::Rst && end - start < cmds.size());
	}

	
	vector<string> CommandList::help(const string& cmdName, HelpFormat format)
	{
		vector<Command*> cmds = CommandList::byName(cmdName);

		// Get help string of each command without data types.
		vector<string> resultsNoDT;
		vector<string> emptyVector;
		for (size_t n = 0; n < cmds.size(); n++)
		{
			resultsNoDT.push_back(cmds[n]->helpString(format, &emptyVector));
		}

		vector<string> results;

		// Group commands that have the same help string without data types.
		for (size_t n = 0; n < cmds.size();)
		{
			size_t m;
			for (m = n + 1; m < cmds.size(); m++)
			{
				if (resultsNoDT[n] != resultsNoDT[m])
					break;
			}

			// Commands [n..m-1] have the same help text without command data types.
			// Generate help with data types.
			results.push_back(getCombinedHelp(format, cmds, n, m));

			n = m;
		}

		// This is the simple version that does not group commands based on data type.
		//vector<string> results;
		//for (size_t n = 0; n < cmds.size(); n++)
		//{
		//	results.push_back(cmds[n]->helpString());
		//}

		return results;
	}

	vector<string> CommandList::apropos(const string& str)
	{
		vector<string> names;

		const auto& cmds = CommandList::all();

		string loStr = str;
		toLower(loStr);

		for (size_t n = 0; n < cmds.size(); n++)
		{
			if (!cmds[n]->isInternal())
			{
				string loName = cmds[n]->name();
				toLower(loName);
				string loHelp = cmds[n]->helpString(HelpFormat::Text);
				toLower(loHelp);
				if (contains(loName, loStr) || contains(loHelp, loStr))
				{
					names.push_back(cmds[n]->toSimpleString());
				}
			}
		}

		removeDuplicates(names);
		return names;
	}

	string CommandList::list(bool forUser)
	{
		vector<tuple<string, bool> > names;
		const vector<unique_ptr<Command> >& commands = CommandList::all();
		for (size_t n = 0; n < commands.size(); n++)
		{
			if(!commands[n]->isInternal())
				names.push_back(make_tuple(commands[n]->toSimpleString(), commands[n]->canDistribute()));
		}
		removeDuplicates(names);
		sort(names.begin(), names.end());

		stringstream msg;
		if (forUser)
			msg << "Letter D placed before command name indicates a command that can be used in distributed computing mode. See help(distribute) for more information." << endl << endl;
		for (size_t n = 0; n < names.size(); n++)
		{
			if (forUser)
				msg << (std::get<1>(names[n]) ? "D " : "  ");
			msg << std::get<0>(names[n]) << endl;
		}
		return msg.str();
	}

}