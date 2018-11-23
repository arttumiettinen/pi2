#pragma once

#include <string>
#include <vector>

#include "commandargument.h"

using std::string;
using std::vector;

namespace pilib
{

	class PISystem;

	/**
	Get an item from the beginning of the arguments list and remove it from the list.
	*/
	template<typename out_t> out_t pop(vector<ParamVariant>& args)
	{
		if (args.size() <= 0)
			throw ITLException("Not enough arguments. (this is a bug: count of parameter descriptors and count of parameters requested do not match.");

		out_t result = get<out_t>(args[0]);
		args.erase(args.begin());
		return result;
	}

	/**
	Base class for commands.
	*/
	class Command
	{
	private:
		/**
		Name of the command.
		*/
		string commandName;

		/**
		Help text for the command.
		*/
		string help;

		/**
		Command argument definitions.
		*/
		vector<CommandArgumentBase> arguments;

	public:

		/**
		Constructor
		*/
		Command(const string& name, const string& help, const vector<CommandArgumentBase>& args = {}) :
			commandName(name),
			help(help)
		{
			arguments.insert(arguments.end(), args.begin(), args.end());
		}

		/**
		Destructor
		*/
		virtual ~Command()
		{

		}

		/**
		Gets the name of this command.
		*/
		const string& name() const
		{
			return commandName;
		}

		/**
		Gets argument definitions.
		*/
		const vector<CommandArgumentBase>& args() const
		{
			return arguments;
		}

		/**
		Make a string representation of this command.
		*/
		string toString() const;

		/**
		Make a simple string representation of this command.
		*/
		string toSimpleString() const;

		/**
		Construct help text for this command.
		*/
		string helpString() const;

		/**
		Gets a value indicating whether this command can be used in distributed mode.
		*/
		bool canDistribute() const;

		/**
		Run method that calls the pure run method or is overridden in special commands.
		There are two run methods so that the most used one is as simple as possible
		(as there will be many derived classes!)
		*/
		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const;

		/**
		Run this command
		*/
		virtual void run(vector<ParamVariant>& args) const = 0;

	};

}
