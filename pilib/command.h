#pragma once

#include <string>
#include <vector>

#include "commandargument.h"
#include "helpformat.h"

using std::string;
using std::vector;

namespace pilib
{

	class PISystem;

	/**
	Get an item from the beginning of the arguments list and remove it from the list.
	*/
	template<typename out_t> out_t pop(std::vector<ParamVariant>& args)
	{
		if (args.size() <= 0)
			throw ITLException("Not enough arguments. (this is a bug: count of parameter descriptors and count of parameters requested do not match)");

		out_t result = std::get<out_t>(args[0]);
		args.erase(args.begin());
		return result;
	}

	/**
	Base class for all commands.
	*/
	class Command
	{
	private:
		/**
		Name of the command.
		*/
		std::string commandName;

		/**
		Help text for the command.
		*/
		std::string help;

		/**
		Commands that are related to this command.
		*/
		std::string seeAlso;

		/**
		Notes about usage of the command in specific UIs.
		*/
		std::string uiNotes;

		/**
		Command argument definitions.
		*/
		std::vector<CommandArgumentBase> arguments;

	protected:
		/**
		Constructor
		*/
		Command(const std::string& name, const std::string& help, const std::vector<CommandArgumentBase>& args = {}, const std::string& seeAlso = "", const std::string& uiNotes = "") :
			commandName(name),
			help(help),
			seeAlso(seeAlso),
			uiNotes(uiNotes)
		{
			arguments.insert(arguments.end(), args.begin(), args.end());
		}

	public:

		/**
		Destructor
		*/
		virtual ~Command()
		{

		}

		Command(Command const&) = delete;
		void operator=(Command const&) = delete;

		/**
		Gets the name of this command.
		*/
		const std::string& name() const
		{
			return commandName;
		}

		/**
		Gets argument definitions.
		*/
		const std::vector<CommandArgumentBase>& args() const
		{
			return arguments;
		}

		/**
		Make a string representation of this command.
		*/
		std::string toString() const;

		/**
		Make a simple string representation of this command.
		*/
		std::string toSimpleString() const;

		/**
		Construct help text for this command.
		@param argumentDataTypeOverride Pass nullptr to print real data types of arguments.
		Pass a pointer to an empty vector to not print argument data types at all.
		Pass a pointer to a vector containing string for each argument to override data type
		of each argument with the corresponding string.
		*/
		std::string helpString(HelpFormat format, const std::vector<std::string>* argumentDataTypeOverride = nullptr, bool syntaxAsTitle = false) const;

		/**
		Gets a value indicating whether this command can be used in distributed mode.
		*/
		bool canDistribute() const;

		/**
		Gets a value indicating whether this command is for internal use only.
		Internal commands should not be listed in help and command reference.
		Returns false by default.
		*/
		virtual bool isInternal() const
		{
			return false;
		}

		/**
		Run method that calls the pure run method or is overridden in special commands.
		There are two run methods so that the most used one is as simple as possible
		(as there will be many derived classes!)
		*/
		virtual void runInternal(PISystem* system, std::vector<ParamVariant>& args) const;

		/**
		Run this command
		*/
		virtual void run(std::vector<ParamVariant>& args) const = 0;

	};

}
