#pragma once

#include "command.h"

namespace pilib
{

	class CommandList
	{
	private:
		/**
		List of available commands.
		*/
		static std::vector<std::unique_ptr<Command> > commands;

		/**
		Static constructor hack for filling commands list automatically.
		*/
		friend class constructor;
		struct constructor {
			constructor();
		};
		static constructor cons;

	public:

		/**
		Gets list of all commands, including internal commands.
		*/
		static const std::vector<std::unique_ptr<Command> >& all()
		{
			return commands;
		}

		/**
		Add a command to the system.
		*/
		template<class command_t> static void add()
		{
			// make_unique cannot be used as we are a friend class of the command class and calling a private constructor.
			//commands.push_back(make_unique<command_t>());
			commands.push_back(std::unique_ptr<command_t>(new command_t()));
		}

		/**
		Get all commands whose name is the given one.
		*/
		static std::vector<Command*> byName(const string& name)
		{
			std::vector<Command*> cmds;
			for (size_t n = 0; n < commands.size(); n++)
			{
				if (commands[n]->name() == name)
					cmds.push_back(commands[n].get());
			}
			return cmds;
		}

		/**
		Finds command based on its type.
		*/
		template<class command_t> static command_t& get()
		{
			for (std::unique_ptr<Command>& cmd : commands)
			{
				command_t* c = dynamic_cast<command_t*>(cmd.get());
				if (c)
					return *c;
			}

			throw std::logic_error(string("Unknown command type: ") + typeid(command_t).name());
		}

		/**
		Gets help topics for given command name.
		*/
		static std::vector<string> help(const string& cmdName, HelpFormat format);

		/**
		Gets command names that are related to the given string.
		*/
		static std::vector<string> apropos(const string& str);

		/**
		Get a list of commands and some basic information about each one, each command separated by newline.
		Does not include internal commands in the returned list.
		*/
		static string list(bool forUser);
	};

}
