#pragma once

#include <string>
#include <sstream>

#include "argumentdatatype.h"

using namespace std;
using itl2::Image;

namespace pilib
{
	
	class CommandArgumentBase
	{
	private:

		ArgumentDataType type;

		ParameterDirection dir;

		string def;

		bool defAllowed;

		string cname;

		string help;

	protected:
		CommandArgumentBase(ParameterDirection direction, ArgumentDataType type, string name, string help, string defaultValue, bool defaultValueAllowed) :
			type(type),
			dir(direction),
			def(defaultValue),
			defAllowed(defaultValueAllowed),
			cname(name),
			help(help)
		{
		}

	public:



		///**
		//Constructor for arguments with default value.
		//*/
		//CommandArgumentBase(ParameterDirection direction, ArgumentDataType type, string name, string help, string defaultValue) :
		//	CommandArgumentBase(direction, type, name, help, defaultValue, true)
		//{
		//}

		///**
		//Constructor for arguments without default value.
		//*/
		//CommandArgumentBase(ParameterDirection direction, ArgumentDataType type, string name, string help) :
		//	CommandArgumentBase(direction, type, name, help, "", false)
		//{
		//}

		ArgumentDataType dataType() const
		{
			return type;
		}

		ParameterDirection direction() const
		{
			return dir;
		}

		const string& defaultValue() const
		{
			return def;
		}

		bool defaultAllowed() const
		{
			return defAllowed;
		}

		const string& name() const
		{
			return cname;
		}

		/*
		Makes a string representation of this argument.
		*/
		string toString() const
		{
			stringstream msg;
			msg << "[" << pilib::toString(direction()) << ", " << pilib::toString(dataType()) << "] " << name();
			if (defaultAllowed())
			{
				msg << " = \"" << defaultValue() << "\"";
			}
			return msg.str();
		}

		/*
		Makes simple string representation of this argument.
		*/
		string toSimpleString() const
		{
			stringstream msg;
			msg << name();
			return msg.str();
		}

		/*
		Creates string containing help for usage of this parameter.
		*/
		string helpString() const
		{
			stringstream msg;
			msg << name() << " [";
			if (direction() == ParameterDirection::In)
				msg << "input";
			else if (direction() == ParameterDirection::Out)
				msg << "output";
			else
				msg << "input & output";
			msg << ", " << pilib::toString(dataType());

			if (direction() == ParameterDirection::In && defAllowed)
				msg << ", default value = " << defaultValue();

			msg << "]";
			msg << ": ";
			msg << help << endl;
			msg << endl;
			return msg.str();
		}
	};

	/*
	Encapsulates information about command argument.
	*/
	template<typename param_t> class CommandArgument : public CommandArgumentBase
	{
	public:
		/**
		Constructor for arguments with default value.
		*/
		CommandArgument(ParameterDirection direction, string name, string help, param_t defaultValue) :
			CommandArgumentBase(direction, parameterType<param_t>(), name, help, itl2::toString(defaultValue), true)
		{
		}

		/**
		Constructor for arguments without default value.
		*/
		CommandArgument(ParameterDirection direction, string name, string help) :
			CommandArgumentBase(direction, parameterType<param_t>(), name, help, "", false)
		{
		}
	};

}
