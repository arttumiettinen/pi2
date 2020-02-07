#pragma once

#include <string>
#include <sstream>

#include "argumentdatatype.h"
#include "helpformat.h"

using itl2::Image;

namespace pilib
{
	
	class CommandArgumentBase
	{
	private:

		ArgumentDataType type;

		ParameterDirection dir;

		std::string def;

		bool defAllowed;

		std::string cname;

		std::string help;

	protected:
		CommandArgumentBase(ParameterDirection direction, ArgumentDataType type, const std::string& name, const std::string& help, const std::string& defaultValue, bool defaultValueAllowed) :
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
		//CommandArgumentBase(ParameterDirection direction, ArgumentDataType type, std::string name, std::string help, std::string defaultValue) :
		//	CommandArgumentBase(direction, type, name, help, defaultValue, true)
		//{
		//}

		///**
		//Constructor for arguments without default value.
		//*/
		//CommandArgumentBase(ParameterDirection direction, ArgumentDataType type, std::string name, std::string help) :
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

		const std::string& defaultValue() const
		{
			return def;
		}

		bool defaultAllowed() const
		{
			return defAllowed;
		}

		const std::string& name() const
		{
			return cname;
		}

		/**
		Gets default value, quoted if it is empty or contains spaces.
		Returns empty string if there is no default value allowed.
		*/
		std::string defaultValueQuoted() const;

		/*
		Makes a string representation of this argument.
		*/
		std::string toString() const;

		/*
		Makes simple string representation of this argument.
		*/
		std::string toSimpleString() const;

		/*
		Makes a string representation of the argument and its default value.
		*/
		std::string dataTypeHelpString(const std::string* dataTypeOverride = nullptr) const;

		/*
		Makes a string representation of the argument but not its default value.
		*/
		std::string dataTypeHelpStringNoDefault(const std::string* dataTypeOverride = nullptr) const;

		/*
		Creates string containing help for usage of this parameter.
		@param dataTypeOverride Pointer to string that contains the string to be printed as data type for this argument.
		Pass nullptr to print the real data type.
		If data type override length is zero, no data type will be shown.
		*/
		std::string helpString(HelpFormat format, const std::string* dataTypeOverride = nullptr) const;
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
		CommandArgument(ParameterDirection direction, const std::string& name, const std::string& help, param_t defaultValue) :
			CommandArgumentBase(direction, parameterType<param_t>(), name, help, itl2::toString(defaultValue), true)
		{
		}

		/**
		Constructor for arguments without default value.
		*/
		CommandArgument(ParameterDirection direction, const std::string& name, const std::string& help) :
			CommandArgumentBase(direction, parameterType<param_t>(), name, help, "", false)
		{
		}
	};

}
