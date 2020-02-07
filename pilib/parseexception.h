#pragma once

#include <string>

namespace pilib
{
	/**
	Parse exception
	*/
	class ParseException : public ITLException
	{
	public:

		/**
		Default constructor
		*/
		ParseException() : ParseException("")
		{
		}

		/**
		Constructor
		*/
		ParseException(const std::string& msg) :
			ITLException(msg)
		{
		}
	};
}


