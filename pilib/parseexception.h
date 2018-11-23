#pragma once

#include <string>

using std::string;

namespace pilib
{
	/*
	Parse exception
	*/
	class ParseException : public ITLException
	{
	public:

		/*
		Default constructor
		*/
		ParseException() : ParseException("")
		{

		}

		/*
		Constructor
		*/
		ParseException(const string& msg) :
			ITLException(msg)
		{
		}
	};
}


