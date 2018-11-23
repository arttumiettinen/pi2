#pragma once

#include <string>
#include <stdexcept>

using std::string;
using std::runtime_error;

namespace pi2
{
	/**
	Exception class for command line argument exceptions.
	*/
	class ArgumentException : public runtime_error
	{
	public:
		ArgumentException(const string& msg) : runtime_error(msg)
		{
		}
	};
}