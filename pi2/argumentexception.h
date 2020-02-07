#pragma once

#include <string>
#include <stdexcept>

namespace pi2
{
	/**
	Exception class for command line argument exceptions.
	*/
	class ArgumentException : public std::runtime_error
	{
	public:
		ArgumentException(const std::string& msg) : std::runtime_error(msg)
		{
		}
	};
}