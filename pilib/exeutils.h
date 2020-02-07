#pragma once

#include <string>

namespace pilib
{

	/**
	Execute a program and return all output.
	@param showOutput Set to true to echo all output to cout.
	*/
	std::string execute(const std::string& cmd, bool showOutput = false);

	/**
	Execute the given program with given arguments and return all output that was generated.
	@param showOutput Set to true to echo all output to cout.
	*/
	std::string execute(const std::string& cmd, const std::string& args, bool showOutput = false);

	/**
	Gets path of this program.
	*/
	std::string getExecutablePath();

}


