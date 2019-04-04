#pragma once

#include <string>

using std::string;

namespace pilib
{

	/**
	Execute a program and return all output.
	@param showOutput Set to true to echo all output to cout.
	*/
	string execute(const string& cmd, bool showOutput = false);

	/**
	Execute the given program with given arguments and return all output that was generated.
	@param showOutput Set to true to echo all output to cout.
	*/
	string execute(const string& cmd, const string& args, bool showOutput = false);

	/**
	Gets path of this program.
	*/
	string getExecutablePath();

}


