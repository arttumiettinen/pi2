#pragma once

#include <string>

using std::string;

namespace pilib
{

	/**
	Execute a program and return all output.
	*/
	string execute(const string& cmd);

	/**
	Execute the given program with given arguments and return all output that was generated.
	*/
	string execute(const string& cmd, const string& args);

	/**
	Gets path of this program.
	*/
	string getExecutablePath();

}


