
#include <iostream>
#include <fstream>
#include <memory>

#include "pilib.h"

#include "argumentexception.h"

#if defined(__linux__)
#include <signal.h>
#include <cstring>
#endif

using namespace std;
using namespace pi2;

/*
Read contents of text file and return true if succesful.
*/
bool readFile(const string& filename, string& code)
{
	ifstream in;
	in.open(filename);
	if (in.good())
	{
		while (in.good())
		{
			string s;
			getline(in, s);
			code += s + "\n";
		}
		return true;
	}
	return false;
}

#if defined(__linux__)

void sigbushandler(int signo, siginfo_t* si, void* data) {

	if (signo == SIGBUS)
	{
		cout << endl << endl;
		cout << "Received sigbus." << endl;
		cout << "signo = " << signo << endl;
		cout << "si->si_signo = " << si->si_signo << endl;
		cout << "si->si_code = " << si->si_code << endl;
		cout << "si->si_errno = " << si->si_errno << endl;
		exit(0);
	}
}
#endif

/*
Main entry point
*/
int main(int argc, char** argv)
{

#if defined(__linux__)
	struct sigaction sa;
	memset(&sa, 0, sizeof(sa));
	sa.sa_flags = SA_SIGINFO;
	sa.sa_sigaction = sigbushandler;
	if (sigaction(SIGBUS, &sa, 0) == -1)
		cout << "Unable to catch SIGBUS, " << strerror(errno) << endl;
#endif
	
	auto handle = unique_ptr<void, decltype(destroyPI)*>(createPI(), destroyPI);
	try
	{
		string code;
		bool file = false;
		
		if (argc == 1)
		{
			// Default input file
			if (!readFile("pi.txt", code))
				if (!readFile("pi2.txt", code))
					throw ArgumentException("No input supplied.");

			file = true;
		}
		else if (argc == 2)
		{
			// Input file or single command

			string fname = argv[1];
			file = true;

			if (!readFile(fname, code))
			{
				file = false;
				code = fname;
			}
		}
		else
		{
			// Arguments as code

			for (size_t n = 1; n < argc; n++)
				code += string(argv[n]) + " ";
		}

		run(handle.get(), "echo(true, false)");

		if (!run(handle.get(), code.c_str()))
		{
			if (file)
				cout << "Error(line " << lastErrorLine(handle.get()) << "): ";
			else
				cout << "Error: ";
			cout << lastErrorMessage(handle.get()) << endl;
			clearLastError(handle.get());
			return 1;
		}

		return 0;
	}
	catch (ArgumentException& e)
	{
		cout << e.what() << endl;
		//cout << "Invalid command line arguments." << endl;
		cout << endl;
		cout << "Usage:" << endl;
		cout << "pi2 script_file_name" << endl;
		cout << "pi2 command; command; command;" << endl;
		cout << endl;
		cout << "If no file name or commands are specified and the current directory contains file named pi.txt, commands in that file are run." << endl;
		cout << endl;
		cout << "Things to try:" << endl;
		cout << "pi2 info" << endl;
		cout << "pi2 help" << endl;

		return 1;
	}
	catch (std::bad_alloc& e)
	{
		cout << "Error: Out of memory (" << e.what() << ")" << endl;
		
		return 2;
	}
	catch (std::exception& e)
	{
		cout << "Error: " << e.what() << endl;

		return 3;
	}
	catch (...)
	{
		cout << "Error: Unknown error" << endl;

		return 4;
	}
}
