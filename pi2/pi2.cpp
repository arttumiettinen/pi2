
#include <iostream>
#include <fstream>

#include "pilib.h"

#include "argumentexception.h"

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

/*
Main entry point
*/
int main(int argc, char** argv)
{
	
	void* handle = createPI();
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

		run(handle, "echo(true, false)");

		if (!run(handle, code.c_str()))
		{
			if (file)
				cout << "Error(line " << lastErrorLine(handle) << "): ";
			else
				cout << "Error: ";
			cout << lastErrorMessage(handle) << endl;
			clearLastError(handle);
		}

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
	}
	catch (std::bad_alloc& e)
	{
		cout << "Error: Out of memory (" << e.what() << ")" << endl;
	}
	catch (std::exception& e)
	{
		cout << "Error: " << e.what() << endl;
	}
	catch (...)
	{
		cout << "Error: Unknown error" << endl;
	}
	destroyPI(handle);
}
