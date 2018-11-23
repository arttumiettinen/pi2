
#include "exeutils.h"

#include <experimental/filesystem>

namespace fs = std::experimental::filesystem;

#if defined(__linux__)

	#include "pstream.h"

#elif defined(_WIN32)

	#define WIN32_LEAN_AND_MEAN
	#define NOMINMAX
	#include <Windows.h>

#else

	#error exeutils.cpp not configured for this platform.

#endif

#include "itlexception.h"

using itl2::ITLException;
using namespace std;

namespace pilib
{
#if defined(__linux__)

	/**
	Execute a program and return all output.
	*/
	string execute(const string& cmd)
	{
		redi::ipstream proc(cmd, redi::pstreams::pstdout | redi::pstreams::pstderr);
		string line;
		stringstream output;

        bool hasOut, hasErr;
        do
        {
            hasOut = false;
            hasErr = false;
        
            if(std::getline(proc.err(), line))
            {
                output << line << endl;
                hasErr = true;
            }
                
            proc.clear();
                
            if(std::getline(proc.out(), line))
            {
                output << line << endl;
                hasOut = true;
            }
                
            proc.clear();
        }
        while(hasOut || hasErr);

        /*
        while (std::getline(proc.err(), line))
			output << line << endl;
			
		proc.clear();

		while (std::getline(proc.out(), line))
			output << line << endl;
        */
        
		return output.str();
	}

	/**
	Gets path of this program.
	*/
	string getExecutablePath()
	{
		if (!fs::exists("/proc/self/exe"))
			throw ITLException("Unable to read executable path.");

		return fs::read_symlink("/proc/self/exe").parent_path().string() + "/";
	}


#elif defined(_WIN32)
	/**
	Execute a program and return all output.
	This code is mostly from
	https://stackoverflow.com/questions/478898/how-to-execute-a-command-and-get-output-of-command-within-c-using-posix
	*/
	string execute(const string& cmd)
	{
		HANDLE hPipeRead, hPipeWrite;

		SECURITY_ATTRIBUTES saAttr = { sizeof(SECURITY_ATTRIBUTES) };
		saAttr.bInheritHandle = TRUE;   //Pipe handles are inherited by child process.
		saAttr.lpSecurityDescriptor = NULL;

		// Create a pipe to get results from child's stdout.
		if (!CreatePipe(&hPipeRead, &hPipeWrite, &saAttr, 0))
			throw ITLException("Unable to create pipes.");

		STARTUPINFOA si = { sizeof(STARTUPINFOA) };
		si.dwFlags = STARTF_USESHOWWINDOW | STARTF_USESTDHANDLES;
		si.hStdOutput = hPipeWrite;
		si.hStdError = hPipeWrite;
		si.wShowWindow = SW_HIDE;       // Prevents cmd window from flashing. Requires STARTF_USESHOWWINDOW in dwFlags.

		PROCESS_INFORMATION pi = { 0 };

		BOOL fSuccess = CreateProcessA(NULL, (LPSTR)cmd.c_str(), NULL, NULL, TRUE, CREATE_NEW_CONSOLE, NULL, NULL, &si, &pi);
		if (!fSuccess)
		{
			CloseHandle(hPipeWrite);
			CloseHandle(hPipeRead);
			throw ITLException(string("Unable to create process to run ") + cmd);
		}

		bool bProcessEnded = false;
		string strResult;
		for (; !bProcessEnded;)
		{
			// Give some timeslice (50ms), so we won't waste 100% cpu.
			bProcessEnded = WaitForSingleObject(pi.hProcess, 50) == WAIT_OBJECT_0;

			// Even if process exited - we continue reading, if there is some data available over pipe.
			for (;;)
			{
				char buf[1024];
				DWORD dwRead = 0;
				DWORD dwAvail = 0;

				if (!::PeekNamedPipe(hPipeRead, NULL, 0, NULL, &dwAvail, NULL))
					break;

				if (!dwAvail) // no data available, return
					break;

				if (!::ReadFile(hPipeRead, buf, std::min((DWORD)(sizeof(buf) - 1), dwAvail), &dwRead, NULL) || !dwRead)
					// error, the child process might ended
					break;

				buf[dwRead] = 0;
				strResult += buf;
			}
		}

		CloseHandle(hPipeWrite);
		CloseHandle(hPipeRead);
		CloseHandle(pi.hProcess);
		CloseHandle(pi.hThread);
		return strResult;
	}

	/**
	Gets path of this program.
	This code is mostly from
	https://stackoverflow.com/questions/1528298/get-path-of-executable
	*/
	string getExecutablePath()
	{
		char ownPth[MAX_PATH];

		HMODULE hModule = GetModuleHandle(NULL);
		if (hModule != NULL)
		{
			GetModuleFileNameA(hModule, ownPth, (sizeof(ownPth)));
			return ownPth;
		}
		else
		{
			throw ITLException("Unable to get module handle.");
		}
	}

#else
#error LocalDistributor.h not configured for this platform.
#endif

	/**
	Execute the given program with given arguments and return all output that was generated.
	*/
	string execute(const string& cmd, const string& args)
	{
		return execute(cmd + " " + args);
	}
}
