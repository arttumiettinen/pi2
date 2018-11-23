#pragma once


#if defined(__linux__)

	#include <unistd.h>
	#include <sys/sysinfo.h>

	inline bool _isatty(int fileHandle)
	{
		return isatty(fileHandle);
	}

	inline bool _fileno(FILE* stream)
	{
		return fileno(stream);
	}

#elif defined(_WIN32)

	#include <stdio.h>
	#include <io.h>
	#define WIN32_LEAN_AND_MEAN
	#define NOMINMAX
	#include <Windows.h>

#else
	#error utilities.h not configured for this platform.
#endif


#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include "math/mathutils.h"
#include "math/vec3.h"
#include "datatypes.h"

using namespace std;

namespace itl2
{
	/**
	Waits given number of milliseconds.
	*/
	inline void sleep(unsigned int ms)
	{
#if defined(__linux__)
		usleep(ms);
#elif defined(_WIN32)
		Sleep(ms);
#else
	#error Sleep not configured for this platform.
#endif
	}

	/**
	Returns true if the standard output is a terminal.
	*/
	inline bool isTerminal()
	{
		return _isatty(_fileno(stdout));
	}

    /**
    Returns amount of usable system RAM in bytes.
    */
	inline size_t memorySize()
	{
#if defined(__linux__)

        struct sysinfo info;
        if(sysinfo(&info) != 0)
            throw ITLException("sysinfo call failed.");
        return info.totalram;

#elif defined(_WIN32)
		MEMORYSTATUSEX status;
		status.dwLength = sizeof(status);
		GlobalMemoryStatusEx(&status);
		return (size_t)status.ullTotalPhys;
#endif
	}

	/**
	Shows progress information in single-threaded process.
	Shows at most 100 steps of progress (i.e. skips printing if the result wouldn't change).
	Selects most suitable progress indicator for terminal and file output.
	@param n The zero-based index of last iteration that has been completed.
	@param max Number of iterations.
	*/
	inline void showProgress(size_t n, size_t max, bool show = true)
	{
		if (!show)
			return;

		if (n > 0)
		{
			if (isTerminal())
			{
				coord_t prevProgress = math::round((double)(n - 1) / (double)(max - 1) * 100);
				coord_t currProgress = math::round((double)(n) / (double)(max - 1) * 100);
				if (currProgress > prevProgress)
				{
					if(currProgress < 100)
						cout << currProgress << " %              \r" << flush;
					else
						cout << "                   \r" << flush;
				}
			}
			else
			{
				coord_t prevProgress = math::round((double)(n - 1) / (double)(max - 1) * 100 / 10);
				coord_t currProgress = math::round((double)(n) / (double)(max - 1) * 100 / 10);
				for(coord_t n = prevProgress; n < currProgress; n++)
				{
					cout << "=" << flush;
				}

				if (prevProgress < 10  && currProgress == 10)
				{
					cout << endl;
				}
			}
		}

		//if (max > 1)
		//{
		//	size_t dispStep = max / 100;
		//	if (dispStep <= 0)
		//		dispStep = 1;
		//	
		//	if (n % dispStep == 0)
		//		cout << ::round((double)n / (double)(max - 1) * 100) << " %              \r" << flush;
		//}
	}

	/**
	Showw progress information for multithreaded processes.
	Call from each thread once in each iteration.
	@param counter Reference to counter variable common to all threads.
	@param max Total number of iterations all threads are going to make.
	*/
	inline void showThreadProgress(size_t& counter, size_t max, bool show = true)
	{
		if (!show)
			return;

		#pragma omp critical(show_progress)
		{
			showProgress(counter, max);
			counter++;
		}
	}

	/**
	Returns size of given file.
	*/
	inline ifstream::pos_type filesize(const string& filename)
	{
		ifstream in(filename.c_str(), ifstream::ate | ifstream::binary);
		return in.tellg();
	}





	/**
	Tests if the given character is whitespace.
	*/
	inline bool isWhiteSpace(const char s)
	{
		return s == ' ' || s == '\t' || s == '\r' || s == '\n';
	}

	/**
	Removes whitespace from the beginning and the end of the string.
	*/
	inline void trim(string& s)
	{
		while (s.length() > 0 && isWhiteSpace(s[0]))
			s.erase(0, 1);

		while (s.length() > 0 && isWhiteSpace(s[s.length() - 1]))
			s.erase(s.length() - 1, 1);
	}



	/**
	Convert from string to type T.
	*/
	template<typename T>
	T fromString(const string& str)
	{
		istringstream stream(str);
		T t;
		stream >> t;
		if(stream.fail())
		{
			ostringstream errstr;
			errstr << "Invalid conversion from string '" << str << "'";
			throw ITLException(errstr.str());
		}
		return t;
	}

	/**
	Dummy string to string converter.
	*/
	template<> inline
	string fromString(const string& str)
	{
		return str;
	}

	/**
	Convert from string to double.
	Supports nan, inf, and -inf.
	*/
	template<>
	inline double fromString(const string& str)
	{
		if(str == "nan" || str == "NAN" || str == "NaN" ||
			str == "-nan" || str == "-NAN" || str == "-NaN")
			return numeric_limits<double>::quiet_NaN();

		if(str == "inf" || str == "INF" || str == "Inf" || str == "Infinity" || str == "infinity")
			return numeric_limits<double>::infinity();

		if(str == "-inf" || str == "-INF" || str == "-Inf" || str == "-Infinity" || str == "-infinity")
			return -numeric_limits<double>::infinity();

		istringstream stream(str);
		double t;
		stream >> t;
		if(stream.fail())
		{
			ostringstream errstr;
			errstr << "Invalid conversion from string '" << str << "'";
			throw ITLException(errstr.str());
		}
		return t;
	}

	/**
	Convert from string to float32_t.
	Supports nan, inf, and -inf.
	*/
	template<>
	inline float32_t fromString(const string& str)
	{
		return (float32_t)fromString<double>(str);
	}

	/**
	Converts string to bool.
	*/
	template<>
	inline bool fromString(const string& value)
	{
		if (value == "true" || value == "TRUE" || value == "True" || value == "1" || value == "on" || value == "ON" || value == "On")
		{
			return true;
		}

		if (value == "false" || value == "FALSE" || value == "False" || value == "0" || value == "off" || value == "OFF" || value == "Off")
		{
			return false;
		}

		ostringstream errstr;
		errstr << "Invalid conversion from string '" << value << "' to boolean.";
		throw ITLException(errstr.str());
	}

	/**
	Convert from string to vector
	*/
	template<>
	inline math::Vec3d fromString(const string& value)
	{
		try
		{
			if (value.length() < 1)
				throw ITLException("Empty string.");

			if (value[0] == '[')
			{
				// Vector notation [1, 2, 3]
				std::stringstream parts;
				parts << value;
				std::string sx, sy, sz;
				std::getline(parts, sx, '[');
				std::getline(parts, sx, ',');
				std::getline(parts, sy, ',');
				std::getline(parts, sz, ']');

				trim(sx);
				trim(sy);
				trim(sz);

				math::Vec3d result;
				result.x = itl2::fromString<double>(sx);
				result.y = itl2::fromString<double>(sy);
				result.z = itl2::fromString<double>(sz);
				return result;
			}
			else
			{
				// Single number
			
				double val = itl2::fromString<double>(value);

				math::Vec3d result;
				result.x = val;
				result.y = val;
				result.z = val;
				return result;
			}
		}
		catch (ITLException)
		{
			ostringstream errstr;
			errstr << "Value '" << value << "' is not a valid 3-component vector.";
			throw ITLException(errstr.str());
		}
	}

	/**
	Convert from string to vector.
	TODO: This is repeat from above!
	*/
	template<>
	inline math::Vec3<coord_t> fromString(const string& value)
	{
		try
		{
			if (value.length() < 1)
				throw ITLException("Empty string.");

			if (value[0] == '[')
			{
				// Vector notation [1, 2, 3]
				std::stringstream parts;
				parts << value;
				std::string sx, sy, sz;
				std::getline(parts, sx, '[');
				std::getline(parts, sx, ',');
				std::getline(parts, sy, ',');
				std::getline(parts, sz, ']');

				trim(sx);
				trim(sy);
				trim(sz);

				math::Vec3<coord_t> result;
				result.x = itl2::fromString<coord_t>(sx);
				result.y = itl2::fromString<coord_t>(sy);
				result.z = itl2::fromString<coord_t>(sz);
				return result;
			}
			else
			{
				// Single number

				coord_t val = itl2::fromString<coord_t>(value);

				math::Vec3<coord_t> result;
				result.x = val;
				result.y = val;
				result.z = val;
				return result;
			}
		}
		catch (ITLException)
		{
			ostringstream errstr;
			errstr << "Value '" << value << "' is not a valid 3-component integer vector.";
			throw ITLException(errstr.str());
		}
	}

	/**
	Converts strings in the given list to type T.
	*/
	template<typename T> vector<T> fromString(const vector<string>& values)
	{
		vector<T> result;
		result.reserve(values.size());

		for (size_t i = 0; i < values.size(); i++)
		{
			result.push_back(fromString<T>(values[i]));
		}

		return result;
	}

	/**
	Converts to string from type T.
	*/
	template<typename T>
	string toString(const T& x)
	{
		stringstream s;
		s << x;
		return s.str();
	}

	template<>
	inline string toString(const double& x)
	{
		if (math::isnan(x))
			return "nan";
		if (math::isinf(x) && x > 0)
			return "inf";
		if (math::isinf(x) && x < 0)
			return "-inf";

		stringstream s;
		s << x;
		return s.str();
	}

	template<>
	inline string toString(const float& x)
	{
		if (math::isnan(x))
			return "nan";
		if (math::isinf(x) && x > 0)
			return "inf";
		if (math::isinf(x) && x < 0)
			return "-inf";

		stringstream s;
		s << x;
		return s.str();
	}

	template<>
	inline string toString(const bool& x)
	{
		if (x)
			return "True";
		return "False";
	}
}
