#pragma once


#if defined(__linux__) || defined(__APPLE__)

	#include <unistd.h>
	#include <stdio.h>

	#if defined(__linux__)	
		#include <sys/sysinfo.h>
	#else
		#include <sys/types.h>
		#include <sys/sysctl.h>
	#endif

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
#include <iomanip>
#include <list>
#include <algorithm>
#include "math/mathutils.h"
#include "math/vec3.h"
#include "math/vec2.h"
#include "datatypes.h"
#include "stringutils.h"

namespace itl2
{
	/**
	Waits given number of milliseconds.
	*/
	inline void sleep(unsigned int ms)
	{
#if defined(__linux__)  || defined(__APPLE__)
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
            throw ITLException("Sysinfo call failed.");
        return info.totalram;
#elif defined(__APPLE__)
		int mib[] = { CTL_HW, HW_MEMSIZE };
		int64_t value = 0;
		size_t length = sizeof(value);

		if (-1 == sysctl(mib, 2, &value, &length, NULL, 0))
		{
			throw ITLException("Sysctl call failed when requesting system RAM size.");
		}

		return value;

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
	NOTE: Might not show the progress indicator if the function is not called for every possible n between [0, max-1].
	@param n The zero-based index of last iteration that has been completed.
	@param max Number of iterations.
	*/
	inline void showProgress(size_t n, size_t max, bool show = true)
	{
		if (!show)
			return;

		if (n > 0)
		{
			static int count = 0;
			count++;

			if (isTerminal())
			{
				coord_t prevProgress = round((double)(n - 1) / (double)(max - 1) * 100);
				coord_t currProgress = round((double)(n) / (double)(max - 1) * 100);
				if (currProgress > prevProgress)
				{
					if(currProgress < 100)
						::std::cout << currProgress << " %              \r" << ::std::flush;
					else
						::std::cout << "                   \r" << ::std::flush;
				}
			}
			else
			{
				coord_t prevProgress = round((double)(n - 1) / (double)(max - 1) * 100 / 10);
				coord_t currProgress = round((double)(n) / (double)(max - 1) * 100 / 10);
				for(coord_t n = prevProgress; n < currProgress; n++)
				{
					::std::cout << "=" << ::std::flush;
				}

				if (prevProgress < 10  && currProgress == 10)
				{
					::std::cout << ::std::endl;
				}
			}
		}
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

		size_t localCounter;

#if defined(_MSC_VER)
		#pragma omp atomic
		counter++;

		localCounter = counter;
#else
		#pragma omp atomic capture
		localCounter = ++counter;
#endif

		if (localCounter > 0)
		{
			if (isTerminal())
			{
				coord_t prevProgress = round((double)(localCounter - 1) / (double)(max - 1) * 100);
				coord_t currProgress = round((double)(localCounter    ) / (double)(max - 1) * 100);
				if (currProgress > prevProgress)
				{
					#pragma omp critical
					{
						if (currProgress < 100)
							::std::cout << currProgress << " %              \r" << ::std::flush;
						else
							::std::cout << "                   \r" << ::std::flush;
					}
				}
			}
			else
			{
				coord_t prevProgress = round((double)(localCounter - 1) / (double)(max - 1) * 100 / 10);
				coord_t currProgress = round((double)(localCounter    ) / (double)(max - 1) * 100 / 10);
				if (currProgress > prevProgress)
				{
					#pragma omp critical
					{
						for (coord_t n = prevProgress; n < currProgress; n++)
						{
							::std::cout << "=" << ::std::flush;
						}

						if (prevProgress < 10 && currProgress == 10)
						{
							::std::cout << ::std::endl;
						}
					}
				}
			}
		}

		/*#pragma omp critical(show_progress)
		{
			showProgress(counter, max);
			counter++;
		}*/
	}







	/**
	Convert from string to type T.
	*/
	template<typename T>
	T fromString(const ::std::string& str)
	{
		::std::istringstream stream(str);
		T t;
		stream >> t;
		if(stream.fail())
		{
			::std::ostringstream errstr;
			errstr << "Invalid conversion from string '" << str << "'";
			throw ITLException(errstr.str());
		}
		return t;
	}

	/**
	Dummy string to string converter.
	*/
	template<> inline
		::std::string fromString(const ::std::string& str)
	{
		return str;
	}

	/**
	Convert from string to double.
	Supports nan, inf, and -inf.
	*/
	template<>
	inline double fromString(const ::std::string& str)
	{
		string lostr = str;
		toLower(lostr);
		trim(lostr);

		// Possible formats that we want to accept are e.g.
		// -1#QNAN
		// -1.#QNAN
		// -1.#IND
		// -nan(ind)
		// -nan
		// -qnan
		// but not e.g.
		// diiba daa nannannaa foo bar bar baa
		if(contains(lostr, "nan") && lostr.length() < 10)
			return std::numeric_limits<double>::quiet_NaN();

		// Old more restrictive version
		//if(str == "nan" || str == "NAN" || str == "NaN" ||
		//	str == "-nan" || str == "-NAN" || str == "-NaN")
		//	return std::numeric_limits<double>::quiet_NaN();

		// Formats that we want to accept:
		// -inf
		// -infinity
		// -1.#INF
		if (contains(lostr, "inf") && lostr.length() < 10)
		{
			if(!startsWith(lostr, "-"))
				return std::numeric_limits<double>::infinity();
			else
				return -std::numeric_limits<double>::infinity();
		}

		// Old more restrictive versions
		//if(str == "inf" || str == "INF" || str == "Inf" || str == "Infinity" || str == "infinity")
		//	return std::numeric_limits<double>::infinity();

		//if(str == "-inf" || str == "-INF" || str == "-Inf" || str == "-Infinity" || str == "-infinity")
		//	return -std::numeric_limits<double>::infinity();

		::std::istringstream stream(str);
		double t;
		stream >> t;
		if(stream.fail())
		{
			::std::ostringstream errstr;
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
	inline float32_t fromString(const ::std::string& str)
	{
		return (float32_t)fromString<double>(str);
	}

	/**
	Converts string to bool.
	*/
	template<>
	inline bool fromString(const ::std::string& value)
	{
		if (value == "true" || value == "TRUE" || value == "True" || value == "1" || value == "on" || value == "ON" || value == "On")
		{
			return true;
		}

		if (value == "false" || value == "FALSE" || value == "False" || value == "0" || value == "off" || value == "OFF" || value == "Off")
		{
			return false;
		}

		::std::ostringstream errstr;
		errstr << "Invalid conversion from string '" << value << "' to boolean.";
		throw ITLException(errstr.str());
	}

	///**
	//Convert from string to vector
	//*/
	//template<>
	//inline Vec2f fromString(const string& value)
	//{
	//	try
	//	{
	//		if (value.length() < 1)
	//			throw ITLException("Empty string.");

	//		if (value[0] == '[')
	//		{
	//			// Vector notation [1, 2]
	//			std::stringstream parts;
	//			parts << value;
	//			std::string sx, sy;
	//			std::getline(parts, sx, '[');
	//			std::getline(parts, sx, ',');
	//			std::getline(parts, sy, ']');

	//			trim(sx);
	//			trim(sy);

	//			Vec2f result;
	//			result.x = itl2::fromString<float32_t>(sx);
	//			result.y = itl2::fromString<float32_t>(sy);
	//			return result;
	//		}
	//		else
	//		{
	//			// Single number

	//			float32_t val = itl2::fromString<float32_t>(value);

	//			Vec2f result;
	//			result.x = val;
	//			result.y = val;
	//			return result;
	//		}
	//	}
	//	catch (ITLException)
	//	{
	//		ostringstream errstr;
	//		errstr << "Value '" << value << "' is not a valid 2-component vector.";
	//		throw ITLException(errstr.str());
	//	}
	//}

	namespace internals
	{
		/**
		Helps to parse vec2 string to elements of generic type.
		*/
		template<typename T> void vec2ParseHelper(::std::string value, T& x, T& y)
		{
			try
			{
				if (value.length() < 1)
					throw ITLException("Empty string.");

				trim(value);
				trimStart(value, " [");
				trimEnd(value, " ]");

				if (contains(value, ",") || contains(value, " "))
				{
					// Vector notation "[1, 2]", "[1 2]", "1,2", etc.

					std::vector<::std::string> parts = split(value, false, ',', true);
					if (parts.size() != 2)
						parts = split(value, false, ' ', true);

					if (parts.size() != 2)
						throw ITLException("Invalid string.");


					x = itl2::fromString<T>(parts[0]);
					y = itl2::fromString<T>(parts[1]);
				}
				else
				{
					// Single number, either 1 or [1]

					T val = itl2::fromString<T>(value);

					x = val;
					y = val;
				}
			}
			catch (ITLException)
			{
				std::ostringstream errstr;
				errstr << "Value '" << value << "' is not a valid 2-component vector.";
				throw ITLException(errstr.str());
			}
		}

		/**
		Helps to parse vec3 string to elements of generic type.
		*/
		template<typename T> void vec3ParseHelper(::std::string value, T& x, T& y, T& z)
		{
			try
			{
				if (value.length() < 1)
					throw ITLException("Empty string.");

				trimStart(value, " [");
				trimEnd(value, " ]");

				if (contains(value, ",") || contains(value, " "))
				{
					// Vector notation "[1, 2, 3]", "[1 2 3]", "1,2,3", etc.

					std::vector<::std::string> parts = split(value, false, ',', true);
					if (parts.size() != 3)
						parts = split(value, false, ' ', true);

					if (parts.size() != 3)
						throw ITLException("Invalid string.");

					
					x = itl2::fromString<T>(parts[0]);
					y = itl2::fromString<T>(parts[1]);
					z = itl2::fromString<T>(parts[2]);
				}
				else
				{
					// Single number, either 1 or [1]

					T val = itl2::fromString<T>(value);

					x = val;
					y = val;
					z = val;
				}
			}
			catch (ITLException)
			{
				std::ostringstream errstr;
				errstr << "Value '" << value << "' is not a valid 3-component vector.";
				throw ITLException(errstr.str());
			}
		}
	}

	/**
	Convert from string to vector
	*/
	template<>
	inline Vec2d fromString(const ::std::string& value)
	{
		Vec2d v;
		internals::vec2ParseHelper(value, v.x, v.y);
		return v;
	}

	/**
	Convert from string to vector
	*/
	template<>
	inline Vec2f fromString(const ::std::string& value)
	{
		Vec2f v;
		internals::vec2ParseHelper(value, v.x, v.y);
		return v;
	}

	/**
	Convert from string to vector.
	*/
	template<>
	inline Vec2<coord_t> fromString(const ::std::string& value)
	{
		Vec2c v;
		internals::vec2ParseHelper(value, v.x, v.y);
		return v;
	}

	/**
	Convert from string to vector
	*/
	template<>
	inline Vec3d fromString(const ::std::string& value)
	{
		Vec3d v;
		internals::vec3ParseHelper(value, v.x, v.y, v.z);
		return v;
	}

	/**
	Convert from string to vector
	*/
	template<>
	inline Vec3f fromString(const ::std::string& value)
	{
		Vec3f v;
		internals::vec3ParseHelper(value, v.x, v.y, v.z);
		return v;
	}

	/**
	Convert from string to vector.
	*/
	template<>
	inline Vec3<coord_t> fromString(const ::std::string& value)
	{
		Vec3c v;
		internals::vec3ParseHelper(value, v.x, v.y, v.z);
		return v;
	}

	/**
	Converts strings in the given list to type T.
	*/
	template<typename T> ::std::vector<T> fromString(const ::std::vector<::std::string>& values)
	{
		::std::vector<T> result;
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
	::std::string toString(const T& x)
	{
		std::stringstream s;
		s << x;
		return s.str();
	}

	template<>
	inline ::std::string toString(const uint8_t& x)
	{
		return toString((int)x);
	}

	template<>
	inline ::std::string toString(const double& x)
	{
		if (::std::isnan(x))
			return "nan";
		if (::std::isinf(x) && x > 0)
			return "inf";
		if (::std::isinf(x) && x < 0)
			return "-inf";

		::std::stringstream s;
		s << std::setprecision(std::numeric_limits<double>::digits10 + 1) << x;
		return s.str();
	}

	template<>
	inline ::std::string toString(const float& x)
	{
		if (::std::isnan(x))
			return "nan";
		if (::std::isinf(x) && x > 0)
			return "inf";
		if (::std::isinf(x) && x < 0)
			return "-inf";

		::std::stringstream s;
		s << std::setprecision(std::numeric_limits<float>::digits10 + 1) << x;
		return s.str();
	}

	template<>
	inline ::std::string toString(const bool& x)
	{
		if (x)
			return "True";
		return "False";
	}


	inline double sizeRound(double size)
	{
		return round(size * 100.0) / 100.0;
	}

	/**
	Converts value in bytes to nice string.
	*/
	inline ::std::string bytesToString(double size)
	{
		double sizeBytes = sizeRound(size);
		double sizeKilos = sizeRound(size / 1024.0);
		double sizeMegas = sizeRound(size / (1024.0 * 1024.0));
		double sizeGigas = sizeRound(size / (1024.0 * 1024.0 * 1024.0));
		double sizeTeras = sizeRound(size / (1024.0 * 1024.0 * 1024.0 * 1024.0));

		if (sizeKilos < 0.5)
			return itl2::toString(sizeBytes) + " bytes";

		if (sizeMegas < 0.5)
			return itl2::toString(sizeKilos) + " kiB";

		if (sizeGigas < 0.5)
			return itl2::toString(sizeMegas) + " MiB";

		if (sizeTeras < 0.5)
			return itl2::toString(sizeGigas) + " GiB";

		return itl2::toString(sizeTeras) + " TiB";
	}

	template<typename T>
	bool listContains(const std::list<T>& l, const T& elem){
		return std::find(l.begin(), l.end(), elem) != l.end();
	}
}
