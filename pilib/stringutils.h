#pragma once

#include <string>
#include <sstream>
#include <algorithm>

#include "utilities.h"
#include "io/fileutils.h"

using std::string;

namespace pilib
{
	/**
	Tests if a string ends with another string.
	*/
	inline bool endsWith(const string &fullString, const string &ending)
	{
		if (fullString.length() >= ending.length())
		{
			return fullString.compare(fullString.length() - ending.length(), ending.length(), ending) == 0;
		}
		else
		{
			return false;
		}
	}

	/**
	Tests if a string starts with another string.
	*/
	inline bool startsWith(const string& fullString, const string& start)
	{
		if (fullString.length() >= start.length())
		{
			return fullString.compare(0, start.length(), start) == 0;
		}
		else
		{
			return false;
		}
	}

	/**
	Gets last non-empty line in a multi-line string.
	*/
	inline string lastLine(string lines)
	{
		while (lines.length() > 0 && (*lines.rbegin() == '\n' || *lines.rbegin() == '\r'))
			lines = lines.substr(0, lines.length() - 1);

		size_t start = lines.rfind('\n');
		if (start == string::npos)
			return lines;

		return lines.substr(start + 1);
	}

	/**
	Splits string at given delimiter.
	*/
	inline vector<string> split(const string& lines, bool includeEmptyItems = true, const char delimiter = '\n', bool trimItems = true)
	{
		vector<string> result;
		std::stringstream ss(lines);
		std::string item;
		while (std::getline(ss, item, delimiter))
		{
			if (trimItems)
				itl2::trim(item);

			if(includeEmptyItems || item.length() > 0)
				result.push_back(item);
		}
		return result;
	}

	/**
	Reads text file.
	Optionally throws exception on failure; otherwise returns empty string.
	*/
	inline string readText(const string& filename, bool throwOnFail = false)
	{
		std::ifstream in(filename);
		if (!in.is_open() && throwOnFail)
			throw std::runtime_error(string("File not found: ") + filename);

		string contents((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
		return contents;
	}

	/**
	Convert all characters in the string to lower case.
	*/
	inline void toLower(string& str)
	{
		std::transform(str.begin(), str.end(), str.begin(), ::tolower);
	}
	
	/**
	Write string to text file.
	*/
	inline void writeText(const string& filename, const string& contents)
	{
	    createFoldersFor(filename);
	    std::ofstream out(filename);
	    if(!out.is_open())
	        throw std::runtime_error(string("Unable to open file for writing: ") + filename);
	    out << contents;
	}

	inline double sizeRound(double size)
	{
		return math::round(size * 100.0) / 100.0;
	}

	/**
	Converts value in bytes to nice string.
	*/
	inline string bytesToString(double size)
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
}
