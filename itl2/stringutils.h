#pragma once

#include <string>
#include <sstream>
#include <algorithm>
#include <vector>

#include "io/fileutils.h"

using std::string;

namespace itl2
{
	/**
	Default set of white-space characters.
	*/
	#define WHITE_SPACE_CHARS " \t\r\n" 

	/**
	Tests if the given character is whitespace.
	*/
	inline bool isWhiteSpace(const char s, const string& whiteSpace = WHITE_SPACE_CHARS)
	{
		//return s == ' ' || s == '\t' || s == '\r' || s == '\n';
		for (size_t n = 0; n < whiteSpace.length(); n++)
			if (s == whiteSpace[n])
				return true;

		return false;
	}

	/**
	Removes whitespace from the beginning of the string.
	*/
	inline void trimStart(string& s, const string& whiteSpace = WHITE_SPACE_CHARS)
	{
		while (s.length() > 0 && isWhiteSpace(s[0], whiteSpace))
			s.erase(0, 1);
	}

	/**
	Removes whitespace from the end of the string.
	*/
	inline void trimEnd(string& s, const string& whiteSpace = WHITE_SPACE_CHARS)
	{
		while (s.length() > 0 && isWhiteSpace(s[s.length() - 1], whiteSpace))
			s.erase(s.length() - 1, 1);
	}

	/**
	Removes whitespace from the beginning and the end of the string.
	*/
	inline void trim(string& s, const string& whiteSpace = WHITE_SPACE_CHARS)
	{
		trimStart(s, whiteSpace);
		trimEnd(s, whiteSpace);
	}

	/**
	Tests if string str contains string part.
	*/
	inline bool contains(const string& str, const string& part)
	{
		return str.find(part) != string::npos;
	}

	/**
	Convert all characters in the string to lower case.
	*/
	inline void toLower(string& str)
	{
		std::transform(str.begin(), str.end(), str.begin(), ::tolower);
	}

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
	Tests if a string starts with another string, performs case-insensitive comparison.
	*/
	inline bool iStartsWith(const string& fullString, const string& start)
	{
		if (fullString.length() >= start.length())
		{
			string iFullString = fullString;
			toLower(iFullString);
			string iStart = start;
			toLower(iStart);
			return iFullString.compare(0, iStart.length(), iStart) == 0;
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
	inline std::vector<string> split(const string& lines, bool includeEmptyItems = true, const char delimiter = '\n', bool trimItems = true)
	{
		std::vector<string> result;
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

	/**
	Get next token from a string and remove it from the string.
	Token is the part of string from its beginning to the first occurence of any character in the delimiters string.
	If no delimiter is found, the whole string is returned and zero is assigned to foundDelimiter.
	@param foundDelimiter The delimiter that stopped the token is assigned here.
	@return The token.
	*/
	inline string getToken(string& str, const char* delimiters, char& foundDelimiter)
	{
		size_t pos = str.find_first_of(delimiters);
		if (pos == string::npos)
		{
			// No delimiter found
			foundDelimiter = 0;
			string res = str;
			str = "";
			return res;
		}

		foundDelimiter = str[pos];
		string res = str.substr(0, pos);
		str.erase(0, pos + 1);
		return res;
	}
}
