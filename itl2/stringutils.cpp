
#include "stringutils.h"
#include <sstream>
#include <algorithm>

#include "io/fileutils.h"


namespace itl2
{
	bool isWhitespace(const char s, const string& whiteSpace)
	{
		for (size_t n = 0; n < whiteSpace.length(); n++)
			if (s == whiteSpace[n])
				return true;

		return false;
	}

	bool isWhitespace(const string& s, const string& whiteSpace)
	{
		for (size_t n = 0; n < s.length(); n++)
			if (!isWhitespace(s[n], whiteSpace))
				return false;
		return true;
	}

	void trimStart(string& s, const string& whiteSpace)
	{
		while (s.length() > 0 && isWhitespace(s[0], whiteSpace))
			s.erase(0, 1);
	}

	void trimEnd(string& s, const string& whiteSpace)
	{
		while (s.length() > 0 && isWhitespace(s[s.length() - 1], whiteSpace))
			s.erase(s.length() - 1, 1);
	}

	void trim(string& s, const string& whiteSpace)
	{
		trimStart(s, whiteSpace);
		trimEnd(s, whiteSpace);
	}

	bool contains(const string& str, const string& part)
	{
		return str.find(part) != string::npos;
	}

	bool containsIgnoreCase(const string& str, const string& part)
	{
		string istr = str;
		toLower(istr);
		string ipart = part;
		toLower(ipart);
		return contains(istr, ipart);
	}

	void toLower(string& str)
	{
		std::transform(str.begin(), str.end(), str.begin(), ::tolower);
	}

	bool endsWith(const string &fullString, const string &ending)
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

	bool endsWithIgnoreCase(const string& fullString, const string& ending)
	{
		string iFullString = fullString;
		toLower(iFullString);
		string iEnding = ending;
		toLower(iEnding);
		return endsWith(iFullString, iEnding);
	}

	bool startsWith(const string& fullString, const string& start)
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

	bool startsWithIgnoreCase(const string& fullString, const string& start)
	{
		string iFullString = fullString;
		toLower(iFullString);
		string iStart = start;
		toLower(iStart);
		return startsWith(iFullString, iStart);
	}

	string lastLine(string lines)
	{
		while (lines.length() > 0 && (*lines.rbegin() == '\n' || *lines.rbegin() == '\r'))
			lines = lines.substr(0, lines.length() - 1);

		size_t start = lines.rfind('\n');
		if (start == string::npos)
			return lines;

		return lines.substr(start + 1);
	}

	std::vector<string> split(const string& lines, bool includeEmptyItems, const char delimiter, bool trimItems)
	{
		std::vector<string> result;
		std::stringstream ss(lines);
		std::string item;
		while (std::getline(ss, item, delimiter))
		{
			if (trimItems)
				itl2::trim(item);

			if (includeEmptyItems || item.length() > 0)
				result.push_back(item);
		}
		return result;
	}

	string readText(const string& filename, bool throwOnFail)
	{
		std::ifstream in(filename);
		if (!in.is_open() && throwOnFail)
			throw std::runtime_error(string("File not found: ") + filename);

		string contents((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
		return contents;
	}

	void writeText(const string& filename, const string& contents)
	{
		createFoldersFor(filename);
		std::ofstream out(filename);
		if (!out.is_open())
			throw std::runtime_error(string("Unable to open file for writing: ") + filename);
		out << contents;
	}

	string getToken(string& str, const char* delimiters, char& foundDelimiter)
	{
		// Old version that does not support quoted strings
		//size_t pos = str.find_first_of(delimiters);
		//if (pos == string::npos)
		//{
		//	// No delimiter found
		//	foundDelimiter = 0;
		//	string res = str;
		//	str = "";
		//	return res;
		//}

		//foundDelimiter = str[pos];
		//string res = str.substr(0, pos);
		//str.erase(0, pos + 1);
		//return res;

		string delimList = delimiters;

		char terminator = 0;
		size_t pos = 0;
		while (true)
		{
			if (terminator == 0)
			{
				// Not in a quoted string

				if (pos >= str.length())
				{
					// No delimiter found
					foundDelimiter = 0;
					string res = str;
					str = "";
					return res;
				}
				else if (delimList.find(str[pos]) != string::npos)
				{
					// Delimiter found
					foundDelimiter = str[pos];
					string res = str.substr(0, pos);
					str.erase(0, pos + 1);
					return res;
				}
				else if (str[pos] == '\"')
				{
					// Start of double-quote string
					terminator = '\"';
				}
				else if (str[pos] == '\'')
				{
					// Start of single-quote string
					terminator = '\'';
				}
			}
			else
			{
				// In a quoted string

				// Test for escape sequence
				if (str[pos] == '\\' && pos < str.length())
				{
					// The next character is an escape sequence that cannot be a delimiter, so skip one extra character
					pos++;
				}
				else if(str[pos] == terminator)
				{
					// End of the quoted string
					terminator = 0;
				}
			}

			pos++;
		}

	}
}