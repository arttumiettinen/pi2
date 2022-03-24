#pragma once

#include <string>
#include <vector>

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
	bool isWhitespace(const char s, const string& whiteSpace = WHITE_SPACE_CHARS);

	/**
	Tests if the given string is entirely whitespace or empty.
	*/
	bool isWhitespace(const string& s, const string& whiteSpace = WHITE_SPACE_CHARS);

	/**
	Removes whitespace from the beginning of the string.
	*/
	void trimStart(string& s, const string& whiteSpace = WHITE_SPACE_CHARS);

	/**
	Removes whitespace from the end of the string.
	*/
	void trimEnd(string& s, const string& whiteSpace = WHITE_SPACE_CHARS);

	/**
	Removes whitespace from the beginning and the end of the string.
	*/
	void trim(string& s, const string& whiteSpace = WHITE_SPACE_CHARS);

	/**
	Tests if string str contains string part.
	*/
	bool contains(const string& str, const string& part);

	/**
	Tests if string str contains tring part, ignoring case of both strings.
	*/
	bool containsIgnoreCase(const string& str, const string& part);

	/**
	Convert all characters in the string to lower case.
	*/
	void toLower(string& str);

	/**
	Tests if a string ends with another string.
	*/
	bool endsWith(const string &fullString, const string &ending);

	/**
	Tests if a string ends with another string, performs case-insensitive comparison.
	*/
	bool endsWithIgnoreCase(const string& fullString, const string& ending);

	/**
	Tests if a string starts with another string.
	*/
	bool startsWith(const string& fullString, const string& start);

	/**
	Tests if a string starts with another string, performs case-insensitive comparison.
	*/
	bool startsWithIgnoreCase(const string& fullString, const string& start);


	/**
	Gets last non-empty line in a multi-line string.
	*/
	string lastLine(string lines);

	/**
	Splits string at given delimiter.
	*/
	std::vector<string> split(const string& lines, bool includeEmptyItems = true, const char delimiter = '\n', bool trimItems = true);

	/**
	Reads text file.
	Optionally throws exception on failure; otherwise returns empty string.
	*/
	string readText(const string& filename, bool throwOnFail = false);

	/**
	Write string to text file.
	*/
	void writeText(const string& filename, const string& contents);

	/**
	Get next token from a string and remove it from the string.
	Token is the part of string from its beginning to the first occurence of any character in the delimiters string.
	If no delimiter is found, the whole string is returned and zero is assigned to foundDelimiter.
	@param foundDelimiter The delimiter that stopped the token is assigned here.
	@return The token.
	*/
	string getToken(string& str, const char* delimiters, char& foundDelimiter);
}
