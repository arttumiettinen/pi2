#pragma once
// Read an INI file into easy-to-access name/value pairs.

// Based on work in ini.h/INIReader.h
// inih and INIReader are released under the New BSD license (see LICENSE.txt).
// Go to the project home page for more info:
//
// http://code.google.com/p/inih/
// 
// Changed pass-by-value to const reference where possible.
// Added const keywords where possible.

#include <map>
#include <string>

#include "utilities.h"

/**
Read an INI file into easy-to-access name/value pairs.
*/
class INIReader
{
public:
		
	/**
	Construct INIReader and parse given filename.
	*/
	INIReader(const std::string& filename);

	// Return the result of ini_parse(), i.e., 0 on success, line number of
	// first error on parse error, or -1 on file open error.
	int parseError() const
	{
		return _error;
	}

	/**
	Get a value from the ini file.
	*/
	template<typename value_t> value_t get(const std::string& name, const value_t defaultValue) const
	{
		return get<value_t>("", name, defaultValue);
	}

	/**
	Get a value from the ini file.
	*/
	template<typename value_t> value_t get(const std::string& section, const std::string& name, const value_t defaultValue) const
	{
		std::string key = makeKey(section, name);
		std::map<std::string, std::string>::const_iterator it = _values.find(key);
		if (it != _values.end())
			return itl2::fromString<value_t>(it->second);
		else
			return defaultValue;
	}

	/**
	Gets the name of the last read ini file, if any.
	*/
	inline std::string getFilename() const
	{
		return filename;
	}

private:
	std::string filename;
	int _error;
	std::map<std::string, std::string> _values;
	static std::string makeKey(const std::string& section, const std::string& name);
	static int valueHandler(void* user, const char* section, const char* name,
						const char* value);
};
