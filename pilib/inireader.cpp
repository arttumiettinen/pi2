// Read an INI file into easy-to-access name/value pairs.

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include "ini.h"
#include "inireader.h"

using std::string;

INIReader::INIReader(const string& filename) :
	filename(filename)
{
    _error = ini_parse(filename.c_str(), valueHandler, this);
}

string INIReader::makeKey(const string& section, const string& name)
{
    string key = section + "." + name;
    // Convert to lower case to make section/name lookups case-insensitive
    std::transform(key.begin(), key.end(), key.begin(), ::tolower);
    return key;
}

int INIReader::valueHandler(void* user, const char* section, const char* name,
                            const char* value)
{
    INIReader* reader = (INIReader*)user;
    string key = makeKey(section, name);
    if (reader->_values[key].size() > 0)
        reader->_values[key] += "\n";
    reader->_values[key] += value;
    return 1;
}
