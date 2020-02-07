#pragma once

#include <vector>
#include <string>

namespace pilib
{
	/**
	Parses each string in the given string list and finds line that ends with given suffix. Converts characters before the suffix to
	integer and calculates the sum of all of them.
	If elements of 'list' are

	kfkfkf
	5 XYZ
	fkgg

	6 XYZ
	sks
	ggigg

	returns 5 + 6 = 11.
	*/
	size_t parseTotalCount(const std::vector<std::string>& list, std::string suffix);

	/**
	Creates name for temporary file and deletes it if it exists.
	@param purpose This string is added to the file name to identify it among other temporary files.
	*/
	std::string createTempFilename(const std::string& purpose);

	/**
	Removes duplicate elements from the given list.
	*/
	template<typename T>
	static void removeDuplicates(std::vector<T>& v)
	{
		for (size_t n = 0; n < v.size(); n++)
		{
			for (size_t m = n + 1; m < v.size(); m++)
			{
				if (v[m] == v[n])
				{
					v.erase(v.begin() + m);
					m--;
				}
			}
		}
	}

	/**
	Prints underlined string like

	This is str
	-----------

	The parameter level controls the type of underlining.
	*/
	void printTitle(std::ostream& msg, const std::string& str, int level);
}
