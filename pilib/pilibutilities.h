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
}
