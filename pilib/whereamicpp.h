#pragma once

#include "whereami.h"

#include <vector>
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;

namespace pilib
{
	/**
	Returns the path to the current module.
	*/
	inline fs::path getModulePath()
	{
		int length = wai_getModulePath(NULL, 0, NULL);

		std::vector<char> path(length + 1, 0);

		wai_getModulePath(&path[0], length, NULL);

		return fs::path(path.begin(), path.end());
	}

}