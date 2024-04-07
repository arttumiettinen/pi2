#pragma once

#include <cstdio>

namespace itl2
{
#ifndef _WIN32	
	inline int fopen_s(FILE** f, const char* name, const char* mode)
	{
		int ret = 0;
		*f = fopen(name, mode);
		if (!*f)
			ret = errno;
		return ret;
	}
#endif
}