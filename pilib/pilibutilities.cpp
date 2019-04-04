
#include "utilities.h"
#include "itlexception.h"
#include "pilibutilities.h"

using std::string;
using std::vector;

namespace pilib
{
	size_t parseTotalCount(const vector<string>& list, string suffix)
	{
		size_t total = 0;
		for (size_t n = 0; n < list.size(); n++)
		{
			string out = list[n];

			size_t startPos = out.find(suffix);
			if (startPos != string::npos)
			{
				size_t lineStart = out.rfind('\n', startPos);
				if (lineStart == string::npos)
					lineStart = 0;
				string number = out.substr(lineStart + 1, startPos - lineStart - 1);
				size_t val = fromString<size_t>(number);
				total += val;
			}
			else
			{
				throw ITLException(string("Invalid output from job ") + itl2::toString(n) + string(", line with suffix '") + suffix + string("' not found:\n") + out);
			}
		}

		return total;
	}
}
