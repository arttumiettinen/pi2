
#include "utilities.h"
#include "itlexception.h"
#include "pilibutilities.h"

#include <random>
#include <experimental/filesystem>

namespace fs = std::experimental::filesystem;

using namespace itl2;
using namespace std;

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


	namespace
	{
		std::mt19937 gen((unsigned int)std::chrono::system_clock::now().time_since_epoch().count());
	}

	string createTempFilename(const string& purpose)
	{
		//unsigned int seed = (unsigned int)std::chrono::system_clock::now().time_since_epoch().count();
		//std::mt19937 gen(seed);

		string tempFilename = string("./tmp_images/") + purpose + "_" + itl2::toString(gen());
		fs::remove(tempFilename);

		return tempFilename;
	}

	void printTitle(ostream& msg, const string& str, int level)
	{
		char uc;
		if (level <= 0)
			uc = '*';
		else if (level == 1)
			uc = '=';
		else if (level == 2)
			uc = '-';
		else
			uc = '~';
		
		msg << str << endl;
		for (size_t n = 0; n < str.length(); n++)
			msg << uc;
		msg << endl << endl;
	}
}
