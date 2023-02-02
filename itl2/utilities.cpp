#include "utilities.h"
#include "testutils.h"

using namespace std;

namespace itl2
{
	void replace_all(string& str, const string& oldValue, const string& newValue)
	{
		std::string::size_type n = 0;
		while ((n = str.find(oldValue, n)) != std::string::npos)
		{
			str.replace(n, oldValue.size(), newValue);
			n += newValue.size();
		}
	}

	void escape(std::string& value)
	{
		replace_all(value, "\\", "\\\\");
		replace_all(value, "[", "\\[");
		replace_all(value, "]", "\\]");
		replace_all(value, "=", "\\=");
		replace_all(value, ",", "\\,");
		replace_all(value, "\n", "\\n");
		replace_all(value, "\r", "\\r");
	}

	void undoEscape(std::string& value)
	{
		replace_all(value, "\\r", "\r");
		replace_all(value, "\\n", "\n");
		replace_all(value, "\\,", ",");
		replace_all(value, "\\=", "=");
		replace_all(value, "\\]", "]");
		replace_all(value, "\\[", "[");
		replace_all(value, "\\\\", "\\");
	}

	namespace tests
	{
		void escapes()
		{
			string s = "ab\\c\n\rdef=k\nc=[1, 2, 3]";
			string orig = s;

			cout << "Before:" << endl << s << endl;

			escape(s);
			cout << "Escaped:" << endl << s << endl;

			undoEscape(s);
			cout << "Un-escaped:" << endl << s << endl;

			testAssert(s == orig, "escape round robin");
		}
	}
}