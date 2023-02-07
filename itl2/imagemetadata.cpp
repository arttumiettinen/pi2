
#include "imagemetadata.h"
#include "stringutils.h"
#include "itlexception.h"
#include "io/fileutils.h"
#include "test.h"

#include <sstream>

using namespace std;

namespace itl2
{

	
	void ImageMetadata::addEmptyItem(const std::string& key)
	{
		data[key] = vector<vector<string>>();
	}

	vector<string> parseValue(std::string value)
	{
		trim(value);

		vector<string> l;
		if (startsWith(value, "[") && endsWith(value, "]"))
		{
			// Multi-element item
			trimStart(value, "[");
			trimEnd(value, "]");

			
			while (value.length() > 0)
			{
				char dummy;
				string s = getToken(value, ",", dummy);
				trim(s);
				undoEscape(s);
				l.push_back(s);
			}

			//auto list = split(value, true, ',');
			//for (string s : list)
			//{
			//	trim(s);
			//	l.push_back(s);
			//}
		}
		else
		{
			// Single-element item
			undoEscape(value);
			l.push_back(value);
		}

		return l;
	}

	void ImageMetadata::addSingleItem(string key, const string& value)
	{
		trim(key);

		remove(key);
		addEmptyItem(key);
		data[key].push_back(parseValue(value));
	}


	void ImageMetadata::readFromFile(const std::string& filename)
	{
		readFromString(readText(filename, true));
	}

	void ImageMetadata::readFromString(const std::string& datastring)
	{
		std::istringstream input;
		input.str(datastring);

		// Name of active array that we are reading.
		string arrayName = "";
		while (true)
		{
			string line;
			if (!std::getline(input, line))
				break;

			trim(line);

			if (line.length() <= 0)
			{
				arrayName = "";
			}
			else if (itl2::contains(line, "="))
			{
				// This is normal key = value

				arrayName = "";

				//vector<string> parts = itl2::split(line, true, '=', false);
				//if (parts.size() == 1)
				//	addSingleItem(parts[0], "");
				//else if (parts.size() == 2)
				//	addSingleItem(parts[0], parts[1]);
				//else
				//	throw ITLException(string("Line contains multiple values: ") + line);

				size_t pos = line.find('=');
				if (pos == string::npos)
				{
					addSingleItem(line, "");
				}
				else
				{
					string key = line.substr(0, pos - 1);
					string value = line.substr(pos + 1);
					addSingleItem(key, value);
				}
			}
			else
			{
				if (arrayName == "")
				{
					// This is a header.
					// Read until next empty line or line that contains =
					arrayName = line;
					addEmptyItem(arrayName);
				}
				else
				{
					// This is an element for the array
					data[arrayName].push_back(parseValue(line));
				}
			}

		}
	}


	void ImageMetadata::setStr(const string& key, const string& value, coord_t column, coord_t row)
	{
		// Rules
		// setStr("key", "value", 7, -1)		=> "key" contains a list (single row), set element 7 to "value"
		// setStr("key", "value")				=> "key" contains a single value "value"
		// setStr("key", "value", 1, 2)			=> "key" contains a matrix, set element (1, 2) to "value"
		// setStr("key", "[0, 1, 2]", -1, n)	=> n:th row of "key" is set to elements (0, 1, 2)

		if (!contains(key))
			addEmptyItem(key);

		vector<vector<string>>& mat = getStringMatrix(key);

		if (column < 0)
		{
			vector<string> subElements = parseValue(value);
			if (subElements.size() > 1)
			{
				// This is a special case of value like [1, 2, 7, 8]
				// Remove possible existing data, and
				// set all the elements one-by-one.
				if(row >= 0 && (size_t)row < mat.size())
					mat[row].clear();
				for (size_t n = 0; n < subElements.size(); n++)
					setStr(key, subElements[n], n, row);
				return;
			}
		}

		if (row >= 0)
		{
			while (mat.size() <= (size_t)row)
				mat.push_back(vector<string>());

			if (column >= 0)
			{
				while (mat[row].size() <= (size_t)column)
					mat[row].push_back("");
			}
		}
		else
		{
			// row < 0 => assume there is only one row
			if(mat.size() > 1)
				mat.erase(mat.begin() + 1, mat.end());
			else if(mat.size() <= 0)
				mat.push_back(vector<string>());

			if (column >= 0)
			{
				while (mat[0].size() <= (size_t)column)
					mat[0].push_back("");
			}
		}

		if (row >= 0 && column >= 0)
		{
			// Set single item in the matrix
			mat[row][column] = value;
		}
		else if (row >= 0 && column < 0)
		{
			// Set entire single row
			mat[row].clear();
			mat[row].push_back(value);
		}
		else if (row < 0 && column >= 0)
		{
			// Set single column
			// row < 0 so we assume there is only one row.
			mat[0][column] = value;
		}
		else
		{
			// Set whole thing
			mat.clear();
			mat.push_back(vector<string>());
			mat[0].push_back(value);
		}
	}



	void ImageMetadata::appendValueOrVector(ostringstream& str, const vector<string>& value)
	{
		if (value.size() == 1)
		{
			// Single value
			string esc = value[0];
			escape(esc);
			str << esc;
		}
		else if (value.size() > 1)
		{
			// Vector
			str << toString(value);
		}
		else
		{
			throw ITLException("Invalid element in data matrix.");
		}
	}

	string toString(const ImageMetadata& meta)
	{
		std::ostringstream str;

		vector<vector<string>> def;

		for(const string& key : meta.keys())
		{
			const vector<vector<string>> mat = meta.getMatrix<string>(key, def);
			if (mat.size() == 1)
			{
				str << key << " = ";
				ImageMetadata::appendValueOrVector(str, mat[0]);
				str << endl;
			}
			else
			{
				// Matrix
				str << key << endl;
				for(const auto& row : mat)
				{
					ImageMetadata::appendValueOrVector(str, row);
					str << endl;
				}
				str << endl; // Append empty line that terminates the multi-line element.
			}
		}


		return str.str();
	}


	void ImageMetadata::writeToFile(const std::string& path) const
	{
		createFoldersFor(path);
		writeText(path, toString(*this));
	}



	namespace tests
	{
		int matrixElement(int column, int row)
		{
			return (column + 1) * (row + 1) + row;
		}

		int vectorElement(int listIndex, int elementIndex)
		{
			return (listIndex + 1) * (elementIndex + 1) + listIndex;
		}

		void testGetters(const ImageMetadata& data)
		{
			testAssert(data.get<Vec3c>("single_vec3", Vec3c(0, 0, 0)) == Vec3c(1, 2, 3), "single vec3");
			testAssert(data.get<int>("single_vec3", -100, 3) == -100, "column overflow");
			testAssert(data.get<int>("single_vec3", -100, 0, 1) == -100, "row overflow");
			testAssert(data.get<string>("string", "") == "test = test test\n\nsecond line test\nthird line , []]", "string");
			testAssert(data.get<int>("key1", 0) == 1, "int");
			testAssert(data.get<double>("pi", 0) == 3.14159265, "double");
			testAssert(data.get<int>("list", 0, 3) == 4, "list int");
			testAssert(data.get<int>("matrix", 0, 4, 2) == matrixElement(4, 2), "int matrix element");

			// List as individual elements
			testAssert(data.get<int>("vec2_list", -1, 0, 0) == vectorElement(0, 0), "list element");
			testAssert(data.get<int>("vec2_list", -1, 1, 0) == vectorElement(0, 1), "list element");
			testAssert(data.get<int>("vec2_list", -1, 0, 1) == vectorElement(1, 0), "list element");
			testAssert(data.get<int>("vec2_list", -1, 1, 1) == vectorElement(1, 1), "list element");
			testAssert(data.get<int>("vec2_list", -1, 0, 2) == vectorElement(2, 0), "list element");
			testAssert(data.get<int>("vec2_list", -1, 1, 2) == vectorElement(2, 1), "list element");

			// Lists in different formats
			vector<vector<int>> def1 = {};
			vector<vector<int>> def2 = { {10} };
			testAssert(data.getMatrix<int>("vec2_list", def1) == data.getMatrix<int>("vec2_list2", def2), "vec2 lists as matrices");

			// List as vector of Vec2f
			vector<Vec2c> def;
			vector<Vec2c> vlist = data.getList<Vec2c>("vec2_list", def);
			if (testAssert(vlist.size() == 3, "vec2 list size"))
			{
				testAssert(vlist[0] == Vec2c(vectorElement(0, 0), vectorElement(0, 1)), "element");
				testAssert(vlist[1] == Vec2c(vectorElement(1, 0), vectorElement(1, 1)), "element");
				testAssert(vlist[2] == Vec2c(vectorElement(2, 0), vectorElement(2, 1)), "element");
			}
		}

		void imagemetadata()
		{
			ImageMetadata data;
			data.set("string", "test = test test\n\nsecond line test\nthird line , []]");
			data.set("key1", 1);
			data.set("pi", 3.14159265);
			data.set("list", 1, 0);
			data.set("list", 2, 1);
			data.set("list", 3, 2);
			data.set("list", 4, 3);
			data.set("list", 5, 4);
			data.set("single_vec3", Vec3c(1, 2, 3));

			data.set("string list", "two-line\nitem 1", 0);
			data.set("string list", "two-line\nitem 2", 1);

			data.set("two_vec3", 0, 0);
			data.set("two_vec3", 1, 1);
			data.set("two_vec3", 2, 2);
			data.set("two_vec3", 3, 3);
			data.set("two_vec3", Vec3c(1, 2, 3), -1, 0);
			data.set("two_vec3", Vec3c(4, 5, 6), -1, 1);

			for (int column = 0; column < 5; column++)
			{
				for (int row = 0; row < 3; row++)
				{
					data.set(string("matrix"), matrixElement(column, row), column, row);
				}
			}

			for (int listIndex = 0; listIndex < 3; listIndex++)
			{
				for (int element = 0; element < 2; element++)
				{
					data.set(string("vec2_list"), vectorElement(listIndex, element), element, listIndex);
				}

				data.set(string("vec2_list2"), Vec2c(vectorElement(listIndex, 0), vectorElement(listIndex, 1)), -1, listIndex);
			}

			testAssert(data.rowCount("string") == 1, "row count");
			testAssert(data.rowCount("list") == 1, "row count");
			testAssert(data.rowCount("string list") == 1, "row count");
			testAssert(data.rowCount("vec2_list") == 3, "row count");

			testAssert(data.columnCount("string", 0) == 1, "column count");
			testAssert(data.columnCount("list", 0) == 5, "column count");
			testAssert(data.columnCount("string list", 0) == 2, "column count");
			testAssert(data.columnCount("vec2_list", 0) == 2, "column count");
			testAssert(data.columnCount("matrix", 0) == 5, "column count");

			cout << "string directly from metadata:" << endl << data.get<string>("string", "") << endl;
			cout << "string list directly from metadata:" << endl << data.get<string>("string list", "") << endl;
			cout << "string list item 1 directly from metadata:" << endl << data.get<string>("string list", "", 0) << endl;

			string str = toString(data);
			cout << "toString output: " << endl << str << endl;

			testGetters(data);

			ImageMetadata data2 = fromString<ImageMetadata>(str);

			testAssert(str == toString(data2), "ImageMetadata before and after saving");

			// Test getters again after loading
			testGetters(data2);
		}
	}
	
}