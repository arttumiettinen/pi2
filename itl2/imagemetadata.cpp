
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
			auto list = split(value, true, ',');
			for (string s : list)
			{
				trim(s);
				l.push_back(s);
			}
		}
		else
		{
			// Single-element item
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
			// TODO: To support equals signs in values, test here if '=' occurs before '"' or the end of the line.
			// TODO: To support multi-line string values, test if the string starts and ends with '"', and if yes, the line should go to current array.
			else if (itl2::contains(line, "="))
			{
				// This is normal key = value

				arrayName = "";

				vector<string> parts = itl2::split(line, true, '=', false);
				if (parts.size() == 1)
					addSingleItem(parts[0], "");
				else if (parts.size() == 2)
					addSingleItem(parts[0], parts[1]);
				else
					throw ITLException(string("Line contains multiple values: ") + line);
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


	/**
	Finds array indices from keys like name[0] and name[7][3].
	*/
	void ImageMetadata::getIndices(string& key, int& i, int& j)
	{
		vector<string> parts = split(key, true, '[');
		if (parts.size() <= 1)
		{
			i = 0;
			j = 0;
			return;
		}
		else if (parts.size() == 2)
		{
			key = parts[0];
			string ind = parts[1];
			trimEnd(ind, "]");
			i = 0;
			j = (int)fromString<unsigned int>(ind);
			return;
		}
		else if (parts.size() == 3)
		{
			key = parts[0];
			string ind1 = parts[1];
			trimEnd(ind1);
			string ind2 = parts[2];
			trimEnd(ind2, "]");
			i = (int)fromString<unsigned int>(ind1);
			j = (int)fromString<unsigned int>(ind2);
			return;
		}
		throw ITLException("Invalid key: " + key);
	}


	void ImageMetadata::setStr(const string& key, const string& value)
	{
		string k = key;

		int i, j;
		getIndices(k, i, j);

		if (!contains(k))
			addEmptyItem(k);

		vector<vector<string>>& mat = getStringMatrix(k);

		while (mat.size() <= i)
			mat.push_back(vector<string>());

		while (mat[i].size() <= j)
			mat[i].push_back("");

		mat[i][j] = value;
	}



	void ImageMetadata::appendValueOrVector(ostringstream& str, const vector<string>& value)
	{
		if (value.size() == 1)
		{
			// Single value
			str << value[0];
		}
		else if (value.size() > 1)
		{
			// Vector
			str << "[";
			for(size_t n = 0; n < value.size(); n++)
			{
				str << value[n];
				if(n < value.size() - 1)
					str << ", ";
			}
			str << "]";
		}
		else
		{
			throw ITLException("Invalid element in data matrix.");
		}
	}

	string toString(const ImageMetadata& meta)
	{
		std::ostringstream str;

		for(const string& key : meta.keys())
		{
			const vector<vector<string>> mat = meta.getMatrix<string>(key);
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
		void imagemetadata()
		{
			ImageMetadata data;
			data.set("string", "test test test");
			data.set("key1", 1);
			data.set("pi", 3.14159265);
			data.set("list[0]", 1);
			data.set("list[1]", 2);
			data.set("list[2]", 3);
			data.set("list[3]", 4);
			data.set("list[4]", 5);

			for (int n = 0; n < 5; n++)
			{
				for (int m = 0; m < 3; m++)
				{
					data.set(string("matrix[") + toString(n) + "][" + toString(m) + "]", n * m);
				}
			}

			for (int n = 0; n < 3; n++)
			{
				for (int m = 0; m < 2; m++)
				{
					data.set(string("vec_list[") + toString(n) + "][" + toString(m) + "]", n + m);
				}
			}

			// getters
			testAssert(data.get<string>("string", "") == "test test test", "string");
			testAssert(data.get<int>("key1", 0) == 1, "int");
			testAssert(data.get<double>("pi", 0) == 3.14159265, "double");
			testAssert(data.get<int>("list[3]", 0) == 4, "list int");
			testAssert(data.get<int>("matrix[4][2]", 0) == 4*2, "matrix int");

			vector<Vec2f> vlist = data.getList<Vec2f>("vec_list");
			testAssert(vlist[0] == Vec2f(0 + 0, 0 + 1), "element");
			testAssert(vlist[1] == Vec2f(1 + 0, 1 + 1), "element");
			testAssert(vlist[2] == Vec2f(2 + 0, 2 + 1), "element");


			string str = toString(data);
			cout << str << endl;

			ImageMetadata data2 = fromString<ImageMetadata>(str);

			testAssert(str == toString(data2), "ImageMetadata before and after saving");

		}
	}
	
}