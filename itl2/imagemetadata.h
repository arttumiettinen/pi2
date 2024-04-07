#pragma once

#include <string>
#include <map>
#include <vector>

#include "utilities.h"

namespace itl2
{

	/**
	Stores image metadata.
	*/
	class ImageMetadata
	{
	private:
		std::map<std::string, std::vector<std::vector<std::string>>> data;

		void addEmptyItem(const std::string& key);
		void addSingleItem(std::string key, const std::string& value);

		/**
		Get reference to specific data element as a matrix of strings.
		*/
		std::vector<std::vector<std::string>>& getStringMatrix(const std::string& key)
		{
			return data.at(key);
		}

		/**
		Get const reference to specific data element as a matrix of strings.
		Throws std::out_of_range if the key is invalid.
		*/
		const std::vector<std::vector<std::string>>& getStringMatrix(const std::string& key) const
		{
			return data.at(key);
		}

		/**
		Converts list of strings to given data type.
		*/
		template<typename T> static std::vector<T> lineTo(const std::vector<std::string>& line)
		{
			std::vector<T> dl;
			dl.reserve(line.size());

			for (const std::string& s : line)
			{
				dl.push_back(fromString<T>(s));
			}

			return dl;
		}

		void setStr(const std::string& key, const std::string& value, coord_t column, coord_t row);

		static std::string listToString(const std::vector<std::string>& list)
		{
			std::ostringstream str;
			appendValueOrVector(str, list);
			return str.str();
		}

		std::string getStr(const std::string& k, coord_t column, coord_t row) const
		{
			const auto& mat = getStringMatrix(k);
			if (row >= 0)
			{
				if (column >= 0)
				{
					// Retrieve single element
					return mat[row][column];
				}
				else
				{
					// Retrieve whole row as a string
					return listToString(mat[row]);
				}
			}
			else
			{
				if (column >= 0)
				{
					// Retrieve one column
					// row < 0 => Assume there is only one row.
					return mat[0][column];
				}
				else if (mat.size() == 1 && mat[0].size() == 1)
				{
					// Retrieve the single item
					return mat[0][0];
				}
				else
				{
					// Retrieve whole thing as a parsable string.
					std::ostringstream str;
					for (size_t n = 0; n < mat.size(); n++)
					{
						appendValueOrVector(str, mat[n]);
						if(n < mat.size() - 1)
							str << std::endl;
					}
					return str.str();
				}
			}
			
		}

	public:
		static void appendValueOrVector(std::ostringstream& str, const std::vector<std::string>& value);

	public:
		/**
		Constructor
		*/
		ImageMetadata()
		{

		}


		/**
		Tests if the data contains the given key.
		*/
		bool contains(const std::string& key) const
		{
			return data.find(key) != data.end();
		}


		/**
		Removes an element from the data.
		*/
		void remove(const std::string& key)
		{
			auto iter = data.find(key);
			if (iter != data.end())
				data.erase(iter);
		}

		/**
		Gets a list of all keys.
		*/
		std::vector<std::string> keys() const
		{
			std::vector<std::string> keys;
			keys.reserve(data.size());
			for (auto it = data.begin(); it != data.end(); ++it)
				keys.push_back(it->first);
			return keys;
		}

		/**
		Erases all data elements.
		*/
		void clear()
		{
			data.clear();
		}

		/**
		Sets value of an item corresponding to a key.
		If row and column are negative, assumes the item is a single item (not a list or a matrix).
		If row is negative and column is non-negative, assumes the item is a list (not a matrix) and sets one list element.
		If both row and column are given, assumes the item is a matrix, and sets one matrix element.
		If column is -1, and value is in list format "[a, b, c, ...]", replaces the entire column with the given data.
		*/
		template<typename T> void set(const std::string& key, T value, coord_t column = -1, coord_t row = -1)
		{
			setStr(key, toString(value), column, row);
		}

		/**
		Retrieves count of rows in a metadata item.
		Returns -1 if the key is not found.
		*/
		coord_t rowCount(const std::string& key) const
		{
			if (!contains(key))
				return -1;
			const auto& mat = getStringMatrix(key);
			return (coord_t)mat.size();
		}

		/**
		Retrieves count of columns in a specific row of a metadata item.
		Returns -1 if the key is not found or if row value is invalid.
		Row value -1 refers to the first row.
		*/
		coord_t columnCount(const std::string& key, coord_t row) const
		{
			if (!contains(key))
				return -1;

			const auto& mat = getStringMatrix(key);

			// If row is given, test that it is less than matrix size.
			if (row >= 0 && (size_t)row >= mat.size())
				return -1;
			
			// If row is 'nothing', return first row info.
			if (row < 0)
				row = 0;

			return mat[row].size();
		}
		
		/**
		Retrieves data with given key.
		If the key does not exist, returns the given default value.
		@param row, column Coordinates of the element to get in the data matrix. Specify -1 to retrieve whole row or column.
		*/
		template<typename T> T get(const std::string& key, T def, coord_t column = -1, coord_t row = -1) const
		{
			coord_t cc = columnCount(key, row);
			if (cc < 0)
				return def; // Here, data does not contain key, or row is invalid.

			if (column >= 0 && // Column specified?
				column >= cc)  // Column invalid?
				return def;

			return fromString<T>(getStr(key, column, row));
		}

		/**
		Gets item as a list of values.
		@param squeeze If set to true, returns list of all elements of the value matrix irrespective of its shape; use, e.g., to always return a list of N numbers despite it being stored either as 1xN or Nx1 matrix. If set to false, the elements in the returned list correspond to the first dimension of the value matrix.
		*/
		template<typename T> std::vector<T> getList(const std::string& key, std::vector<T>& def, bool squeeze) const
		{
			if (!contains(key))
				return def;

			auto& mat = getStringMatrix(key);

			std::vector<T> lst;
			lst.reserve(mat.size());

			for(const std::vector<string>& row : mat)
			{
				if (row.size() == 1)
				{
					lst.push_back(fromString<T>(row[0]));
				}
				else
				{
					if (squeeze)
					{
						// Convert each element of the row to the output type.
						for (const string& elem : row)
						{
							lst.push_back(fromString<T>(elem));
						}
					}
					else
					{
						// Convert row of values to one value of output type.
						std::ostringstream str;
						appendValueOrVector(str, row);
						lst.push_back(fromString<T>(str.str()));
					}
				}
			}

			return lst;
		}

		

		/**
		Gets item as matrix.
		*/
		template<typename T> std::vector<std::vector<T>> getMatrix(const std::string& key, std::vector<std::vector<T>>& def) const
		{
			if (!contains(key))
				return def;

			auto& mat = getStringMatrix(key);

			std::vector<std::vector<T>> cMat;
			cMat.reserve(mat.size());

			for(const std::vector<string>& line : mat)
			{
				cMat.push_back(lineTo<T>(line));
			}

			return cMat;
		}

		/**
		Write all data elements to given file.
		*/
		void writeToFile(const std::string& path) const;

		/**
		Read data elements from given file.
		*/
		void readFromFile(const std::string& filename);

		/**
		Read data elements from given data string.
		*/
		void readFromString(const std::string& data);
	};

	/**
	Convert all data elements to a string.
	*/
	std::string toString(const ImageMetadata& meta);

	template<>
	inline ImageMetadata fromString(const std::string& data)
	{
		ImageMetadata md;
		md.readFromString(data);
		return md;
	}

	namespace tests
	{
		void imagemetadata();
	}
}