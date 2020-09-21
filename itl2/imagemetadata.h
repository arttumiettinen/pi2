#pragma once

#include <string>
#include <map>
#include <vector>

#include "utilities.h"

namespace itl2
{
	// TODO: Add support for multi-line string values and values that contain equals signs.

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

		void setStr(const std::string& key, const std::string& value, coord_t i, coord_t j);

		static std::string listToString(const std::vector<std::string>& list)
		{
			std::ostringstream str;
			appendValueOrVector(str, list);
			return str.str();
		}

		std::string getStr(const std::string& k, coord_t i, coord_t j) const
		{
			const auto& mat = getStringMatrix(k);
			if (i >= 0)
			{
				if (j >= 0)
				{
					// Retrieve single element
					return mat[i][j];
				}
				else
				{
					// Retrieve whole row i as string
					return listToString(mat[i]);
				}
			}
			else
			{
				if (j >= 0)
				{
					// Retrieve one column
					throw ITLException("Unimplemented: Retrieval of a column of data from image metadata object.");
				}
				else
				{
					// Retrieve whole thing
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
		Sets value corresponding to a key.
		*/
		template<typename T> void set(const std::string& key, T value, coord_t i = -1, coord_t j = -1)
		{
			setStr(key, toString(value), i, j);
		}

		
		/**
		Retrieves data with given key.
		If the key does not exist, returns the given default value.
		@param i, j Coordinates of the element to get in the data matrix. Specify -1 to retrieve whole row or column.
		*/
		template<typename T> T get(const std::string& key, T def, coord_t i = -1, coord_t j = -1) const
		{
			if (!contains(key))
				return def;
			return fromString<T>(getStr(key, i, j));
		}

		/**
		Gets item as list of double values.
		*/
		template<typename T> std::vector<T> getList(const std::string& key) const
		{
			auto& mat = getStringMatrix(key);

			std::vector<T> lst;
			lst.reserve(mat.size());

			for(const std::vector<string>& row : mat)
			{
				if (row.size() == 1)
					lst.push_back(fromString<T>(row[0]));
				else
				{
					std::ostringstream str;
					appendValueOrVector(str, row);
					lst.push_back(fromString<T>(str.str()));
				}
			}

			return lst;
		}

		

		/**
		Gets item as matrix.
		*/
		template<typename T> std::vector<std::vector<T>> getMatrix(const std::string& key) const
		{
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