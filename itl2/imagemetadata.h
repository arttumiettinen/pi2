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

		void setStr(const std::string& key, const std::string& value);

		std::string getStr(const std::string& k, int i, int j) const
		{
			const auto& mat = getStringMatrix(k);
			return mat[i][j];
		}

		static void getIndices(string& key, int& i, int& j);

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
		Key string can be followed by at most two indices, e.g. key[0][1] would set element (0, 1) at given key.
		*/
		template<typename T> void set(const std::string& key, T value)
		{
			setStr(key, toString(value));
		}


		

		/**
		Retrieves data with given key.
		If the key does not exist, returns the given default value.
		*/
		template<typename T> T get(const std::string& key, T def) const
		{
			string k = key;
			int i, j;
			getIndices(k, i, j);
			if (!contains(k))
				return def;
			return fromString<T>(getStr(k, i, j));
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