#pragma once
#include <type_traits>
#include <fstream>
#include "itlexception.h"

namespace itl2
{
	
	/**
	Writes any trivially copyable value to a stream.
	*/
	template<typename T>
	typename ::std::enable_if<std::is_trivially_copyable_v<T> >::type writeItem(std::ofstream& out, const T& item)
	{
		out.write((char*)&item, sizeof(T));
	}

	/**
	Reads any trivially copyable value from a stream.
	*/
	template<typename T>
	typename ::std::enable_if<std::is_trivially_copyable_v<T> >::type readItem(std::ifstream& in, T& item)
	{
		in.read((char*)&item, sizeof(T));
	}

	/**
	Reads a list of values from a stream, given a function that can read one item.
	*/
	template<typename item_t, typename ReadItem = decltype(readItem<item_t>)> void readList(std::ifstream& in, std::vector<item_t>& v, ReadItem readItem = itl2::readItem<item_t>)
	{
		size_t s = 0;
		in.read((char*)&s, sizeof(size_t));

		v.reserve(v.size() + s);

		for (size_t n = 0; n < s; n++)
		{
			item_t val;
			readItem(in, val);

			if (in.eof() || !in.good())
				throw ITLException("Invalid input file or unable to read.");

			v.push_back(val);
		}
	}

	/**
	Reads a list of values from a file, given a function that can read one item.
	*/
	template<typename item_t, typename ReadItem = decltype(readItem<item_t>)> void readListFile(const std::string& filename, std::vector<item_t>& v, ReadItem readItem = itl2::readItem<item_t>)
	{
		std::ifstream in(filename, std::ios_base::in | std::ios_base::binary);
		if (!in)
			throw ITLException(string("Unable to open file: ") + filename);

		try
		{
			readList<item_t, ReadItem>(in, v, readItem);
		}
		catch (const ITLException& e)
		{
			throw ITLException(e.message() + " (" + filename + ")");
		}
	}

	/**
	Writes a list of values to a stream, given a function that can write one item.
	*/
	template<typename item_t, typename WriteItem = decltype(writeItem<item_t>)> void writeList(std::ofstream& out, const std::vector<item_t>& v, WriteItem writeItem = itl2::writeItem<item_t>)
	{
		size_t s = v.size();
		out.write((const char*)&s, sizeof(size_t));

		for (size_t n = 0; n < v.size(); n++)
		{
			writeItem(out, v[n]);
		}
	}

	/**
	Writes a list of values to a file, given a function that can write one item.
	*/
	template<typename item_t, typename WriteItem = decltype(writeItem<item_t>)> void writeListFile(const std::string& filename, const std::vector<item_t>& v, WriteItem writeItem = itl2::writeItem<item_t>)
	{
		createFoldersFor(filename);

		std::ofstream out(filename, std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
		if (!out)
			throw ITLException(string("Unable to write to: ") + filename);

		writeList<item_t, WriteItem>(out, v, writeItem);
	}

}