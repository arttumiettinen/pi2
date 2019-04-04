#pragma once

#include <vector>
#include <string>

#include "stringutils.h"
#include "itlexception.h"
#include "image.h"

namespace itl2
{
	/**
	Encapsulates headers of (particle analysis) result table.
	*/
	class Headers : public std::vector<std::string>
	{
	public:

		/**
		Find in index of element name in headers array.
		*/
		size_t getColumnIndex(const std::string& name) const
		{
			for (size_t n = 0; n < size(); n++)
			{
				if (iStartsWith((*this)[n], name))
					return n;
			}

			throw ITLException(std::string("Column ") + name + std::string(" not found."));
		}

		/**
		Converts this object to std::string in csv format.
		*/
		friend ostream& operator<<(ostream& stream, const Headers& v)
		{
			for (size_t n = 0; n < v.size(); n++)
			{
				stream << v[n];
				if (n < v.size() - 1)
					stream << ", ";
			}
			return stream;
		}

		/**
		Converts this object to std::string in csv format.
		*/
		std::string str() const
		{
			std::stringstream s;
			s << *this;
			return s.str();
		}
	};

	/**
	Encapsulates results of (particle analysis) as a table.
	*/
	class Results : public std::vector<std::vector<double> >
	{
	private:

		Headers titles;

		/**
		Outputs line n of the table to the given stream (as csv text).
		*/
		void toStream(ostream& stream, size_t n) const;

	public:

		/**
		Gets a reference to column headers.
		*/
		Headers& headers()
		{
			return titles;
		}

		/**
		Gets a reference to column headers.
		*/
		const Headers& headers() const
		{
			return titles;
		}

		/**
		Returns index of column that has the given name, or throws an exception if the column is not found.
		@param name Case-insensitive prefix of the name of the column to retrieve.
		*/
		size_t getColumnIndex(const string& columnName) const
		{
			return headers().getColumnIndex(columnName);
		}

		/**
		Gets result from results table.
		@param name Case-insensitive prefix of the name of the column to retrieve.
		*/
		double get(const std::string& name, size_t row)
		{
			return (*this)[row][headers().getColumnIndex(name)];
		}

		/**
		Removes column from this results table.
		*/
		void removeColumn(const std::string& namePrefix)
		{
			removeColumn(headers().getColumnIndex(namePrefix));
		}

		/**
		Removes column from this results table.
		*/
		void removeColumn(size_t columnIndex);

		/**
		Converts this object to std::string in csv format.
		*/
		friend ostream& operator<<(ostream& stream, const Results& v);

		/**
		Converts this object to std::string in csv format.
		*/
		std::string str() const;

		/**
		Returns one row of results table as string.
		*/
		std::string str(size_t n) const;

		/**
		Reads results table in the given file and appends it to this table.
		*/
		void readText(const string& filename);

		/**
		Convert this results table to image.
		*/
		void toImage(Image<float32_t>& img) const;

		/**
		Read data from the given image to this results table.
		Throws exception if the headers don't match exactly.
		*/
		void fromImage(const Headers& headers, const Image<float32_t>& img);
	};
}
