
#include "resultstable.h"

using namespace std;

namespace itl2
{
	void Results::removeColumn(size_t columnIndex)
	{
		titles.erase(titles.begin() + columnIndex);
		for (size_t row = 0; row < size(); row++)
		{
			std::vector<double>& r = (*this)[row];
			r.erase(r.begin() + columnIndex);
		}
	}

	void Results::toStream(ostream& stream, size_t n) const
	{
		const std::vector<double>& row = (*this)[n];
		for (size_t m = 0; m < row.size(); m++)
		{
			stream << row[m];
			if (m < row.size() - 1)
				stream << ", ";
		}
	}

	ostream& operator<<(ostream& stream, const Results& v)
	{
		stream << v.headers();
		if (v.size() > 0)
			stream << endl;

		for (size_t n = 0; n < v.size(); n++)
		{
			v.toStream(stream, n);
			if (n < v.size() - 1)
				stream << endl;
		}
		return stream;
	}

	std::string Results::str() const
	{
		std::stringstream s;
		s << *this;
		return s.str();
	}

	std::string Results::str(size_t n) const
	{
		if (n >= size())
			throw ITLException("Invalid row number.");

		std::stringstream s;
		toStream(s, n);
		return s.str();
	}

	void Results::readText(const string& filename)
	{
		// TODO: This is not very memory efficient for big tables - the data is read twice to memory (first as text, then as doubles)
		string txt = itl2::readText(filename, true);
		vector<string> data = split(txt);

		for (size_t n = 1; n < data.size(); n++) // Skip first line, it is the header line.
		{
			const string& s = data[n];
			vector<string> strrow = split(s, true, ',');
			vector<double> row;
			row.reserve(strrow.size());
			for (string item : strrow)
			{
				double val = fromString<double>(item);
				row.push_back(val);
			}
			push_back(row);
		}
	}

	void Results::toImage(Image<float32_t>& img) const
	{
		img.metadata.set("column headers", headers().str());
		if (size() > 0)
		{
			size_t cols = headers().size();
			img.ensureSize(cols, size());
			
			for (size_t rowi = 0; rowi < size(); rowi++)
			{
				const vector<double>& row = (*this)[rowi];
				for (size_t coli = 0; coli < cols; coli++)
				{
					img(coli, rowi) = (float32_t)row[coli];
				}
			}
		}
		else
		{
			// We have no results to store in the table.
			// There is no functionality in the Image class to represent images of zero height, so we'll throw an exception.
			throw ITLException("The results table is empty and cannot therefore be expressed as an image.");
		}
	}

	void Results::fromImage(const Headers& headers, const Image<float32_t>& img)
	{
		// TODO: Make a version of this function that reads headers from image metadata, if available.

		if (this->size() > 0)
		{
			// Check that the new data is compatible with the old data.

			if (headers != this->headers())
				throw ITLException("Headers are not equal between two results tables.");

			if (img.width() != headers.size())
				throw ITLException("Headers and image have different column count.");
		}
		else
		{
			this->headers() = headers;
		}

		reserve(size() + img.height());

		for (coord_t rowi = 0; rowi < img.height(); rowi++)
		{
			vector<double> row(img.width(), 0.0);
			for (coord_t coli = 0; coli < img.width(); coli++)
			{
				row[coli] = img(coli, rowi);
			}
			push_back(row);
		}
	}
}