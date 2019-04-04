#pragma once

#include <cstring>

#include "math/numberutils.h"
#include "test.h"
#include "itlexception.h"

namespace math
{

	/**
	Generic m x n matrix.
	*/
	template <typename real_t = double> class Matrix
	{
	private:
		/**
		Elements of the matrix
		*/
		real_t* pElem;

		/**
		Count of rows and count of columns.
		*/
		size_t rows, columns;

		/**
		Calculate index in pElem array from row and column indices.
		*/
		size_t getIndex(const size_t row, const size_t column) const
		{
			if (row >= rows || column >= columns)
				throw itl2::ITLException("Index out of bounds.");

			return column * rows + row;
		}

		/**
		Check that dimensions of this and r matrix are equal. Throw exception on error.
		*/
		void checkDimensions(const Matrix<real_t>& r) const
		{
			if (!isSameSize(r))
				throw itl2::ITLException("Matrix size mismatch.");
		}

	public:
		Matrix(const size_t rows = 1, const size_t columns = 1)
		{
			if (rows == 0 || columns == 0)
				throw itl2::ITLException("Invalid initialization of matrix with zero rows or columns.");

			this->rows = rows;
			this->columns = columns;
			pElem = new real_t[rows * columns];
			memset(pElem, 0, count() * sizeof(real_t));
		}

		~Matrix()
		{
			delete [] pElem;
		}

		Matrix(const Matrix<real_t>& r)
		{
			rows = r.rows;
			columns = r.columns;

			pElem = new real_t[rows * columns];
			operator=(r);
		}

		Matrix<real_t>& operator=(const Matrix<real_t>& r)
		{
			if (this == &r)
				return *this;

			if (!isSameSize(r))
			{
				delete pElem;
				rows = r.rows;
				columns = r.columns;
				pElem = new real_t[rows * columns];
			}

			for (size_t n = 0; n < count(); n++)
				pElem[n] = r.pElem[n];

			return *this;
		}

		operator real_t() const
		{
			if (count() != 1)
				throw itl2::ITLException("Invalid conversion to number from matrix with more than one element.");

			return pElem[0];
		}

		size_t count() const
		{
			return rows * columns;
		}

		size_t rowCount() const
		{
			return rows;
		}

		size_t columnCount() const
		{
			return columns;
		}

		bool isSquare() const
		{
			return rows == columns;
		}

		real_t& operator()(const size_t row, const size_t column = 0)
		{
			return pElem[getIndex(row, column)];
		}

		const real_t& operator()(const size_t row, const size_t column = 0) const
		{
			return pElem[getIndex(row, column)];
		}

		bool isSameSize(const Matrix<real_t>& r) const
		{
			return r.rows == rows && r.columns == columns;
		}

		const bool equals(const Matrix<real_t>& r, real_t tolerance = NumberUtils<real_t>::tolerance()) const
		{
			if (!isSameSize(r))
				return false;

			for (size_t n = 0; n < count(); n++)
			{
				if (!NumberUtils<real_t>::equals(r.pElem[n], pElem[n], tolerance))
					return false;
			}
			return true;
		}

		const bool operator==(const Matrix<real_t>& r) const
		{
			return equals(r);
		}

		const bool operator!=(const Matrix<real_t>& r) const
		{
			return !(*this == r);
		}

		Matrix<real_t>& operator+=(const real_t r)
		{
			for (size_t n = 0; n < count(); n++)
				pElem[n] += r;

			return *this;
		}

		Matrix<real_t>& operator-=(const real_t r)
		{
			for (size_t n = 0; n < count(); n++)
				pElem[n] -= r;

			return *this;
		}

		Matrix<real_t>& operator*=(const real_t r)
		{
			for (size_t n = 0; n < count(); n++)
				pElem[n] *= r;

			return *this;
		}

		Matrix<real_t>& operator/=(const real_t r)
		{
			for (size_t n = 0; n < count(); n++)
				pElem[n] /= r;

			return *this;
		}

		Matrix<real_t>& operator+=(const Matrix<real_t> &r)
		{
			checkDimensions(r);

			for (size_t n = 0; n < count(); n++)
				pElem[n] += r.pElem[n];

			return *this;
		}

		Matrix<real_t>& operator-=(const Matrix<real_t> &r)
		{
			checkDimensions(r);

			for (size_t n = 0; n < count(); n++)
				pElem[n] -= r.pElem[n];

			return *this;
		}

		const Matrix<real_t> operator+(const Matrix<real_t>& r) const
		{
			Matrix<real_t> result = *this;
			result += r;
			return result;
		}

		const Matrix<real_t> operator-(const Matrix<real_t>& r) const
		{
			Matrix<real_t> result = *this;
			result -= r;
			return result;
		}

		const Matrix<real_t> operator+(const real_t r) const
		{
			Matrix<real_t> result = *this;
			result += r;
			return result;
		}

		const Matrix<real_t> operator-(const real_t r) const
		{
			Matrix<real_t> result = *this;
			result -= r;
			return result;
		}

		const Matrix<real_t> operator*(const real_t r) const
		{
			Matrix<real_t> result = *this;
			result *= r;
			return result;
		}

		const Matrix<real_t> operator/(const real_t r) const
		{
			Matrix<real_t> result = *this;
			result /= r;
			return result;
		}
		
		const Matrix<real_t> operator*(const Matrix<real_t>& r) const
		{
			if (columns != r.rows)
				throw itl2::ITLException("Incompatible size of multiplicant matrices.");

			Matrix<real_t> result(rows, r.columns);

			for (size_t i = 0; i < rows; i++)
			{
				for (size_t j = 0; j < r.columns; j++)
				{
					for (size_t k = 0; k < columns; k++)
					{
						result(i, j) += (*this)(i, k) * r(k, j);
					}
				}
			}

			return result;
		}

		Matrix<real_t> transpose() const
		{
			Matrix<real_t> result(columns, rows);

			for (size_t i = 0; i < rows; i++)
			{
				for (size_t j = 0; j < columns; j++)
				{
					result(j, i) = (*this)(i, j);
				}
			}

			return result;
		}

		/**
		Tests if this matrix is symmetric.
		*/
		bool isSymmetric() const
		{
			if (!isSquare())
				return false;

			for (size_t i = 0; i < rows; i++)
			{
				for (size_t j = 0; j < rows; j++)
				{
					if (!NumberUtils<real_t>::equals((*this)(i, j), (*this)(j, i)))
						return false;
				}
			}

			return true;
		}
	};

	template<typename real_t> Matrix<real_t> operator+(const real_t l, const Matrix<real_t>& r)
	{
		return r + l;
	}

	template<typename real_t> Matrix<real_t> operator-(const real_t l, const Matrix<real_t>& r)
	{
		return r + (-l);
	}

	template<typename real_t> Matrix<real_t> operator*(const real_t l, const Matrix<real_t>& r)
	{
		return r * l;
	}


	namespace tests
	{
		inline void matrix()
		{
			Matrix<double> mat(2, 2);

			mat(0, 0) = 1;
			mat(1, 0) = 2;
			mat(0, 1) = 3;
			mat(1, 1) = 4;

			mat(1, 1) = 7;

			Matrix<double> mat2 = mat;

			mat2(1, 1) = 0;

			itl2::testAssert(mat2 != mat, "Matrix inequality");

			mat2(1, 1) = mat(1, 1);

			itl2::testAssert(mat2 == mat, "Matrix equality");

			itl2::testAssert(mat.isSymmetric() == false, "Matrix isSymmetric");

			mat2 * mat2;

			mat + 5.0;
			5.0 + mat;

		}
	}
}

