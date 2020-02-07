#pragma once

#include <cstring>
#include <initializer_list>
#include <algorithm>
#include <iostream>

#include "math/numberutils.h"
#include "itlexception.h"
#include "math/mathutils.h"

namespace itl2
{

	/**
	Generic m x n matrix.
	Supports access as row matrix if column count n is zero.
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

		/**
		Creates identity matrix.
		The matrix has ones on the main diagonal and zeros everywhere else.
		*/
		static Matrix identity(size_t rows, size_t columns)
		{
			Matrix A(rows, columns);

			for (size_t i = 0; i < std::min(rows, columns); i++)
				A(i, i) = 1;

			return A;
		}

		/**
		Creates identity matrix of the same size and type as the given matrix.
		*/
		static Matrix identity(const Matrix& r)
		{
			return identity(r.rowCount(), r.columnCount());
		}

		/**
		Creates zero matrix.
		*/
		static Matrix zero(size_t rows, size_t columns)
		{
			return Matrix(rows, columns);
		}

		/**
		Constructor.
		Initializes the matrix to all zeroes.
		*/
		Matrix(size_t rows = 1, size_t columns = 1) : pElem(0)
		{
			if (rows == 0 || columns == 0)
				throw itl2::ITLException("Invalid initialization of matrix with zero rows or columns.");

			this->rows = rows;
			this->columns = columns;
			pElem = new real_t[rows * columns];
			memset(pElem, 0, count() * sizeof(real_t));
		}

		/**
		Destructor
		*/
		~Matrix()
		{
			delete[] pElem;
		}

		/**
		Copy constructor.
		*/
		Matrix(const Matrix<real_t>& r) : pElem(0)
		{
			rows = r.rows;
			columns = r.columns;

			pElem = new real_t[rows * columns];
			operator=(r);
		}

		/**
		Constructor.
		Create matrix from initializer list using syntax
		Matrix A = {{1, 2, 3},
					{4, 5, 6}}
		*/
		Matrix(const std::initializer_list<std::initializer_list<real_t> > elems) : pElem(0)
		{
			if (elems.size() <= 0)
				throw itl2::ITLException("Invalid initialization of matrix with zero elements.");

			rows = elems.size();
			columns = elems.begin()->size();

			pElem = new real_t[rows * columns];

			int r = 0;
			for (const std::initializer_list<real_t>& row : elems)
			{
				int c = 0;
				for (real_t val : row)
				{
					(*this)(r, c) = val;
					c++;
				}
				r++;
			}
		}

		/**
		Assignment operator.
		Returns reference to this matrix.
		*/
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

		/**
		Assigns the same value to all the elements.
		*/
		Matrix<real_t>& operator=(real_t r)
		{
			for (size_t n = 0; n < count(); n++)
				pElem[n] = r;
		}

		/**
		Converts this matrix to single number.
		Throws exception if the matrix has more than one element.
		*/
		operator real_t() const
		{
			if (count() != 1)
				throw itl2::ITLException("Invalid conversion to number from matrix with more than one element.");

			return pElem[0];
		}

		/**
		Gets count of elements in this matrix.
		*/
		size_t count() const
		{
			return rows * columns;
		}

		/**
		Gets count of rows in this matrix.
		*/
		size_t rowCount() const
		{
			return rows;
		}

		/**
		Gets count of columns in this matrix.
		*/
		size_t columnCount() const
		{
			return columns;
		}

		/**
		Gets a value indicating whether this matrix is square or not.
		*/
		bool isSquare() const
		{
			return rows == columns;
		}

		/**
		Tests if this matrix is symmetric up to given numerical tolerance.
		*/
		bool isSymmetric(real_t tolerance = NumberUtils<real_t>::tolerance()) const
		{
			if (!isSquare())
				return false;

			for (size_t i = 0; i < rows; i++)
			{
				for (size_t j = 0; j < rows; j++)
				{
					if (!NumberUtils<real_t>::equals((*this)(i, j), (*this)(j, i), tolerance))
						return false;
				}
			}

			return true;
		}

		/**
		Gets a reference to an element of this matrix.
		*/
		real_t& operator()(const size_t row, const size_t column = 0)
		{
			return pElem[getIndex(row, column)];
		}

		const real_t& operator()(const size_t row, const size_t column = 0) const
		{
			return pElem[getIndex(row, column)];
		}

		/**
		Gets a value indicating whether the size of this matrix equals the size of the given matrix.
		*/
		bool isSameSize(const Matrix<real_t>& r) const
		{
			return r.rows == rows && r.columns == columns;
		}

		/**
		Tests if this matrix equals the given matrix up to the given tolerance.
		*/
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

		/**
		Multiplies two matrices.
		*/
		Matrix<real_t> operator*(const Matrix<real_t>& r) const
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

		/**
		Creates a new matrix that is transpose of this matrix.
		*/
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
		Solve A*X = B
		@param B Right hand side
		@return Solution if A is square, least squares solution otherwise.
		*/
		Matrix solve(Matrix B) const;

		/**
		Solve X*A = B, which is also A'*X' = B'
		@param B Right hand side
		@return Solution if A is square, least squares solution otherwise.
		*/
		Matrix<real_t> solveTranspose(Matrix B) const
		{
			return transpose().solve(B.transpose());
		}

		/**
		Calculate inverse or pseudoinverse of this matrix.
		@return Inverse of this matrix if this matrix is is square, pseudoinverse otherwise.
		*/
		Matrix<real_t> inverse() const
		{
			return solve(identity(rowCount(), rowCount()));
		}

		/**
		Calculate determinant of this matrix
		*/
		real_t det() const;

		/**
		Get a submatrix.
		@param r Array of row indices.
		@param j0 Initial column index
		@param j1 Final column index
		@return A(r(:), j0:j1)
		*/
		Matrix<real_t> subMatrix(const Matrix<size_t> &rows, size_t j0, size_t j1) const
		{
			size_t rowCount = rows.count();

			Matrix<real_t> X(rowCount, j1 - j0 + 1);

			for (size_t i = 0; i < rowCount; i++)
			{
				for (size_t j = j0; j <= j1; j++)
					X(i, j - j0) = (*this)(rows(i), j);
			}

			return X;
		}

		/**
		Converts this object to string.
		*/
		friend std::ostream& operator<<(std::ostream& stream, const Matrix<real_t>& v)
		{
			stream << "[";
			for (size_t i = 0; i < v.rowCount(); i++)
			{
				if (i > 0)
					stream << " ";
				stream << "[";
				for (size_t j = 0; j < v.columnCount(); j++)
				{
					stream << v(i, j);
					if (j < v.columnCount() - 1)
						stream << ", ";
				}
				stream << "]";
				if (i < v.rowCount() - 1)
					stream << "," << std::endl;
			}
			stream << "]";
			
			return stream;
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



}


// Include matrix algorithms
#include "math/ludecomposition.h"
#include "math/qrdecomposition.h"


// Implementation of methods requiring algorithms
namespace itl2
{

	template<typename real_t> Matrix<real_t> Matrix<real_t>::solve(Matrix B) const
	{
		if (isSquare())
			return LUDecomposition(*this).solve(B);
		else
			return QRDecomposition(*this).solve(B);
	}

	template<typename real_t> real_t Matrix<real_t>::det() const
	{
		return LUDecomposition(this).det();
	}





	namespace tests
	{
		void matrix();
		void solve();
		void leastSquares();
	}
}

