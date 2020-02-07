#pragma once

/*
This code is based on public-domain reference implementation of JAMA : A Java Matrix Package (https://math.nist.gov/javanumerics/jama/)
and public-domain library TNT: Template Numerical Toolkit and its JAMA/C++ part (https://math.nist.gov/tnt/).
*/

#include "math/matrix.h"

namespace itl2
{
	/**
	LU decomposition of a matrix.

	For an m-by-n matrix A with m >= n, the LU decomposition is an m-by-n
	unit lower triangular matrix L, an n-by-n upper triangular matrix U,
	and a permutation vector piv of length m so that A(piv,:) = L*U.
	If m < n, then L is m-by-m and U is m-by-n.
	
	The LU decomposition with pivoting always exists, even if the matrix is
	singular, so the constructor will never fail. The primary use of the
	LU decomposition is in the solution of square systems of simultaneous
	linear equations. This will fail if isNonsingular() returns false.
	*/
	template <class real_t> class LUDecomposition
	{

		/**
		Array for internal storage of decomposition.
		*/
		Matrix<real_t> LU_;

		/**
		Row and column dimensions.
		*/
		size_t m, n; // TODO: These are the same than dimensions of LU_ so these can be eliminated.

		/**
		Pivot sign
		*/
		int pivsign;

		/**
		Internal storage of pivot vector.
		*/
		Matrix<size_t> piv;

	public:

		/**
		Constructor
		@param A Square matrix.
		*/
		LUDecomposition(const Matrix<real_t> &A) :
			LU_(A),
			m(A.rowCount()), n(A.columnCount()),
			piv(A.rowCount())
		{
			if (!A.isSquare())
				throw itl2::ITLException("Cannot calculate LU decomposition for a non-square matrix.");

			// Use a "left-looking", dot-product, Crout/Doolittle algorithm.

			for (size_t i = 0; i < m; i++)
				piv(i) = i;
			
			pivsign = 1;
			Matrix<real_t> LUcolj(m);

			// Outer loop (over columns)
			for (size_t j = 0; j < n; j++)
			{

				// Make a copy of the j-th column to localize references.
				for (size_t i = 0; i < m; i++)
				{
					LUcolj(i) = LU_(i, j);
				}

				// Apply previous transformations.
				for (size_t i = 0; i < m; i++)
				{
					// Most of the time is spent in the following dot product.
					size_t kmax = std::min(i, j);
					real_t s = 0;
					for (size_t k = 0; k < kmax; k++)
					{
						s += LU_(i, k) * LUcolj(k);
					}

					LU_(i, j) = LUcolj(i) -= s;
				}

				// Find pivot and exchange if necessary.
				size_t p = j;
				for (size_t i = j + 1; i < m; i++)
				{
					if (std::abs(LUcolj(i)) > std::abs(LUcolj(p)))
						p = i;
				}

				if (p != j)
				{
					size_t k = 0;
					for (k = 0; k < n; k++)
					{
						real_t t = LU_(p, k);
						LU_(p, k) = LU_(j, k);
						LU_(j, k) = t;
					}
					k = piv(p);
					piv(p) = piv(j);
					piv(j) = k;
					pivsign = -pivsign;
				}

				// Compute multipliers.
				if ((j < m) && (LU_(j, j) != (real_t)0.0)) // TODO: numerical tolerance?
				{
					for (size_t i = j + 1; i < m; i++)
					{
						LU_(i, j) /= LU_(j, j);
					}
				}
			}
		}


		/**
		Gets a value indicating if the matrix is nonsingular?
		@return True if upper triangular factor U (and hence A) is nonsingular, false otherwise.
		*/
		bool isNonsingular() const
		{
			for (size_t j = 0; j < n; j++)
			{
				if (LU_(j, j) == 0) // TODO: numerical tolerance?
					return false;
			}
			return true;
		}

		/**
		Return lower triangular factor.
		*/
		Matrix<real_t> getL() const
		{
			Matrix<real_t> L_(m, n);
			for (size_t i = 0; i < m; i++)
			{
				for (size_t j = 0; j < n; j++)
				{
					if (i > j)
					{
						L_(i, j) = LU_(i, j);
					}
					else if (i == j)
					{
						L_(i, j) = 1;
					}
					else
					{
						L_(i, j) = 0;
					}
				}
			}
			return L_;
		}

		/**
		Return upper triangular factor
		@return U portion of LU factorization.
		*/
		Matrix<real_t> getU() const
		{
			Matrix<real_t> U_(n, n);
			for (size_t i = 0; i < n; i++)
			{
				for (size_t j = 0; j < n; j++)
				{
					if (i <= j)
						U_(i, j) = LU_(i, j);
					else
						U_(i, j) = 0;
				}
			}
			return U_;
		}

		/**
		Return pivot permutation vector
		*/
		Matrix<size_t> getPivot() const
		{
			return piv;
		}


		/**
		Compute determinant using LU factors.
		@return Determinant of A, or 0 if A is not square.
		*/
		real_t det() const
		{
			if (m != n)
			{
				return 0;
			}
			real_t d = (real_t)pivsign;
			for (size_t j = 0; j < n; j++)
			{
				d *= LU_(j, j);
			}
			return d;
		}

		/**
		Solve A*X = B
		@param B Matrix with as many rows as A and any number of columns.
		@return X so that L*U*X = B(piv,:). If B is nonconformant, throws ITLException. The (piv,:) notation means that this method accounts for pivoting internally.
		*/
		Matrix<real_t> solve(const Matrix<real_t> &B) const
		{

			/* Dimensions: A is mxn, X is nxk, B is mxk */

			if (B.rowCount() != m)
				throw itl2::ITLException("Row count of right side matrix B does not equal row count of coefficient matrix A.");
			
			if (!isNonsingular())
				throw itl2::ITLException("Coefficient matrix A is singular.");
			
			// Copy right hand side with pivoting
			size_t nx = B.columnCount();

			Matrix<real_t> X = B.subMatrix(piv, 0, nx - 1);

			// Solve L*Y = B(piv,:)
			for (size_t k = 0; k < n; k++)
			{
				for (size_t i = k + 1; i < n; i++)
				{
					for (size_t j = 0; j < nx; j++)
					{
						X(i, j) -= X(k, j) * LU_(i, k);
					}
				}
			}

			// Solve U*X = Y;
			for (itl2::coord_t k = n - 1; k >= 0; k--)
			{
				for (size_t j = 0; j < nx; j++)
				{
					X(k, j) /= LU_(k, k);
				}
				for (itl2::coord_t i = 0; i < k; i++)
				{
					for (size_t j = 0; j < nx; j++)
					{
						X(i, j) -= X(k, j) * LU_(i, k);
					}
				}
			}

			return X;
		}
	};
}