#pragma once

/*
This code is based on public-domain reference implementation of JAMA : A Java Matrix Package (https://math.nist.gov/javanumerics/jama/)
and public-domain library TNT: Template Numerical Toolkit and its JAMA/C++ part (https://math.nist.gov/tnt/).
*/

#include "math/matrix.h"
#include "math/mathutils.h"

namespace itl2
{

	/**
	Classical QR Decompisition:
	for an m-by-n matrix A with m >= n, the QR decomposition is an m-by-n
	orthogonal matrix Q and an n-by-n upper triangular matrix R so that
	A = Q*R.
	
	The QR decompostion always exists, even if the matrix does not have
	full rank, so the constructor will never fail. The primary use of the
	QR decomposition is in the least squares solution of nonsquare systems
	of simultaneous linear equations. This will fail if isFullRank()
	returns false.

	The Q and R factors can be retrived via the getQ() and getR()
	methods. Furthermore, a solve() method is provided to find the
	least squares solution of Ax=b using the QR factors.
	*/
	template <class real_t> class QRDecomposition
	{
		/**
		Array for internal storage of decomposition.
		*/
		Matrix<real_t> QR_;

		/**
		Row and column dimensions.
		*/
		size_t m, n;

		/**
		Array for internal storage of diagonal of R.
		*/
		Matrix<real_t> Rdiag;


	public:

		/**
		Constructor
		@param A rectangular (m>=n) matrix.
		*/
		QRDecomposition(const Matrix<real_t> &A) :
			QR_(A),
			m(A.rowCount()),
			n(A.columnCount()),
			Rdiag(A.columnCount())
		{
			for (size_t k = 0; k < n; k++)
			{
				// Compute 2-norm of k-th column without under/overflow.
				real_t nrm = 0;
				for (size_t i = k; i < m; i++)
				{
					nrm = hypot(nrm, QR_(i, k));
				}

				if (nrm != 0.0)
				{
					// Form k-th Householder vector.
					if (QR_(k, k) < 0)
					{
						nrm = -nrm;
					}
					for (size_t i = k; i < m; i++)
					{
						QR_(i, k) /= nrm;
					}
					QR_(k, k) += 1.0;

					// Apply transformation to remaining columns.
					for (size_t j = k + 1; j < n; j++)
					{
						real_t s = 0.0;
						for (size_t i = k; i < m; i++)
						{
							s += QR_(i, k) * QR_(i, j);
						}
						s = -s / QR_(k, k);
						for (size_t i = k; i < m; i++) 
						{
							QR_(i, j) += s * QR_(i, k);
						}
					}
				}
				Rdiag(k) = -nrm;
			}
		}

		/**
		Gets a value indicating if the matrix where the decomposition was calculate from is full rank.
		*/
		bool isFullRank() const
		{
			for (size_t j = 0; j < n; j++)
			{
				if (Rdiag(j) == 0) // TODO: numerical tolerance?
					return false;
			}
			return true;
		}

		/**
		Retrieve the Householder vectors from QR factorization.
		@returns lower trapezoidal matrix whose columns define the reflections
		*/
		Matrix<real_t> getHouseholder() const
		{
			Matrix<real_t> H(m, n);

			// NOTE: H is completely filled in by algorithm, so
			// initialization of H is not necessary.
			for (size_t i = 0; i < m; i++)
			{
				for (size_t j = 0; j < n; j++)
				{
					if (i >= j)
						H(i, j) = QR_(i, j);
					else
						H(i, j) = 0;
				}
			}
			return H;
		}

		/**
		Return the upper triangular factor, R, of the QR factorization.
		*/
		Matrix<real_t> getR() const
		{
			Matrix<real_t> R(n, n);
			for (size_t i = 0; i < n; i++)
			{
				for (size_t j = 0; j < n; j++)
				{
					if (i < j)
						R(i, j) = QR_(i, j);
					else if (i == j)
						R(i, j) = Rdiag(i);
					else
						R(i, j) = 0;
				}
			}
			return R;
		}

		/**
		Generate and return the (economy-sized) orthogonal factor.
		@param Q the (ecnomy-sized) orthogonal factor (Q*R=A).
		*/
		Matrix<real_t> getQ() const
		{
			Matrix<real_t> Q(m, n);
			for (coord_t k = n - 1; k >= 0; k--) {
				
				for (size_t i = 0; i < m; i++)
					Q(i, k) = 0;

				Q(k, k) = 1;

				for (size_t j = k; j < n; j++)
				{
					if (QR_(k, k) != 0)
					{
						real_t s = 0;
						for (size_t i = k; i < m; i++)
							s += QR_(i, k) * Q(i, j);

						s = -s / QR_(k, k);

						for (size_t i = k; i < m; i++)
							Q(i, j) += s * QR_(i, k);
					}
				}
			}
			return Q;
		}


		///** Least squares solution of A*x = b
		//@param B     m-length array (vector).
		//@return x    n-length array (vector) that minimizes the two norm of Q*R*X-B.
		//	 If B is non-conformant, or if QR.isFullRank() is false,
		//					 the routine returns a null (0-length) vector.
		//*/
		//TNT::Array1D<real_t> solve(const TNT::Array1D<real_t> &b) const
		//{
		//	if (b.rowCount() != m)		/* arrays must be conformant */
		//		return TNT::Array1D<real_t>();

		//	if (!isFullRank())		/* matrix is rank deficient */
		//	{
		//		return TNT::Array1D<real_t>();
		//	}

		//	TNT::Array1D<real_t> x = b.copy();

		//	// Compute Y = transpose(Q)*b
		//	for (int k = 0; k < n; k++)
		//	{
		//		real_t s = 0.0;
		//		for (int i = k; i < m; i++)
		//		{
		//			s += QR_(i, k) * x(i);
		//		}
		//		s = -s / QR_(k, k);
		//		for (int i = k; i < m; i++)
		//		{
		//			x(i) += s * QR_(i, k);
		//		}
		//	}
		//	// Solve R*X = Y;
		//	for (int k = n - 1; k >= 0; k--)
		//	{
		//		x(k) /= Rdiag(k);
		//		for (int i = 0; i < k; i++) {
		//			x(i) -= x(k) * QR_(i, k);
		//		}
		//	}


		//	/* return n x nx portion of X */
		//	TNT::Array1D<real_t> x_(n);
		//	for (int i = 0; i < n; i++)
		//		x_(i) = x(i);

		//	return x_;
		//}

		/**
		Least squares solution of A*X = B
		@param B m x k array.
		@return X n x k Array that minimizes the two norm of Q*R*X-B. If B is non-conformant, or if QR.isFullRank() is false, the routine throws an exception.
		*/
		Matrix<real_t> solve(const Matrix<real_t> &B) const
		{
			if (B.rowCount() != m)		/* arrays must be conformant */
				throw itl2::ITLException("Row count of right side matrix B does not equal row count of coefficient matrix A.");

			if (!isFullRank())
				throw itl2::ITLException("Matrix A is rank deficient, there is no solution.");


			size_t nx = B.columnCount();
			Matrix<real_t> X = B;

			// Compute Y = transpose(Q)*B
			for (size_t k = 0; k < n; k++)
			{
				for (size_t j = 0; j < nx; j++)
				{
					real_t s = 0;
					for (size_t i = k; i < m; i++)
						s += QR_(i, k) * X(i, j);
					
					s = -s / QR_(k, k);

					for (size_t i = k; i < m; i++)
						X(i, j) += s * QR_(i, k);
				}
			}

			// Solve R*X = Y
			for (itl2::coord_t k = n - 1; k >= 0; k--)
			{
				for (size_t j = 0; j < nx; j++)
					X(k, j) /= Rdiag(k);

				for (itl2::coord_t i = 0; i < k; i++)
				{
					for (size_t j = 0; j < nx; j++)
					{
						X(i, j) -= X(k, j) * QR_(i, k);
					}
				}
			}

			// return n x nx portion of X
			Matrix<real_t> X_(n, nx);
			for (size_t i = 0; i < n; i++)
			{
				for (size_t j = 0; j < nx; j++)
					X_(i, j) = X(i, j);
			}

			return X_;
		}


	};


}