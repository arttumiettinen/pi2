#pragma once

#include <stdexcept>

#include "math/matrix.h"
#include "math/numberutils.h"
#include "test.h"

namespace math
{
	/**
	Solves Ax = b using conjugate gradientMagnitude method.
	*/
	template<typename real_t = double> void conjugateGradient(const Matrix<real_t>& A, const Matrix<real_t>& b, Matrix<real_t>& x, const real_t tolerance = NumberUtils<real_t>::tolerance(), const size_t maxiter = numeric_limits<size_t>::max())
	{
		if (!A.isSymmetric())
		{
			throw std::runtime_error("Matrix A must be symmetric in conjugate gradient algorithm.");
		}

		Matrix<real_t> r = b - A * x;
		Matrix<real_t> p = r;
		real_t k = 0;

		real_t rsqnormold = r.transpose() * r;

		Matrix<real_t> Ap;

		// Perform at most as many iterations as there are elements in b.
		//maxiter = min(b.count(), maxiter);
		for (size_t k = 0; k < maxiter; k++)
		{
			Ap = A * p;

			real_t alpha = rsqnormold / (p.transpose() * Ap);
			x = x + alpha * p;
			r = r - alpha * Ap;

			real_t rsqnormnew = r.transpose() * r;
			if (sqrt(rsqnormnew) < tolerance)
				break;

			real_t beta = rsqnormnew / rsqnormold;
			p = r + beta * p;
			rsqnormold = rsqnormnew;
		}

	}

	/**
	Solves A^T A x = A^T b using conjugate gradientMagnitude method.
	*/
	template<typename real_t = double> void cgne(const Matrix<real_t>& A, const Matrix<real_t>& b, Matrix<real_t>& x, const real_t tolerance = NumberUtils<real_t>::tolerance(), const size_t maxiter = numeric_limits<size_t>::max())
	{
		Matrix<real_t> At = A.transpose();

		Matrix<real_t> r = At * b - At * A * x;
		Matrix<real_t> p = r;
		real_t k = 0;

		real_t rsqnormold = r.transpose() * r;

		Matrix<real_t> Ap;

		// Perform at most as many iterations as there are elements in b.
		//maxiter = min(b.count(), maxiter);
		for (size_t k = 0; k < maxiter; k++)
		{
			Ap = At * A * p;

			real_t alpha = rsqnormold / (p.transpose() * Ap);
			x = x + alpha * p;
			r = r - alpha * Ap;

			real_t rsqnormnew = r.transpose() * r;
			if (sqrt(rsqnormnew) < tolerance)
				break;

			real_t beta = rsqnormnew / rsqnormold;
			p = r + beta * p;
			rsqnormold = rsqnormnew;
		}

	}


	namespace tests
	{
		/**
		Test for conjugate gradientMagnitude algorithm.
		*/
		inline void conjugateGradient()
		{
			Matrix<double> A(2, 2);
			A(0, 0) = 4;
			A(1, 0) = 1;
			A(0, 1) = 1;
			A(1, 1) = 3;

			Matrix<double> b(2, 1);
			b(0) = 1;
			b(1) = 2;

			Matrix<double> x(2, 1);
			x(0) = 0;
			x(1) = 0;

			conjugateGradient(A, b, x);

			Matrix<double> xtrue(2, 1);
			xtrue(0) = 1.0 / 11.0;
			xtrue(1) = 7.0 / 11.0;

			itl2::testAssert(x == xtrue, "Conjugate gradient");

		}

		inline void cgne()
		{
			Matrix<double> A(3, 2);
			A(0, 0) = 1;
			A(0, 1) = 2;
			A(1, 0) = 3;
			A(1, 1) = 4;
			A(2, 0) = 1.1;
			A(2, 1) = 1.9;

			Matrix<double> b(3, 1);
			b(0) = 3;
			b(1) = 4;
			b(2) = 3.1;

			Matrix<double> x(2, 1);
			x(0) = 0;
			x(1) = 0;

			cgne(A, b, x);

			Matrix<double> xtrue(2, 1);
			xtrue(0) = -2.437716262975772;
			xtrue(1) = 2.842560553633213;

			itl2::testAssert(x.equals(xtrue, 1e-5), "CGNE result");
		}
	}

}
