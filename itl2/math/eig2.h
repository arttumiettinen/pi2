#pragma once

#include <iostream>

#include <cmath>
#include "math/mathutils.h"
#include "math/numberutils.h"

namespace math
{
	/**
	Calculate eigendecomposition of 2x2 real matrix A.
	Eigenvector i will be in V[n][i] and eigenvalue in L[i].
	See http://www.math.harvard.edu/archive/21b_fall_04/exhibits/2dmatrices/index.html
	@return true if the eigenvalues are real; false otherwise. If false is returned, the decomposition is not computed.
	*/
	inline bool eigen_decomposition2(double A[2][2], double V[2][2], double L[2])
	{
		double a = A[0][0];
		double b = A[0][1];
		double c = A[1][0];
		double d = A[1][1];

		double T = a + d;
		double D = a * d - b * c;

		double discr2 = T * T / 4.0 - D;
		if(discr2 <= 0.0)
			return false;

		double discr = sqrt(discr2);
		double L1 = T / 2.0 + discr;
		double L2 = T / 2.0 - discr;

		L[0] = L1;
		L[1] = L2;

		if(!NumberUtils<double>::equals(c, 0.0))
		{
			V[0][0] = L1 - d;
			V[1][0] = c;
			V[0][1] = L2 - d;
			V[1][1] = c;
		}
		else if(!NumberUtils<double>::equals(b, 0.0))
		{
			V[0][0] = b;
			V[1][0] = L1 - a;
			V[0][1] = b;
			V[1][1] = L2 - a;
		}
		else
		{
			// Eigenvectors are in coordinate directions.
			// Put the first vector point towards maximum deviation.
			if(a >= d)
			{
				V[0][0] = 1;
				V[1][0] = 0;
				V[0][1] = 0;
				V[1][1] = 1;
			}
			else
			{
				V[0][1] = 1;
				V[1][1] = 0;
				V[0][0] = 0;
				V[1][0] = 1;
			}
		}

		if(L[0] > L[1])
		{
			std::swap(L[0], L[1]);
			std::swap(V[0][0], V[0][1]);
			std::swap(V[1][0], V[1][1]);
		}

		return true;
	}
}
