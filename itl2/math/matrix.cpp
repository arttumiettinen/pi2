
#include "matrix.h"
#include "test.h"

#include <iostream>
#include <iomanip>

using namespace std;

namespace itl2
{
	namespace tests
	{
		void matrix()
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


		void solve()
		{
			Matrix A({ {0.8311,   0.8351,    0.0349,    0.1548,    0.4018},
						{0.1565,    0.8954,    0.8854,    0.1439,    0.4064},
						{0.4573,    0.5825,    0.4077,    0.6060,    0.3862},
						{0.6181,    0.5827,    0.0364,    0.2545,    0.6098},
						{0.9322,    0.8549,    0.7461,    0.3242,    0.1669} });

			Matrix AinvGT = { {-0.9956, -0.9328,		-1.0209,	1.4755,		1.6394},
								{3.2183,	0.8507,		1.0545,		-2.9191,	-1.5939 },
								{-2.1758,	0.3661,		-0.8454,	1.3408,		1.4041},
								{0.4845,	-0.5868,	2.6242,		-1.4007,	-0.6920},
								{-2.1385,	0.3557,		-1.0176,	3.4382,		0.0663} };

			cout << "A = " << endl;
			cout << A << endl;

			cout << "A^-1 ground truth = " << endl;
			cout << AinvGT << endl;

			Matrix Ainv = A.inverse();

			cout << "A^-1 = " << endl;
			cout << Ainv << endl;

			cout << "A * A^-1 = " << endl;
			cout << fixed << setprecision(4) << A * Ainv << endl;

			itl2::testAssert(Ainv.equals(AinvGT, 0.0001), "A inverse vs ground truth"); // NOTE: 0.0001 tolerance as we only store that many decimals in the ground truth above!
			itl2::testAssert(A * Ainv == Matrix<double>::identity(A), "A times inverse of A");
			itl2::testAssert(Ainv * A == Matrix<double>::identity(A), "inverse of A times A");

		}

		void leastSquares()
		{
			Matrix A = { {0.1881,    0.7404,    0.3094},
						{0.0946,    0.6928,    0.5230},
						{0.3232,    0.8241,    0.3253},
						{0.7696,    0.8280,    0.8318},
						{0.2341,    0.2934,    0.8103} };

			Matrix B = { {1.0},
							{2.0},
							{3.0},
							{4.0},
							{5.0} };

			Matrix XGT = { {0.4340},
						{-0.2743},
						{5.3453} };

			Matrix X = A.solve(B);

			cout << "A * X = B" << endl;
			cout << A << endl;
			cout << "*" << endl;
			cout << "[[X1]," << endl;
			cout << " [X2]," << endl;
			cout << " [X3]]" << endl;
			cout << "=" << endl;
			cout << B << endl;
			cout << "X = " << endl;
			cout << X << endl;
			cout << "X ground truth = " << endl;
			cout << XGT << endl;

			itl2::testAssert(X.equals(XGT, 0.0001), "least squares solution vs ground truth"); // NOTE: 0.0001 tolerance as we only store that many decimals in the ground truth above!
		}
	}
}