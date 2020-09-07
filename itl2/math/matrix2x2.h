#pragma once

#include <string>
#include <sstream>
#include <cstring>

#include "math/vec2.h"
#include "math/eig2.h"
#include "math/numberutils.h"


namespace itl2
{

	/**
	2x2 matrix
	@param T data type used to store values.
	*/
	template<typename T> class Matrix2x2
	{
	public:
		union
		{
			struct
			{
				T a00, a01,
					a10, a11;
			};
			T t[2][2];
			T e[4];
		};

		/**
		Constructor.
		*/
		Matrix2x2<T>()
		{
			zero();
		}

		/**
		Constructor
		*/
		Matrix2x2<T>(T a00, T a01, T a10, T a11) :
			a00(a00), a01(a01),
			a10(a10), a11(a11)
		{
		}

		/**
		Constructor that constructs matrix from column vectors.
		*/
		Matrix2x2<T>(const Vec2<T>& c1, const Vec2<T>& c2) :
			a00(c1.x), a01(c2.x),
			a10(c1.y), a11(c2.y)
		{
		}


		/**
		Creates a rotation matrix that rotates around given axis by given amount in radians.
		*/
		static Matrix2x2<T> rotationMatrix(const T angle)
		{
			// See https://en.wikipedia.org/wiki/Rotation_matrix
			T c = cos(angle);
			T s = sin(angle);
			return Matrix2x2<T>(c, -s, s, c);
		}

		/**
		Returns identity matrix.
		*/
		static Matrix2x2<T> identity()
		{
			Matrix2x2<T> mat;
			mat.loadIdentity();
			return mat;
		}

		/**
		Load identity matrix into this matrix.
		*/
		void loadIdentity()
		{
			a00 = 1; a01 = 0;
			a10 = 0; a11 = 1;
		}

		/**
		Sets all elements to zero.
		*/
		void zero()
		{
			a00 = 0; a01 = 0;
			a10 = 0; a11 = 0;
		}

		/**
		Gets total number of elements in this matrix.
		*/
		size_t size() const
		{
			return 2 * 2;
		}

		/**
		Compute determinant of this matrix.
		*/
		T det() const
		{
			return a00 * a11 - a10 * a01;
		}

		/**
		Invert this matrix.
		@return True if this matrix is invertible, false if this matrix is singular. This matrix is not modified if it is singular.
		*/
		bool invert()
		{
			return inverse(*this);
		}

		/**
		Calculate inverse of this matrix.
		@param invTarget The inverse is placed to this matrix.
		@return True if this matrix is invertible, false if this matrix is singular.
		*/
		bool inverse(Matrix2x2<T>& invTarget) const
		{
			T determinant = det();
			if (NumberUtils<T>::equals(determinant, 0))
				return false;

			invTarget = Matrix2x2<T>(a11, -a01, -a10, a00) / determinant;

			return true;
		}

		/**
		Transposes this matrix.
		*/
		void transpose()
		{
			std::swap(a01, a10);
		}




		/**
		Multiply this matrix from right by a 2x1 matrix.
		@param right Right operand matrix
		@return Multiplication result matrix. (2x1)
		*/
		Vec2<T> operator *(const Vec2<T>& r) const
		{
			return Vec2<T>(t[0][0] * r[0] + t[0][1] * r[1],
				t[1][0] * r[0] + t[1][1] * r[1]);
		}





		/**
		Multiplication by constant in-place.
		*/
		Matrix2x2& operator*=(const T r)
		{
			for (size_t n = 0; n < size(); n++)
				e[n] *= r;
			return *this;
		}

		/**
		Multiplication by constant.
		*/
		Matrix2x2 operator*(const T r) const
		{
			Matrix2x2 result(*this);
			return result *= r;
		}

		/**
		Multiplication by constant.
		*/
		friend Matrix2x2 operator*(const T c, const Matrix2x2& r)
		{
			return r * c;
		}




		/**
		Division by constant in-place.
		*/
		Matrix2x2& operator/=(const T r)
		{
			for (size_t n = 0; n < size(); n++)
				e[n] /= r;
			return *this;
		}

		/**
		Division by constant.
		*/
		Matrix2x2 operator/(const T r) const
		{
			Matrix2x2 result(*this);
			return result /= r;
		}





		/**
		Addition in-place.
		*/
		Matrix2x2& operator+=(const Matrix2x2& r)
		{
			for (size_t n = 0; n < size(); n++)
				e[n] += r.e[n];
			return *this;
		}

		/**
		Addition operator.
		*/
		Matrix2x2 operator+(const Matrix2x2& r) const
		{
			Matrix2x2 result(*this);
			return result += r;
		}



		/**
		Subtraction in-place.
		*/
		Matrix2x2& operator-=(const Matrix2x2& r)
		{
			for (size_t n = 0; n < size(); n++)
				e[n] -= r.e[n];
			return *this;
		}

		/**
		Subtract operator.
		*/
		Matrix2x2 operator-(const Matrix2x2& r) const
		{
			Matrix2x2 result(*this);
			return result -= r;
		}




		/**
		Calculates trace of this matrix.
		*/
		T trace() const
		{
			return a00 + a11;
		}


		/**
		Creates a matrix as outer product of the two vectors, i.e. u * v.transpose()
		*/
		static Matrix2x2<T> outer(const Vec2<T>& u, const Vec2<T>& v)
		{
			return Matrix2x2<T>(u.x * v.x, u.x * v.y,
				u.y * v.x, u.y * v.y);
		}

		/**
		Tests if this matrix is symmetric.
		*/
		bool isSymmetric(T tolerance = NumberUtils<T>::tolerance()) const
		{
			for (size_t i = 0; i < 2; i++)
			{
				for (size_t j = 0; j < 2; j++)
				{
					if (!NumberUtils<T>::equals(t[i][j], t[j][i], tolerance))
						return false;
				}
			}

			return true;
		}

		/**
		Calculates eigendecomposition of this matrix, assuming it is symmetric, and stores sorted eigenvalues in
		lambda1 and lambda2, lambda1 corresponding to the largest eigenvalue.
		The corresponding eigenvectors are stored in v1 and v2.
		*/
		void eigsym(Vec2<T>& v1, Vec2<T>& v2, T& lambda1, T& lambda2) const
		{
			double A[2][2];
			double V[2][2];
			double d[2];

			for (size_t i = 0; i < 2; i++)
				for (size_t j = 0; j < 2; j++)
					A[i][j] = t[i][j];

			eigen_decomposition2(A, V, d);

			
			lambda1 = (T)d[1];
			lambda2 = (T)d[0];

			T V1mult = V[0][1] < 0 ? (T)-1 : (T)1;
			T V2mult = V[0][0] < 0 ? (T)-1 : (T)1;

			v1 = Vec2<T>((T)V[0][1], (T)V[1][1]) * V1mult;
			v2 = Vec2<T>((T)V[0][0], (T)V[1][0]) * V2mult;

			// Sort eigenvalues and eigenvectors
			if (lambda1 < lambda2)
			{
				swap(lambda1, lambda2);
				swap(v1, v2);
			}

		}

		/**
		Converts this object to string.
		*/
		friend std::ostream& operator<<(std::ostream& stream, const Matrix2x2<T>& v)
		{
			stream << "[ [" << v.a00 << ", " << v.a01  << "] [" <<
				v.a10 << ", " << v.a11  << "] ]";
			return stream;
		}

	};

	typedef Matrix2x2<float> Matrix2x2f;
	typedef Matrix2x2<double> Matrix2x2d;

}

