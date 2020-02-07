#pragma once

#include <string>
#include <sstream>
#include <cstring>

#include "math/vec3.h"
#include "math/dsyevh3.h"
#include "math/numberutils.h"


namespace itl2
{

	/**
	3x3 matrix
	@param T data type used to store values.
	*/
	template<typename T> class Matrix3x3
	{
	public:
		union
		{
			struct
			{
				T a00, a01, a02,
				  a10, a11, a12,
				  a20, a21, a22;
			};
			T t[3][3];
			T e[9];
		};

		/**
		Constructor. Initializes the matrix to the identity matrix.
		TODO: Identity initialization is probably bad choice.
		*/
		Matrix3x3<T>()
		{
			//loadIdentity();
			zero();
		}

		/**
		Constructor
		*/
		Matrix3x3<T>(T a00, T a01, T a02, T a10, T a11, T a12, T a20, T a21, T a22) :
			a00(a00), a01(a01), a02(a02),
			a10(a10), a11(a11), a12(a12),
			a20(a20), a21(a21), a22(a22)
		{
		}

		/**
		Constructor that constructs matrix from three column vectors.
		*/
		Matrix3x3<T>(const Vec3<T>& c1, const Vec3<T>& c2, const Vec3<T>& c3) :
			a00(c1.x), a01(c2.x), a02(c3.x),
			a10(c1.y), a11(c2.y), a12(c3.y),
			a20(c1.z), a21(c2.z), a22(c3.z)
		{
		}


		/**
		Creates a rotation matrix that rotates around given axis by given amount in radians.
		*/
		static Matrix3x3<T> rotationMatrix(const T angle, Vec3<T> axis)
		{
			// See https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
			axis.normalize();
			T c = cos(angle);
			T s = sin(angle);
			return c * Matrix3x3<T>::identity() + s * Matrix3x3<T>::crossProductMatrix(axis) + (1 - c) * Matrix3x3<T>::outer(axis, axis);
		}

		/**
		Creates a rotation matrix that rotates (1, 0, 0) to u1, (0, 1, 0) to u2, and (0, 0, 1) to u1 x u2.
		If u1 and u2 are not orthogonal, u2 is adjusted such that all three vectors u1, u2, and u1 x u2 are orthogonal.
		The vectors u1 and u2 do not have to be unit vectors (they are normalized automatically).
		*/
		static Matrix3x3<T> rotationMatrix(const Vec3<T>& u1, const Vec3<T>& u2)
		{
			// Normalize the axes to avoid stupid errors with non-unit axis vectors.
			Vec3<T> u1n = u1.normalized();
			Vec3<T> u2n = u2.normalized();

			// Calculate the third axis
			Vec3<T> u3n = u1n.cross(u2n).normalized();

			// This makes sure that all the three axes are orthogonal to each other.
			u2n = u3n.cross(u1n);

			return Matrix3x3<T>(u1n.x, u2n.x, u3n.x,
				u1n.y, u2n.y, u3n.y,
				u1n.z, u2n.z, u3n.z);
		}

		/**
		Returns identity matrix.
		*/
		static Matrix3x3<T> identity()
		{
			Matrix3x3<T> mat;
			mat.loadIdentity();
			return mat;
		}

		/**
		Constructor
		*/
		//Matrix3x3<T>(const Matrix3x3<T>& m)
		//{
		//	set(m);
		//}

		/**
		Set the values of the elements of this matrix to those of the given matrix m.
		*/
		//void set(const Matrix3x3<T>& m)
		//{
		//	memcpy(this->t, m.t, sizeof(t));
		//}

		/**
		Load identity matrix into this matrix.
		*/
		void loadIdentity()
		{
			a00=1; a01=0; a02=0;
			a10=0; a11=1; a12=0;
			a20=0; a21=0; a22=1;
		}

		/**
		Sets all elements to zero.
		*/
		void zero()
		{
			a00 = 0; a01 = 0; a02 = 0;
			a10 = 0; a11 = 0; a12 = 0;
			a20 = 0; a21 = 0; a22 = 0;
		}

		/**
		Gets total number of elements in this matrix.
		*/
		size_t size() const
		{
			return 3 * 3;
		}

		/**
        Calculate skew symmetric cross product matrix from this vector.
        This method creates a matrix so that if this vector is a and
        some other vector is b,
        a.Cross(b) = a.CrossProductMatrix() * b
        */
        static Matrix3x3 crossProductMatrix(const Vec3<T>& v)
        {
            return Matrix3x3(0.0,  -v.z,    v.y,
                            v.z,    0.0,    -v.x,
                            -v.y,   v.x,    0.0);
        }

		/**
		Compute determinant of this matrix.
		*/
		T det() const
		{
			return t[0][0] * (t[1][1] * t[2][2] - t[1][2] * t[2][1])
					- t[0][1] * (t[1][0] * t[2][2] - t[1][2] * t[2][0])
					+ t[0][2] * (t[1][0] * t[2][1] - t[1][1] * t[2][0]);
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
		bool inverse(Matrix3x3<T> &invTarget) const
		{
			T determinant = det();
			if(NumberUtils<T>::equals(determinant, 0))
				return false;

			T invDet = 1 / determinant;

			//Make a temporary matrix so that one can safely invert this matrix.
			Matrix3x3<T> inv;

			inv.t[0][0] = invDet * (t[1][1] * t[2][2] - t[1][2] * t[2][1]);
			inv.t[0][1] = invDet * (t[0][2] * t[2][1] - t[0][1] * t[2][2]);
			inv.t[0][2] = invDet * (t[0][1] * t[1][2] - t[0][2] * t[1][1]);
			inv.t[1][0] = invDet * (t[1][2] * t[2][0] - t[1][0] * t[2][2]);
			inv.t[1][1] = invDet * (t[0][0] * t[2][2] - t[0][2] * t[2][0]);
			inv.t[1][2] = invDet * (t[0][2] * t[1][0] - t[0][0] * t[1][2]);
			inv.t[2][0] = invDet * (t[1][0] * t[2][1] - t[1][1] * t[2][0]);
			inv.t[2][1] = invDet * (t[0][1] * t[2][0] - t[0][0] * t[2][1]);
			inv.t[2][2] = invDet * (t[0][0] * t[1][1] - t[0][1] * t[1][0]);

			invTarget = inv;

			return true;
		}

		/**
		Transposes this matrix.
		*/
		void transpose()
		{
			std::swap(a01, a10);
			std::swap(a02, a20);
			
			std::swap(a12, a21);
		}




		/**
		Multiply this matrix from right by a 3x1 matrix.
		@param right Right operand matrix
		@return Multiplication result matrix. (3x1)
		*/
		Vec3<T> operator *(const Vec3<T>& r) const
		{
			return Vec3<T>(t[0][0] * r[0] + t[0][1] * r[1] + t[0][2] * r[2],
								t[1][0] * r[0] + t[1][1] * r[1] + t[1][2] * r[2],
								t[2][0] * r[0] + t[2][1] * r[1] + t[2][2] * r[2]);
		}





		/**
		Multiplication by constant in-place.
		*/
		Matrix3x3& operator*=(const T r)
		{
			for (size_t n = 0; n < size(); n++)
				e[n] *= r;
			return *this;
		}

		/**
		Multiplication by constant.
		*/
		Matrix3x3 operator*(const T r) const
		{
			Matrix3x3 result(*this);
			return result *= r;
		}

		/**
		Multiplication by constant.
		*/
		friend Matrix3x3 operator*(const T c, const Matrix3x3& r)
		{
			return r * c;
		}




		/**
		Division by constant in-place.
		*/
		Matrix3x3& operator/=(const T r)
		{
			for (size_t n = 0; n < size(); n++)
				e[n] /= r;
			return *this;
		}

		/**
		Division by constant.
		*/
		Matrix3x3 operator/(const T r) const
		{
			Matrix3x3 result(*this);
			return result /= r;
		}





		/**
		Addition in-place.
		*/
		Matrix3x3& operator+=(const Matrix3x3& r)
		{
			for (size_t n = 0; n < size(); n++)
				e[n] += r.e[n];
			return *this;
		}

		/**
		Addition operator.
		*/
		Matrix3x3 operator+(const Matrix3x3& r) const
		{
			Matrix3x3 result(*this);
			return result += r;
		}



		/**
		Subtraction in-place.
		*/
		Matrix3x3& operator-=(const Matrix3x3& r)
		{
			for (size_t n = 0; n < size(); n++)
				e[n] -= r.e[n];
			return *this;
		}

		/**
		Subtract operator.
		*/
		Matrix3x3 operator-(const Matrix3x3& r) const
		{
			Matrix3x3 result(*this);
			return result -= r;
		}




		/**
		Calculates trace of this matrix.
		*/
		T trace() const
		{
			return a00 + a11 + a22;
		}

		/**
		Calculates x^T * M * y,
		where M is this matrix,
		assuming x and y are column vectors
		(or x * M * y^T if x and y are row vectors).
		*/
		T bilinear(const Vec3<T>& x, const Vec3<T>& y) const
		{
			T x0 = x[0];
			T x1 = x[1];
			T x2 = x[2];
			Vec3<T> tmp(x0 * a00 + x1 * a10 + x2 * a20,
				x0 * a01 + x1 * a11 + x2 * a21,
				x0 * a02 + x1 * a12 + x2 * a22);
			return tmp.dot(y);
		}

		/**
		Creates a matrix as outer product of the two vectors, i.e. u * v.transpose()
		*/
		static Matrix3x3<T> outer(const Vec3<T>& u, const Vec3<T>& v)
		{
			return Matrix3x3<T>(u.x * v.x, u.x * v.y, u.x * v.z,
								u.y * v.x, u.y * v.y, u.y * v.z, 
								u.z * v.x, u.z * v.y, u.z * v.z);
		}

		/**
		Tests if this matrix is symmetric.
		*/
		bool isSymmetric(T tolerance = NumberUtils<T>::tolerance()) const
		{
			for (size_t i = 0; i < 3; i++)
			{
				for (size_t j = 0; j < 3; j++)
				{
					if (!NumberUtils<T>::equals(t[i][j], t[j][i], tolerance))
						return false;
				}
			}

			return true;
		}

		/**
		Calculates eigendecomposition of this matrix, assuming it is symmetric, and stores sorted eigenvalues in
		lambda1, lambda2 and lambda3, lambda1 corresponding to the largest eigenvalue.
		The corresponding eigenvectors are stored in v1, v2 and v3.
		*/
		void eigsym(Vec3<T>& v1, Vec3<T>& v2, Vec3<T>& v3, T& lambda1, T& lambda2, T& lambda3) const
		{
			double A[3][3];
			double V[3][3];
			double d[3];

			for (size_t i = 0; i < 3; i++)
				for (size_t j = 0; j < 3; j++)
					A[i][j] = t[i][j];

			dsyevh3(A, V, d);
			
			// Sort eigenvalues and eigenvectors
			constexpr int n = 3;
			for (int i = 0; i < n - 1; i++)
			{
				int k = i;
				double p = d[i];
				for (int j = i + 1; j < n; j++)
				{
					if (d[j] < p)
					{
						k = j;
						p = d[j];
					}
				}
				if (k != i)
				{
					d[k] = d[i];
					d[i] = p;
					for (int j = 0; j < n; j++)
					{
						p = V[j][i];
						V[j][i] = V[j][k];
						V[j][k] = p;
					}
				}
			}

			lambda1 = (T)d[2];
			lambda2 = (T)d[1];
			lambda3 = (T)d[0];

			// Sanity check
			if (lambda1 < lambda2 || lambda2 < lambda3 || lambda1 < lambda3)
				throw itl2::ITLException("Eigenvalues are unsorted.");

			T V1mult = V[0][2] < 0 ? (T)-1 : (T)1;
			T V2mult = V[0][1] < 0 ? (T)-1 : (T)1;
			T V3mult = V[0][0] < 0 ? (T)-1 : (T)1;
			
			v1 = Vec3<T>((T)V[0][2], (T)V[1][2], (T)V[2][2]) * V1mult;
			v2 = Vec3<T>((T)V[0][1], (T)V[1][1], (T)V[2][1]) * V2mult;
			v3 = Vec3<T>((T)V[0][0], (T)V[1][0], (T)V[2][0]) * V3mult;
		}

		/**
        Converts this object to string.
        */
        friend std::ostream& operator<<(std::ostream& stream, const Matrix3x3<T>& v)
        {
            stream << "[ [" << v.a00 << ", " << v.a01 << ", " << v.a02 << "] [" <<
				               v.a10 << ", " << v.a11 << ", " << v.a12 << "] [" <<
				               v.a20 << ", " << v.a21 << ", " << v.a22 << "] ]";
            return stream;
        }

	};

	typedef Matrix3x3<float> Matrix3x3f;
	typedef Matrix3x3<double> Matrix3x3d;

	namespace tests
	{
		void matrix3x3();
	}
}

