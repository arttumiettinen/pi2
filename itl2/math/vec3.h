#pragma once

#include <sstream>
#include <vector>
#include "datatypes.h"
#include "itlexception.h"
#include "math/numberutils.h"
#include "math/mathutils.h"

namespace itl2
{
    /**
    Three-component vector.
    */
    template <typename T> class Vec3
    {
        public:
            /**
            Components
            */
			union
			{
				struct
				{
					T x;
					T y;
					T z;
				};
				T components[3];
			};

            /**
            Default constructor, initializes the vector to zero.
            */
            Vec3() :
                x(0),
                y(0),
                z(0)
            {
            }

            /**
            Constructor
            */
            Vec3(T x, T y, T z)
            {
                this->x = x;
                this->y = y;
                this->z = z;
            }

			/**
			Constructor that makes it possible to cast/initialize vector from vector of another type.
			*/
			//template<typename Tother> explicit Vec3(const Vec3<Tother>& other) :
			//	x(pixelRound<T, Tother>(other.x)),
			//	y(pixelRound<T, Tother>(other.y)),
			//	z(pixelRound<T, Tother>(other.z))
			//{
			//}

			/**
			Constructs Vec3 from the first three elements of an array type. If there are not enough elements in the array, the remaining elements are set to zero.
			*/
			template<typename VEC> explicit Vec3(const VEC& v)
			{
				x = 0;
				y = 0;
				z = 0;
				for (size_t n = 0; n < std::min<size_t>(3, v.size()); n++)
					components[n] = pixelRound<T>(v[n]);
			}

			/**
			Array access operator
			*/
			T& operator[] (size_t index)
			{
				return components[index];
			}

			/**
			Array access operator
			*/
			const T& operator[] (size_t index) const
			{
				return components[index];
			}

			/**
			Get count of components.
			*/
			size_t size() const
			{
				return 3;
			}

            /**
            Addition in-place
            */
            Vec3& operator+=(const Vec3& r)
            {
                x += r.x;
                y += r.y;
                z += r.z;
                return *this;
            }

            /**
            Subtraction in-place
            */
            Vec3& operator-=(const Vec3& r)
            {
                x -= r.x;
                y -= r.y;
                z -= r.z;
                return *this;
            }

			/**
            Multiplication in-place
            */
			Vec3& operator*=(const T r)
            {
                x *= r;
                y *= r;
				z *= r;
                return *this;
            }

            /**
            Division in-place
            */
            Vec3& operator/=(const T r)
            {
                x /= r;
                y /= r;
				z /= r;
                return *this;
            }

			/**
			Comparison of greater or equal in all components
			*/
			bool operator>=(const Vec3& r)
			{
				return x >= r.x && y >= r.y && z >= r.z;
			}


            /**
            Addition of constant in-place
            */
            /*
            Vec3& operator+=(const T c)
            {
                x += c;
                y += c;
                z += c;
                return *this;
            }
            */


            /**
            Subtraction of constant in-place
            */
            /*
            Vec3& operator-=(const T c)
            {
                x -= c;
                y -= c;
                z -= c;
                return *this;
            }
            */

            /**
            Unary minus
            */
            Vec3 operator-() const
            {
                return Vec3(-x, -y, -z);
            }

            /**
            Binary addition
            */
            Vec3 operator+(const Vec3 &r) const
            {
                return Vec3(*this) += r;
            }

            /**
            Binary subtraction
            */
            Vec3 operator-(const Vec3 &r) const
            {
                return Vec3(*this) -= r;
            }

            /**
            Multiplication by constant.
            */
            Vec3 operator*(const T c) const
            {
                return Vec3(x * c, y * c, z * c);
            }

            /**
            Multiplication by constant.
            */
            friend Vec3 operator*(const T c, const Vec3& r)
            {
                return r * c;
            }

            /**
            Division by constant.
            */
            Vec3 operator/(const T c) const
            {
                return Vec3(x / c, y / c, z / c);
            }

            /**
            Equality, tests for strict equality even for numeric storage types.
            */
            bool operator==(const Vec3& r) const
            {
                return x == r.x && y == r.y && z == r.z;
            }

            /**
            Inequality, tests for strict inequality even for numeric storage types.
            */
            bool operator!=(const Vec3& r) const
            {
                return !(*this == r);
            }

			/**
			Equality, uses NumberUtils<T>::equals(...) to test for equality of elements.
			*/
			bool equals(const Vec3& r, T tol = NumberUtils<T>::tolerance()) const
			{
				return NumberUtils<T>::equals(x, r.x, tol) && NumberUtils<T>::equals(y, r.y, tol) && NumberUtils<T>::equals(z, r.z, tol);
			}

			/**
			Calculates dot product between this vector and the given vector.
			*/
			template<typename real_t = typename NumberUtils<T>::FloatType>
			real_t dot(const Vec3& right) const
			{
				return ((real_t)x * (real_t)right.x) + ((real_t)y * (real_t)right.y) + ((real_t)z * (real_t)right.z);
			}

			/**
			Calculates the squared Euclidean norm of this vector.
			*/
			template<typename real_t = typename NumberUtils<T>::FloatType>
			real_t normSquared() const
			{
				return this->dot<real_t>(*this);
			}

			/**
			Calculates the Euclidean norm of this vector.
			*/
			template<typename real_t = typename NumberUtils<T>::FloatType>
			real_t norm() const
			{
				return pixelRound<real_t>(sqrt(normSquared<real_t>()));
			}

			/**
			Calculates cross product between this vector and the given vector.
			*/
			Vec3 cross(const Vec3& right) const
			{
				return Vec3(y * right.z - z * right.y,
							z * right.x - x * right.z,
							x * right.y - y * right.x);
			}

			/**
			Rotates this vector counterclockwise around the given axis by the given angle (in radians).
			Uses Rodrigues' rotation formula.
			*/
			template<typename real_t = typename NumberUtils<T>::FloatType>
			Vec3 rotate(const Vec3& axis, real_t angle) const
			{
				real_t c = cos(angle);
				real_t s = sin(angle);
				Vec3<real_t> ra(axis);
				Vec3<real_t> rv(*this);
				rv = rv * c + (ra.cross(rv)) * s + (ra.dot(rv)) * ((real_t)1.0 - c) * ra;
				return Vec3(rv);
			}

			/**
			Returns normalized version of this vector and calculates its original length.
			*/
			template<typename real_t = typename NumberUtils<T>::FloatType>
			Vec3 normalized(real_t& length) const
			{
				length = norm<real_t>();
				if (!NumberUtils<real_t>::equals(length, 0.0))
				{
					real_t m = 1 / length;
					return Vec3(pixelRound<T, real_t>(x * m), pixelRound<T, real_t>(y * m), pixelRound<T, real_t>(z * m));
				}
				else
					return *this;
			}

			/**
			Returns normalized version of this vector.
			*/
			template<typename real_t = typename NumberUtils<T>::FloatType>
			Vec3 normalized() const
			{
				real_t dummy;
				return normalized<real_t>(dummy);
			}

			/**
			Calculates angle between this vector and the given vector.
			Returns result in range [0, PI].
			*/
			template<typename real_t = typename NumberUtils<T>::FloatType>
			real_t angleTo(const Vec3<T>& r)
			{
				real_t dp = this->dot<real_t>(r) / (this->norm<real_t>() * r.norm<real_t>());
				if (dp <= -1)
					return (real_t)PI;
				else if (dp >= 1)
					return 0;
				return acos(dp);
			}

			/**
			Calculates sharp angle between this vector and the given vector.
			Returns result in range [0, PI/2].
			*/
			template<typename real_t = typename NumberUtils<T>::FloatType>
			real_t sharpAngleTo(const Vec3<T>& r)
			{
				real_t theta = angleTo<real_t>(r);
				if (theta < (real_t)(PI / 2))
					return theta;
				
				return (real_t)PI - theta;
			}

			/**
			Normalizes this vector.
			*/
			void normalize()
			{
				*this = normalized();
			}

			/**
			Normalizes this vector and calculates its original length.
			*/
			template<typename real_t = typename NumberUtils<T>::FloatType>
			void normalize(real_t& length)
			{
				*this = normalized<real_t>(length);
			}

			/**
			Returns a new vector whose components are absolute values of components of this vector.
			*/
			Vec3<T> abs() const
			{
				return Vec3<T>(::abs(x), ::abs(y), ::abs(z));
			}

			/**
			Returns largest component of this vector.
			*/
			T max() const
			{
				return std::max(x, std::max(y, z));
			}

			/**
			Returns smallest component of this vector.
			*/
			T min() const
			{
				return std::min(x, std::min(y, z));
			}

			/**
			Returns sum of all elements.
			*/
			T sum() const
			{
				return x + y + z;
			}

			/**
			Returns product of all elements.
			*/
			T product() const
			{
				return x * y * z;
			}

			Vec3<T> componentwiseMultiply(const Vec3<T>& r) const
			{
				return Vec3<T>(x * r.x, y * r.y, z * r.z);
			}

			Vec3<T> componentwiseDivide(const Vec3<T>& r) const
			{
				return Vec3<T>(x / r.x, y / r.y, z / r.z);
			}

			/**
			Returns if is a permutation of [0, 1, 2] and therefore a valid order for transposing
			*/
			const bool isPermutation() const
			{
				if(max()>2 || min()<0)
					return false;
				Vec3<T> counter(0,0,0);
				counter[x]++;
				counter[y]++;
				counter[z]++;
				return counter==Vec3<T>(1,1,1);
			}

			Vec3<T> transposed(Vec3<T> order)
			{
				if(!order.isPermutation()) throw ITLException("invalid order: "+ toString(order) + "expected a permutation of [0, 1, 2]");
				Vec3<T> transposed(0, 0, 0);
				transposed.x = operator[](order.x);
				transposed.y = operator[](order.y);
				transposed.z = operator[](order.z);
				return transposed;
			}

			const Vec3<T> inverseOrder()const
			{
				if(!isPermutation()) throw ITLException("invalid order: "+ toString(this) + "expected a permutation of [0, 1, 2]");
				Vec3<T> inverse(0,0,0);
				inverse[x] = 0;
				inverse[y] = 1;
				inverse[z] = 2;
				return inverse;
			}

            /**
            Converts this object to string.
            */
            friend std::ostream& operator<<(std::ostream& stream, const Vec3<T>& v)
            {
                stream << "[" << v.x << ", " << v.y << ", " << v.z << "]";
                return stream;
            }
    };

	template<typename T>
	std::istream& operator>>(std::istream& i, Vec3<T>& out)
	{
		i.ignore(1, '[');
		i >> out.x;
		i.ignore(2, ',');
		i >> out.y;
		i.ignore(2, ',');
		i >> out.z;
		i.ignore(1, ']');
		return i;
	}

	extern template class Vec3<itl2::float32_t>;
	extern template class Vec3<double>;
	extern template class Vec3<itl2::coord_t>;
	extern template class Vec3<int32_t>;

	typedef Vec3<itl2::float32_t> Vec3f;
    typedef Vec3<double> Vec3d;
	typedef Vec3<itl2::coord_t> Vec3c;
	typedef Vec3<int32_t> Vec3sc; // This is for use in situations where signed coordinate vector is needed and Vec3<coord_t> requires too much memory.

	/**
	Clamps the given value to range [lower, upper].
	*/
	template<typename T> void clamp(Vec3<T>& value, const Vec3<T>& lower, const Vec3<T>& upper)
	{
		clamp(value.x, lower.x, upper.x);
		clamp(value.y, lower.y, upper.y);
		clamp(value.z, lower.z, upper.z);
	}

	///**
	//Calculates componentwise minimum of a and b.
	//*/
	//template<typename T> Vec3<T> componentwiseMin(const Vec3<T>& a, const Vec3<T>& b)
	//{
	//	Vec3<T> res = a;
	//	if (b.x < res.x)
	//		res.x = b.x;
	//	if (b.y < res.y)
	//		res.y = b.y;
	//	if (b.z < res.z)
	//		res.z = b.z;
	//	return res;
	//}

	///**
	//Calculates componentwise maximum of a and b.
	//*/
	//template<typename T> Vec3<T> componentwiseMax(const Vec3<T>& a, const Vec3<T>& b)
	//{
	//	Vec3<T> res = a;
	//	if (b.x > res.x)
	//		res.x = b.x;
	//	if (b.y > res.y)
	//		res.y = b.y;
	//	if (b.z > res.z)
	//		res.z = b.z;
	//	return res;
	//}

	/*
	Calculates componentwise ceiling of a.
	*/
	template<typename T> Vec3c ceil(const Vec3<T>& a)
	{
		return Vec3c(itl2::ceil(a.x), itl2::ceil(a.y), itl2::ceil(a.z));
	}

	/*
	Calculates componentwise floor of a.
	*/
	template<typename T> Vec3c floor(const Vec3<T>& a)
	{
		return Vec3c(itl2::floor(a.x), itl2::floor(a.y), itl2::floor(a.z));
	}

	/**
	Rounds floating point vector to coordinate vector.
	*/
	template<typename T> Vec3c round(const Vec3<T>& value)
	{
		return Vec3c(itl2::round(value.x), itl2::round(value.y), itl2::round(value.z));
	}

	/*
	Trivial pixel rounding for vector value.
	*/
	template<> inline Vec3d pixelRound(Vec3d value)
	{
		return value;
	}

	/*
	Trivial pixel rounding for vector value.
	*/
	template<> inline Vec3<float> pixelRound(Vec3<float> value)
	{
		return value;
	}

	/*
	Trivial pixel rounding for vector value.
	*/
	template<> inline Vec3<float> pixelRound(Vec3d value)
	{
		return Vec3<float>((float)value.x, (float)value.y, (float)value.z);
	}

	/**
	Calculates minimum of two vectors elementwise.
	*/
	template<typename T> Vec3<T> min(const Vec3<T>& a, const Vec3<T>& b)
	{
		return Vec3<T>(std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z));
	}

	/**
	Calculates maximum of two vectors elementwise.
	*/
	template<typename T> Vec3<T> max(const Vec3<T>& a, const Vec3<T>& b)
	{
		return Vec3<T>(std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z));
	}

// *** Projections etc

	/**
	Project vector v to plane defined by plane normal vector.
	*/
	template<typename T> Vec3<T> projectToPlane(const Vec3<T>& v, const Vec3<T>& planeNormal)
	{
		return v - planeNormal * (v.dot(planeNormal) / planeNormal.dot(planeNormal));
	}

	/**
	Converts from spherical coordinates (r, azimuthal, polar) to cartesian coordinates (x, y, z).
	*/
	inline Vec3d sphericalToCartesian(double r, double azimuthal, double polar)
	{
		double x, y, z;
		sphericalToCartesian(r, azimuthal, polar, x, y, z);
		return Vec3d(x, y, z);
	}

	/**
	Converts from cartesian coordinates (x, y, z) to spherical coordinates (r, azimuthal, polar).
	*/
	inline void cartesianToSpherical(const Vec3d& v, double& r, double& azimuthal, double& polar)
	{
		cartesianToSpherical(v.x, v.y, v.z, r, azimuthal, polar);
	}

	/**
	Convert from cartesian coordinates (x, y, z) to cylindrical coordinates (r, azimuthal, z).
	*/
	inline void cartesianToCylindrical(const Vec3d& v, double& r, double& azimuthal, double& z)
	{
		z = v.z;
		cartesianToPolar(v.x, v.y, r, azimuthal);
	}

	/**
	Convert from cylindrical coordinates (r, azimuthal, z) to cartesian coordinates (x, y, z).
	*/
	inline Vec3d cylindricalToCartesian(double r, double azimuthal, double z)
	{
		double x, y;
		polarToCartesian(r, azimuthal, x, y);
		return Vec3d(x, y, z);
	}

	/**
	Compares two vectors to establish order so that containers containing vectors can be sorted.
	Uses exact comparisons even for floating point types.
	*/
	template<typename T> bool vecComparer(const Vec3<T>& a, const Vec3<T>& b)
	{
		if (a.z != b.z)
			return a.z < b.z;
		if (a.y != b.y)
			return a.y < b.y;
		return a.x < b.x;
	}

	namespace tests
	{
		void vectorAngles();
	}
}

