#pragma once

#include <sstream>
#include <vector>
#include "datatypes.h"
#include "itlexception.h"
#include "math/numberutils.h"
#include "math/mathutils.h"

using namespace itl2;

namespace math
{
    /**
    Three-component vector.
    */
    template <typename T, typename real_t = typename NumberUtils<T>::FloatType> class Vec3
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
            Copy constructor
            */
            Vec3(const Vec3& other) :
                x(other.x),
                y(other.y),
                z(other.z)
            {
            }

			/**
			Constructor that makes it possible to cast/initialize vector from vector of another type.
			*/
			template<typename Tother, typename other_real_t> explicit Vec3(const Vec3<Tother, other_real_t>& other) :
				x(pixelRound<T, Tother>(other.x)),
				y(pixelRound<T, Tother>(other.y)),
				z(pixelRound<T, Tother>(other.z))
			{
			}

			/**
			Constructs Vec3 from vector<T>. Throws exception if the vector does not contain 3 elements.
			*/
			Vec3(const vector<T>& other)
			{
				if(other.size() != 3)
					throw itl2::ITLException("Invalid vector in Vec3 constructor.");

				x = other[0];
				y = other[1];
				z = other[2];
			}

            /**
            Assignment
            */
            Vec3& operator=(const Vec3& other)
            {
                x = other.x;
                y = other.y;
                z = other.z;
                return *this;
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
			Equality, uses equals(T, T) to test for equality of elements.
			*/
			bool equals(const Vec3& r) const
			{
				return NumberUtils<T>::equals(x, r.x) && NumberUtils<T>::equals(y, r.y) && NumberUtils<T>::equals(z, r.z);
			}

			/**
			Calculates dot product between this vector and the given vector.
			*/
			real_t dot(const Vec3& right) const
			{
				return ((real_t)x * (real_t)right.x) + ((real_t)y * (real_t)right.y) + ((real_t)z * (real_t)right.z);
			}

			/**
			Calculates the squared Euclidean norm of this vector.
			*/
			real_t normSquared() const
			{
				return this->dot(*this);
			}

			/**
			Calculates the Euclidean norm of this vector.
			*/
			real_t norm() const
			{
				return sqrt(normSquared());
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
			Vec3 normalized(real_t& length) const
			{
				length = normSquared();
				if (!NumberUtils<real_t>::equals(length, 0.0))
				{
					real_t m = 1 / sqrt(length);
					return Vec3(pixelRound<T, real_t>(x * m), pixelRound<T, real_t>(y * m), pixelRound<T, real_t>(z * m));
				}
				else
					return *this;
			}

			/**
			Returns normalized version of this vector.
			*/
			Vec3 normalized() const
			{
				real_t dummy;
				return normalized(dummy);
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
			void normalize(real_t& length)
			{
				*this = normalized(length);
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
				return math::max(x, math::max(y, z));
			}

			/**
			Returns smallest component of this vector.
			*/
			T min() const
			{
				return math::min(x, math::min(y, z));
			}

            /**
            Converts this object to string.
            */
            friend ostream& operator<<(ostream& stream, const Vec3<T, real_t>& v)
            {
                stream << "[" << v.x << ", " << v.y << ", " << v.z << "]";
                return stream;
            }
    };

	typedef Vec3<float> Vec3f;
    typedef Vec3<double> Vec3d;
    typedef Vec3<int> Vec3i;
	typedef Vec3<coord_t> Vec3c;

	/**
	Clamps the given value to range [lower, upper].
	*/
	template<typename T> void clamp(Vec3<T>& value, const Vec3<T>& lower, const Vec3<T>& upper)
	{
		clamp(value.x, lower.x, upper.x);
		clamp(value.y, lower.y, upper.y);
		clamp(value.z, lower.z, upper.z);
	}

	// TODO: We have componentWiseMax/Min and max/min that do the same thing!

	/**
	Calculates componentwise minimum of a and b.
	*/
	template<typename T> Vec3<T> componentwiseMin(const Vec3<T>& a, const Vec3<T>& b)
	{
		Vec3<T> res = a;
		if (b.x < res.x)
			res.x = b.x;
		if (b.y < res.y)
			res.y = b.y;
		if (b.z < res.z)
			res.z = b.z;
		return res;
	}

	/**
	Calculates componentwise maximum of a and b.
	*/
	template<typename T> Vec3<T> componentwiseMax(const Vec3<T>& a, const Vec3<T>& b)
	{
		Vec3<T> res = a;
		if (b.x > res.x)
			res.x = b.x;
		if (b.y > res.y)
			res.y = b.y;
		if (b.z > res.z)
			res.z = b.z;
		return res;
	}

	/*
	Calculates componentwise ceiling of a.
	*/
	template<typename T> Vec3c componentwiseCeil(const Vec3<T>& a)
	{
		return Vec3c((coord_t)ceil(a.x), (coord_t)ceil(a.y), (coord_t)ceil(a.z));
	}

	/*
	Calculates componentwise floor of a.
	*/
	template<typename T> Vec3c componentwiseFloor(const Vec3<T>& a)
	{
		return Vec3c((coord_t)floor(a.x), (coord_t)floor(a.y), (coord_t)floor(a.z));
	}

	/**
	Rounds double vector to coordinate vector.
	*/
	template<typename T> Vec3<coord_t> round(const Vec3<T>& value)
	{
		return Vec3<coord_t>(round(value.x), round(value.y), round(value.z));
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
		return Vec3<T>(min(a.x, b.x), min(a.y, b.y), min(a.z, b.z));
	}

	/**
	Calculates maximum of two vectors elementwise.
	*/
	template<typename T> Vec3<T> max(const Vec3<T>& a, const Vec3<T>& b)
	{
		return Vec3<T>(max(a.x, b.x), max(a.y, b.y), max(a.z, b.z));
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
	See toCartesian overloads.
	*/
	inline Vec3d toCartesian(double r, double azimuthal, double polar)
	{
		double x, y, z;
		toCartesian(r, azimuthal, polar, x, y, z);
		return Vec3d(x, y, z);
	}

	/**
	See toSpherical overloads.
	*/
	inline void toSpherical(const Vec3d& v, double& r, double& azimuthal, double& polar)
	{
		toSpherical(v.x, v.y, v.z, r, azimuthal, polar);
	}
}

