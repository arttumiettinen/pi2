#pragma once

#include <sstream>
#include "datatypes.h"
#include "math/numberutils.h"
#include "math/mathutils.h"

using namespace itl2;

namespace math
{
    /**
    Three-component vector.
    */
    template <typename T, typename real_t = typename NumberUtils<T>::FloatType> class Vec2
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
				};
				T components[2];
				struct // If the vector is used as range.
				{
					T min;
					T max;
				};
			};

            /**
            Default constructor, initializes the vector to zero.
            */
            Vec2() :
                x(0),
                y(0)
            {
            }

            /**
            Constructor
            */
            Vec2(T x, T y)
            {
                this->x = x;
                this->y = y;
            }

            /**
            Copy constructor
            */
            Vec2(const Vec2& other) :
                x(other.x),
                y(other.y)
            {
            }

			/**
			Constructor that makes it possible to initialize vector from vector of another type.
			Compiler warning identifies the cases where the conversion could lose precision.
			*/
			template<typename Tother> Vec2(const Vec2<Tother>& other) :
				x(other.x),
				y(other.y)
			{
			}

            /**
            Assignment
            */
            Vec2& operator=(const Vec2& other)
            {
                x = other.x;
                y = other.y;
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
			Gets count of components in this vector.
			*/
			size_t size() const
			{
				return 2;
			}

            /**
            Addition in-place
            */
            Vec2& operator+=(const Vec2& r)
            {
                x += r.x;
                y += r.y;
                return *this;
            }

            /**
            Subtraction in-place
            */
            Vec2& operator-=(const Vec2& r)
            {
                x -= r.x;
                y -= r.y;
                return *this;
            }

			/**
            Multiplication in-place
            */
			Vec2& operator*=(const T r)
            {
                x *= r;
                y *= r;
                return *this;
            }

            /**
            Division in-place
            */
            Vec2& operator/=(const T r)
            {
                x /= r;
                y /= r;
                return *this;
            }

            /**
            Addition of constant in-place
            */
            /*
            Vec2& operator+=(const T c)
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
            Vec2& operator-=(const T c)
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
            Vec2 operator-() const
            {
                return Vec2(-x, -y);
            }

            /**
            Binary addition
            */
            const Vec2 operator+(const Vec2 &r) const
            {
                return Vec2(*this) += r;
            }

            /**
            Binary subtraction
            */
            const Vec2 operator-(const Vec2 &r) const
            {
                return Vec2(*this) -= r;
            }

            /**
            Multiplication by constant.
            */
            const Vec2 operator*(const T c) const
            {
                return Vec2(x * c, y * c);
            }

            /**
            Multiplication by constant.
            */
            friend const Vec2 operator*(const T c, const Vec2& r)
            {
                return r * c;
            }

            /**
            Division by constant.
            */
            const Vec2 operator/(const T c) const
            {
                return Vec2(x / c, y / c);
            }

			/**
			Calculates dot product between this vector and the given vector.
			*/
			real_t dot(const Vec2& right) const
			{
				return (real_t)x * (real_t)right.x + (real_t)y * (real_t)right.y;
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
			Returns normalized version of this vector.
			*/
			Vec2 normalized() const
			{
				real_t l = normSquared();
				if (!NumberUtils<real_t>::equals(l, 0.0))
					return *this / sqrt(l);
				else
					return *this;
			}

			/**
			Returns normalized version of this vector and calculates its original length.
			*/
			Vec2 normalized(real_t& length) const
			{
				length = norm();
				if (!NumberUtils<real_t>::equals(length, 0.0))
					return *this / length;
				else
					return *this;
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
            Converts this object to string.
            */
            friend ostream& operator<<(ostream& stream, const Vec2<T>& v)
            {
                stream << "[" << v.x << ", " << v.y << "]";
                return stream;
            }

			/**
            Equality, tests for strict equality even for numeric storage types.
            */
            bool operator==(const Vec2& r) const
            {
                return x == r.x && y == r.y;
            }

            /**
            Inequality, tests for strict inequality even for numeric storage types.
            */
            bool operator!=(const Vec2& r) const
            {
                return !(*this == r);
            }

			/**
			Equality, uses equals(T, T) to test for equality of elements.
			*/
			bool equals(const Vec2& r) const
			{
				return NumberUtils<T>::equals(x, r.x) && NumberUtils<T>::equals(y, r.y);
			}

			bool operator< (const Vec2<T>& right) const
			{
				if(y < right.y)
					return true;
				else if(y > right.y)
					return false;
				else
					return x < right.x;
			}

			bool operator> (const Vec2<T>& rhs) const
			{
				return rhs < *this;
			}

			bool operator<=(const Vec2<T>& rhs) const
			{
				return !(*this > rhs);
			}

			bool operator>=(const Vec2<T>& rhs) const
			{
				return !(*this < rhs);
			}
    };

    typedef Vec2<double> Vec2d;
	typedef Vec2<float> Vec2f;
    typedef Vec2<int> Vec2i;
	typedef Vec2<coord_t> Vec2c;

	/**
	Rounds double vector to coordinate vector.
	*/
	template<typename T> Vec2<coord_t> round(Vec2<T> value)
	{
		return Vec2<coord_t>(round(value.x), round(value.y));
	}

	/**
	See tocartesian overloads.
	*/
	inline Vec2d toCartesian(double r, double azimuthal)
	{
		double x, y;
		toCartesian(r, azimuthal, x, y);
		return Vec2d(x, y);
	}

	/**
	Converts v=(x, y) in cartesian coordinates to v=(r, atzimuthal)
	*/
	inline void toPolar(Vec2d& v)
	{
		double r, atz;
		toPolar(v.x, v.y, r, atz);
		v.x = r;
		v.y = atz;
	}

	/**
	Converts v=(r, atzimuthal) in cartesian coordinates to v=(x, y)
	*/
	inline void toCartesian(Vec2d& v)
	{
		double x, y;
		toCartesian(v.x, v.y, x, y);
		v.x = x;
		v.y = y;
	}
}

