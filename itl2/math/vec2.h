#pragma once

#include <sstream>
#include "datatypes.h"
#include "math/numberutils.h"
#include "math/mathutils.h"

namespace itl2
{
    /**
    Three-component vector.
    */
    template <typename T> class Vec2
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
			Constructor that makes it possible to cast/initialize vector from vector of another type.
			*/
			template<typename Tother> explicit Vec2(const Vec2<Tother>& other) :
				x(pixelRound<T, Tother>(other.x)),
				y(pixelRound<T, Tother>(other.y))
			{
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
			template<typename real_t = typename NumberUtils<T>::FloatType>
			real_t dot(const Vec2& right) const
			{
				return (real_t)x * (real_t)right.x + (real_t)y * (real_t)right.y;
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
			Returns normalized version of this vector.
			*/
			template<typename real_t = typename NumberUtils<T>::FloatType>
			Vec2 normalized() const
			{
				real_t l = normSquared<real_t>();
				if (!NumberUtils<real_t>::equals(l, 0.0))
					return *this / sqrt(l);
				else
					return *this;
			}

			/**
			Returns normalized version of this vector and calculates its original length.
			*/
			template<typename real_t = typename NumberUtils<T>::FloatType>
			Vec2 normalized(real_t& length) const
			{
				length = norm<real_t>();
				if (!NumberUtils<real_t>::equals(length, 0.0))
					return *this / length;
				else
					return *this;
			}

			/**
			Normalizes this vector.
			*/
			template<typename real_t = typename NumberUtils<T>::FloatType>
			void normalize()
			{
				*this = normalized<real_t>();
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
            Converts this object to string.
            */
            friend std::ostream& operator<<(std::ostream& stream, const Vec2<T>& v)
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

			/**
			Returns largest component of this vector.
			*/
			T max() const
			{
				return std::max(x, y);
			}

			/**
			Returns smallest component of this vector.
			*/
			T min() const
			{
				return std::min(x, y);
			}
    };

    typedef Vec2<double> Vec2d;
	typedef Vec2<itl2::float32_t> Vec2f;
	typedef Vec2<itl2::coord_t> Vec2c;

	/**
	Rounds floating point vector to coordinate vector.
	*/
	template<typename T> Vec2<itl2::coord_t> round(Vec2<T> value)
	{
		return Vec2<itl2::coord_t>(round(value.x), round(value.y));
	}

	/**
	Converts polar coordinates (r, azimuthal) to cartesian coordinates (x, y).
	*/
	inline Vec2d polarToCartesian(double r, double azimuthal)
	{
		double x, y;
		polarToCartesian(r, azimuthal, x, y);
		return Vec2d(x, y);
	}

	/**
	Converts (x, y) in cartesian coordinates to polar coordinates (r, atzimuthal)
	*/
	inline void cartesianToPolar(Vec2d& v, double& r, double& azimuthal)
	{
		cartesianToPolar(v.x, v.y, r, azimuthal);
	}

	/*
	Calculates componentwise ceiling of a.
	*/
	template<typename T> Vec2c ceil(const Vec2<T>& a)
	{
		return Vec2c(itl2::ceil(a.x), itl2::ceil(a.y));
	}

	/*
	Calculates componentwise floor of a.
	*/
	template<typename T> Vec2c floor(const Vec2<T>& a)
	{
		return Vec2c(itl2::floor(a.x), itl2::floor(a.y));
	}

	/**
	Calculates componentwise minimum of a and b.
	*/
	template<typename T> Vec2<T> min(const Vec2<T>& a, const Vec2<T>& b)
	{
		return Vec2<T>(std::min(a.x, b.x), std::min(a.y, b.y));
	}

	/**
	Calculates componentwise maximum of a and b.
	*/
	template<typename T> Vec2<T> max(const Vec2<T>& a, const Vec2<T>& b)
	{
		return Vec2<T>(std::max(a.x, b.x), std::max(a.y, b.y));
	}
}

