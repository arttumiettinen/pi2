#pragma once

#include <sstream>
#include <vector>
#include "itlexception.h"
#include "math/numberutils.h"

namespace itl2
{
    /**
    Four-component vector.
    @param T Storage type of elements.
    */
    template <typename T> class Vec4
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
					T w;
				};
				T components[4];
			};

            /**
            Default constructor, initializes the vector to zero.
            */
            Vec4() :
                x(0),
                y(0),
                z(0),
                w(0)
            {
            }

            /**
            Constructor
            */
            Vec4(T x, T y, T z, T w)
            {
                this->x = x;
                this->y = y;
                this->z = z;
                this->w = w;
            }

			/**
			Constructs Vec4 from the first four elements of an array type. If there are not enough elements in the array, the remaining elements are set to zero.
			*/
			//Vec4(const std::vector<T>& other)
			//{
			//	if(other.size() != 4)
			//		throw itl2::ITLException("Invalid vector in Vec4 constructor.");

			//	x = other[0];
			//	y = other[1];
			//	z = other[2];
			//	w = other[3];
			//}
            template<typename VEC> explicit Vec4(const VEC& v)
            {
                x = 0;
                y = 0;
                z = 0;
                w = 0;
                for (size_t n = 0; n < std::min<size_t>(4, v.size()); n++)
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
				return 4;
			}

            /**
            Addition in-place
            */
            Vec4& operator+=(const Vec4& r)
            {
                x += r.x;
                y += r.y;
                z += r.z;
                w += r.w;
                return *this;
            }

            /**
            Subtraction in-place
            */
            Vec4& operator-=(const Vec4& r)
            {
                x -= r.x;
                y -= r.y;
                z -= r.z;
                w -= r.w;
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
            Vec4 operator-() const
            {
                return Vec4(-x, -y, -z, -w);
            }

            /**
            Binary addition
            */
            Vec4 operator+(const Vec4 &r) const
            {
                return Vec4(*this) += r;
            }

            /**
            Binary subtraction
            */
            Vec4 operator-(const Vec4 &r) const
            {
                return Vec4(*this) -= r;
            }

            /**
            Multiplication by constant.
            */
            Vec4<T> operator*(const T c) const
            {
                return Vec4<T>(x * c, y * c, z * c, w * c);
            }

            /**
            Multiplication by constant.
            */
            friend Vec4 operator*(const T c, const Vec4& r)
            {
                return r * c;
            }

            /**
            Division by constant.
            */
            Vec4 operator/(const T c) const
            {
                return Vec4(x / c, y / c, z / c, w / c);
            }

            /**
            Equality, tests for strict equality even for numeric storage types.
            */
            bool operator==(const Vec4& r) const
            {
                return x == r.x && y == r.y && z == r.z && w == r.w;
            }

            /**
            Inequality, tests for strict inequality even for numeric storage types.
            */
            bool operator!=(const Vec4& r) const
            {
                return !(*this == r);
            }

			/**
			Calculates dot product between this vector and the given vector.
			*/
			template<typename real_t = typename NumberUtils<T>::FloatType>
			real_t dot(const Vec4& right) const
			{
				return ((real_t)x * (real_t)right.x) + ((real_t)y * (real_t)right.y) + ((real_t)z * (real_t)right.z) + ((real_t)w * (real_t)right.w);
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
			Returns normalized version of this vector and calculates its original length.
			*/
			template<typename real_t = typename NumberUtils<T>::FloatType>
			Vec4 normalized(real_t& length) const
			{
				length = normSquared<real_t>();
				if (!NumberUtils<real_t>::equals(length, 0.0))
				{
					real_t m = 1 / sqrt(length);
					return Vec4(pixelRound<T, real_t>(x * m), pixelRound<T, real_t>(y * m), pixelRound<T, real_t>(z * m), pixelRound<T, real_t>(w * m));
				}
				else
					return *this;
			}

			/**
			Returns normalized version of this vector.
			*/
			template<typename real_t = typename NumberUtils<T>::FloatType>
			Vec4 normalized() const
			{
				real_t dummy;
				return normalized<real_t>(dummy);
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
            friend std::ostream& operator<<(std::ostream& stream, const Vec4<T>& v)
            {
                stream << "[" << v.x << ", " << v.y << ", " << v.z << ", " << v.w << "]";
                return stream;
            }
    };

	typedef Vec4<float32_t> Vec4f;
    typedef Vec4<double> Vec4d;

}

