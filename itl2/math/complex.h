#pragma once

#include "math/equals.h"

namespace math
{

	/**
	Complex number
	*/
	template<typename T> struct Complex
	{
	public:

		union
		{
			T components[2];
			struct
			{
				T re;
				T im;
			};
		};


		Complex() : re(0), im(0)
		{
		}

		Complex(const Complex& right) : re(right.re), im(right.im)
		{
		}

		Complex(T re, T im) : re(re), im(im)
		{
		}

		
		explicit Complex(T reim[2]) : re(reim[0]), im(reim[1])
		{
		}
		
		Complex(T re) : re(re), im(0)
		{
		}

		/**
        Assignment
        */
        Complex& operator=(const Complex<T>& other)
        {
            re = other.re;
            im = other.im;
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
		double size() const
		{
			return 2;
		}

		/**
		Return complex conjugate of this number.
		*/
		Complex<T> conjugate() const
		{
			return Complex<T>(re, -im);
		}

		/**
        Addition in-place
        */
        Complex<T>& operator+=(const Complex<T>& r)
        {
            re += r.re;
            im += r.im;
            return *this;
        }

        /**
        Subtraction in-place
        */
        Complex<T>& operator-=(const Complex<T>& r)
        {
            re -= r.re;
            im -= r.im;
            return *this;
        }

		/**
		Multiplication in-place
		*/
		Complex<T>& operator*=(const Complex<T>& r)
		{
			T a = re;
			T b = im;
			T c = r.re;
			T d = r.im;

			re = a*c-b*d;
			im = a*d+b*c;
			return *this;
		}

		Complex<T>& operator*=(const T x)
		{
			re *= x;
			im *= x;
			return *this;
		}

		Complex<T>& operator/=(const T x)
		{
			re /= x;
			im /= x;
			return *this;
		}

		/**
        Unary minus
        */
        Complex<T> operator-() const
        {
            return Complex<T>(-re, -im);
        }

		/**
        Binary addition
        */
        const Complex<T> operator+(const Complex<T> &r) const
        {
            return Complex<T>(*this) += r;
        }

        /**
        Binary subtraction
        */
        const Complex<T> operator-(const Complex<T> &r) const
        {
            return Complex<T>(*this) -= r;
        }

		/**
        Binary multiplication
        */
        const Complex<T> operator*(const Complex<T> &r) const
        {
            return Complex<T>(*this) *= r;
        }

		/**
        Binary multiplication
        */
        const Complex<T> operator*(const T r) const
        {
            return Complex<T>(*this) *= r;
        }

		/**
        Binary division
        */
        const Complex<T> operator/(const T r) const
        {
            return Complex<T>(*this) /= r;
        }

		/**
		Calculates the squared Euclidean norm of this vector.
		*/
		T normSquared() const
		{
			return re * re + im * im;
		}

		/**
		Calculates the Euclidean norm of this vector.
		*/
		T norm() const
		{
			return (T)sqrt(normSquared());
		}

		/**
		Returns normalized version of this vector.
		*/
		Complex<T> normalized() const
		{
			double l = normSquared();
			if (!equals(l, 0.0))
				return *this / (T)sqrt(l);
			else
				return *this;
		}

		/**
		Returns normalized version of this vector and calculates its original length.
		*/
		Complex<T> normalized(double& length) const
		{
			length = norm();
			if (!equals(length, 0.0))
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
		void normalize(double& length)
		{
			*this = normalized(length);
		}
	};

	typedef Complex<float> Complexf;
	typedef Complex<double> Complexd;

}
