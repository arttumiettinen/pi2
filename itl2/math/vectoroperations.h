#pragma once

#include <vector>
#include <map>
#include "math/vec2.h"

using std::vector;
using std::map;
using std::logic_error;

namespace math
{

	inline void toPolar(vector<Vec2d>& v)
	{
		for(size_t n = 0; n < v.size(); n++)
		{
			toPolar(v[n]);
		}
	}

	inline void toCartesian(vector<Vec2d>& v)
	{
		for(size_t n = 0; n < v.size(); n++)
		{
			toCartesian(v[n]);
		}
	}

	/**
	* Calculates product of elements of the vector.
	*/
	template<typename T> T product(const vector<T>& items)
	{
		T result = 1;
		for(size_t n = 0; n < items.size(); n++)
			result *= items[n];
		return result;
	}

	/**
	* Create vector containing cumulative product of the given vector.
	* I.e. result[0] = items[0], result[1] = items[0] * items[1] etc.
	*/
	template<typename T> vector<T> cumulativeProduct(const vector<T>& items)
	{
		vector<T> result;
		T currVal = 1.0;
		for(size_t n = 0; n < items.size(); n++)
		{
			currVal *= items[n];
			result.push_back(currVal);
		}
		return result;
	}

	/**
	* Calculates maximum element in the vector.
	* Returns numeric_limits<T>::min() for empty list.
	*/
	template<typename T> T max(const vector<T>& items)
	{
		T result = numeric_limits<T>::lowest();
		for(size_t n = 0; n < items.size(); n++)
		{
			if(items[n] > result)
				result = items[n];
		}
		return result;
	}

	/**
	* Calculates minimum element in the vector.
	* Returns numeric_limits<T>::max() for empty list.
	*/
	template<typename T> T min(const vector<T>& items)
	{
		T result = numeric_limits<T>::max();
		for(size_t n = 0; n < items.size(); n++)
		{
			if(items[n] < result)
				result = items[n];
		}
		return result;
	}

	/**
	* Calculate min(a, b) itemwise.
	*/
	template<typename T> vector<T> min(const vector<T>& a, const vector<T>& b)
	{
		if(a.size() != b.size())
			throw ITLException("Sizes of vectors do not match.");

		vector<T> result;
		result.reserve(a.size());
		for(size_t n = 0; n < a.size(); n++)
		{
			result.push_back(min<T>(a[n], b[n]));
		}
		return result;
	}

	/**
	* Calculate max(a, b) itemwise.
	*/
	template<typename T> vector<T> max(const vector<T>& a, const vector<T>& b)
	{
		if(a.size() != b.size())
			throw ITLException("Sizes of vectors do not match.");

		vector<T> result;
		result.reserve(a.size());
		for(size_t n = 0; n < a.size(); n++)
		{
			result.push_back(max<T>(a[n], b[n]));
		}
		return result;
	}

	/**
	* Calculate a+b itemwise.
	*/
	template<typename T> vector<T> operator+(const vector<T>& a, const vector<T>& b)
	{
		if(a.size() != b.size())
			throw ITLException("Sizes of vectors do not match.");

		vector<T> result;
		result.reserve(a.size());
		for(size_t n = 0; n < a.size(); n++)
		{
			result.push_back(a[n] + b[n]);
		}
		return result;
	}

	/**
	Calculate a += const itemwise.
	*/
	template<typename T> vector<T>& operator+=(vector<T>& a, T b)
	{
		for(size_t n = 0; n < a.size(); n++)
		{
			a[n] += b;
		}
		return a;
	}

	/**
	* Calculate a+const itemwise.
	*/
	template<typename T> vector<T> operator+(vector<T> a, T b)
	{
		return a += b;
	}

	/**
	* Calculate a-b itemwise.
	*/
	template<typename T> vector<T> operator-(const vector<T>& a, const vector<T>& b)
	{
		if(a.size() != b.size())
			throw ITLException("Sizes of vectors do not match.");

		vector<T> result;
		result.reserve(a.size());
		for(size_t n = 0; n < a.size(); n++)
		{
			result.push_back(a[n] - b[n]);
		}
		return result;
	}

	/**
	Calculate a -= const itemwise.
	*/
	template<typename T> vector<T>& operator-=(vector<T>& a, T b)
	{
		for(size_t n = 0; n < a.size(); n++)
		{
			a[n] -= b;
		}
		return a;
	}

	/**
	* Calculate a-const itemwise.
	*/
	template<typename T> vector<T> operator-(vector<T> a, T b)
	{
		return a -= b;
	}

	/**
	 * Calculate -a itemwise.
	 */
	template<typename T> vector<T> operator-(const vector<T>& a)
	{
		vector<T> result;
		result.reserve(a.size());
		for(size_t n = 0; n < a.size(); n++)
		{
			result.push_back(-a[n]);
		}
		return result;
	}

	/**
	 * Calculates absolute value itemwise.
	 */
	template<typename T> vector<T> abs(const vector<T>& a)
	{
		vector<T> result;
		result.reserve(a.size());
		for(size_t n = 0; n < a.size(); n++)
		{
			result.push_back(::abs(a[n]));
		}
		return result;
	}

	/**
	* Calculate a*const itemwise.
	*/
	template<typename T> vector<double> operator*(const vector<T>& a, double b)
	{
		vector<double> result;
		result.reserve(a.size());
		for(size_t n = 0; n < a.size(); n++)
		{
			result.push_back(a[n] * b);
		}
		return result;
	}

	/**
	* Calculate a*b itemwise.
	*/
	template<typename T1, typename T2, typename Tresult> vector<Tresult> multiply(const vector<T1>& a, const vector<T2>& b)
	{
		if(a.size() != b.size())
			throw ITLException("Sizes of vectors do not match.");

		vector<Tresult> result;
		result.reserve(a.size());
		for(size_t n = 0; n < a.size(); n++)
		{
			result.push_back(pixelRound<Tresult>((double)a[n] * (double)b[n]));
		}
		return result;
	}

	/**
	Calculates sum of values of the vector.
	*/
	template<typename T> T sum(const vector<T>& a)
	{
		T result = T();
		for(size_t n = 0; n < a.size(); n++)
		{
			result += a[n];
		}
		return result;
	}

	/**
	Calculate squared norm of the given vector.
	*/
	template<typename T, typename Tresult> Tresult norm2(const vector<T>& a)
	{
		return sum<Tresult>(multiply<T, T, Tresult>(a, a));
	}

	/**
	Calculate norm of the given vector.
	*/
	template<typename T, typename Tresult> Tresult norm(const vector<T>& a)
	{
		return sqrt(norm2<T, Tresult>(a));
	}


	/**
	* Calculate a*b itemwise.
	*/
	template<typename T> vector<double> operator*(double b, const vector<T>& a)
	{
		return a * b;
	}

	/**
	* Calculate a/const itemwise.
	*/
	template<typename T> vector<double> operator/(const vector<T>& a, double b)
	{
		vector<double> result;
		result.reserve(a.size());
		for(size_t n = 0; n < a.size(); n++)
		{
			result.push_back(a[n] / b);
		}
		return result;
	}

	/**
	* Calculate a/b itemwise.
	*/
	template<typename T1, typename T2, typename Tresult> vector<Tresult> divide(const vector<T1>& a, const vector<T2>& b)
	{
		if(a.size() != b.size())
			throw ITLException("Sizes of vectors do not match.");

		vector<Tresult> result;
		result.reserve(a.size());
		for(size_t n = 0; n < a.size(); n++)
		{
			result.push_back(pixelRound<Tresult>((double)a[n] / (double)b[n]));
		}
		return result;
	}



	/**
	Calculates mean of values of the vector.
	*/
	template<typename T> T mean(const vector<T>& a)
	{
		return sum(a) / (double)a.size();
	}

	/**
	Calculates mean of values of the vector using given weights.
	*/
	template<typename T> T mean(const vector<T>& a, const vector<T>& weights)
	{
		return sum(multiply<T, T, T>(a, weights)) / sum(weights);
	}

	/**
	Compares two map entries x and y. Returns true if x.second < y.second.
	*/
	template<class T> bool pairCompare(const T & x, const T & y)
	{
		return x.second < y.second;
	}

	/**
	Finds the element of a map A whose .second is the largest.
	*/
	template<class T> typename T::const_iterator mapMaxElement(const T &A)
	{
		typedef typename T::value_type pair_type;
		return max_element(A.begin(), A.end(), pairCompare<typename T::value_type>);
	}

	/**
	Calculates mode of samples.
	*/
	template<typename T> T mode(const vector<T>& data)
	{
		if (data.size() <= 0)
			throw logic_error("Mode of no data");

		map<T, size_t> counts;
		for (const T& x : data)
			counts[x]++;

		return mapMaxElement(counts)->first;
	}


	/**
	 * Cast all items in a vector to some other type.
	 */
	template<typename Tdest, typename Tsrc> vector<Tdest> cast(const vector<Tsrc>& v)
	{
		vector<Tdest> result;
		result.reserve(v.size());
		for(size_t n = 0; n < v.size(); n++)
		{
			result.push_back((Tdest)(v[n]));
		}
		return result;
	}

	/**
	Round all values in a vector.
	*/
	template<typename Tout, typename Tin> vector<Tout> roundv(const vector<Tin>& a)
	{
		vector<Tout> result;
		result.reserve(a.size());
		for(size_t n = 0; n < a.size(); n++)
		{
			result.push_back(round(a[n]));
		}
		return result;
	}

	/**
	Ceil all values in a vector.
	*/
	template<typename Tout, typename Tin> vector<Tout> ceilv(const vector<Tin>& a)
	{
		vector<Tout> result;
		result.reserve(a.size());
		for(size_t n = 0; n < a.size(); n++)
		{
			result.push_back(ceil(a[n]));
		}
		return result;
	}
	
	/**
	Create vector containing three values.
	*/
	template<typename T> vector<T> createv(const T& a, const T& b, const T& c)
	{
		vector<T> result;
		result.reserve(3);
		result.push_back(a);
		result.push_back(b);
		result.push_back(c);
		return result;
	}

}
