#pragma once

#include <vector>
#include <map>
#include "math/vec2.h"

namespace itl2
{

	//inline void toPolar(std::vector<Vec2d>& v)
	//{
	//	for(size_t n = 0; n < v.size(); n++)
	//	{
	//		cartesianToPolar(v[n]);
	//	}
	//}

	//inline void toCartesian(std::vector<Vec2d>& v)
	//{
	//	for(size_t n = 0; n < v.size(); n++)
	//	{
	//		toCartesian(v[n]);
	//	}
	//}
	
	/**
	Cast all items in a std::vector to some other type.
	*/
	template<typename Tdest, typename Tsrc> std::vector<Tdest> cast(const std::vector<Tsrc>& v)
	{
		std::vector<Tdest> result;
		result.reserve(v.size());
		for(size_t n = 0; n < v.size(); n++)
		{
			result.push_back((Tdest)(v[n]));
		}
		return result;
	}

	/**
	Calculates product of elements of the std::vector.
	*/
	template<typename T> T product(const std::vector<T>& items)
	{
		T result = 1;
		for(size_t n = 0; n < items.size(); n++)
			result *= items[n];
		return result;
	}

	/**
	Create std::vector containing cumulative product of the given std::vector.
	I.e. result[0] = items[0], result[1] = items[0] * items[1] etc.
	*/
	template<typename T> std::vector<T> cumulativeProduct(const std::vector<T>& items)
	{
		std::vector<T> result;
		T currVal = 1.0;
		for(size_t n = 0; n < items.size(); n++)
		{
			currVal *= items[n];
			result.push_back(currVal);
		}
		return result;
	}

	/**
	Calculates maximum element in the std::vector.
	Returns std::numeric_limits<T>::min() for empty list.
	*/
	template<typename T> T max(const std::vector<T>& items)
	{
		T result = std::numeric_limits<T>::lowest();
		for(size_t n = 0; n < items.size(); n++)
		{
			if(items[n] > result)
				result = items[n];
		}
		return result;
	}

	/**
	Calculates minimum element in the std::vector.
	Returns std::numeric_limits<T>::max() for empty list.
	*/
	template<typename T> T min(const std::vector<T>& items)
	{
		T result = std::numeric_limits<T>::max();
		for(size_t n = 0; n < items.size(); n++)
		{
			if(items[n] < result)
				result = items[n];
		}
		return result;
	}

	/**
	Calculate min(a, b) itemwise.
	*/
	template<typename T> std::vector<T> min(const std::vector<T>& a, const std::vector<T>& b)
	{
		if(a.size() != b.size())
			throw ITLException("Sizes of std::vectors do not match.");

		std::vector<T> result;
		result.reserve(a.size());
		for(size_t n = 0; n < a.size(); n++)
		{
			result.push_back(min<T>(a[n], b[n]));
		}
		return result;
	}

	/**
	Calculate max(a, b) itemwise.
	*/
	template<typename T> std::vector<T> max(const std::vector<T>& a, const std::vector<T>& b)
	{
		if(a.size() != b.size())
			throw ITLException("Sizes of std::vectors do not match.");

		std::vector<T> result;
		result.reserve(a.size());
		for(size_t n = 0; n < a.size(); n++)
		{
			result.push_back(max<T>(a[n], b[n]));
		}
		return result;
	}

	/**
	Calculate a+b itemwise.
	*/
	template<typename T> std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
	{
		if(a.size() != b.size())
			throw ITLException("Sizes of std::vectors do not match.");

		std::vector<T> result;
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
	template<typename T> std::vector<T>& operator+=(std::vector<T>& a, T b)
	{
		for(size_t n = 0; n < a.size(); n++)
		{
			a[n] += b;
		}
		return a;
	}

	/**
	Calculate a+const itemwise.
	*/
	template<typename T> std::vector<T> operator+(std::vector<T> a, T b)
	{
		return a += b;
	}

	/**
	Calculate a-b itemwise.
	*/
	template<typename T> std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b)
	{
		if(a.size() != b.size())
			throw itl2::ITLException("Sizes of std::vectors do not match.");

		std::vector<T> result;
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
	template<typename T> std::vector<T>& operator-=(std::vector<T>& a, T b)
	{
		for(size_t n = 0; n < a.size(); n++)
		{
			a[n] -= b;
		}
		return a;
	}

	/**
	Calculate a-const itemwise.
	*/
	template<typename T> std::vector<T> operator-(std::vector<T> a, T b)
	{
		return a -= b;
	}

	/**
	Calculate -a itemwise.
	*/
	template<typename T> std::vector<T> operator-(const std::vector<T>& a)
	{
		std::vector<T> result;
		result.reserve(a.size());
		for(size_t n = 0; n < a.size(); n++)
		{
			result.push_back(-a[n]);
		}
		return result;
	}

	/**
	Calculates absolute value itemwise.
	*/
	template<typename T> std::vector<T> abs(const std::vector<T>& a)
	{
		std::vector<T> result;
		result.reserve(a.size());
		for(size_t n = 0; n < a.size(); n++)
		{
			result.push_back(::abs(a[n]));
		}
		return result;
	}

	/**
	Calculate a*const itemwise.
	*/
	template<typename T> std::vector<double> operator*(const std::vector<T>& a, double b)
	{
		std::vector<double> result;
		result.reserve(a.size());
		for(size_t n = 0; n < a.size(); n++)
		{
			result.push_back(a[n] * b);
		}
		return result;
	}

	/**
	Calculate a*b itemwise.
	*/
	template<typename T1, typename T2, typename Tresult> std::vector<Tresult> multiply(const std::vector<T1>& a, const std::vector<T2>& b)
	{
		if(a.size() != b.size())
			throw ITLException("Sizes of std::vectors do not match.");

		std::vector<Tresult> result;
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
	template<typename T, typename Tresult> Tresult sum(const std::vector<T>& a)
	{
		Tresult result = Tresult();
		for(size_t n = 0; n < a.size(); n++)
		{
			result += (Tresult)a[n];
		}
		return result;
	}

	/**
	Calculate squared norm of the given vector.
	*/
	template<typename T, typename Tresult> Tresult norm2(const std::vector<T>& a)
	{
		return sum<Tresult>(multiply<T, T, Tresult>(a, a));
	}

	/**
	Calculate norm of the given vector.
	*/
	template<typename T, typename Tresult> Tresult norm(const std::vector<T>& a)
	{
		return sqrt(norm2<T, Tresult>(a));
	}


	/**
	Calculate a*b itemwise.
	*/
	template<typename T> std::vector<double> operator*(double b, const std::vector<T>& a)
	{
		return a * b;
	}

	/**
	Calculate a/const itemwise.
	*/
	template<typename T> std::vector<double> operator/(const std::vector<T>& a, double b)
	{
		std::vector<double> result;
		result.reserve(a.size());
		for(size_t n = 0; n < a.size(); n++)
		{
			result.push_back(a[n] / b);
		}
		return result;
	}

	/**
	Calculate a/b itemwise.
	*/
	template<typename T1, typename T2, typename Tresult> std::vector<Tresult> divide(const std::vector<T1>& a, const std::vector<T2>& b)
	{
		if(a.size() != b.size())
			throw ITLException("Sizes of std::vectors do not match.");

		std::vector<Tresult> result;
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
	template<typename T, typename Tresult = T, typename Treal = double> Tresult mean(const std::vector<T>& a)
	{
		return pixelRound<Tresult>(sum<T, Tresult>(a) / (Treal)a.size());
	}

	/**
	Calculates mean of values of the vector using given weights.
	*/
	template<typename T1, typename T2> T1 mean(const std::vector<T1>& a, const std::vector<T2>& weights)
	{
		return pixelRound<T1>(sum<double, double>(multiply<T1, T2, double>(a, weights)) / sum<T2, double>(weights));
	}

	template<typename T, typename Tresult = double> Tresult variance(const std::vector<T>& a)
	{
		double mu = mean<T, double>(a);
		auto diff = cast<double, T>(a) - mu;
		return pixelRound<Tresult>(mean<double, double>(multiply<double, double, double>(diff, diff)));
	}

	template<typename T> T stddev(const std::vector<T>& a)
	{
		return sqrt(variance(a));
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
	template<typename T> T mode(const std::vector<T>& data)
	{
		if (data.size() <= 0)
			throw std::logic_error("Mode of no data");

		std::map<T, size_t> counts;
		for (const T& x : data)
			counts[x]++;

		return mapMaxElement(counts)->first;
	}


	

	/**
	Round all values in a vector.
	*/
	template<typename Tout, typename Tin> std::vector<Tout> roundv(const std::vector<Tin>& a)
	{
		std::vector<Tout> result;
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
	template<typename Tout, typename Tin> std::vector<Tout> ceilv(const std::vector<Tin>& a)
	{
		std::vector<Tout> result;
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
	template<typename T> std::vector<T> createv(const T& a, const T& b, const T& c)
	{
		std::vector<T> result;
		result.reserve(3);
		result.push_back(a);
		result.push_back(b);
		result.push_back(c);
		return result;
	}

}
