#pragma once

#include "datatypes.h"
#include "utilities.h"

namespace itl2
{
	/**
	Test if the program is running on a little-endian architecture.
	*/
	inline bool isLittleEndian()
	{
		uint16_t num = 1;
		return (*(uint8_t*)&num == 1);
	}

	/**
	Test if the program is running on a big-endian architecture.
	*/
	inline bool isBigEndian()
	{
		return !isLittleEndian();
	}

	enum class Endianness
	{
		Little,
		Big
	};

	template<>
	inline std::string toString(const Endianness& x)
	{
		switch (x)
		{
		case Endianness::Little: return "Little endian";
		case Endianness::Big: return "Big endian";
		}
		throw ITLException("Invalid endianness.");
	}

	template<>
	inline Endianness fromString(const string& str0)
	{
		string str = str0;
		toLower(str);
		if (str == "little endian" || str == "little")
			return Endianness::Little;
		else if (str == "big endian" || str == "big")
			return Endianness::Big;
		
		throw ITLException(string("Invalid endianness: ") + str);
	}

	/**
	Returns the endianness of the architecture the program is currently running on.
	*/
	inline Endianness nativeByteOrder()
	{
		if (isLittleEndian())
			return Endianness::Little;
		return Endianness::Big;
	}

	/**
	Swaps byte order of the given value.
	This is from StackOverflow: https://stackoverflow.com/questions/105252/how-do-i-convert-between-big-endian-and-little-endian-values-in-c
	NOTE: The use of union here is undefined behaviour (or implementation defined?) but... it seems to work anyway.
	*/
	template <typename T> T swapByteOrder(T u)
	{
		// This does not work in gcc
		//static_assert (CHAR_BIT == 8, "CHAR_BIT != 8");

		union
		{
			T u;
			unsigned char u8[sizeof(T)];
		} source, dest;

		source.u = u;

		for (size_t k = 0; k < sizeof(T); k++)
			dest.u8[k] = source.u8[sizeof(T) - k - 1];

		return dest.u;
	}

	template<> inline complex32_t swapByteOrder(complex32_t u)
	{
		float32_t r = u.real();
		float32_t i = u.imag();
		r = swapByteOrder(r);
		i = swapByteOrder(i);
		return complex32_t(r, i);
	}

}

