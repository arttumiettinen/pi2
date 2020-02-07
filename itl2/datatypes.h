#pragma once

#include <stdlib.h>
#include <cstddef>
#include <stdint.h>
#include <complex>
#include <utility>

using namespace std::literals;
using namespace std::complex_literals;
using namespace std::literals::complex_literals;

namespace itl2
{

    /**
    Type for 32bit float images.
    */
    typedef float float32_t;

	/**
	Complex number consisting of two 32-bit floats.
	*/
	typedef std::complex<float32_t> complex32_t;
	
	/**
    Defines type for image coordinates.
    This type is signed integer.
    */
    typedef int64_t coord_t;




// Byte order conversion methods for all data types.

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


	//inline uint8_t swapByteOrder(uint8_t value)
	//{
	//	return value;
	//}

	///**
 //    * Swaps byte order of 16-bit value.
 //    */
 //   inline uint16_t swapByteOrder(uint16_t value)
 //   {
 //       union Be
 //       {
 //           struct
 //           {
 //               uint8_t b0;
 //               uint8_t b1;
 //           };
 //           uint16_t word;
 //       };

 //       Be conv;
 //       conv.word = value;
 //       std::swap(conv.b0, conv.b1);
 //       return conv.word;
 //   }

	///**
 //    * Swaps byte order of 32-bit value.
 //    */
 //   inline uint32_t swapByteOrder(uint32_t value)
 //   {
 //       union Be
 //       {
 //           struct
 //           {
 //               uint8_t b0;
 //               uint8_t b1;
 //               uint8_t b2;
 //               uint8_t b3;
 //           };
 //           uint32_t flo;
 //       };

 //       Be conv;
 //       conv.flo = value;
	//	std::swap(conv.b0, conv.b3);
	//	std::swap(conv.b1, conv.b2);
 //       return conv.flo;
 //   }

 //   /**
 //    * Swaps byte order of 32-bit value.
 //    */
 //   inline float32_t swapByteOrder(float32_t value)
 //   {
 //       union Be
 //       {
 //           struct
 //           {
 //               uint8_t b0;
 //               uint8_t b1;
 //               uint8_t b2;
 //               uint8_t b3;
 //           };
 //           float32_t flo;
 //       };

 //       Be conv;
 //       conv.flo = value;
	//	std::swap(conv.b0, conv.b3);
	//	std::swap(conv.b1, conv.b2);
 //       return conv.flo;
 //   }

 //   // NOTE: These assume the system is little endian!

 //   inline uint8_t toLittleEndian(uint8_t value)
 //   {
 //       return value;
 //   }

 //   inline uint8_t toBigEndian(uint8_t value)
 //   {
 //       return value;
 //   }

 //   inline uint16_t toLittleEndian(uint16_t bigEndianValue)
 //   {
 //       return swapByteOrder(bigEndianValue);
 //   }

	//inline uint32_t toLittleEndian(uint32_t bigEndianValue)
 //   {
 //       return swapByteOrder(bigEndianValue);
 //   }

 //   inline uint16_t toBigEndian(uint16_t littleEndianValue)
 //   {
 //       return swapByteOrder(littleEndianValue);
 //   }

	//inline uint32_t toBigEndian(uint32_t littleEndianValue)
 //   {
 //       return swapByteOrder(littleEndianValue);
 //   }

 //   inline float32_t toLittleEndian(float32_t bigEndianValue)
 //   {
 //       return swapByteOrder(bigEndianValue);
 //   }

 //   inline float32_t toBigEndian(float32_t littleEndianValue)
 //   {
 //       return swapByteOrder(littleEndianValue);
 //   }
}
