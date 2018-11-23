#pragma once

#include <stdlib.h>
#include <cstddef>
#include <stdint.h>
#include <complex>

using namespace std::literals;
using namespace std::complex_literals;
using namespace std::literals::complex_literals;

namespace itl2
{

    /**
    Defines type for image coordinates.
    This type is signed integer.
    */
    typedef long long coord_t;

    ///**
    //Type for image size and pixel index.
    //*/
    // typedef size_t size_t

    /**
    Type for 32bit float images.
    */
    typedef float float32_t;

	/**
	8-bit unsigned integer.
	This is defined in stdint.h.
	*/
	typedef unsigned char uint8_t;

	/**
	16-bit unsigned integer.
	This is defined in stdint.h.
	*/
	typedef unsigned short uint16_t;

	/**
	32-bit unsigned integer.
	This is defined in stdint.h.
	*/
	typedef unsigned int uint32_t;

	/**
	Complex number consisting of two 32-bit floats.
	*/
	typedef std::complex<float32_t> complex32_t;



// Byte order conversion methods for all data types.

    /**
    Swap the two given values.
    */
    template<typename T> void swap(T& a, T& b)
    {
        T temp = a;
        a = b;
        b = temp;
    }

	/**
     * Swaps byte order of 16-bit value.
     */
    inline uint16_t swapByteOrder(uint16_t value)
    {
        union Be
        {
            struct
            {
                uint8_t b0;
                uint8_t b1;
            };
            uint16_t word;
        };

        Be conv;
        conv.word = value;
        swap(conv.b0, conv.b1);
        return conv.word;
    }

	/**
     * Swaps byte order of 32-bit value.
     */
    inline uint32_t swapByteOrder(uint32_t value)
    {
        union Be
        {
            struct
            {
                uint8_t b0;
                uint8_t b1;
                uint8_t b2;
                uint8_t b3;
            };
            uint32_t flo;
        };

        Be conv;
        conv.flo = value;
        swap(conv.b0, conv.b3);
        swap(conv.b1, conv.b2);
        return conv.flo;
    }

    /**
     * Swaps byte order of 32-bit value.
     */
    inline float32_t swapByteOrder(float32_t value)
    {
        union Be
        {
            struct
            {
                uint8_t b0;
                uint8_t b1;
                uint8_t b2;
                uint8_t b3;
            };
            float32_t flo;
        };

        Be conv;
        conv.flo = value;
        swap(conv.b0, conv.b3);
        swap(conv.b1, conv.b2);
        return conv.flo;
    }

    // NOTE: These assume the system is little endian!

    inline uint8_t toLittleEndian(uint8_t value)
    {
        return value;
    }

    inline uint8_t toBigEndian(uint8_t value)
    {
        return value;
    }

    inline uint16_t toLittleEndian(uint16_t bigEndianValue)
    {
        return swapByteOrder(bigEndianValue);
    }

	inline uint32_t toLittleEndian(uint32_t bigEndianValue)
    {
        return swapByteOrder(bigEndianValue);
    }

    inline uint16_t toBigEndian(uint16_t littleEndianValue)
    {
        return swapByteOrder(littleEndianValue);
    }

	inline uint32_t toBigEndian(uint32_t littleEndianValue)
    {
        return swapByteOrder(littleEndianValue);
    }

    inline float32_t toLittleEndian(float32_t bigEndianValue)
    {
        return swapByteOrder(bigEndianValue);
    }

    inline float32_t toBigEndian(float32_t littleEndianValue)
    {
        return swapByteOrder(littleEndianValue);
    }
}
