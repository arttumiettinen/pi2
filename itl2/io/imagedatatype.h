#pragma once

#include "datatypes.h"

namespace itl2
{
	/*
	Enumerates supported pixel data types.
	*/
	enum ImageDataType
	{
		/**
		Unknown pixel data type.
		*/
		Unknown = 0,
		/**
		Unsigned 8-bit integer pixel data type.
		*/
		UInt8 = 1,
		/**
		Unsigned 16-bit integer pixel data type.
		*/
		UInt16 = 2,
		/**
		Unsigned 32-bit integer pixel data type.
		*/
		UInt32 = 3,
		/**
		Unsigned 64-bit integer pixel data type.
		*/
		UInt64 = 4,
		/**
		32-bit floating point pixel data type.
		*/
		Float32 = 5,
		/**
		Complex value consisting of two 32-bit floating point values.
		*/
		Complex32 = 6
	};

	/*
	Use to convert data type to ImageDataType value.
	*/
	template<typename pixel_t> ImageDataType imageDataType()
	{
		//pixel_t::error_not_implemented_parameter_data_type;
		return Unknown;
	}

	template<> inline ImageDataType imageDataType<uint8_t>()
	{
		return UInt8;
	}

	template<> inline ImageDataType imageDataType<uint16_t>()
	{
		return UInt16;
	}

	template<> inline ImageDataType imageDataType<uint32_t>()
	{
		return UInt32;
	}

	template<> inline ImageDataType imageDataType<uint64_t>()
	{
		return UInt64;
	}

	template<> inline ImageDataType imageDataType<float32_t>()
	{
		return Float32;
	}

	template<> inline ImageDataType imageDataType<complex32_t>()
	{
		return Complex32;
	}

	inline ImageDataType fromString(const string& dt)
	{
		if (dt == "uint8")
			return UInt8;
		if (dt == "uint16")
			return UInt16;
		if (dt == "uint32")
			return UInt32;
		if (dt == "uint64")
			return UInt64;
		if (dt == "float32")
			return Float32;
		if (dt == "complex32")
			return Complex32;
		return Unknown;
	}

	inline string toString(ImageDataType dt)
	{
		switch (dt)
		{
		case UInt8: return "uint8";
		case UInt16: return "uint16";
		case UInt32: return "uint32";
		case UInt64: return "uint64";
		case Float32: return "float32";
		case Complex32: return "complex32";
		default: return "Unknown";
		}
	}
}
