#pragma once

#include "datatypes.h"
#include "utilities.h"

namespace itl2
{
	/*
	Enumerates supported pixel data types.
	*/
	enum class ImageDataType
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
		return ImageDataType::Unknown;
	}

	template<> inline ImageDataType imageDataType<uint8_t>()
	{
		return ImageDataType::UInt8;
	}

	template<> inline ImageDataType imageDataType<uint16_t>()
	{
		return ImageDataType::UInt16;
	}

	template<> inline ImageDataType imageDataType<uint32_t>()
	{
		return ImageDataType::UInt32;
	}

	template<> inline ImageDataType imageDataType<uint64_t>()
	{
		return ImageDataType::UInt64;
	}

	template<> inline ImageDataType imageDataType<float32_t>()
	{
		return ImageDataType::Float32;
	}

	template<> inline ImageDataType imageDataType<complex32_t>()
	{
		return ImageDataType::Complex32;
	}

	template<>
	inline ImageDataType fromString(const string& dt)
	{
		string dt2 = dt;
		trim(dt2);
		toLower(dt2);
		if (dt2 == "uint8")
			return ImageDataType::UInt8;
		if (dt2 == "uint16")
			return ImageDataType::UInt16;
		if (dt2 == "uint32")
			return ImageDataType::UInt32;
		if (dt2 == "uint64")
			return ImageDataType::UInt64;
		if (dt2 == "float32")
			return ImageDataType::Float32;
		if (dt2 == "complex32")
			return ImageDataType::Complex32;
		return ImageDataType::Unknown;
	}

	inline string toString(ImageDataType dt)
	{
		switch (dt)
		{
		case ImageDataType::UInt8: return "uint8";
		case ImageDataType::UInt16: return "uint16";
		case ImageDataType::UInt32: return "uint32";
		case ImageDataType::UInt64: return "uint64";
		case ImageDataType::Float32: return "float32";
		case ImageDataType::Complex32: return "complex32";
		default: return "Unknown";
		}
	}
}
