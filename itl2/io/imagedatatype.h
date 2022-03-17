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
		Complex32 = 6,
		/**
		Signed 8-bit integer pixel data type.
		*/
		Int8 = 7,
		/**
		Signed 16-bit integer pixel data type.
		*/
		Int16 = 8,
		/**
		Signed 32-bit integer pixel data type.
		*/
		Int32 = 9,
		/**
		Signed 64-bit integer pixel data type.
		*/
		Int64 = 10,
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

	template<> inline ImageDataType imageDataType<int8_t>()
	{
		return ImageDataType::Int8;
	}

	template<> inline ImageDataType imageDataType<int16_t>()
	{
		return ImageDataType::Int16;
	}

	template<> inline ImageDataType imageDataType<int32_t>()
	{
		return ImageDataType::Int32;
	}

	template<> inline ImageDataType imageDataType<int64_t>()
	{
		return ImageDataType::Int64;
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
		if (dt2 == "int8")
			return ImageDataType::Int8;
		if (dt2 == "int16")
			return ImageDataType::Int16;
		if (dt2 == "int32")
			return ImageDataType::Int32;
		if (dt2 == "int64")
			return ImageDataType::Int64;
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
		case ImageDataType::Int8: return "int8";
		case ImageDataType::Int16: return "int16";
		case ImageDataType::Int32: return "int32";
		case ImageDataType::Int64: return "int64";
		case ImageDataType::Float32: return "float32";
		case ImageDataType::Complex32: return "complex32";
		default: return "Unknown";
		}
	}

	inline size_t pixelSize(ImageDataType dt)
	{
		switch (dt)
		{
		case ImageDataType::UInt8: return sizeof(uint8_t);
		case ImageDataType::UInt16: return sizeof(uint16_t);
		case ImageDataType::UInt32: return sizeof(uint32_t);
		case ImageDataType::UInt64: return sizeof(uint64_t);
		case ImageDataType::Int8: return sizeof(int8_t);
		case ImageDataType::Int16: return sizeof(int16_t);
		case ImageDataType::Int32: return sizeof(int32_t);
		case ImageDataType::Int64: return sizeof(int64_t);
		case ImageDataType::Float32: return sizeof(float32_t);
		case ImageDataType::Complex32: return sizeof(complex32_t);
		default: return 0;
		}
	}

	/**
	Calls F<data_type>::run(args...), where data_type is data type corresponding to dt enumeration value.
	Use instead of long if/switch chains where ImageDataType must be converted to real type.
	The run method must be static.
	*/
	template<template<class> class F, class... Args> void pick(ImageDataType dt, Args&&... args)
	{
		switch (dt)
		{
		case ImageDataType::UInt8: F<uint8_t>::run(args...); break;
		case ImageDataType::UInt16: F<uint16_t>::run(args...); break;
		case ImageDataType::UInt32: F<uint32_t>::run(args...); break;
		case ImageDataType::UInt64: F<uint64_t>::run(args...); break;
		case ImageDataType::Int8: F<int8_t>::run(args...); break;
		case ImageDataType::Int16: F<int16_t>::run(args...); break;
		case ImageDataType::Int32: F<int32_t>::run(args...); break;
		case ImageDataType::Int64: F<int64_t>::run(args...); break;
		case ImageDataType::Float32: F<float32_t>::run(args...); break;
		case ImageDataType::Complex32: F<complex32_t>::run(args...); break;
		default: throw ITLException(string("Unsupported data type: ") + toString(dt));
		}
	}

	/**
	Parses data type string dts and calls F<parsed_data_type>::run(args...).
	Use instead of long if/switch chains where ImageDataType must be converted to real type.
	The run method must be static.
	*/
	template<template<class> class F, class... Args> void pick(const string& dts, Args&&... args)
	{
		ImageDataType dt = fromString<ImageDataType>(dts);

		if (dt == ImageDataType::Unknown)
			throw ITLException(string("Invalid data type: ") + dts);

		pick<F, Args...>(dt, std::forward<Args>(args)...);
	}
}
