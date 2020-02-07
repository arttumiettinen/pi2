#pragma once

#include "io/imagedatatype.h"
#include "parseexception.h"

namespace pilib
{
	/**
	Calls F<data_type>::run(args...), where data_type is data type corresponding to dt enumeration value.
	Use instead of long if/switch chains where ImageDataType must be converted to real type.
	The run method must be static.
	*/
	template<template<class> class F, class... Args> void pick(ImageDataType dt, Args&&... args)
	{
		switch (dt)
		{
		case itl2::ImageDataType::UInt8: F<uint8_t>::run(args...); break;
		case itl2::ImageDataType::UInt16: F<uint16_t>::run(args...); break;
		case itl2::ImageDataType::UInt32: F<uint32_t>::run(args...); break;
		case itl2::ImageDataType::UInt64: F<uint64_t>::run(args...); break;
		case itl2::ImageDataType::Int8: F<int8_t>::run(args...); break;
		case itl2::ImageDataType::Int16: F<int16_t>::run(args...); break;
		case itl2::ImageDataType::Int32: F<int32_t>::run(args...); break;
		case itl2::ImageDataType::Int64: F<int64_t>::run(args...); break;
		case itl2::ImageDataType::Float32: F<float32_t>::run(args...); break;
		case itl2::ImageDataType::Complex32: F<complex32_t>::run(args...); break;
		default: throw ParseException(string("Unsupported data type: ") + toString(dt));
		}
	}

	/**
	Parses data type string dts and calls F<parsed_data_type>::run(args...).
	Use instead of long if/switch chains where ImageDataType must be converted to real type.
	The run method must be static.
	*/
	template<template<class> class F, class... Args> void pick(const string& dts, Args&&... args)
	{
		itl2::ImageDataType dt = itl2::fromString<itl2::ImageDataType>(dts);

		if (dt == itl2::ImageDataType::Unknown)
			throw ParseException(string("Invalid data type: ") + dts);

		pick<F, Args...>(dt, std::forward<Args>(args)...);
	}
}
