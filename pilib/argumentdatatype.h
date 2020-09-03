#pragma once

#include <stdexcept>
#include <string>

#include <variant>

#include "datatypes.h"
#include "neighbourhoodtype.h"
#include "boundarycondition.h"
#include "connectivity.h"
#include "interpolationmode.h"
#include "image.h"
#include "distributedimage.h"
#include "math/vec3.h"

using std::string;
using itl2::float32_t;
using itl2::complex32_t;
using itl2::coord_t;
using itl2::NeighbourhoodType;
using itl2::BoundaryCondition;
using itl2::Connectivity;
using itl2::Image;
using itl2::InterpolationMode;
using itl2::Vec3d;
using itl2::Vec3c;



namespace pilib
{
	/**
	Command parameter data types.
	If you add a data type, change also
	toString, ParamVariant, and PISystem::tryConvert.
	If you add an image data type, change also parameterType specializations.
	*/
	enum class ArgumentDataType
	{
		String,
		Double,
		Int,
		Size,
		Bool,
		NBType,
		BoundaryCond,
		Connectiv,
		InterpolationMode,
		ImageUInt8,
		ImageUInt16,
		ImageUInt32,
		ImageUInt64,
		ImageInt8,
		ImageInt16,
		ImageInt32,
		ImageInt64,
		ImageFloat32,
		ImageComplex32,
		Vect3d,
		Vect3c,
	};

	template<typename T> ArgumentDataType parameterType()
	{
		T::error_unimplemented_parameter_data_type;
	}

	template<> inline ArgumentDataType parameterType<string>()
	{
		return ArgumentDataType::String;
	}

	template<> inline ArgumentDataType parameterType<double>()
	{
		return ArgumentDataType::Double;
	}

	template<> inline ArgumentDataType parameterType<coord_t>()
	{
		return ArgumentDataType::Int;
	}

	template<> inline ArgumentDataType parameterType<size_t>()
	{
		return ArgumentDataType::Size;
	}

	template<> inline ArgumentDataType parameterType<NeighbourhoodType>()
	{
		return ArgumentDataType::NBType;
	}

	template<> inline ArgumentDataType parameterType<BoundaryCondition>()
	{
		return ArgumentDataType::BoundaryCond;
	}

	template<> inline ArgumentDataType parameterType<Connectivity>()
	{
		return ArgumentDataType::Connectiv;
	}

	template<> inline ArgumentDataType parameterType<InterpolationMode>()
	{
		return ArgumentDataType::InterpolationMode;
	}

	template<> inline ArgumentDataType parameterType<bool>()
	{
		return ArgumentDataType::Bool;
	}

	template<> inline ArgumentDataType parameterType<Image<uint8_t> >()
	{
		return ArgumentDataType::ImageUInt8;
	}

	template<> inline ArgumentDataType parameterType<Image<uint16_t> >()
	{
		return ArgumentDataType::ImageUInt16;
	}

	template<> inline ArgumentDataType parameterType<Image<uint32_t> >()
	{
		return ArgumentDataType::ImageUInt32;
	}

	template<> inline ArgumentDataType parameterType<Image<uint64_t> >()
	{
		return ArgumentDataType::ImageUInt64;
	}

	template<> inline ArgumentDataType parameterType<Image<int8_t> >()
	{
		return ArgumentDataType::ImageInt8;
	}

	template<> inline ArgumentDataType parameterType<Image<int16_t> >()
	{
		return ArgumentDataType::ImageInt16;
	}

	template<> inline ArgumentDataType parameterType<Image<int32_t> >()
	{
		return ArgumentDataType::ImageInt32;
	}

	template<> inline ArgumentDataType parameterType<Image<int64_t> >()
	{
		return ArgumentDataType::ImageInt64;
	}

	template<> inline ArgumentDataType parameterType<Image<float32_t> >()
	{
		return ArgumentDataType::ImageFloat32;
	}

	template<> inline ArgumentDataType parameterType<Image<complex32_t> >()
	{
		return ArgumentDataType::ImageComplex32;
	}

	template<> inline ArgumentDataType parameterType<Vec3d>()
	{
		return ArgumentDataType::Vect3d;
	}

	template<> inline ArgumentDataType parameterType<Vec3c>()
	{
		return ArgumentDataType::Vect3c;
	}

	/*
	Convert parameter data type to string.
	*/
	static inline string toString(ArgumentDataType t)
	{
		if (t == ArgumentDataType::String)
			return "string";
		if (t == ArgumentDataType::Double)
			return "real";
		if (t == ArgumentDataType::Int)
			return "integer";
		if (t == ArgumentDataType::Size)
			return "positive integer";
		if (t == ArgumentDataType::NBType)
			return "neighbourhood type";
		if (t == ArgumentDataType::BoundaryCond)
			return "boundary condition";
		if (t == ArgumentDataType::Connectiv)
			return "connectivity";
		if (t == ArgumentDataType::InterpolationMode)
			return "interpolation mode";
		if (t == ArgumentDataType::Bool)
			return "boolean";
		if (t == ArgumentDataType::ImageUInt8)
			return "uint8 image";
		if (t == ArgumentDataType::ImageUInt16)
			return "uint16 image";
		if (t == ArgumentDataType::ImageUInt32)
			return "uint32 image";
		if (t == ArgumentDataType::ImageUInt64)
			return "uint64 image";
		if (t == ArgumentDataType::ImageInt8)
			return "int8 image";
		if (t == ArgumentDataType::ImageInt16)
			return "int16 image";
		if (t == ArgumentDataType::ImageInt32)
			return "int32 image";
		if (t == ArgumentDataType::ImageInt64)
			return "int64 image";
		if (t == ArgumentDataType::ImageFloat32)
			return "float32 image";
		if (t == ArgumentDataType::ImageComplex32)
			return "complex32 image";
		if (t == ArgumentDataType::Vect3d)
			return "3-component real vector";
		if (t == ArgumentDataType::Vect3c)
			return "3-component integer vector";
		throw std::runtime_error("Not implemented");
	}

	/**
	String that can be shown to the user to list supported image data types.
	*/
	inline string listSupportedImageDataTypes()
	{
		return "uint8, uint16, uint32, uint64, int8, int16, int32, int64, float32, or complex32";
	}

	/**
	Converts argument data type to image data type.
	*/
	static inline ImageDataType argumentDataTypeToImageDataType(ArgumentDataType t)
	{
		switch (t)
		{
		case ArgumentDataType::ImageUInt8: return ImageDataType::UInt8;
		case ArgumentDataType::ImageUInt16: return ImageDataType::UInt16;
		case ArgumentDataType::ImageUInt32: return ImageDataType::UInt32;
		case ArgumentDataType::ImageUInt64: return ImageDataType::UInt64;
		case ArgumentDataType::ImageInt8: return ImageDataType::Int8;
		case ArgumentDataType::ImageInt16: return ImageDataType::Int16;
		case ArgumentDataType::ImageInt32: return ImageDataType::Int32;
		case ArgumentDataType::ImageInt64: return ImageDataType::Int64;
		case ArgumentDataType::ImageFloat32: return ImageDataType::Float32;
		case ArgumentDataType::ImageComplex32: return ImageDataType::Complex32;
		}

		return ImageDataType::Unknown;
	}

	static inline ArgumentDataType imageDataTypeToArgumentDataType(ImageDataType t)
	{
		switch (t)
		{
		case ImageDataType::UInt8: return ArgumentDataType::ImageUInt8;
		case ImageDataType::UInt16: return ArgumentDataType::ImageUInt16;
		case ImageDataType::UInt32: return ArgumentDataType::ImageUInt32;
		case ImageDataType::UInt64: return ArgumentDataType::ImageUInt64;
		case ImageDataType::Int8: return ArgumentDataType::ImageInt8;
		case ImageDataType::Int16: return ArgumentDataType::ImageInt16;
		case ImageDataType::Int32: return ArgumentDataType::ImageInt32;
		case ImageDataType::Int64: return ArgumentDataType::ImageInt64;
		case ImageDataType::Float32: return ArgumentDataType::ImageFloat32;
		case ImageDataType::Complex32: return ArgumentDataType::ImageComplex32;
		}

		throw ITLException("Unsupported image to argument type conversion.");
		//return ArgumentDataType::Int;
	}

	/**
	Gets pixel size of image of given data type.
	Returns 0 if the data type does not describe image.
	*/
	static inline size_t pixelSize(ArgumentDataType t)
	{
		return pixelSize(argumentDataTypeToImageDataType(t));
	}

	/**
	Gets a value indicating whether given data type describes an image.
	*/
	static inline bool isImage(ArgumentDataType t)
	{
		return pixelSize(t) > 0;
	}

	typedef std::variant<coord_t,
		size_t,
		double,
		bool,
		string,
		string*,
		NeighbourhoodType,
		BoundaryCondition,
		Connectivity,
		InterpolationMode,
		Vec3d,
		Vec3c,
		Image<uint8_t>*, Image<uint16_t>*, Image<uint32_t>*, Image<uint64_t>*,
		Image<int8_t>*, Image<int16_t>*, Image<int32_t>*, Image<int64_t>*,
		Image<float32_t>*,
		Image<complex32_t>*,
		DistributedImage<uint8_t>*, DistributedImage<uint16_t>*, DistributedImage<uint32_t>*, DistributedImage<uint64_t>*,
		DistributedImage<int8_t>*, DistributedImage<int16_t>*, DistributedImage<int32_t>*, DistributedImage<int64_t>*,
		DistributedImage<float32_t>*,
		DistributedImage<complex32_t>*>
		ParamVariant;

	/**
	Gets ImageBase* from ParamVariant. Returns 0 if the variant does not contain ImageBase*.
	*/
	static inline ImageBase* getImageNoThrow(ParamVariant& v)
	{
		ImageBase* p = 0;

		std::visit(
			[&p](auto& item)
				{
					using T = std::decay_t<decltype(item)>;
					if constexpr (std::is_convertible_v<T, ImageBase*>)
					{
						p = (ImageBase*)item;
					}
				},
			v);

		return p;
	}


	/**
	Gets ImageBase* from ParamVariant. Throws ITLException if the variant does not contain any ImageBase*.
	*/
	static inline ImageBase* getImage(ParamVariant& v)
	{
		ImageBase* p = getImageNoThrow(v);

		if (p == 0)
			throw ITLException("No Image found in variant.");

		return p;
	}


	/**
	Gets DistributedImageBase* from ParamVariant. Returns 0 if the variant does not contain DistributedImageBase*.
	*/
	static inline DistributedImageBase* getDistributedImageNoThrow(ParamVariant& v)
	{
		DistributedImageBase* p = 0;

		std::visit(
			[&p](auto& item)
				{
					using T = std::decay_t<decltype(item)>;
					if constexpr (std::is_convertible_v<T, DistributedImageBase*>)
					{
						p = (DistributedImageBase*)item;
					}
				},
			v);

		return p;
	}


	/**
	Gets DistributedImage* from ParamVariant. Returns 0 if the variant does not contain DistributedImage*.
	*/
	static inline const DistributedImageBase* getDistributedImageNoThrow(const ParamVariant& v)
	{
		const DistributedImageBase* p = 0;

		std::visit(
			[&p](auto& item)
				{
					using T = std::decay_t<decltype(item)>;
					if constexpr (std::is_convertible_v<T, const DistributedImageBase*>)
					{
						p = (const DistributedImageBase*)item;
					}
				},
			v);

		return p;
	}

	/**
	Gets DistributedImage* from ParamVariant. Throws ITLException if the variant does not contain any DistributedImage*.
	*/
	static inline DistributedImageBase* getDistributedImage(ParamVariant& v)
	{
		DistributedImageBase* p = getDistributedImageNoThrow(v);

		if (p == 0)
			throw ITLException("No DistributedImage found in variant.");
		
		return p;
	}

	/**
	Gets DistributedImage* from ParamVariant. Throws ITLException if the variant does not contain any DistributedImage*.
	*/
	static inline const DistributedImageBase* getDistributedImage(const ParamVariant& v)
	{
		const DistributedImageBase* p = getDistributedImageNoThrow(v);

		if (p == 0)
			throw ITLException("No DistributedImage found in variant.");

		return p;
	}

	/*
	Determines whether parameter is input or output parameter.
	Input images must exist before call to the function.
	If output image does not exist, it is created before call to the function (with size 1 x 1 x 1).
	*/
	enum class ParameterDirection
	{
		/**
		The value of the corresponding argument is used as input data for the command.
		*/
		In,
		/**
		The value of the corresponding argument is set by when processing the command.
		*/
		Out,
		/**
		The value of the corresponding argument is used as input data and changed when processing the command.
		*/
		InOut
	};

	static inline string toString(ParameterDirection dir)
	{
		if (dir == ParameterDirection::In)
			return "in";
		else if (dir == ParameterDirection::Out)
			return "out";
		else
			return "in & out";
	}
}
