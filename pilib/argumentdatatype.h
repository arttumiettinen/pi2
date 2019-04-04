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
using std::runtime_error;
using itl2::uint8_t;
using itl2::uint16_t;
using itl2::float32_t;
using itl2::coord_t;
using itl2::NeighbourhoodType;
using itl2::BoundaryCondition;
using itl2::Connectivity;
using itl2::Image;
using math::Vec3c;
using math::Vec3d;

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
		ImageFloat32,
		ImageComplex32,
		Vect3d,
		Vect3c,
		DImageUInt8,
		DImageUInt16,
		DImageUInt32,
		DImageUInt64,
		DImageFloat32,
		DImageComplex32,
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

	template<> inline ArgumentDataType parameterType<DistributedImage<uint8_t> >()
	{
		return ArgumentDataType::DImageUInt8;
	}

	template<> inline ArgumentDataType parameterType<DistributedImage<uint16_t> >()
	{
		return ArgumentDataType::DImageUInt16;
	}

	template<> inline ArgumentDataType parameterType<DistributedImage<uint32_t> >()
	{
		return ArgumentDataType::DImageUInt32;
	}

	template<> inline ArgumentDataType parameterType<DistributedImage<uint64_t> >()
	{
		return ArgumentDataType::DImageUInt64;
	}

	template<> inline ArgumentDataType parameterType<DistributedImage<float32_t> >()
	{
		return ArgumentDataType::DImageFloat32;
	}

	template<> inline ArgumentDataType parameterType<DistributedImage<complex32_t> >()
	{
		return ArgumentDataType::DImageComplex32;
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
		if (t == ArgumentDataType::ImageFloat32)
			return "float32 image";
		if (t == ArgumentDataType::ImageComplex32)
			return "complex32 image";
		if (t == ArgumentDataType::Vect3d)
			return "3-component real vector";
		if (t == ArgumentDataType::Vect3c)
			return "3-component integer vector";
		if (t == ArgumentDataType::DImageUInt8)
			return "distributed uint8 image";
		if (t == ArgumentDataType::DImageUInt16)
			return "distributed uint16 image";
		if (t == ArgumentDataType::DImageUInt32)
			return "distributed uint32 image";
		if (t == ArgumentDataType::DImageUInt64)
			return "distributed uint64 image";
		if (t == ArgumentDataType::DImageFloat32)
			return "distributed float32 image";
		if (t == ArgumentDataType::DImageComplex32)
			return "distributed complex32 image";
		throw runtime_error("Not implemented");
	}

	/**
	Gets pixel size of image of given data type.
	Returns 0 if the data type does not describe image.
	*/
	static inline size_t pixelSize(ArgumentDataType t)
	{
		if (t == ArgumentDataType::ImageUInt8 || t == ArgumentDataType::DImageUInt8)
			return 1;
		if (t == ArgumentDataType::ImageUInt16 || t == ArgumentDataType::DImageUInt16)
			return 2;
		if (t == ArgumentDataType::ImageUInt32)
			return 4;
		if (t == ArgumentDataType::ImageUInt64)
			return 8;
		if (t == ArgumentDataType::ImageFloat32 || t == ArgumentDataType::DImageFloat32)
			return 4;
		if (t == ArgumentDataType::ImageComplex32 || t == ArgumentDataType::DImageComplex32)
			return 2*4;

		return 0;
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
		NeighbourhoodType,
		BoundaryCondition,
		Connectivity,
		InterpolationMode,
		Vec3d,
		Vec3c,
		Image<uint8_t>*, Image<uint16_t>*, Image<uint32_t>*, Image<uint64_t>*, Image<float32_t>*, Image<complex32_t>*,
		DistributedImage<uint8_t>*, DistributedImage<uint16_t>*, DistributedImage<uint32_t>*, DistributedImage<uint64_t>*, DistributedImage<float32_t>*, DistributedImage<complex32_t>*>
		ParamVariant;

	/**
	Gets DistributedImage* from ParamVariant. Throws ITLException if the variant does not contain any DistributedImage*.
	*/
	static inline DistributedImageBase* getDistributedImage(ParamVariant& v)
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
