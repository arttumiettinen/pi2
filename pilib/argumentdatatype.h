#pragma once

#include <stdexcept>
#include <string>

#include "datatypes.h"
#include "neighbourhoodtype.h"
#include "boundarycondition.h"
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
	enum ArgumentDataType
	{
		String,
		Double,
		Int,
		Size,
		NBType,
		BoundaryCond,
		Bool,
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
		return String;
	}

	template<> inline ArgumentDataType parameterType<double>()
	{
		return Double;
	}

	template<> inline ArgumentDataType parameterType<coord_t>()
	{
		return Int;
	}

	template<> inline ArgumentDataType parameterType<size_t>()
	{
		return Size;
	}

	template<> inline ArgumentDataType parameterType<NeighbourhoodType>()
	{
		return NBType;
	}

	template<> inline ArgumentDataType parameterType<BoundaryCondition>()
	{
		return BoundaryCond;
	}

	template<> inline ArgumentDataType parameterType<bool>()
	{
		return Bool;
	}

	template<> inline ArgumentDataType parameterType<Image<uint8_t> >()
	{
		return ImageUInt8;
	}

	template<> inline ArgumentDataType parameterType<Image<uint16_t> >()
	{
		return ImageUInt16;
	}

	template<> inline ArgumentDataType parameterType<Image<uint32_t> >()
	{
		return ImageUInt32;
	}

	template<> inline ArgumentDataType parameterType<Image<uint64_t> >()
	{
		return ImageUInt64;
	}

	template<> inline ArgumentDataType parameterType<Image<float32_t> >()
	{
		return ImageFloat32;
	}

	template<> inline ArgumentDataType parameterType<Image<complex32_t> >()
	{
		return ImageComplex32;
	}

	template<> inline ArgumentDataType parameterType<Vec3d>()
	{
		return Vect3d;
	}

	template<> inline ArgumentDataType parameterType<Vec3c>()
	{
		return Vect3c;
	}

	template<> inline ArgumentDataType parameterType<DistributedImage<uint8_t> >()
	{
		return DImageUInt8;
	}

	template<> inline ArgumentDataType parameterType<DistributedImage<uint16_t> >()
	{
		return DImageUInt16;
	}

	template<> inline ArgumentDataType parameterType<DistributedImage<uint32_t> >()
	{
		return DImageUInt32;
	}

	template<> inline ArgumentDataType parameterType<DistributedImage<uint64_t> >()
	{
		return DImageUInt64;
	}

	template<> inline ArgumentDataType parameterType<DistributedImage<float32_t> >()
	{
		return DImageFloat32;
	}

	template<> inline ArgumentDataType parameterType<DistributedImage<complex32_t> >()
	{
		return DImageComplex32;
	}

	/*
	Convert parameter data type to string.
	*/
	static inline string toString(ArgumentDataType t)
	{
		if (t == String)
			return "string";
		if (t == Double)
			return "real";
		if (t == Int)
			return "integer";
		if (t == Size)
			return "positive integer";
		if (t == NBType)
			return "neighbourhood type";
		if (t == BoundaryCond)
			return "boundary condition";
		if (t == Bool)
			return "boolean";
		if (t == ImageUInt8)
			return "uint8 image";
		if (t == ImageUInt16)
			return "uint16 image";
		if (t == ImageUInt32)
			return "uint32 image";
		if (t == ImageUInt64)
			return "uint64 image";
		if (t == ImageFloat32)
			return "float32 image";
		if (t == ImageComplex32)
			return "complex32 image";
		if (t == Vect3d)
			return "3-component real vector";
		if (t == Vect3c)
			return "3-component integer vector";
		if (t == DImageUInt8)
			return "distributed uint8 image";
		if (t == DImageUInt16)
			return "distributed uint16 image";
		if (t == DImageUInt32)
			return "distributed uint32 image";
		if (t == DImageUInt64)
			return "distributed uint64 image";
		if (t == DImageFloat32)
			return "distributed float32 image";
		if (t == DImageComplex32)
			return "distributed complex32 image";
		throw runtime_error("Not implemented");
	}

	/**
	Gets pixel size of image of given data type.
	Returns 0 if the data type does not describe image.
	*/
	static inline size_t pixelSize(ArgumentDataType t)
	{
		if (t == ImageUInt8 || t == DImageUInt8)
			return 1;
		if (t == ImageUInt16 || t == DImageUInt16)
			return 2;
		if (t == ImageUInt32)
			return 4;
		if (t == ImageUInt64)
			return 8;
		if (t == ImageFloat32 || t == DImageFloat32)
			return 4;
		if (t == ImageComplex32 || t == DImageComplex32)
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

	/*
	Encapsulates value of an argument.
	*/
	struct ParamVariant
	{
		string userPassedString;
		union
		{
			string* sval;
			double dval;
			coord_t ival;
			bool bval;
			ImageBase* imgval;
			Image<uint8_t>* img8;
			Image<uint16_t>* img16;
			Image<uint32_t>* img32;
			Image<uint64_t>* img64;
			Image<float32_t>* imgf32;
			Image<complex32_t>* imgc32;
			DistributedImageBase* dimgval;
			DistributedImage<uint8_t>* dimg8;
			DistributedImage<uint16_t>* dimg16;
			DistributedImage<uint32_t>* dimg32;
			DistributedImage<uint64_t>* dimg64;
			DistributedImage<float32_t>* dimgf32;
			DistributedImage<complex32_t>* dimgc32;
			NeighbourhoodType nbtval;
			BoundaryCondition bcval;
			struct
			{
				double vx, vy, vz;
			};
			struct
			{
				coord_t vix, viy, viz;
			};
		};

	};

	template<typename out_t> out_t get(ParamVariant& v)
	{
		out_t::this_data_type_is_not_implemented_in_ParamVariant_get;
	}

	template<typename out_t> const out_t get(const ParamVariant& v)
	{
		out_t::this_data_type_is_not_implemented_in_ParamVariant_get;
	}

	template<> inline string get(ParamVariant& v)
	{
		return *v.sval;
	}

	template<> inline const string get(const ParamVariant& v)
	{
		return *v.sval;
	}

	template<> inline coord_t get(ParamVariant& v)
	{
		return v.ival;
	}

	template<> inline const coord_t get(const ParamVariant& v)
	{
		return v.ival;
	}

	template<> inline double get(ParamVariant& v)
	{
		return v.dval;
	}

	template<> inline const double get(const ParamVariant& v)
	{
		return v.dval;
	}

	template<> inline bool get(ParamVariant& v)
	{
		return v.bval;
	}

	template<> inline const bool get(const ParamVariant& v)
	{
		return v.bval;
	}

	template<> inline NeighbourhoodType get(ParamVariant& v)
	{
		return v.nbtval;
	}

	template<> inline const NeighbourhoodType get(const ParamVariant& v)
	{
		return v.nbtval;
	}

	template<> inline BoundaryCondition get(ParamVariant& v)
	{
		return v.bcval;
	}

	template<> inline const BoundaryCondition get(const ParamVariant& v)
	{
		return v.bcval;
	}

	template<> inline Image<uint8_t>* get(ParamVariant& v)
	{
		return v.img8;
	}

	template<> inline Image<uint8_t>* const get(const ParamVariant& v)
	{
		return v.img8;
	}

	template<> inline Image<uint16_t>* get(ParamVariant& v)
	{
		return v.img16;
	}

	template<> inline Image<uint16_t>* const get(const ParamVariant& v)
	{
		return v.img16;
	}

	template<> inline Image<uint32_t>* get(ParamVariant& v)
	{
		return v.img32;
	}

	template<> inline Image<uint32_t>* const get(const ParamVariant& v)
	{
		return v.img32;
	}

	template<> inline Image<uint64_t>* get(ParamVariant& v)
	{
		return v.img64;
	}

	template<> inline Image<uint64_t>* const get(const ParamVariant& v)
	{
		return v.img64;
	}

	template<> inline Image<float32_t>* get(ParamVariant& v)
	{
		return v.imgf32;
	}
	
	template<> inline Image<float32_t>* const get(const ParamVariant& v)
	{
		return v.imgf32;
	}

	template<> inline Image<complex32_t>* get(ParamVariant& v)
	{
		return v.imgc32;
	}

	template<> inline Image<complex32_t>* const get(const ParamVariant& v)
	{
		return v.imgc32;
	}

	template<> inline DistributedImage<uint8_t>* get(ParamVariant& v)
	{
		return v.dimg8;
	}

	template<> inline DistributedImage<uint8_t>* const get(const ParamVariant& v)
	{
		return v.dimg8;
	}

	template<> inline DistributedImage<uint16_t>* get(ParamVariant& v)
	{
		return v.dimg16;
	}

	template<> inline DistributedImage<uint16_t>* const get(const ParamVariant& v)
	{
		return v.dimg16;
	}

	template<> inline DistributedImage<uint32_t>* get(ParamVariant& v)
	{
		return v.dimg32;
	}

	template<> inline DistributedImage<uint32_t>* const get(const ParamVariant& v)
	{
		return v.dimg32;
	}

	template<> inline DistributedImage<uint64_t>* get(ParamVariant& v)
	{
		return v.dimg64;
	}

	template<> inline DistributedImage<uint64_t>* const get(const ParamVariant& v)
	{
		return v.dimg64;
	}

	template<> inline DistributedImage<float32_t>* get(ParamVariant& v)
	{
		return v.dimgf32;
	}

	template<> inline DistributedImage<float32_t>* const get(const ParamVariant& v)
	{
		return v.dimgf32;
	}

	template<> inline DistributedImage<complex32_t>* get(ParamVariant& v)
	{
		return v.dimgc32;
	}

	template<> inline DistributedImage<complex32_t>* const get(const ParamVariant& v)
	{
		return v.dimgc32;
	}

	template<> inline Vec3d get(ParamVariant& v)
	{
		return Vec3d(v.vx, v.vy, v.vz);
	}

	template<> inline const Vec3d get(const ParamVariant& v)
	{
		return Vec3d(v.vx, v.vy, v.vz);
	}

	template<> inline Vec3c get(ParamVariant& v)
	{
		return Vec3c(v.vix, v.viy, v.viz);
	}

	template<> inline const Vec3c get(const ParamVariant& v)
	{
		return Vec3c(v.vix, v.viy, v.viz);
	}

	// Extra getters
	template<> inline size_t get(ParamVariant& v)
	{
		if (v.ival < 0)
			return 0;
		return (size_t)v.ival;
	}
	template<> inline const size_t get(const ParamVariant& v)
	{
		if (v.ival < 0)
			return 0;
		return (size_t)v.ival;
	}

	/*
	Determines whether parameter is input or output parameter.
	Input images must exist before call to the function.
	If output image does not exist, it is created before call to the function (with size 1 x 1 x 1).
	*/
	enum ParameterDirection
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
		if (dir == In)
			return "in";
		else if (dir == Out)
			return "out";
		else
			return "in & out";
	}
}
