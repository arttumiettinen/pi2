#pragma once

#include "datatypes.h"

namespace pilib
{
	enum class ValueType
	{
		String,
		Int,
		Real,
		Bool
	};

	inline string toString(const ValueType type)
	{
		switch (type)
		{
		case ValueType::String: return "string";
		case ValueType::Int: return "int";
		case ValueType::Real: return "real";
		case ValueType::Bool: return "bool";
		default: throw ITLException("Invalid value class.");
		}
	}

	/**
	Changeable single value (non-image) in the pilib.
	*/
	class Value
	{
	private:
		ValueType currentType;

	public:

		std::string stringValue;
		int64_t intValue;
		double realValue;
		bool boolValue;

		Value(ValueType type) :
			currentType(type),
			stringValue(""),
			intValue(0),
			realValue(0.0),
			boolValue(false)
		{

		}

		ValueType getType() const
		{
			return currentType;
		}
	};

}