#pragma once

namespace pilib
{
	enum class ParsedArgType
	{
		/**
		Parsed argument is a name of a named value.
		*/
		Name,
		/**
		Parsed argument is a string.
		*/
		String,
		/**
		Parsed argument is a value.
		*/
		Value
	};
}