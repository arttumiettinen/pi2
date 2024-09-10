
#include "pisystem.h"
#include "commandmacros.h"
#include "io/io.h"
#include "commandlist.h"
#include "pilibutilities.h"
#include "slurmdistributor.h"
#include "localdistributor.h"
#include "lsfdistributor.h"
#include "timing.h"

using namespace std;

namespace pilib
{

	void addSpecialSystemCommands()
	{
		ADD_ALL(NewLikeCommand);
		ADD_ALL(NewLike2Command);
		CommandList::add<NewLikeFileCommand>();
		CommandList::add<NewLikeFile2Command>();
	}

	/**
	Parses string from argument section.
	Reads characters until terminator character is found.
	Supports \[terminator] escape sequence.
	*/
	string parseString(string& argSection, char terminator)
	{
		string result = "";
		while (argSection.length() > 0)
		{
			char s = argSection[0];
			argSection.erase(argSection.begin());
			if (s == terminator)
				break;

			if (s == '\\' && argSection.length() > 0)
			{
				if (argSection[0] == terminator)
				{
					s = terminator;
					argSection.erase(argSection.begin());
				}
			}

			result += s;
		}

		return result;
	}

	/**
	Test if the string represents a valid boolean value.
	*/
	bool isValidBoolean(const string& value)
	{
		try
		{
			bool val = fromString<bool>(value);
			return true;
		}
		catch (ITLException)
		{
			return false;
		}
	}

	/**
	Parse line expected to contain function call
	funcname(param1, param2, param3, ...)
	*/
	void PISystem::parseFunctionCall(const string& line, string& name, vector<tuple<ParsedArgType, string>>& args)
	{
		args.clear();

		stringstream s;
		s << line;

		getline(s, name, '(');
		trim(name);

		if (s.good())
		{
			// There are also parameters

			string rem;
			getline(s, rem);

			if (rem[rem.length() - 1] != ')')
				throw ParseException("Expected closing bracket");

			string argSection;
			argSection = rem.substr(0, rem.length() - 1);
			trim(argSection);

			if (argSection.length() > 0)
			{
				while (argSection.length() > 0)
				{
					if (argSection[0] == '[')
					{
						// This is start of vector argument
						// Get tokens until ending ]
						char delim;
						string arg = getToken(argSection, "]", delim);

						string mid = getToken(argSection, ",)", delim);
						trim(mid);
						if (mid.length() > 0)
							throw ParseException("Trailing characters after vector value.");

						arg = arg + ']'; // There is already [ in the beginning of the arg.
						args.push_back(make_tuple(ParsedArgType::Value, arg));
					}
					else if (argSection[0] == '"')
					{
						// Start of double-quote string
						argSection.erase(argSection.begin());
						string arg = parseString(argSection, '"');
						args.push_back(make_tuple(ParsedArgType::String, arg));
						trim(argSection);
						if (argSection.length() > 0 && argSection[0] != ',')
							throw ParseException("Trailing characters ("+argSection+") after arg: "+ arg);
						if (argSection.length() > 0)
    						argSection.erase(argSection.begin());
					}
					else if (argSection[0] == '\'')
					{
						// Start of single-quote string
						argSection.erase(argSection.begin());
						string arg = parseString(argSection, '\'');
						args.push_back(make_tuple(ParsedArgType::String, arg));
						trim(argSection);
						if (argSection.length() > 0 && argSection[0] != ',')
							throw ParseException("Trailing characters ("+argSection+") after arg: "+ arg);
						if (argSection.length() > 0)
    						argSection.erase(argSection.begin());
					}
					else
					{
						// This is standard argument
						// Get tokens until ending ,
						char delim;
						string arg = getToken(argSection, ",", delim);
						trim(arg);

						// Is this a boolean value?
						ParsedArgType type = ParsedArgType::Value;

						string dummy;
						if (isValidBoolean(arg))
							type = ParsedArgType::Value;
						else if (isValidImageName(arg, dummy))
							type = ParsedArgType::Name;
						else
							type = ParsedArgType::Value;

						args.push_back(make_tuple(type, arg));
					}

					trim(argSection);
				}
			}
		}
	}



	/**
	Tests if the given value is valid image name.
	*/
	bool PISystem::isValidImageName(const string& value, string& reason)
	{
		if (value.length() <= 0)
		{
			reason = "Image name must not be empty.";
			return false;
		}

		if (!isalpha(value[0]))
		{
			reason = "Image name must start with letter. Wrong name: " + value;
			return false;
		}

		return true;
	}

	/**
	Converts image to variable name.
	*/
	string PISystem::imageName(const ImageBase* img) const
	{
		if (img)
		{
			for (auto it = images.begin(); it != images.end(); it++)
			{
				//ImageBase* ptr = pilib::getImageNoThrow(*it->second.second.get());
				ImageBase* ptr = it->second.get();
				if (ptr == img)
					return it->first;
			}
		}
		throw ITLException("Unable to find name for image.");
	}

	/**
	Converts distributed image to variable name.
	*/
	//string PISystem::distributedImageName(const DistributedImageBase* img) const
	//{
	//	return img->varName();
	//}

	/**
	Functor that extracts Vec3d from image.
	*/
	template<typename pixel_t> struct ConvertToVec3
	{
		static void run(const ImageBase* p, Vec3d& vec)
		{
			const Image<pixel_t>* img = (const Image<pixel_t>*)p;
			vec.x = (double)(*img)(0);
			vec.y = (double)(*img)(1);
			vec.z = (double)(*img)(2);
		}
	};

	template<> struct ConvertToVec3<complex32_t>
	{
		static void run(const ImageBase* p, Vec3d& vec)
		{
			const Image<complex32_t>* img = (const Image<complex32_t>*)p;
			vec.x = (double)(*img)(0).real();
			vec.y = (double)(*img)(1).real();
			vec.z = (double)(*img)(2).real();
		}
	};

	Vec3d toVec3(const ImageBase* p)
	{
		//const Image<uint8_t>* p8 = dynamic_cast<const Image<uint8_t>*>(p);
		//const Image<uint16_t>* p16 = dynamic_cast<const Image<uint16_t>*>(p);
		//const Image<uint32_t>* p32 = dynamic_cast<const Image<uint32_t>*>(p);
		//const Image<uint64_t>* p64 = dynamic_cast<const Image<uint64_t>*>(p);
		//const Image<float32_t>* pf32 = dynamic_cast<const Image<float32_t>*>(p);

		//if (p8)
		//	return Vec3d((*p8)(0), (*p8)(1), (*p8)(2));
		//if (p16)
		//	return Vec3d((*p16)(0), (*p16)(1), (*p16)(2));
		//if (p32)
		//	return Vec3d((*p32)(0), (*p32)(1), (*p32)(2));
		//if (p64)
		//	return Vec3d((double)(*p64)(0), (double)(*p64)(1), (double)(*p64)(2));
		//if (pf32)
		//	return Vec3d((*pf32)(0), (*pf32)(1), (*pf32)(2));

		//return Vec3d(0, 0, 0);
		
		Vec3d vec;
		pick<ConvertToVec3>(p->dataType(), p, vec);
		return vec;
	}

	/**
	Finds image with given name, checks that it has three pixels, and converts it to Vec3.
	*/
	bool PISystem::getImageAsVec3(const string& name, Vec3d& v)
	{
		if (images.find(name) != images.end() || distributedImgs.find(name) != distributedImgs.end())
		{
			// First get info as getImage reads the whole image into memory.
			coord_t w = 0, h = 0, d = 0;
			ImageDataType dt;
			getImageInfoNoThrow(name, w, h, d, dt);

			coord_t pixelCount = w * h * d;
			if (pixelCount == 3 && dt != ImageDataType::Complex32 && dt != ImageDataType::Unknown)
			{
				ImageBase* pValueImage = getImage(name);

				v = toVec3(pValueImage);
				return true;
			}
		}
		return false;
	}

	/**
	Tries to convert string to given data type using fromString implementation.
	Returns true if succesfull (and if dt matches target_t).
	*/
	template<typename target_t> bool trySimpleConversion(ArgumentDataType dt, const string& value, bool doConversion, ParamVariant& result, string& reason, bool allowLiterals)
	{
		if (parameterType<target_t>() != dt)
			return false;

		if (!allowLiterals)
		{
			reason = "Literal value is not allowed here.";
			return false;
		}

		try
		{
			target_t val = itl2::fromString<target_t>(value);
			if (doConversion)
				result = val;
			return true;
		}
		catch (ITLException& e)
		{
			reason = e.message();
			return false;
		}
	}

	/**
	Functor that casts a pointer to an image to another pixel type and assigns it to a ParamVariant.
	Typically used to cast from ImageBase* to Image<pixel_t>*.
	*/
	template<typename pixel_t> struct CastImage
	{
		static void run(ImageBase* p, ParamVariant& target)
		{
			target = dynamic_cast<Image<pixel_t>*>(p);
		}
	};

	/**
	Functor that casts distributed image base to another (derived) type and assigns it to ParamVariant.
	*/
	template<typename pixel_t> struct CastDistributedImage
	{
		static void run(DistributedImageBase* p, ParamVariant& target)
		{
			target = dynamic_cast<DistributedImage<pixel_t>*>(p);
		}
	};

	/**
	Try to convert the given string to match the given argument.
	The conversion is done and, e.g., images are created if doConversion flag is set to true.
	The result of the conversion is assigned to the result variable.
	If conversion is not possible, reason string is assigned an explanation of the error.
	@return 0 if there is no match, 1 if the argument type and parameter type match, and 2 if they match after creation of new images.
	*/
	int PISystem::tryConvert(tuple<ParsedArgType, string>& value, const CommandArgumentBase& type, bool doConversion, ParamVariant& result, string& reason)
	{
		ArgumentDataType targetDt = type.dataType();
		ParameterDirection direction = type.direction();
		ParsedArgType parsedType = get<0>(value);
		string parsedValue = get<1>(value);

		bool variableValuesAsPointers = false;
		bool variablesMustExist = false;
		bool allowLiterals = false;
		if (direction == ParameterDirection::In)
		{
			// Accept literal values.
			// Named variables must exist and are converted to their values.
			// Images must exist.

			allowLiterals = true;
			variableValuesAsPointers = false;
			variablesMustExist = true;
		}
		else if (direction == ParameterDirection::InOut)
		{
			// No literal values allowed.
			// Named variables must exist and are converted to pointers.
			// Images must exist.

			allowLiterals = false;
			variableValuesAsPointers = true;
			variablesMustExist = true;
		}
		else if (direction == ParameterDirection::Out)
		{
			// No literal values allowed.
			// Named variables can be created and are converted to pointers.
			// Images can be created and are converted to pointers.

			allowLiterals = false;
			variableValuesAsPointers = true;
			variablesMustExist = false;
		}
		else
		{
			throw new ITLException("Unsupported parameter direction in TryConvert.");
		}


		if (targetDt == ArgumentDataType::String)
		{
			// Decide if the value is literal or name of a variable.

			// The value is a variable name, if variable with such a name exists.
			auto it = namedValues.find(parsedValue);
			bool isExistingVarName = it != namedValues.end();

			// If the value is parsed from a quoted string, it is never a variable name.
			bool isStringLiteral = false;
			if (parsedType == ParsedArgType::String)
				isStringLiteral = true;

			string dummyReason;
			bool isNewVariableName = isValidImageName(parsedValue, dummyReason);

			if (isExistingVarName)
			{
				// The value is the name of an existing variable.
				// Such names cannot be literal string values.

				if (doConversion)
				{
					if (variableValuesAsPointers)
						result = &(it->second->stringValue);
					else
						result = it->second->stringValue;
					namedValueStore.push_back(it->second);
				}
				return 1;
			}
			else
			{
				// The value is not an existing variable name, but could be either
				// string literal or name of a new variable.

				if (isStringLiteral)
				{
					// The value has been parsed as a quoted string, so it is not a variable name.

					if (!allowLiterals)
					{
						reason = "A variable name is required here.";
						return 0;
					}

					if (doConversion)
						result = parsedValue;
					return 1;
				}
				else
				{
					// The value might be name of a new variable or string literal.
					if (isNewVariableName && !variablesMustExist)
					{
						// Name is OK to be a variable name, and new variables are allowed.

						// Create new string variable.
						if (doConversion)
						{
							shared_ptr<Value> ptr = make_shared<Value>(ValueType::String);
							if (variableValuesAsPointers)
								result = &(ptr.get()->stringValue);
							else
								result = ptr.get()->stringValue;
							namedValues[parsedValue] = ptr;
							namedValueStore.push_back(ptr);
						}
						return 2;
					}
					else
					{
						// Name is not a variable name, or creation of new variables is not allowed.
						if (!allowLiterals)
						{
							reason = "A variable name is required here.";
							return 0;
						}

						if (doConversion)
							result = parsedValue;
						return 1;
					}

				}
			}

			throw logic_error("This line of code should never be reached (String data type handling in TryConvert).");
		}
		else if (targetDt == ArgumentDataType::JSON)
		{
			try
			{
				result = nlohmann::json::parse(parsedValue);
				return 1;
			}
			catch (const nlohmann::json::parse_error& e)
			{
				reason = e.what();
				return 0;
			}
		}
		else if(targetDt == ArgumentDataType::Int ||
			targetDt == ArgumentDataType::Double ||
			targetDt == ArgumentDataType::Bool)
		{
			// Decide if the value is literal or name of a variable.

			// For these data types, literal values are either numerical or "true" or "false".
			// It should be enough to try conversion and assume variable name if conversion does not work.
			if (trySimpleConversion<bool>(targetDt, parsedValue, doConversion, result, reason, allowLiterals))
				return 1;
			if (trySimpleConversion<double>(targetDt, parsedValue, doConversion, result, reason, allowLiterals))
				return 1;
			if (trySimpleConversion<coord_t>(targetDt, parsedValue, doConversion, result, reason, allowLiterals))
				return 1;


			// At this point, the value is not literal. Find named variable.
			auto it = namedValues.find(parsedValue);
			if (it != namedValues.end())
			{
				if (doConversion)
				{
					if (variableValuesAsPointers)
					{
						if (targetDt == ArgumentDataType::Int)
							result = &(it->second->intValue);
						else if (targetDt == ArgumentDataType::Double)
							result = &(it->second->realValue);
						else if (targetDt == ArgumentDataType::Bool)
							result = &(it->second->boolValue);
						else
							throw new ITLException("Variable type not correctly implemented (1).");
					}
					else
					{
						if (targetDt == ArgumentDataType::Int)
							result = it->second->intValue;
						else if (targetDt == ArgumentDataType::Double)
							result = it->second->realValue;
						else if (targetDt == ArgumentDataType::Bool)
							result = it->second->boolValue;
						else
							throw new ITLException("Variable type not correctly implemented (2).");
					}
					namedValueStore.push_back(it->second);
				}
				return 1;
			}
			else
			{
				if (!variablesMustExist)
				{
					// Create a new variable
					if (doConversion)
					{
						shared_ptr<Value> ptr;
						if (targetDt == ArgumentDataType::Int)
						{
							ptr = make_shared<Value>(ValueType::Int);
							if(variableValuesAsPointers)
								result = &(ptr.get()->intValue);
							else
								result = ptr.get()->intValue;
						}
						else if (targetDt == ArgumentDataType::Double)
						{
							ptr = make_shared<Value>(ValueType::Real);
							if(variableValuesAsPointers)
								result = &(ptr.get()->realValue);
							else
								result = ptr.get()->realValue;
						}
						else if (targetDt == ArgumentDataType::Bool)
						{
							ptr = make_shared<Value>(ValueType::Bool);
							if(variableValuesAsPointers)
								result = &(ptr.get()->boolValue);
							else
								result = ptr.get()->boolValue;
						}
						else
							throw new ITLException("Variable type not correctly implemented (3).");
						namedValues[parsedValue] = ptr;
						namedValueStore.push_back(ptr);
					}
					return 2;
				}
			}

			reason = string("Value '") + parsedValue + "' is not a valid literal or a variable name.";
			return 0;
		}
		else if (isImage(targetDt))
		{
			// Find image with given name

			ImageDataType idt = argumentDataTypeToImageDataType(targetDt);
			Vec3c oldSize(1, 1, 1);
			if (!isDistributed())
			{
				// Non-distributed case

				auto it = images.find(parsedValue);
				if (it != images.end())
				{
					// Check that image data type is correct
					ImageBase* p = it->second.get();
					oldSize = p->dimensions();

					if (idt == p->dataType())
					{
						if (doConversion)
						{
							imageStore.push_back(it->second);
							pick<CastImage>(idt, p, result);
						}

						return 1;
					}
					else
					{
						reason = string("Expected ") + toString(idt) + string(" image, but ") + parsedValue + string(" is ") + toString(p->dataType()) + string(" image.");
					}
				}
				else
				{
					reason = string("Image ") + parsedValue + string(" does not exist.");
				}

				// Do not return 0 here, as we might be able to create a suitable image below.
				//return 0;
			}
			else
			{
				auto it = distributedImgs.find(parsedValue);
				if (it != distributedImgs.end())
				{
					// Check that image data type is correct
					DistributedImageBase* p = it->second.get();
					oldSize = p->dimensions();

					if (idt == p->dataType())
					{
						if (doConversion)
						{
							distributedImageStore.push_back(it->second);
							pick<CastDistributedImage>(idt, p, result);
						}

						return 1;
					}
					else
					{
						reason = string("Expected ") + toString(idt) + string(" image, but ") + parsedValue + string(" is ") + toString(p->dataType()) + string(" image.");
					}
				}
				else
				{
					reason = string("Image ") + parsedValue + string(" does not exist.");
				}

				// Do not return 0 here, as we might be able to create a suitable image below.
				//return 0;
			}

			if(!variablesMustExist)
			{
				// Image can be created.
				// No match found from existing images.
				if (doConversion)
				{
					if (!isDistributed())
					{
						// Create normal image
						// Retain size of old image of the same name
						pick<CreateImage>(idt, oldSize, parsedValue, this);
						auto ptr = images.at(parsedValue);
						imageStore.push_back(ptr);
						ImageBase* p = ptr.get();
						pick<CastImage>(idt, p, result);
					}
					else
					{
						// Create distributed image
						// Retain size of old image of the same name
						pick<CreateEmptyDistributedImage>(idt, parsedValue, oldSize, this);
						auto ptr = distributedImgs.at(parsedValue);
						distributedImageStore.push_back(ptr);
						DistributedImageBase* p = ptr.get();
						pick<CastDistributedImage>(idt, p, result);
					}
				}

				return 2;
			}

			// reason should be already set here.
			return 0;
		}
		else if (targetDt == ArgumentDataType::Vect3d)
		{
			// Literal value or image name

			// Here we assume that image names are effectively literals.
			if (!allowLiterals)
			{
				reason = "Literal value is not allowed here.";
				return 0;
			}

			// Literals are of form "[1.0, 2.0, 3.0]" or number "1.0" => variable name cannot be of same form than literal.
			if (!isValidImageName(parsedValue, reason))
			{
				// Must be literal
				if (trySimpleConversion<Vec3d>(targetDt, parsedValue, doConversion, result, reason, allowLiterals))
					return 1;
				return 0;
			}
			else
			{
				Vec3d v;
				if (getImageAsVec3(parsedValue, v))
				{
					if (doConversion)
						result = v;
					return 2;
				}

				reason = string("Value '") + parsedValue + "' is not a valid vector literal or name of a 3-pixel image.";
				return 0;
			}
		}
		else if (targetDt == ArgumentDataType::Vect3c)
		{
			// Literal value or image name

			// Here we assume that image names are effectively literals.
			if (!allowLiterals)
			{
				reason = "Literal value is not allowed here.";
				return 0;
			}


			// Literals are of form "[1, 2, 3]" or number "1" => variable name cannot be of same form than literal.
			if (!isValidImageName(parsedValue, reason))
			{
				// Must be literal
				if (trySimpleConversion<Vec3c>(targetDt, parsedValue, doConversion, result, reason, allowLiterals))
					return 1;
				return 0;
			}
			else
			{
				Vec3d v;
				if (getImageAsVec3(parsedValue, v))
				{
					if (doConversion)
						result = round(v);
					return 2;
				}

				reason = string("Value '") + parsedValue + "' is not a valid integer vector literal or name of a 3-pixel image.";
				return 0;
			}
		}
		else
		{
			// All other data types must be literal values.

			if (trySimpleConversion<size_t>(targetDt, parsedValue, doConversion, result, reason, allowLiterals))
				return 1;
			if (trySimpleConversion<NeighbourhoodType>(targetDt, parsedValue, doConversion, result, reason, allowLiterals))
				return 1;
			if (trySimpleConversion<BoundaryCondition>(targetDt, parsedValue, doConversion, result, reason, allowLiterals))
				return 1;
			if (trySimpleConversion<Connectivity>(targetDt, parsedValue, doConversion, result, reason, allowLiterals))
				return 1;
			if (trySimpleConversion<InterpolationMode>(targetDt, parsedValue, doConversion, result, reason, allowLiterals))
				return 1;
			return 0;
		}

		throw logic_error("This line of code should be impossible to reach. This means a data type is not properly configured in TryParse.");





		//ArgumentDataType dt = type.dataType();

		//ParsedArgType parsedType = get<0>(value);
		//string parsedValue = get<1>(value);

		//if (type.direction() == ParameterDirection::In || type.direction() == ParameterDirection::InOut)
		//{
		//	// Input parameter.
		//	// Strings must be convertible to primitives, images and variables must exist.

		//	if (parsedType == ParsedArgType::Name /* && isValidImageName(parsedValue, reason) */)
		//	{
		//		if (dt == ArgumentDataType::String ||
		//			dt == ArgumentDataType::Int ||
		//			dt == ArgumentDataType::Double ||
		//			dt == ArgumentDataType::Bool)
		//		{
		//			// Find named value

		//			auto it = namedValues.find(parsedValue);
		//			if (it != namedValues.end())
		//			{
		//				if (doConversion)
		//				{
		//					if (dt == ArgumentDataType::String)
		//						result = &(it->second->stringValue);
		//					else if (dt == ArgumentDataType::Int)
		//						result = &(it->second->intValue);
		//					else if (dt == ArgumentDataType::Double)
		//						result = &(it->second->realValue);
		//					else if (dt == ArgumentDataType::Bool)
		//						result = &(it->second->boolValue);
		//					else
		//						throw new ITLException("In or InOut value type not implemented.");
		//					namedValueStore.push_back(it->second);
		//				}
		//				return 1;
		//			}
		//			reason = string("No variable '") + parsedValue + string("' found.");
		//		}
		//		else if (isImage(dt))
		//		{
		//			// Find image with given name

		//			ImageDataType idt = argumentDataTypeToImageDataType(dt);

		//			if (!isDistributed())
		//			{
		//				// Non-distributed case

		//				auto it = images.find(parsedValue);
		//				if (it == images.end())
		//				{
		//					reason = string("Image ") + parsedValue + string(" does not exist.");
		//					return 0;
		//				}

		//				// Check that image data type is correct
		//				ImageBase* p = it->second.get();

		//				if (idt != p->dataType())
		//				{
		//					reason = string("Expected ") + toString(idt) + string(" image, but ") + parsedValue + string(" is ") + toString(p->dataType()) + string(" image.");
		//					return 0;
		//				}

		//				if (doConversion)
		//				{
		//					imageStore.push_back(it->second);
		//					pick<CastImage>(idt, p, result);
		//				}

		//				return 1;
		//			}
		//			else
		//			{
		//				auto it = distributedImgs.find(parsedValue);
		//				if (it == distributedImgs.end())
		//				{
		//					reason = string("Image ") + parsedValue + string(" does not exist.");
		//					return 0;
		//				}

		//				// Check that image data type is correct
		//				DistributedImageBase* p = it->second.get();

		//				if (idt != p->dataType())
		//				{
		//					reason = string("Expected ") + toString(idt) + string(" image, but ") + parsedValue + string(" is ") + toString(p->dataType()) + string(" image.");
		//					return 0;
		//				}

		//				if (doConversion)
		//				{
		//					distributedImageStore.push_back(it->second);
		//					pick<CastDistributedImage>(idt, p, result);
		//				}

		//				return 1;
		//			}
		//		}
		//		else if (dt == ArgumentDataType::Vect3d)
		//		{
		//			Vec3d v;
		//			if (getImageAsVec3(parsedValue, v))
		//			{
		//				if (doConversion)
		//					result = v;
		//				return 2;
		//			}
		//		}
		//		else if (dt == ArgumentDataType::Vect3c)
		//		{
		//			Vec3d v;
		//			if (getImageAsVec3(parsedValue, v))
		//			{
		//				if (doConversion)
		//					result = round(v);
		//				return 2;
		//			}
		//		}
		//	}

		//	// We accept quoted strings and, e.g., numbers as strings.
		//	if ((parsedType == ParsedArgType::String || parsedType == ParsedArgType::Value) &&
		//		trySimpleConversion<string>(dt, parsedValue, doConversion, result, reason))
		//		return 1;

		//	if (parsedType == ParsedArgType::Value)
		//	{
		//		if (trySimpleConversion<bool>(dt, parsedValue, doConversion, result, reason))
		//			return 1;
		//		if (trySimpleConversion<double>(dt, parsedValue, doConversion, result, reason))
		//			return 1;
		//		if (trySimpleConversion<Vec3d>(dt, parsedValue, doConversion, result, reason))
		//			return 1;
		//		if (trySimpleConversion<Vec3c>(dt, parsedValue, doConversion, result, reason))
		//			return 1;
		//		if (trySimpleConversion<coord_t>(dt, parsedValue, doConversion, result, reason))
		//			return 1;
		//		if (trySimpleConversion<size_t>(dt, parsedValue, doConversion, result, reason))
		//			return 1;
		//		if (trySimpleConversion<NeighbourhoodType>(dt, parsedValue, doConversion, result, reason))
		//			return 1;
		//		if (trySimpleConversion<BoundaryCondition>(dt, parsedValue, doConversion, result, reason))
		//			return 1;
		//		if (trySimpleConversion<Connectivity>(dt, parsedValue, doConversion, result, reason))
		//			return 1;
		//		if (trySimpleConversion<InterpolationMode>(dt, parsedValue, doConversion, result, reason))
		//			return 1;
		//	}
		//	return 0;
		//}
		//else
		//{
		//	// Output parameter

		//	if (parsedType != ParsedArgType::Name /*!isValidImageName(parsedValue, reason)*/)
		//		return 0;

		//	if (isImage(dt))
		//	{
		//		ImageDataType idt = argumentDataTypeToImageDataType(dt);

		//		// Any image type can be output automatically, but existing images get higher priority

		//		// Search match from existing images.
		//		Vec3c oldSize(1, 1, 1);
		//		if (!isDistributed())
		//		{
		//			auto it = images.find(parsedValue);
		//			if (it != images.end())
		//			{
		//				ImageBase* p = it->second.get();
		//				oldSize = p->dimensions();

		//				if (p->dataType() == idt)
		//				{
		//					// Match
		//					if (doConversion)
		//					{
		//						imageStore.push_back(it->second);
		//						pick<CastImage>(idt, p, result);
		//					}
		//					return 1;
		//				}
		//			}
		//		}
		//		else
		//		{
		//			auto it = distributedImgs.find(parsedValue);
		//			if (it != distributedImgs.end())
		//			{
		//				DistributedImageBase* p = it->second.get();
		//				oldSize = p->dimensions();

		//				if (p->dataType() == idt)
		//				{
		//					// Match
		//					if (doConversion)
		//					{
		//						distributedImageStore.push_back(it->second);
		//						pick<CastDistributedImage>(idt, p, result);
		//					}
		//					return 1;
		//				}
		//			}
		//		}

		//		// No match found from existing images.
		//		// Create image if it does not exist.
		//		if (doConversion)
		//		{
		//			if (!isDistributed())
		//			{
		//				// Create normal image
		//				// Retain size of old image of the same name
		//				pick<CreateImage>(idt, oldSize, parsedValue, this);
		//				auto ptr = images.at(parsedValue);
		//				imageStore.push_back(ptr);
		//				ImageBase* p = ptr.get();
		//				pick<CastImage>(idt, p, result);
		//			}
		//			else
		//			{
		//				// Create distributed image
		//				// Retain size of old image of the same name
		//				pick<CreateEmptyDistributedImage>(idt, parsedValue, oldSize, this);
		//				auto ptr = distributedImgs.at(parsedValue);
		//				distributedImageStore.push_back(ptr);
		//				DistributedImageBase* p = ptr.get();
		//				pick<CastDistributedImage>(idt, p, result);
		//			}
		//		}

		//		return 2;
		//	}
		//	else if (dt == ArgumentDataType::String ||
		//			dt == ArgumentDataType::Int ||
		//			dt == ArgumentDataType::Double ||
		//			dt == ArgumentDataType::Bool)
		//	{
		//		// Search match from existing values.
		//		auto it = namedValues.find(parsedValue);
		//		if(it != namedValues.end())
		//		{
		//			if (doConversion)
		//			{
		//				if (dt == ArgumentDataType::String)
		//					result = &(it->second->stringValue);
		//				else if (dt == ArgumentDataType::Int)
		//					result = &(it->second->intValue);
		//				else if (dt == ArgumentDataType::Double)
		//					result = &(it->second->realValue);
		//				else if (dt == ArgumentDataType::Bool)
		//					result = &(it->second->boolValue);
		//				else
		//					throw new ITLException("In or InOut value type not implemented.");
		//				namedValueStore.push_back(it->second);
		//			}
		//			return 1;
		//		}

		//		// No match found from existing values. Create new value.
		//		if (doConversion)
		//		{
		//			//shared_ptr<string> ptr = make_shared<string>("");
		//			//stringStore.push_back(ptr);
		//			//strings[value] = ptr;
		//			//result = ptr.get();
		//			shared_ptr<Value> ptr;
		//			if (dt == ArgumentDataType::String)
		//			{
		//				ptr = make_shared<Value>(ValueType::String);
		//				result = &(ptr.get()->stringValue);
		//			}
		//			else if (dt == ArgumentDataType::Int)
		//			{
		//				ptr = make_shared<Value>(ValueType::Int);
		//				result = &(ptr.get()->intValue);
		//			}
		//			else if (dt == ArgumentDataType::Double)
		//			{
		//				ptr = make_shared<Value>(ValueType::Real);
		//				result = &(ptr.get()->realValue);
		//			}
		//			else if (dt == ArgumentDataType::Bool)
		//			{
		//				ptr = make_shared<Value>(ValueType::Bool);
		//				result = &(ptr.get()->boolValue);
		//			}
		//			else
		//				throw new ITLException("In or InOut value type not implemented.");
		//			namedValues[parsedValue] = ptr;
		//			namedValueStore.push_back(ptr);
		//		}

		//		return 2;
		//	}
		//	else
		//	{
		//		// Any other output data type is not supported...
		//		throw runtime_error("Output parameter data type not implemented.");
		//	}
		//}


	}

	/**
	Checks whether supplied string values can be converted to supplied types.
	Values are not changed but the reference is non-const as required by tryConvert.
	If 0 is returned, reason parameter is assigned an explanation why this match does not succeed.
	@return 0 if there is no match, 1 if there is match, and 2 if there is match after creation of new images.
	*/
	int PISystem::matchParameterTypes(const vector<CommandArgumentBase>& types, vector<tuple<ParsedArgType, string>>& values, string& reason)
	{
		// No match if there are more values than arguments.
		if (values.size() > types.size())
		{
			reason = "Supplied more parameters than there are arguments.";
			return 0;
		}

		// Check that values are convertible to types
		int matchPriority = 1;
		for (size_t n = 0; n < values.size(); n++)
		{
			ParamVariant dummy;
			string convertReason;
			int result = tryConvert(values[n], types[n], false, dummy, convertReason);
			if (result == 0) // No match
			{
				reason = types[n].name() + ": " + convertReason;
				return 0;
			}

			// If some parameter matched with 2 and some with 1, the overall result is 2.
			if (result > matchPriority)
				matchPriority = result;
		}

		// Check that arguments that are not supplied have default values.
		for (size_t n = values.size(); n < types.size(); n++)
		{
			if (!types[n].defaultAllowed())
			{
				reason = types[n].name() + string(": No value supplied and no default value available.");
				return 0;
			}
		}

		return matchPriority;
	}

	

	/**
	Finds some command of given priority from the given list, and returns count of items with given priority.
	*/
	size_t PISystem::getByPriority(const vector<tuple<int, Command*> >& candidates, int allowedPriority, Command*& command)
	{
		size_t count = 0;

		for (size_t n = 0; n < candidates.size(); n++)
		{
			int priority = get<0>(candidates[n]);
			if (priority == allowedPriority)
			{
				command = get<1>(candidates[n]);
				count++;
			}
		}

		return count;
	}

	/**
	Execute given command with given arguments.
	Searches the command based on name and arguments, and if found, executes it.
	If not found, throws ParseException.
	*/
	void PISystem::executeCommand(const string& name, vector<tuple<ParsedArgType, string>>& args)
	{
		// Match name
		vector<Command*> candidates = CommandList::byName(name);

		if (candidates.size() <= 0)
			throw ParseException(string("Unknown command name: ") + name);

		// Match parameter types
		vector<tuple<int, Command*> > candidates2; // Priority, command
		vector<tuple<Command*, string> > nonCandidates2;
		for (size_t n = 0; n < candidates.size(); n++)
		{
			Command* c = candidates[n];

			string reason;
			int result = matchParameterTypes(c->args(), args, reason);
			if (result != 0)
				candidates2.push_back(make_tuple(result, c));
			else
				nonCandidates2.push_back(make_tuple(c, reason));
		}

		// If we have only one candidate with priority 1, we run that.
		// If we have no candidates of priority 1 and one candidate of priority 2 we run that.
		// Otherwise we have multiple command candidates and are unable to decide which one should be run.
		Command* cmd;
		size_t count = getByPriority(candidates2, 1, cmd);
		if (count <= 0)
			count = getByPriority(candidates2, 2, cmd);

		if (count != 1)
		{
			// Error: none or multiple candidates

			if (candidates2.size() <= 0)
			{
				stringstream msg;
				msg << "No command overload found that could accept the given parameters. Could be" << endl;

				// Build reason messages but remove duplicates
				vector<string> messages;
				for (size_t n = 0; n < nonCandidates2.size(); n++)
				{
					stringstream submsg;
					submsg << get<0>(nonCandidates2[n])->toSimpleString();
					submsg << " -- " << get<1>(nonCandidates2[n]);

					messages.push_back(submsg.str());
				}

				removeDuplicates(messages);

				for (size_t n = 0; n < messages.size(); n++)
				{
					msg << messages[n];
					if (n < messages.size() - 1)
						msg << endl;
				}

				throw ParseException(msg.str());
			}

			if (candidates2.size() > 1)
			{
				stringstream msg;
				msg << "Multiple command candidates. Could be " << endl;
				for (size_t n = 0; n < candidates2.size(); n++)
				{
					msg << get<1>(candidates2[n])->toString();
					if (n < candidates2.size() - 1)
						msg << endl;
				}

				throw ParseException(msg.str());
			}
		}
		//Command* cmd = get<1>(candidates2[0]);

		// Sanity check
		if (!cmd)
			throw logic_error("null command");

		// Add defaults to the parameter array
		vector<tuple<ParsedArgType, string>> realArgs;
		realArgs.reserve(cmd->args().size());
		size_t N0 = args.size();
		for (size_t n = 0; n < N0; n++)
			realArgs.push_back(args[n]);
		for (size_t n = N0; n < cmd->args().size(); n++)
		{
			// Note: default argument parsed type is never Name.
			ParsedArgType defType;
			if (cmd->args()[n].dataType() == ArgumentDataType::String && cmd->args()[n].direction() == ParameterDirection::In)
				defType = ParsedArgType::String;
			else
				defType = ParsedArgType::Value;
			realArgs.push_back(make_tuple(defType, cmd->args()[n].defaultValue()));
		}

		// Commands that should not be echoed to screen
		bool isNoShow = cmd->name() == "help" || cmd->name() == "info" || cmd->name() == "license";

		// Show command
		if (showRunCommands && !isNoShow)
		{
			cout << cmd->name() << "(";
			for (size_t n = 0; n < realArgs.size(); n++)
			{
				ParsedArgType type = get<0>(realArgs[n]);
				string val = get<1>(realArgs[n]);
				if (type == ParsedArgType::String) // Only strings should be quoted.
					cout << "\"" << val << "\"";
				else
					cout << val;
				if (n < realArgs.size() - 1)
					cout << ", ";
			}
			cout << ")" << endl;
		}

		// Convert string parameters to values. This must succeed as the process was tested above.
		vector<ParamVariant> convertedArgs;
		convertedArgs.reserve(realArgs.size());

		//vector<shared_ptr<ImageBase> > imageStore;
		
		for (size_t n = 0; n < realArgs.size(); n++)
		{
			ParamVariant res;
			string dummy;
			tryConvert(realArgs[n], cmd->args()[n], true, res, dummy);

			//// Store shared_ptrs to image arguments so that they don't get deleted if the next parameters override them.
			//// NOTE: We just hold the pointers but don't do anything with them.
			//DistributedImageBase* db = getDistributedImageNoThrow(res);
			//if (db)
			//{
			//	distributedImageStore.push_back(getDistributedImagePointer(db));
			//}

			//ImageBase* ib = pilib::getImageNoThrow(res);
			//if (ib)
			//{
			//	imageStore.push_back(getImagePointer(ib));
			//}

			//if (resPtr)
			//{
			//	imageStore.push_back(resPtr);
			//}


			convertedArgs.push_back(res);
		}

		// Run command with timing
		Timer timer;
		timer.start();
		if (!isDistributed())
		{
			// Normal processing without distribution or anything fancy
			TimingFlag flag(TimeClass::Computation);
			cmd->runInternal(this, convertedArgs);
		}
		else
		{
			// Distributed processing
			Distributable* dist = dynamic_cast<Distributable*>(cmd);
			if (dist)
			{
				// NOTE: GlobalTimer continues in Overhead mode.
				dist->runDistributed(*distributor, convertedArgs);
			}
			else
			{
				throw ITLException(string("Command ") + cmd->toSimpleString() + " has not been setup to run in a distributed computing environment.");
			}
		}

		timer.stop();

		if (showTiming && !isNoShow)
			cout << "Operation took " << setprecision(3) << timer.getSeconds() << " s" << endl;

		imageStore.clear();
		namedValueStore.clear();
		distributedImageStore.clear();
	}

	/**
	Parse one statement of input code.
	*/
	void PISystem::parseStatement(string& statement)
	{
		trim(statement);

		string cmd;
		vector<tuple<ParsedArgType, string>> args;
		parseFunctionCall(statement, cmd, args);

		//cout << "Name: " << cmd << endl;
		//for (size_t n = 0; n < args.size(); n++)
		//	cout << "Argument " << (n + 1) << ": " << args[n] << endl;

		executeCommand(cmd, args);
	}

	/**
	Parse one line of input code.
	*/
	void PISystem::parseLine(string& line)
	{
		trim(line);

		//cout << "Parsing line: " << line << endl;

		if (line.length() > 0)
		{
			if (line[0] != '%'
				&& line[0] != '#'
				&& line.find("//") != 0)
			{
				// Not a comment. Parse statements one by one.

				while (line.length() > 0)
				{
					char delim = 0;
					string token = getToken(line, ";\n", delim);

					parseStatement(token);
				}
			}

		}
	}


	/**
	Gets a string representing data type of given image.
	*/
	string PISystem::getDataTypeString(const ImageBase* p)
	{
		return toString(p->dataType()) + " image";
	}


	/**
	Gets a string representing data type of given image.
	*/
	string PISystem::getDataTypeString(const DistributedImageBase* p)
	{
		return string("distributed ") + toString(p->dataType()) + " image";
	}






	PISystem::PISystem()
	{
		clearLastError();
	}

	PISystem::~PISystem()
	{
	}

	/**
	Returns names of images in the system.
	*/
	vector<string> PISystem::getImageNames() const
	{
		std::vector<string> keys;
		keys.reserve(images.size());
		for (auto const& item : images)
		{
			//if(pilib::isImage(item.second.first))
				keys.push_back(item.first);
		}
		return keys;
	}

	std::vector<std::string> PISystem::getValueNames() const
	{
		std::vector<string> keys;
		keys.reserve(namedValues.size());
		for (auto const& item : namedValues)
		{
			keys.push_back(item.first);
		}
		return keys;
	}
	

	/**
	Returns names of distributed images in the system.
	*/
	vector<string> PISystem::getDistributedImageNames() const
	{
		std::vector<string> keys;
		keys.reserve(distributedImgs.size());
		for (auto const& item : distributedImgs)
			keys.push_back(item.first);
		return keys;
	}

	/**
	Retrieve image having given name.
	*/
	ImageBase* PISystem::getImage(const string& name)
	{
		if (isDistributed())
		{
			// Convert the distributed image to normal image and return that.
			unique_ptr<ImageBase> ptr = distributedImgs.at(name)->toNormalImage();

			//std::shared_ptr<ParamVariant> sptr = std::make_shared<ParamVariant>((coord_t)0);
			//ImageDataType dt = ptr->dataType();
			//pick<CastImage>(dt, ptr.release(), *sptr);

			//replaceImage(name, sptr);

			replaceImage(name, std::move(ptr));
		}

		//return imgs[name].get();
		//return pilib::getImage(*namedValues.at(name).second);
		return images.at(name).get();
	}

	Value* PISystem::getValue(const string& name)
	{
		return namedValues.at(name).get();
	}

	Value* PISystem::getValueNoThrow(const string& name)
	{
		return noThrow([&]
			{
				return getValue(name);
			},
			"Value with given name not found.");
	}

	/**
	Flush local changes to distributed image.
	Assumes that the distributed image has been retrieved using getImage.
	Deletes the local image.
	*/
	void PISystem::flushIfDistributed(const string& imgName)
	{
		if (isDistributed())
		{
			distributedImgs.at(imgName)->setData(images.at(imgName).get());
			//distributedImgs[imgName]->setData(pilib::getImage(*namedValues[imgName].second));
			//replaceNamedValue(imgName, ArgumentDataType::Int, 0);
		}
	}

	/**
	Same than flushIfDistributed but sets lastException instead of throwing it.
	@return false if an error occurs.
	*/
	bool PISystem::flushIfDistributedNoThrow(const string& name)
	{
		return noThrow([&]
			{
				flushIfDistributed(name);
				return true;
			});
	}

	/**
	Retrieve image having given name, do not throw exception but set last error and return 0 if exception would be thrown.
	*/
	void PISystem::getImageInfoNoThrow(const string& name, coord_t& width, coord_t& height, coord_t& depth, ImageDataType& dt)
	{
		noThrow([&]
			{
				// NOTE: Do not use getImage as that results in reading distributed images to RAM.
				if (!isDistributed())
				{
					shared_ptr<ImageBase> p = images.at(name);
					width = p->width();
					height = p->height();
					depth = p->depth();
					dt = p->dataType();
				}
				else
				{
					auto p = distributedImgs.at(name);
					width = p->width();
					height = p->height();
					depth = p->depth();
					dt = p->dataType();
				}
				return true;
			},
			"Image not found.");
	}

	/**
	Retrieve image having given name, do not throw exception but set last error and return 0 if exception would be thrown.
	*/
	ImageBase* PISystem::getImageNoThrow(const string& name)
	{
		ImageBase* img = noThrow([&]
			{
				return getImage(name);
			});

		if (!img)
			lastException = "Image not found.";

		return img;
	}



	/**
	Replace image with given image.
	Replace image by null pointer to remove it from the system.
	*/
	void PISystem::replaceImage(const string& name, shared_ptr<ImageBase> newValue)
	{
		if (images.find(name) != images.end())
			images.erase(name);

		if (newValue)
			images[name] = newValue;

		//if (img)
		//{
		//	ArgumentDataType dt = imageDataTypeToArgumentDataType(pilib::getImage(*img)->dataType());
		//	replaceNamedValue(name, dt, img);
		//}
		//else
		//{
		//	replaceNamedValue(name, ArgumentDataType::Int, img);
		//}
	}

	void PISystem::replaceValue(const std::string& name, std::shared_ptr<Value> newValue)
	{
		if (namedValues.find(name) != namedValues.end())
			namedValues.erase(name);

		if (newValue)
			namedValues[name] = newValue;
	}

	/**
	Get distributed image having given name.
	*/
	DistributedImageBase* PISystem::getDistributedImage(const string& name)
	{
		return distributedImgs[name].get();
	}

	/**
	Replace distributed image with a new one.
	Set newImg to zero to remove the image with given name.
	*/
	void PISystem::replaceDistributedImage(const string& name, shared_ptr<DistributedImageBase> newImg)
	{
		if (distributedImgs.find(name) != distributedImgs.end())
			distributedImgs.erase(name);
		
		if (newImg)
			distributedImgs[name] = newImg;
	}

	/**
	Gets error message for last exception that occured.
	*/
	const char* PISystem::getLastErrorMessage() const
	{
		return lastException.c_str();
	}

	/**
	Gets line of code that caused last error.
	*/
	int PISystem::getLastErrorLine() const
	{
		return lastExceptionLine;
	}

	/**
	Clear last error message.
	*/
	void PISystem::clearLastError()
	{
		lastException = "";
		lastExceptionLine = 0;
	}

	/**
	Sets show commands flag.
	*/
	void PISystem::showCommands(bool echo, bool timing)
	{
		showRunCommands = echo;
		showTiming = timing;
	}

	/**
	Gets a value indicating whethe distributed processing mode is active.
	*/
	bool PISystem::isDistributed() const
	{
		return distributor != 0;
	}

	/**
	Enables or disables distributed processing.
	@param provider Provider name. Pass empty string to disable distributed processing.
	*/
	void PISystem::distributedProcessing(string provider)
	{
		if (distributor)
		{
			distributor->flush();
			delete distributor;
			distributor = 0;
			cout << "Distributed computing mode is disabled." << endl;
		}

		toLower(provider);

		if (provider.length() > 0 && provider != "none" && provider != "disable" && provider != "false" && provider != "0" && provider != "off" && provider != "no")
		{
			if (provider == "slurm")
			{
				cout << "Enabling distributed computing mode using SLURM workload manager." << endl;
				distributor = new SLURMDistributor(this);
			}
			else if (provider == "local")
			{
				cout << "Enabling distributed computing mode using local sequential processing." << endl;
				distributor = new LocalDistributor(this);
			}
			else if (provider == "lsf")
			{
				cout << "Enabling distributed computing mode using LSF workload manager." << endl;
				distributor = new LSFDistributor(this);
			}
			else
				throw ITLException(string("Invalid distributed computing system name: ") + provider + ". Valid names are SLURM or Local.");
		}
	}

	/**
	Parses commands in the given string and runs them.
	*/
	bool PISystem::run(const string& commands)
	{
		if (!running)
			commandsWaiting.clear();

		commandsWaiting.push_back(commands);

		if (!running)
		{
			running = true;

			// Reset global timer by calling start explicitly.
			GlobalTimer::start();
			TimingFlag flag(TimeClass::Overhead);

			while (commandsWaiting.size() > 0)
			{
				string nextItem = commandsWaiting[0];
				commandsWaiting.erase(commandsWaiting.begin());

				try
				{
					string rest = nextItem;
					lastExceptionLine = 1;
					while (rest.length() > 0)
					{
						char delim = 0;
						string token = getToken(rest, "\n", delim);

						parseLine(token);

						if (delim == '\n')
							lastExceptionLine++;
					}
				}
				catch (ITLException& e)
				{
					lastException = e.message();
					running = false;
					return false;
				}
				catch (exception& e)
				{
					lastException = e.what();
					running = false;
					return false;
				}
			}

			lastExceptionLine = 0;
			running = false;
		}

		return true;
	}


	shared_ptr<DistributedImageBase> PISystem::getDistributedImagePointer(DistributedImageBase* img) const
	{
		for (auto& item : distributedImgs)
		{
			if (item.second.get() == img)
				return item.second;
		}

		for (auto& item : distributedImageStore)
		{
			if (item.get() == img)
				return item;
		}

		throw ITLException("Unknown distributed image pointer.");
	}

	bool PISystem::isDistributedImage(DistributedImageBase* img) const
	{
		for (auto& item : distributedImgs)
		{
			if (item.second.get() == img)
				return true;
		}

		return false;
	}



	void NewLikeFileCommand::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string filename = pop<string>(args);
		string dts = pop<string>(args);
		coord_t w = pop<coord_t>(args);
		coord_t h = pop<coord_t>(args);
		coord_t d = pop<coord_t>(args);

		ImageDataType templDT;
		Vec3c templDims;
		string reason;
		if (!itl2::io::getInfo(filename, templDims, templDT, reason))
			throw ITLException(string("Unable to find dimensions and data type of file ") + filename + ". " + reason);

		createImage(name, dts, w, h, d, templDT, templDims.x, templDims.y, templDims.z, system);
	}

	vector<string> NewLikeFileCommand::runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string filename = pop<string>(args);
		string dts = pop<string>(args);
		coord_t w = pop<coord_t>(args);
		coord_t h = pop<coord_t>(args);
		coord_t d = pop<coord_t>(args);

		ImageDataType templDT;
		Vec3c templDims;
		string reason;
		if (!itl2::io::getInfo(filename, templDims, templDT, reason))
			throw ITLException(string("Unable to find dimensions and data type of file ") + filename + ". " + reason);

		PISystem* system = distributor.getSystem();

		createDistributedImage(name, dts, w, h, d, templDT, templDims.x, templDims.y, templDims.z, system);

		return vector<string>();
	}

	void NewLikeFile2Command::runInternal(PISystem* system, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string filename = pop<string>(args);
		string dts = pop<string>(args);
		Vec3c dim = pop<Vec3c>(args);

		ImageDataType templDT;
		Vec3c templDims;
		string reason;
		if (!itl2::io::getInfo(filename, templDims, templDT, reason))
			throw ITLException(string("Unable to find dimensions and data type of file ") + filename + ". " + reason);

		createImage(name, dts, dim.x, dim.y, dim.z, templDT, templDims.x, templDims.y, templDims.z, system);
	}

	vector<string> NewLikeFile2Command::runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
	{
		string name = pop<string>(args);
		string filename = pop<string>(args);
		string dts = pop<string>(args);
		Vec3c dim = pop<Vec3c>(args);

		ImageDataType templDT;
		Vec3c templDims;
		string reason;
		if (!itl2::io::getInfo(filename, templDims, templDT, reason))
			throw ITLException(string("Unable to find dimensions and data type of file ") + filename + ". " + reason);

		PISystem* system = distributor.getSystem();

		createDistributedImage(name, dts, dim.x, dim.y, dim.z, templDT, templDims.x, templDims.y, templDims.z, system);

		return vector<string>();
	}

}
