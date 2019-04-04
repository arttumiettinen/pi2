
#include "pisystem.h"
#include "commandmacros.h"
#include "io/io.h"

namespace pilib
{
	/**
	Creates all commands and adds them to the list.
	Used in the constructor of PISystem and defined elsewhere to keep this file shorter and cleaner,
	and to cut dependencies between PISystem header and command definitions.
	*/
	void addFilterCommands(vector<Command*>& commands);
	void addSpecialCommands(vector<Command*>& commands);
	void addIOCommands(vector<Command*>& commands);
	void addOtherCommands(vector<Command*>& commands);
	void addPointProcessCommands(vector<Command*>& commands);
	void addThinAndSkeletonCommands(vector<Command*>& commands);
	void addSpecialSystemCommands(vector<Command*>& commands);
	void addTransformCommands(vector<Command*>& commands);
	void addParticlesCommands(vector<Command*>& commands);
	void addGenerationCommands(vector<Command*>& commands);
	void addProjectionCommands(vector<Command*>& commands);


	void addSpecialSystemCommands(vector<Command*>& commands)
	{
		commands.insert(commands.end(),
			{
			ADD_REAL(ConvertCommand),
			ADD_ALL(NewLikeCommand),
			new NewLikeFileCommand()
			});
	}



	/**
	Parse line expected to contain function call
	funcname(param1, param2, param3, ...)
	*/
	void PISystem::parseFunctionCall(const string& line, string& name, vector<string>& args)
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
					char delim;
					string arg = getToken(argSection, ",[", delim);

					if (delim == ',' || delim == 0)
					{
						// This is standard argument
						trim(arg);

						args.push_back(arg);
					}
					else if (delim == '[')
					{
						// This is start of vector argument
						// Get tokens until ending ]
						arg = getToken(argSection, "]", delim);
						getToken(argSection, ",", delim);
						arg = '[' + arg + ']';
						args.push_back(arg);
					}
				}

				//stringstream ss;
				//ss << argSection;
				//while (ss.good())
				//{
				//	string arg;
				//	getline(ss, arg, ',');
				//	trim(arg);
				//	args.push_back(arg);
				//}
			}
		}

	}

	/**
	Get all commands whose name is the given one.
	*/
	vector<Command*> PISystem::getCommands(const string& name) const
	{
		vector<Command*> cmds;
		for (size_t n = 0; n < commands.size(); n++)
		{
			if (commands[n]->name() == name)
				cmds.push_back(commands[n]);
		}
		return cmds;
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
			reason = "Image name must start with letter.";
			return false;
		}

		return true;
	}

	/**
	Converts image to variable name.
	*/
	string PISystem::imageName(const ImageBase* img) const
	{
		for (map<string, ImageBase*>::const_iterator it = imgs.begin(); it != imgs.end(); it++)
		{
			if (it->second == img)
				return it->first;
		}

		throw ITLException("Unable to find name for image.");
	}

	/**
	Converts distributed image to variable name.
	*/
	string PISystem::distributedImageName(const DistributedImageBase* img) const
	{
		return img->varName();
	}

	Vec3d toVec3(const ImageBase* p)
	{
		const Image<uint8_t>* p8 = dynamic_cast<const Image<uint8_t>*>(p);
		const Image<uint16_t>* p16 = dynamic_cast<const Image<uint16_t>*>(p);
		const Image<uint32_t>* p32 = dynamic_cast<const Image<uint32_t>*>(p);
		const Image<uint64_t>* p64 = dynamic_cast<const Image<uint64_t>*>(p);
		const Image<float32_t>* pf32 = dynamic_cast<const Image<float32_t>*>(p);
		//const Image<complex32_t >* pc32 = dynamic_cast<const Image<complex32_t >*>(p);

		if (p8)
			return Vec3d((*p8)(0), (*p8)(1), (*p8)(2));
		if (p16)
			return Vec3d((*p16)(0), (*p16)(1), (*p16)(2));
		if (p32)
			return Vec3d((*p32)(0), (*p32)(1), (*p32)(2));
		if (p64)
			return Vec3d((double)(*p64)(0), (double)(*p64)(1), (double)(*p64)(2));
		if (pf32)
			return Vec3d((*pf32)(0), (*pf32)(1), (*pf32)(2));
		

		return Vec3d(0, 0, 0);
	}

	/**
	Finds image with given name, checks that it has three pixels, and converts it to Vec3.
	*/
	bool PISystem::getImageAsVec3(const string& name, Vec3d& v)
	{
		if (imgs.find(name) != imgs.end() || distributedImgs.find(name) != distributedImgs.end())
		{
			ImageBase* pValueImage = getImage(name);

			if (pValueImage->pixelCount() == 3)
			{
				v = toVec3(pValueImage);
				return true;
			}
		}
		return false;
	}

	/**
	Try to convert the given string to match the given argument.
	The conversion is done and, e.g., images are created if doConversion flag is set to true.
	The result of the conversion is assigned to the result variable.
	If conversion is not possible, reason string is assigned an explanation of the error.
	@return 0 if there is no match, 1 if the argument type and parameter type match, and 2 if they match after creation of new images.
	*/
	int PISystem::tryConvert(string& value, const CommandArgumentBase& type, bool doConversion, ParamVariant& result, string& reason)
	{
		// TODO: Clean this up. This function has become a mess.

		if (type.direction() == ParameterDirection::In || type.direction() == ParameterDirection::InOut)
		{
			// Input parameter.
			// Strings must be convertible to primitives, images must exist.

			ArgumentDataType dt = type.dataType();

			if (dt == ArgumentDataType::String)
			{
				if (doConversion)
					result = value;
				return 1;
			}
			else if (dt == ArgumentDataType::Bool)
			{
				try
				{
					bool val = itl2::fromString<bool>(value);
					if (doConversion)
						result = val;
					return 1;
				}
				catch (ITLException& e)
				{
					reason = e.message();
					return 0;
				}
			}
			else if (dt == ArgumentDataType::Double)
			{
				try
				{
					double val = itl2::fromString<double>(value);
					if (doConversion)
						result = val;
					return 1;
				}
				catch (ITLException& e)
				{
					reason = e.message();
					return 0;
				}
			}
			else if (dt == ArgumentDataType::Vect3d)
			{
				try
				{
					Vec3d val = itl2::fromString<Vec3d>(value);
					if (doConversion)
					{
						result = val;
					}
					return 1;
				}
				catch (ITLException& e)
				{
					reason = e.message();
				}

				Vec3d v;
				if (getImageAsVec3(value, v))
				{
					if (doConversion)
						result = v;
					return 1;
				}

				return 0;
			}
			else if (dt == ArgumentDataType::Vect3c)
			{
				try
				{
					Vec3c val = itl2::fromString<Vec3c>(value);
					if (doConversion)
					{
						result = val;
					}
					return 1;
				}
				catch (ITLException& e)
				{
					reason = e.message();
				}

				Vec3d v;
				if (getImageAsVec3(value, v))
				{
					if (doConversion)
						result = math::round(v);
					return 1;
				}

				return 0;
			}
			else if (dt == ArgumentDataType::Int || dt == ArgumentDataType::Size)
			{
				char* p;
				coord_t val = strtol(value.c_str(), &p, 10);
				if (*p)
				{
					reason = string("Value ") + value + string(" is not a valid integer number.");
					return 0;
				}
				if (doConversion)
				{
					if (dt == ArgumentDataType::Int)
						result = val;
					else
						result = (size_t)math::max((size_t)0, (size_t)val);
				}
				return 1;
			}
			else if (dt == ArgumentDataType::NBType)
			{
				string lv = value;
				toLower(lv);

				if (lv == "rect" || lv == "rectangular" || lv == "box")
				{
					if (doConversion)
						result = NeighbourhoodType::Rectangular;
					return 1;
				}

				if (lv == "ell" || lv == "ellipsoidal" || lv == "ellipsoid" || lv == "spherical" || lv == "sphere")
				{
					if (doConversion)
						result = NeighbourhoodType::Ellipsoidal;
					return 1;
				}

				reason = string("Unknown neighbourhood type: ") + value;
				return 0;
			}
			else if (dt == ArgumentDataType::BoundaryCond)
			{
				string lv = value;
				toLower(lv);

				if (lv == "zero" || lv == "0")
				{
					if (doConversion)
						result = BoundaryCondition::Zero;
					return 1;
				}

				if (lv == "nearest" || lv == "1")
				{
					if (doConversion)
						result = BoundaryCondition::Nearest;
					return 1;
				}

				reason = string("Unknown boundary condition: ") + value;
				return 0;
			}
			else if (dt == ArgumentDataType::Connectiv)
			{
				string lv = value;
				toLower(lv);

				if (lv == "nearest" || lv == "nearestneigbours" || lv == "nearest_neighbours" ||
					lv == "nearestneigbors" || lv == "nearest_neighbors" ||
					value == "4" || value == "6" || value == "0")
				{
					if (doConversion)
						result = Connectivity::NearestNeighbours;
					return 1;
				}

				if (lv == "all" || lv == "allneigbours" || lv == "all_neighbours" ||
					lv == "allneigbors" || lv == "all_neighbors" ||
					value == "8" || value == "27" || value == "1")
				{
					if (doConversion)
						result = Connectivity::AllNeighbours;
					return 1;
				}

				reason = string("Unknown connectivity value: ") + value;
				return 0;
			}
			else if (dt == ArgumentDataType::InterpolationMode)
			{
				string lv = value;
				toLower(lv);

				if (startsWith(lv, "nearest")
					|| lv == "0" || lv == "no" || lv == "none" || lv == "off")
				{
					if (doConversion)
						result = InterpolationMode::Nearest;
					return 1;
				}

				if (startsWith(lv, "linear")
					|| lv == "1")
				{
					if (doConversion)
						result = InterpolationMode::Linear;
					return 1;
				}

				if (startsWith(lv, "cubic")
					|| lv == "2")
				{
					if (doConversion)
						result = InterpolationMode::Cubic;
					return 1;
				}

				reason = string("Unknown interpolation mode: ") + value;
				return 0;
			}
			else if (dt == ArgumentDataType::ImageUInt8 ||
				dt == ArgumentDataType::ImageUInt16 ||
				dt == ArgumentDataType::ImageUInt32 ||
				dt == ArgumentDataType::ImageUInt64 ||
				dt == ArgumentDataType::ImageFloat32 ||
				dt == ArgumentDataType::ImageComplex32)
			{
				// Find image with given name

				if (!isValidImageName(value, reason))
					return 0;

				if (!isDistributed())
				{
					// Non-distributed case

					map<string, ImageBase*>::const_iterator it = imgs.find(value);
					if (it == imgs.end())
					{
						reason = string("Image ") + value + string(" does not exist.");
						return 0;
					}

					ImageBase* p = it->second;

					// Check image data type
					if (dt == ArgumentDataType::ImageUInt8)
					{
						Image<uint8_t>* p2 = dynamic_cast<Image<uint8_t>*>(p);
						if (!p2)
						{
							reason = string("Expected ") + toString(ArgumentDataType::ImageUInt8) + string(", but ") + value + string(" is ") + getDataTypeString(p);
							return 0;
						}

						if (doConversion)
							result = p2;

						return 1;
					}

					if (dt == ArgumentDataType::ImageUInt16)
					{
						Image<uint16_t>* p2 = dynamic_cast<Image<uint16_t>*>(p);
						if (!p2)
						{
							reason = string("Expected ") + toString(ArgumentDataType::ImageUInt16) + string(", but ") + value + string(" is ") + getDataTypeString(p);
							return 0;
						}

						if (doConversion)
							result = p2;

						return 1;
					}

					if (dt == ArgumentDataType::ImageUInt32)
					{
						Image<uint32_t>* p2 = dynamic_cast<Image<uint32_t>*>(p);
						if (!p2)
						{
							reason = string("Expected ") + toString(ArgumentDataType::ImageUInt32) + string(", but ") + value + string(" is ") + getDataTypeString(p);
							return 0;
						}

						if (doConversion)
							result = p2;

						return 1;
					}

					if (dt == ArgumentDataType::ImageUInt64)
					{
						Image<uint64_t>* p2 = dynamic_cast<Image<uint64_t>*>(p);
						if (!p2)
						{
							reason = string("Expected ") + toString(ArgumentDataType::ImageUInt64) + string(", but ") + value + string(" is ") + getDataTypeString(p);
							return 0;
						}

						if (doConversion)
							result = p2;

						return 1;
					}

					if (dt == ArgumentDataType::ImageFloat32)
					{
						Image<float32_t>* p2 = dynamic_cast<Image<float32_t>*>(p);
						if (!p2)
						{
							reason = string("Expected ") + toString(ArgumentDataType::ImageFloat32) + string(", but ") + value + string(" is ") + getDataTypeString(p);
							return 0;
						}

						if (doConversion)
							result = p2;

						return 1;
					}

					if (dt == ArgumentDataType::ImageComplex32)
					{
						Image<complex32_t >* p2 = dynamic_cast<Image<complex32_t >*>(p);
						if (!p2)
						{
							reason = string("Expected ") + toString(ArgumentDataType::ImageComplex32) + string(", but ") + value + string(" is ") + getDataTypeString(p);
							return 0;
						}

						if (doConversion)
							result = p2;

						return 1;
					}

					throw runtime_error("tryConvert not implemented for given input image data type.");
				}
				else
				{
					// Distributed case

					map<string, DistributedImageBase*>::const_iterator it = distributedImgs.find(value);
					if (it == distributedImgs.end())
					{
						reason = string("Image ") + value + string(" does not exist.");
						return 0;
					}

					DistributedImageBase* p = it->second;

					// Check image data type
					if (dt == ArgumentDataType::ImageUInt8)
					{
						DistributedImage<uint8_t>* p2 = dynamic_cast<DistributedImage<uint8_t>*>(p);
						if (!p2)
						{
							reason = string("Expected ") + toString(ArgumentDataType::ImageUInt8) + string(", but ") + value + string(" is ") + getDataTypeString(p);
							return 0;
						}

						if (doConversion)
							result = p2;

						return 1;
					}

					if (dt == ArgumentDataType::ImageUInt16)
					{
						DistributedImage<uint16_t>* p2 = dynamic_cast<DistributedImage<uint16_t>*>(p);
						if (!p2)
						{
							reason = string("Expected ") + toString(ArgumentDataType::ImageUInt16) + string(", but ") + value + string(" is ") + getDataTypeString(p);
							return 0;
						}

						if (doConversion)
							result = p2;

						return 1;
					}

					if (dt == ArgumentDataType::ImageUInt32)
					{
						DistributedImage<uint32_t>* p2 = dynamic_cast<DistributedImage<uint32_t>*>(p);
						if (!p2)
						{
							reason = string("Expected ") + toString(ArgumentDataType::ImageUInt32) + string(", but ") + value + string(" is ") + getDataTypeString(p);
							return 0;
						}

						if (doConversion)
							result = p2;

						return 1;
					}

					if (dt == ArgumentDataType::ImageUInt64)
					{
						DistributedImage<uint64_t>* p2 = dynamic_cast<DistributedImage<uint64_t>*>(p);
						if (!p2)
						{
							reason = string("Expected ") + toString(ArgumentDataType::ImageUInt64) + string(", but ") + value + string(" is ") + getDataTypeString(p);
							return 0;
						}

						if (doConversion)
							result = p2;

						return 1;
					}

					if (dt == ArgumentDataType::ImageFloat32)
					{
						DistributedImage<float32_t>* p2 = dynamic_cast<DistributedImage<float32_t>*>(p);
						if (!p2)
						{
							reason = string("Expected ") + toString(ArgumentDataType::ImageFloat32) + string(", but ") + value + string(" is ") + getDataTypeString(p);
							return 0;
						}

						if (doConversion)
							result = p2;

						return 1;
					}

					if (dt == ArgumentDataType::ImageComplex32)
					{
						DistributedImage<complex32_t >* p2 = dynamic_cast<DistributedImage<complex32_t >*>(p);
						if (!p2)
						{
							reason = string("Expected ") + toString(ArgumentDataType::ImageComplex32) + string(", but ") + value + string(" is ") + getDataTypeString(p);
							return 0;
						}

						if (doConversion)
							result = p2;

						return 1;
					}

					throw runtime_error("tryConvert not implemented for given input image data type.");
				}
			}
			else
			{
				throw runtime_error("Input parameter data type not implemented.");
			}
		}
		else
		{
			// Output parameter

			ArgumentDataType dt = type.dataType();

			if (dt == ArgumentDataType::ImageUInt8 ||
				dt == ArgumentDataType::ImageUInt16 ||
				dt == ArgumentDataType::ImageUInt32 ||
				dt == ArgumentDataType::ImageUInt64 ||
				dt == ArgumentDataType::ImageFloat32 ||
				dt == ArgumentDataType::ImageComplex32)
			{
				if (!isValidImageName(value, reason))
					return 0;

				// Any image type can be output automatically, but existing images get higher priority

				// Search match from existing images.
				if (!isDistributed())
				{
					map<string, ImageBase*>::const_iterator it = imgs.find(value);
					if (it != imgs.end())
					{
						ImageBase* p = it->second;
						if ((dt == ArgumentDataType::ImageUInt8 && dynamic_cast<Image<uint8_t>*>(p)) ||
							(dt == ArgumentDataType::ImageUInt16 && dynamic_cast<Image<uint16_t>*>(p)) ||
							(dt == ArgumentDataType::ImageUInt32 && dynamic_cast<Image<uint32_t>*>(p)) ||
							(dt == ArgumentDataType::ImageUInt64 && dynamic_cast<Image<uint64_t>*>(p)) ||
							(dt == ArgumentDataType::ImageFloat32 && dynamic_cast<Image<float32_t>*>(p)) ||
							(dt == ArgumentDataType::ImageComplex32 && dynamic_cast<Image<complex32_t>*>(p)))
						{
							// Match
							if (doConversion)
							{
								switch (dt)
								{
								case ArgumentDataType::ImageUInt8: result = dynamic_cast<Image<uint8_t>*>(p); break;
								case ArgumentDataType::ImageUInt16: result = dynamic_cast<Image<uint16_t>*>(p); break;
								case ArgumentDataType::ImageUInt32: result = dynamic_cast<Image<uint32_t>*>(p); break;
								case ArgumentDataType::ImageUInt64: result = dynamic_cast<Image<uint64_t>*>(p); break;
								case ArgumentDataType::ImageFloat32: result = dynamic_cast<Image<float32_t>*>(p); break;
								case ArgumentDataType::ImageComplex32: result = dynamic_cast<Image<complex32_t>*>(p); break;
								default: throw logic_error("Image data type not configured.");
								}
							}
							return 1;
						}
					}
				}
				else
				{
					map<string, DistributedImageBase*>::const_iterator it = distributedImgs.find(value);
					if (it != distributedImgs.end())
					{
						DistributedImageBase* p = it->second;
						if ((dt == ArgumentDataType::ImageUInt8 && dynamic_cast<DistributedImage<uint8_t>*>(p)) ||
							(dt == ArgumentDataType::ImageUInt16 && dynamic_cast<DistributedImage<uint16_t>*>(p)) ||
							(dt == ArgumentDataType::ImageUInt32 && dynamic_cast<DistributedImage<uint32_t>*>(p)) ||
							(dt == ArgumentDataType::ImageUInt64 && dynamic_cast<DistributedImage<uint64_t>*>(p)) ||
							(dt == ArgumentDataType::ImageFloat32 && dynamic_cast<DistributedImage<float32_t>*>(p)) ||
							(dt == ArgumentDataType::ImageComplex32 && dynamic_cast<DistributedImage<complex32_t>*>(p)))
						{
							// Match
							if (doConversion)
							{
								switch (dt)
								{
								case ArgumentDataType::ImageUInt8: result = dynamic_cast<DistributedImage<uint8_t>*>(p); break;
								case ArgumentDataType::ImageUInt16: result = dynamic_cast<DistributedImage<uint16_t>*>(p); break;
								case ArgumentDataType::ImageUInt32: result = dynamic_cast<DistributedImage<uint32_t>*>(p); break;
								case ArgumentDataType::ImageUInt64: result = dynamic_cast<DistributedImage<uint64_t>*>(p); break;
								case ArgumentDataType::ImageFloat32: result = dynamic_cast<DistributedImage<float32_t>*>(p); break;
								case ArgumentDataType::ImageComplex32: result = dynamic_cast<DistributedImage<complex32_t>*>(p); break;
								default: throw logic_error("Image data type not configured.");
								}
							}
							return 1;
						}
					}
				}

				// Create image if it does not exist.
				if (doConversion)
				{
					if (!isDistributed())
					{
						// Create normal image
						if (imgs[value] != 0)
							delete imgs[value];

						if (dt == ArgumentDataType::ImageUInt8)
						{
							result = new Image<uint8_t>();
							imgs[value] = get<Image<uint8_t>*>(result);
						}
						else if (dt == ArgumentDataType::ImageUInt16)
						{
							result = new Image<uint16_t>();
							imgs[value] = get<Image<uint16_t>*>(result);
						}
						else if (dt == ArgumentDataType::ImageUInt32)
						{
							result = new Image<uint32_t>();
							imgs[value] = get<Image<uint32_t>*>(result);
						}
						else if (dt == ArgumentDataType::ImageUInt64)
						{
							result = new Image<uint64_t>();
							imgs[value] = get<Image<uint64_t>*>(result);
						}
						else if (dt == ArgumentDataType::ImageFloat32)
						{
							result = new Image<float32_t>();
							imgs[value] = get<Image<float32_t>*>(result);
						}
						else if (dt == ArgumentDataType::ImageComplex32)
						{
							result = new Image<complex32_t>();
							imgs[value] = get<Image<complex32_t>*>(result);
						}
						else
						{
							throw runtime_error("Output image data type not implemented.");
						}
					}
					else
					{
						// Create distributed image
						if (distributedImgs[value] != 0)
							delete distributedImgs[value];


						if (dt == ArgumentDataType::ImageUInt8)
						{
							result = new DistributedImage<uint8_t>(value);
							distributedImgs[value] = get<DistributedImage<uint8_t>*>(result);
						}
						else if (dt == ArgumentDataType::ImageUInt16)
						{
							result = new DistributedImage<uint16_t>(value);
							distributedImgs[value] = get<DistributedImage<uint16_t>*>(result);
						}
						else if (dt == ArgumentDataType::ImageUInt32)
						{
							result = new DistributedImage<uint32_t>(value);
							distributedImgs[value] = get<DistributedImage<uint32_t>*>(result);
						}
						else if (dt == ArgumentDataType::ImageUInt64)
						{
							result = new DistributedImage<uint64_t>(value);
							distributedImgs[value] = get<DistributedImage<uint64_t>*>(result);
						}
						else if (dt == ArgumentDataType::ImageFloat32)
						{
							result = new DistributedImage<float32_t>(value);
							distributedImgs[value] = get<DistributedImage<float32_t>*>(result);
						}
						else if (dt == ArgumentDataType::ImageComplex32)
						{
							result = new DistributedImage<complex32_t>(value);
							distributedImgs[value] = get<DistributedImage<complex32_t>*>(result);
						}
						else
						{
							throw runtime_error("Output image data type not implemented.");
						}
					}
				}

				return 2;
			}
			else
			{
				// Any other output data type is not supported...
				throw runtime_error("Output parameter data type not implemented.");
			}
		}


	}

	/**
	Checks whether supplied string values can be converted to supplied types.
	Values are not changed but the reference is non-const as required by tryConvert.
	If 0 is returned, reason parameter is assigned an explanation why this match does not succeed.
	@return 0 if there is no match, 1 if there is match, and 2 if there is match after creation of new images.
	*/
	int PISystem::matchParameterTypes(const vector<CommandArgumentBase>& types, vector<string>& values, string& reason)
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
	void PISystem::executeCommand(const string& name, vector<string>& args)
	{
		// Match name
		vector<Command*> candidates = getCommands(name);

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

		// Add defaults to the parameter array
		vector<string> realArgs;
		realArgs.reserve(cmd->args().size());
		size_t N0 = args.size();
		for (size_t n = 0; n < N0; n++)
			realArgs.push_back(args[n]);
		for (size_t n = N0; n < cmd->args().size(); n++)
			realArgs.push_back(cmd->args()[n].defaultValue());

		// Commands that should not be echoed to screen
		bool isNoShow = cmd->name() == "help" || cmd->name() == "info" || cmd->name() == "license";

		// Show command
		if (showRunCommands && !isNoShow)
		{
			cout << cmd->name() << "(";
			for (size_t n = 0; n < realArgs.size(); n++)
			{
				cout << "\"" << realArgs[n] << "\"";
				if (n < realArgs.size() - 1)
					cout << ", ";
			}
			cout << ")" << endl;
		}

		// Convert string parameters to values. This must succeed as the process was tested above.
		vector<ParamVariant> convertedArgs;
		convertedArgs.reserve(realArgs.size());
		for (size_t n = 0; n < realArgs.size(); n++)
		{
			ParamVariant res;
			string dummy;
			tryConvert(realArgs[n], cmd->args()[n], true, res, dummy);
			convertedArgs.push_back(res);
		}

		// Run command with timing
		Timer timer;
		timer.start();
		if (!isDistributed())
		{
			// Normal processing without distribution or anything fancy
			cmd->runInternal(this, convertedArgs);
		}
		else
		{
			// Distributed processing
			Distributable* dist = dynamic_cast<Distributable*>(cmd);
			if (dist)
			{
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
	}

	/**
	Parse one line (one command) of input code.
	*/
	void PISystem::parseLine(string& line)
	{
		trim(line);

		//cout << "Parsing line: " << line << endl;

		if (line.length() > 0)
		{
			if (line[0] != '%' && line[0] != '#' && line.find("//") != 0)
			{
				// Not a comment.
				string cmd;
				vector<string> args;
				parseFunctionCall(line, cmd, args);

				//cout << "Name: " << cmd << endl;
				//for (size_t n = 0; n < args.size(); n++)
				//	cout << "Argument " << (n + 1) << ": " << args[n] << endl;

				executeCommand(cmd, args);
			}

		}
	}


	/**
	Gets a string representing data type of given image.
	*/
	string PISystem::getDataTypeString(const ImageBase* p)
	{
		const Image<uint8_t>* p8 = dynamic_cast<const Image<uint8_t>*>(p);
		const Image<uint16_t>* p16 = dynamic_cast<const Image<uint16_t>*>(p);
		const Image<uint32_t>* p32 = dynamic_cast<const Image<uint32_t>*>(p);
		const Image<uint64_t>* p64 = dynamic_cast<const Image<uint64_t>*>(p);
		const Image<float32_t>* pf32 = dynamic_cast<const Image<float32_t>*>(p);
		const Image<complex32_t >* pc32 = dynamic_cast<const Image<complex32_t >*>(p);

		if (p8)
			return toString(ArgumentDataType::ImageUInt8);
		if (p16)
			return toString(ArgumentDataType::ImageUInt16);
		if (p32)
			return toString(ArgumentDataType::ImageUInt32);
		if (p64)
			return toString(ArgumentDataType::ImageUInt64);
		if (pf32)
			return toString(ArgumentDataType::ImageFloat32);
		if (pc32)
			return toString(ArgumentDataType::ImageComplex32);

		throw runtime_error("getDataTypeString not implemented for given data type.");
	}


	/**
	Gets a string representing data type of given image.
	*/
	string PISystem::getDataTypeString(const DistributedImageBase* p)
	{
		const DistributedImage<uint8_t>* p8 = dynamic_cast<const DistributedImage<uint8_t>*>(p);
		const DistributedImage<uint16_t>* p16 = dynamic_cast<const DistributedImage<uint16_t>*>(p);
		const DistributedImage<uint32_t>* p32 = dynamic_cast<const DistributedImage<uint32_t>*>(p);
		const DistributedImage<uint64_t>* p64 = dynamic_cast<const DistributedImage<uint64_t>*>(p);
		const DistributedImage<float32_t>* pf32 = dynamic_cast<const DistributedImage<float32_t>*>(p);
		const DistributedImage<complex32_t >* pc32 = dynamic_cast<const DistributedImage<complex32_t >*>(p);

		if (p8)
			return toString(ArgumentDataType::ImageUInt8);
		if (p16)
			return toString(ArgumentDataType::ImageUInt16);
		if (p32)
			return toString(ArgumentDataType::ImageUInt32);
		if (p64)
			return toString(ArgumentDataType::ImageUInt64);
		if (pf32)
			return toString(ArgumentDataType::ImageFloat32);
		if (pc32)
			return toString(ArgumentDataType::ImageComplex32);

		throw runtime_error("getDataTypeString not implemented for given data type.");
	}






	PISystem::PISystem()
	{
		clearLastError();

		addFilterCommands(commands);
		addSpecialCommands(commands);
		addIOCommands(commands);
		addOtherCommands(commands);
		addPointProcessCommands(commands);
		addThinAndSkeletonCommands(commands);
		addSpecialSystemCommands(commands);
		addTransformCommands(commands);
		addParticlesCommands(commands);
		addGenerationCommands(commands);
		addProjectionCommands(commands);
	}

	PISystem::~PISystem()
	{
		// Delete all images
		for (map<string, ImageBase*>::iterator it = imgs.begin(); it != imgs.end(); it++)
			delete it->second;
		imgs.clear();

		for (map<string, DistributedImageBase*>::iterator it = distributedImgs.begin(); it != distributedImgs.end(); it++)
			delete it->second;
		distributedImgs.clear();

		// Delete all commands
		for (size_t n = 0; n < commands.size(); n++)
		{
			delete commands[n];
		}
	}

	/**
	Returns names of images in the system.
	*/
	vector<string> PISystem::getImageNames() const
	{
		std::vector<string> keys;
		keys.reserve(imgs.size());
		for (auto const& item : imgs)
			keys.push_back(item.first);
		return keys;
	}

	/**
	Returns name of the given image.
	*/
	string PISystem::getImageName(const ImageBase* img) const
	{
		for (auto const& item : imgs)
			if (item.second == img)
				return item.first;
		throw ITLException("Invalid image name query.");
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
			replaceImage(name, distributedImgs[name]->toNormalImage());
		}

		return imgs[name];
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
			distributedImgs[imgName]->setData(imgs[imgName]);
			replaceImage(imgName, 0);
		}
	}

	/**
	Same than flushIfDistributed but sets lastException instead of throwing it.
	@return false if an error occurs.
	*/
	bool PISystem::flushIfDistributedNoThrow(const string& name)
	{
		try
		{
			flushIfDistributed(name);
			return true;
		}
		catch (ITLException& e)
		{
			lastException = e.message();
		}
		catch (exception& e)
		{
			lastException = e.what();
		}
		return false;
	}

	/**
	Retrieve image having given name, do not throw exception but set last error and return 0 if exception would be thrown.
	*/
	ImageBase* PISystem::getImageNoThrow(const string& name)
	{
		try
		{
			return getImage(name);
		}
		catch (ITLException& e)
		{
			lastException = e.message();
		}
		catch (exception& e)
		{
			lastException = e.what();
		}
		if (lastException == "")
			lastException = "Image not found.";
		return 0;
	}

	/**
	Replace image with given image.
	Replace image by null pointer to remove it from the system.
	*/
	void PISystem::replaceImage(const string& name, ImageBase* newImg)
	{
		if (imgs.find(name) != imgs.end())
		{
			if (imgs[name] != newImg)
				delete imgs[name];
			imgs.erase(name);
		}
		if (newImg)
			imgs[name] = newImg;
	}

	/**
	Get distributed image having given name.
	*/
	DistributedImageBase* PISystem::getDistributedImage(const string& name)
	{
		return distributedImgs[name];
	}

	/**
	Replace distributed image with a new one.
	Set newImg to zero to remove the image with given name.
	*/
	void PISystem::replaceDistributedImage(const string& name, DistributedImageBase* newImg, bool free)
	{
		if (distributedImgs.find(name) != distributedImgs.end())
		{
			if (distributedImgs[name] != newImg && free)
				delete distributedImgs[name];
			distributedImgs.erase(name);
		}
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
			delete distributor;
			distributor = 0;
			cout << "Disabled distributed computing mode." << endl;
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
			else
				throw ITLException(string("Invalid distributed computing system name: ") + provider + ". Valid names are SLURM or Local.");
		}
	}

	/**
	Parses commands in the given string and runs them.
	*/
	bool PISystem::run(const string& commands)
	{
		try
		{
			string rest = commands;
			lastExceptionLine = 1;
			while (rest.length() > 0)
			{
				char delim = 0;
				string token = getToken(rest, ";\n", delim);

				parseLine(token);

				if (delim == '\n')
					lastExceptionLine++;
			}

			lastExceptionLine = 0;
			return true;
		}
		catch (ITLException& e)
		{
			lastException = e.message();
			return false;
		}
		catch (exception& e)
		{
			lastException = e.what();
			return false;
		}
	}

	/**
	Gets help topics for given command name.
	*/
	vector<string> PISystem::getHelp(const string& cmdName) const
	{
		vector<string> results;

		vector<Command*> cmds = getCommands(cmdName);

		for (size_t n = 0; n < cmds.size(); n++)
		{
			results.push_back(cmds[n]->helpString());
		}

		return results;
	}

	/**
	Get a list of commands and some basic information about each one, each command separated by newline.
	*/
	string PISystem::commandList(bool forUser) const
	{
		vector<tuple<string, bool> > names;
		for (size_t n = 0; n < commands.size(); n++)
		{
			names.push_back(make_tuple(commands[n]->toSimpleString(), commands[n]->canDistribute()));
		}
		removeDuplicates(names);
		sort(names.begin(), names.end());
		stringstream msg;
		if (forUser)
			msg << "Letter D placed before command name indicates a command that can be used in distributed computing mode. See help(distribute) for more information." << endl << endl;
		for (size_t n = 0; n < names.size(); n++)
		{
			if (forUser)
				msg << (get<1>(names[n]) ? "D " : "  ");
			msg << get<0>(names[n]) << endl;
		}
		return msg.str();
	}

	/**
	Convert argument to string.
	@param isDistributed Set to true to convert arguments of type Image* to DistributedImage*
	*/
	string PISystem::argumentToString(const CommandArgumentBase& argument, const ParamVariant& value, bool isDistributed) const
	{
		ArgumentDataType dt = argument.dataType();

		if (isDistributed)
		{
			switch (dt)
			{
			case ArgumentDataType::ImageUInt8: dt = ArgumentDataType::DImageUInt8; break;
			case ArgumentDataType::ImageUInt16: dt = ArgumentDataType::DImageUInt16; break;
			case ArgumentDataType::ImageUInt32: dt = ArgumentDataType::DImageUInt32; break;
			case ArgumentDataType::ImageUInt64: dt = ArgumentDataType::DImageUInt64; break;
			case ArgumentDataType::ImageFloat32: dt = ArgumentDataType::DImageFloat32; break;
			case ArgumentDataType::ImageComplex32: dt = ArgumentDataType::DImageComplex32; break;
			}
		}

		string s;
		switch (dt)
		{
		case ArgumentDataType::String: s = get<string>(value); return s;
		case ArgumentDataType::Double: return itl2::toString(get<double>(value));
		case ArgumentDataType::Int: return itl2::toString(get<coord_t>(value));
		case ArgumentDataType::Size: return itl2::toString(get<size_t>(value));
		case ArgumentDataType::NBType: return itl2::toString(get<NeighbourhoodType>(value));
		case ArgumentDataType::BoundaryCond: return itl2::toString(get<BoundaryCondition>(value));
		case ArgumentDataType::Connectiv: return itl2::toString(get<Connectivity>(value));
		case ArgumentDataType::InterpolationMode: return itl2::toString(get<InterpolationMode>(value));
		case ArgumentDataType::Bool: return itl2::toString(get<bool>(value));
		case ArgumentDataType::Vect3d: return itl2::toString(get<Vec3d>(value));
		case ArgumentDataType::Vect3c: return itl2::toString(get<Vec3c>(value));
		case ArgumentDataType::ImageUInt8: return imageName(get<Image<uint8_t>* >(value));
		case ArgumentDataType::ImageUInt16: return imageName(get<Image<uint16_t>* >(value));
		case ArgumentDataType::ImageUInt32: return imageName(get<Image<uint32_t>* >(value));
		case ArgumentDataType::ImageUInt64: return imageName(get<Image<uint64_t>* >(value));
		case ArgumentDataType::ImageFloat32: return imageName(get<Image<float32_t>* >(value));
		case ArgumentDataType::ImageComplex32: return imageName(get<Image<complex32_t>* >(value));
		case ArgumentDataType::DImageUInt8: return distributedImageName(get<DistributedImage<uint8_t>* >(value));
		case ArgumentDataType::DImageUInt16: return distributedImageName(get<DistributedImage<uint16_t>* >(value));
		case ArgumentDataType::DImageUInt32: return distributedImageName(get<DistributedImage<uint32_t>* >(value));
		case ArgumentDataType::DImageFloat32: return distributedImageName(get<DistributedImage<float32_t>* >(value));
		case ArgumentDataType::DImageComplex32: return distributedImageName(get<DistributedImage<complex32_t>* >(value));
		default: throw ITLException("Data type not configured.");
		}
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
		math::Vec3c templDims;
		if (!itl2::io::getInfo(filename, templDims, templDT))
			throw ITLException(string("Unable to find dimensions and data type of file ") + filename);

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
		math::Vec3c templDims;
		if (!itl2::io::getInfo(filename, templDims, templDT))
			throw ITLException(string("Unable to find dimensions and data type of file ") + filename);

		PISystem* system = distributor.getSystem();

		createDistributedImage(name, dts, w, h, d, templDT, templDims.x, templDims.y, templDims.z, system);

		return vector<string>();
	}

}
