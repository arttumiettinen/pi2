
#include "pisystem.h"
#include "commandmacros.h"
#include "io/io.h"
#include "commandlist.h"
#include "pilibutilities.h"

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
						args.push_back(arg);
					}
					else if (argSection[0] == '"')
					{
						argSection.erase(argSection.begin());
						string arg = parseString(argSection, '"');
						args.push_back(arg);
						trim(argSection);
						if (argSection.length() > 0 && argSection[0] != ',')
							throw ParseException("Trailing characters after string value.");
						if (argSection.length() > 0)
    						argSection.erase(argSection.begin());
					}
					else if (argSection[0] == '\'')
					{
						argSection.erase(argSection.begin());
						string arg = parseString(argSection, '\'');
						args.push_back(arg);
						trim(argSection);
						if (argSection.length() > 0 && argSection[0] != ',')
							throw ParseException("Trailing characters after string value.");
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
						args.push_back(arg);
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
		for (map<string, shared_ptr<ImageBase> >::const_iterator it = imgs.begin(); it != imgs.end(); it++)
		{
			if (it->second.get() == img)
				return it->first;
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
		if (imgs.find(name) != imgs.end() || distributedImgs.find(name) != distributedImgs.end())
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
	template<typename target_t> bool trySimpleConversion(ArgumentDataType dt, const string& value, bool doConversion, ParamVariant& result, string& reason)
	{
		if (parameterType<target_t>() != dt)
			return false;

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
	Functor that casts image to another type and assigns it to ParamVariant.
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
	int PISystem::tryConvert(string& value, const CommandArgumentBase& type, bool doConversion, ParamVariant& result, string& reason)
	{
		ArgumentDataType dt = type.dataType();

		if (type.direction() == ParameterDirection::In || type.direction() == ParameterDirection::InOut)
		{
			// Input parameter.
			// Strings must be convertible to primitives, images must exist.

			if (trySimpleConversion<string>(dt, value, doConversion, result, reason))
				return 1;
			if (trySimpleConversion<bool>(dt, value, doConversion, result, reason))
				return 1;
			if (trySimpleConversion<double>(dt, value, doConversion, result, reason))
				return 1;
			if (trySimpleConversion<Vec3d>(dt, value, doConversion, result, reason))
				return 1;
			if (dt == ArgumentDataType::Vect3d)
			{
				Vec3d v;
				if (getImageAsVec3(value, v))
				{
					if (doConversion)
						result = v;
					return 1;
				}
			}
			if (trySimpleConversion<Vec3c>(dt, value, doConversion, result, reason))
				return 1;
			if (dt == ArgumentDataType::Vect3c)
			{
				Vec3d v;
				if (getImageAsVec3(value, v))
				{
					if (doConversion)
						result = round(v);
					return 1;
				}
			}
			if (trySimpleConversion<coord_t>(dt, value, doConversion, result, reason))
				return 1;
			if (trySimpleConversion<size_t>(dt, value, doConversion, result, reason))
				return 1;
			if (trySimpleConversion<NeighbourhoodType>(dt, value, doConversion, result, reason))
				return 1;
			if (trySimpleConversion<BoundaryCondition>(dt, value, doConversion, result, reason))
				return 1;
			if (trySimpleConversion<Connectivity>(dt, value, doConversion, result, reason))
				return 1;
			if (trySimpleConversion<InterpolationMode>(dt, value, doConversion, result, reason))
				return 1;

			if (isImage(dt))
			{
				// Find image with given name

				if (!isValidImageName(value, reason))
					return 0;

				ImageDataType idt = argumentDataTypeToImageDataType(dt);

				if (!isDistributed())
				{
					// Non-distributed case

					auto it = imgs.find(value);
					if (it == imgs.end())
					{
						reason = string("Image ") + value + string(" does not exist.");
						return 0;
					}

					// Check that image data type is correct
					ImageBase* p = it->second.get();

					if (idt != p->dataType())
					{
						reason = string("Expected ") + toString(idt) + string(" image, but ") + value + string(" is ") + toString(p->dataType()) + string(" image.");
						return 0;
					}
					
					if (doConversion)
						pick<CastImage>(idt, p, result);

					return 1;
				}
				else
				{
					auto it = distributedImgs.find(value);
					if (it == distributedImgs.end())
					{
						reason = string("Image ") + value + string(" does not exist.");
						return 0;
					}

					// Check that image data type is correct
					DistributedImageBase* p = it->second.get();

					if (idt != p->dataType())
					{
						reason = string("Expected ") + toString(idt) + string(" image, but ") + value + string(" is ") + toString(p->dataType()) + string(" image.");
						return 0;
					}

					if (doConversion)
						pick<CastDistributedImage>(idt, p, result);

					return 1;
				}


			}

			return 0;
		}
		else
		{
			// Output parameter


			if (isImage(dt))
			{
				if (!isValidImageName(value, reason))
					return 0;

				ImageDataType idt = argumentDataTypeToImageDataType(dt);

				// Any image type can be output automatically, but existing images get higher priority

				// Search match from existing images.
				Vec3c oldSize(1, 1, 1);
				if (!isDistributed())
				{
					auto it = imgs.find(value);
					if (it != imgs.end())
					{
						ImageBase* p = it->second.get();
						oldSize = p->dimensions();

						if (p->dataType() == idt)
						{
							// Match
							if (doConversion)
								pick<CastImage>(idt, p, result);
							return 1;
						}
					}
				}
				else
				{
					auto it = distributedImgs.find(value);
					if (it != distributedImgs.end())
					{
						DistributedImageBase* p = it->second.get();
						oldSize = p->dimensions();

						if (p->dataType() == idt)
						{
							// Match
							if (doConversion)
								pick<CastDistributedImage>(idt, p, result);
							return 1;
						}
					}
				}

				// No match found from existing images.
				// Create image if it does not exist.
				if (doConversion)
				{
					if (!isDistributed())
					{
						// Create normal image
						// Retain size of old image of the same name
						pick<CreateImage>(idt, oldSize, value, this);
						ImageBase* p = imgs[value].get();
						pick<CastImage>(idt, p, result);
					}
					else
					{
						// Create distributed image
						// Retain size of old image of the same name
						pick<CreateEmptyDistributedImage>(idt, value, oldSize, this);
						DistributedImageBase* p = distributedImgs[value].get();
						pick<CastDistributedImage>(idt, p, result);
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

		vector<shared_ptr<ImageBase> > imageStore;
		
		for (size_t n = 0; n < realArgs.size(); n++)
		{
			ParamVariant res;
			string dummy;
			tryConvert(realArgs[n], cmd->args()[n], true, res, dummy);

			// Store shared_ptrs to image arguments so that they don't get deleted if the next parameters override them.
			// NOTE: We just hold the pointers but don't do anything with them.
			DistributedImageBase* db = getDistributedImageNoThrow(res);
			if (db)
			{
				distributedImageStore.push_back(getDistributedImagePointer(db));
			}

			ImageBase* ib = pilib::getImageNoThrow(res);
			if (ib)
			{
				imageStore.push_back(getImagePointer(ib));
			}


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

		distributedImageStore.clear();
	}

	/**
	Parse one statement of input code.
	*/
	void PISystem::parseStatement(string& statement)
	{
		trim(statement);

		string cmd;
		vector<string> args;
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
			if (item.second.get() == img)
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

		return imgs[name].get();
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
			distributedImgs[imgName]->setData(imgs[imgName].get());
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
	void PISystem::getImageInfoNoThrow(const string& name, coord_t& width, coord_t& height, coord_t& depth, ImageDataType& dt)
	{
		try
		{
			// NOTE: Do not use getImage as that results in reading distributed images to RAM.
			if (!isDistributed())
			{
				auto p = imgs[name];
				width = p->width();
				height = p->height();
				depth = p->depth();
				dt = p->dataType();
			}
			else
			{
				auto p = distributedImgs[name];
				width = p->width();
				height = p->height();
				depth = p->depth();
				dt = p->dataType();
			}
			return;
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
	void PISystem::replaceImage(const string& name, shared_ptr<ImageBase> newImg)
	{
		if (imgs.find(name) != imgs.end())
			imgs.erase(name);
		
		if (newImg)
			imgs[name] = newImg;
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
				string token = getToken(rest, "\n", delim);

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

	shared_ptr<ImageBase> PISystem::getImagePointer(ImageBase* img) const
	{
		for (auto& item : imgs)
		{
			if (item.second.get() == img)
				return item.second;
		}

		throw ITLException("Unknown image pointer.");
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
