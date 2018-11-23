#pragma once

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <sstream>
#include <iomanip>
#include <algorithm>

#include "itl2.h"

#include "command.h"
#include "parseexception.h"
#include "timer.h"
#include "distributor.h"
#include "distributable.h"
#include "distributedimage.h"
#include "slurmdistributor.h"
#include "localdistributor.h"
#include "pointprocesscommands.h"

using namespace std;
using namespace itl2;

namespace pilib
{
	/**
	Creates all commands and adds them to the list.
	Used in the constructor of PISystem and defined elsewhere to keep this file shorter and cleaner.
	*/
	void addFilterCommands(vector<Command*>& commands);
	void addSpecialCommands(vector<Command*>& commands);
	void addIOCommands(vector<Command*>& commands);
	void addOtherCommands(vector<Command*>& commands);
	void addPointProcessCommands(vector<Command*>& commands);
	void addThinAndSkeletonCommands(vector<Command*>& commands);
	void addSpecialSystemCommands(vector<Command*>& commands);
	void addTransformCommands(vector<Command*>& commands);

	/**
	Allows usage of itl2 functionality by simple commands.
	*/
	class PISystem
	{
	private:
		/**
		Maps image variable names to image objects.
		*/
		map<string, ImageBase*> imgs;

		/**
		Maps image variable names to distributed image objects.
		*/
		map<string, DistributedImageBase*> distributedImgs;

		/**
		List of available commands.
		*/
		vector<Command*> commands;

		/**
		Last exception that has occured.
		*/
		string lastException;
		int lastExceptionLine;

		/**
		Set to true to show all executed commands with parameters.
		*/
		bool showRunCommands = false;
		bool showTiming = false;

		/**
		Set to Distributor object to enable distributed processing of commands.
		*/
		Distributor* distributor = 0;

		/**
		Get next token from a string and remove it from the string.
		*/
		static string getToken(string& str, const char* delimiters, char& foundDelimiter)
		{
			size_t pos = str.find_first_of(delimiters);
			if (pos == string::npos)
			{
				// No delimiter found
				foundDelimiter = 0;
				string res = str;
				str = "";
				return res;
			}

			foundDelimiter = str[pos];
			string res = str.substr(0, pos);
			str.erase(0, pos + 1);
			return res;
		}

		/**
		Parse line expected to contain function call
		funcname(param1, param2, param3, ...)
		*/
		static void parseFunctionCall(const string& line, string& name, vector<string>& args)
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
		vector<Command*> getCommands(const string& name) const
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
		static bool isValidImageName(const string& value, string& reason)
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
		string imageName(const ImageBase* img) const
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
		string distributedImageName(const DistributedImageBase* img) const
		{
			return img->varName();
		}

		/**
		Try to convert the given string to match the given argument.
		The conversion is done and, e.g., images are created if doConversion flag is set to true.
		The result of the conversion is assigned to the result variable.
		If conversion is not possible, reason string is assigned an explanation of the error.
		@return 0 if there is no match, 1 if the argument type and parameter type match, and 2 if they match after creation of new images.
		*/
		int tryConvert(string& value, const CommandArgumentBase& type, bool doConversion, ParamVariant& result, string& reason)
		{
			if (doConversion)
				result.userPassedString = value;

			if (type.direction() == In || type.direction() == InOut)
			{
				// Input parameter.
				// Strings must be convertible to primitives, images must exist.

				ArgumentDataType dt = type.dataType();

				if (dt == String)
				{
					if (doConversion)
						result.sval = &value;
					return 1;
				}
				else if (dt == Bool)
				{
					try
					{
						bool val = itl2::fromString<bool>(value);
						if (doConversion)
							result.bval = val;
						return 1;
					}
					catch (ITLException& e)
					{
						reason = e.message();
						return 0;
					}
				}
				else if (dt == Double)
				{
					try
					{
						double val = itl2::fromString<double>(value);
						if (doConversion)
							result.dval = val;
						return 1;
					}
					catch (ITLException& e)
					{
						reason = e.message();
						return 0;
					}
				}
				else if (dt == Vect3d)
				{
					try
					{
						Vec3d val = itl2::fromString<Vec3d>(value);
						if (doConversion)
						{
							result.vx = val.x;
							result.vy = val.y;
							result.vz = val.z;
						}
						return 1;
					}
					catch (ITLException& e)
					{
						reason = e.message();
						return 0;
					}
				}
				else if (dt == Vect3c)
				{
					try
					{
						Vec3c val = itl2::fromString<Vec3c>(value);
						if (doConversion)
						{
							result.vix = val.x;
							result.viy = val.y;
							result.viz = val.z;
						}
						return 1;
					}
					catch (ITLException& e)
					{
						reason = e.message();
						return 0;
					}
				}
				else if (dt == Int || dt == Size)
				{
					char* p;
					int val = strtol(value.c_str(), &p, 10);
					if (*p)
					{
						reason = string("Value ") + value + string(" is not a valid integer number.");
						return 0;
					}
					if (doConversion)
					{
						if (dt == Int)
							result.ival = val;
						else
							result.ival = math::max(0, val);
					}
					return 1;
				}
				else if (dt == NBType)
				{
					if (value == "rect" || value == "rectangular" || value == "Rectangular")
					{
						if (doConversion)
							result.nbtval = Rectangular;
						return 1;
					}

					if (value == "ell" || value == "ellipsoidal" || value == "Ellipsoidal")
					{
						if (doConversion)
							result.nbtval = Ellipsoidal;
						return 1;
					}

					reason = string("Unknown neighbourhood type: ") + value;
					return 0;
				}
				else if (dt == BoundaryCond)
				{
					if (value == "zero" || value == "Zero" || value == "ZERO" || value == "0")
					{
						if (doConversion)
							result.bcval = Zero;
						return 1;
					}

					if (value == "nearest" || value == "Nearest" || value == "NEAREST" || value == "1")
					{
						if (doConversion)
							result.bcval = Nearest;
						return 1;
					}

					reason = string("Unknown boundary condition: ") + value;
					return 0;
				}
				else if (dt == ImageUInt8 ||
					dt == ImageUInt16 ||
					dt == ImageUInt32 ||
					dt == ImageUInt64 ||
					dt == ImageFloat32 ||
					dt == ImageComplex32)
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
						if (dt == ImageUInt8)
						{
							Image<uint8_t>* p2 = dynamic_cast<Image<uint8_t>*>(p);
							if (!p2)
							{
								reason = string("Expected ") + toString(ImageUInt8) + string(", but ") + value + string(" is ") + getDataTypeString(p);
								return 0;
							}

							if (doConversion)
								result.img8 = p2;

							return 1;
						}

						if (dt == ImageUInt16)
						{
							Image<uint16_t>* p2 = dynamic_cast<Image<uint16_t>*>(p);
							if (!p2)
							{
								reason = string("Expected ") + toString(ImageUInt16) + string(", but ") + value + string(" is ") + getDataTypeString(p);
								return 0;
							}

							if (doConversion)
								result.img16 = p2;

							return 1;
						}

						if (dt == ImageUInt32)
						{
							Image<uint32_t>* p2 = dynamic_cast<Image<uint32_t>*>(p);
							if (!p2)
							{
								reason = string("Expected ") + toString(ImageUInt32) + string(", but ") + value + string(" is ") + getDataTypeString(p);
								return 0;
							}

							if (doConversion)
								result.img32 = p2;

							return 1;
						}

						if (dt == ImageUInt64)
						{
							Image<uint64_t>* p2 = dynamic_cast<Image<uint64_t>*>(p);
							if (!p2)
							{
								reason = string("Expected ") + toString(ImageUInt64) + string(", but ") + value + string(" is ") + getDataTypeString(p);
								return 0;
							}

							if (doConversion)
								result.img64 = p2;

							return 1;
						}

						if (dt == ImageFloat32)
						{
							Image<float32_t>* p2 = dynamic_cast<Image<float32_t>*>(p);
							if (!p2)
							{
								reason = string("Expected ") + toString(ImageFloat32) + string(", but ") + value + string(" is ") + getDataTypeString(p);
								return 0;
							}

							if (doConversion)
								result.imgf32 = p2;

							return 1;
						}

						if (dt == ImageComplex32)
						{
							Image<complex32_t >* p2 = dynamic_cast<Image<complex32_t >*>(p);
							if (!p2)
							{
								reason = string("Expected ") + toString(ImageComplex32) + string(", but ") + value + string(" is ") + getDataTypeString(p);
								return 0;
							}

							if (doConversion)
								result.imgc32 = p2;

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
						if (dt == ImageUInt8)
						{
							DistributedImage<uint8_t>* p2 = dynamic_cast<DistributedImage<uint8_t>*>(p);
							if (!p2)
							{
								reason = string("Expected ") + toString(ImageUInt8) + string(", but ") + value + string(" is ") + getDataTypeString(p);
								return 0;
							}

							if (doConversion)
								result.dimg8 = p2;

							return 1;
						}

						if (dt == ImageUInt16)
						{
							DistributedImage<uint16_t>* p2 = dynamic_cast<DistributedImage<uint16_t>*>(p);
							if (!p2)
							{
								reason = string("Expected ") + toString(ImageUInt16) + string(", but ") + value + string(" is ") + getDataTypeString(p);
								return 0;
							}

							if (doConversion)
								result.dimg16 = p2;

							return 1;
						}

						if (dt == ImageUInt32)
						{
							DistributedImage<uint32_t>* p2 = dynamic_cast<DistributedImage<uint32_t>*>(p);
							if (!p2)
							{
								reason = string("Expected ") + toString(ImageUInt32) + string(", but ") + value + string(" is ") + getDataTypeString(p);
								return 0;
							}

							if (doConversion)
								result.dimg32 = p2;

							return 1;
						}

						if (dt == ImageUInt64)
						{
							DistributedImage<uint64_t>* p2 = dynamic_cast<DistributedImage<uint64_t>*>(p);
							if (!p2)
							{
								reason = string("Expected ") + toString(ImageUInt64) + string(", but ") + value + string(" is ") + getDataTypeString(p);
								return 0;
							}

							if (doConversion)
								result.dimg64 = p2;

							return 1;
						}

						if (dt == ImageFloat32)
						{
							DistributedImage<float32_t>* p2 = dynamic_cast<DistributedImage<float32_t>*>(p);
							if (!p2)
							{
								reason = string("Expected ") + toString(ImageFloat32) + string(", but ") + value + string(" is ") + getDataTypeString(p);
								return 0;
							}

							if (doConversion)
								result.dimgf32 = p2;

							return 1;
						}

						if (dt == ImageComplex32)
						{
							DistributedImage<complex32_t >* p2 = dynamic_cast<DistributedImage<complex32_t >*>(p);
							if (!p2)
							{
								reason = string("Expected ") + toString(ImageComplex32) + string(", but ") + value + string(" is ") + getDataTypeString(p);
								return 0;
							}

							if (doConversion)
								result.dimgc32 = p2;

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

				if (dt == ImageUInt8 ||
					dt == ImageUInt16 ||
					dt == ImageUInt32 ||
					dt == ImageUInt64 ||
					dt == ImageFloat32 ||
					dt == ImageComplex32)
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
							if ((dt == ImageUInt8 && dynamic_cast<Image<uint8_t>*>(p)) ||
								(dt == ImageUInt16 && dynamic_cast<Image<uint16_t>*>(p)) ||
								(dt == ImageUInt32 && dynamic_cast<Image<uint32_t>*>(p)) ||
								(dt == ImageUInt64 && dynamic_cast<Image<uint64_t>*>(p)) ||
								(dt == ImageFloat32 && dynamic_cast<Image<float32_t>*>(p)) ||
								(dt == ImageComplex32 && dynamic_cast<Image<complex32_t>*>(p)))
							{
								// Match
								if (doConversion)
									result.imgval = p;
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
							if ((dt == ImageUInt8 && dynamic_cast<DistributedImage<uint8_t>*>(p)) ||
								(dt == ImageUInt16 && dynamic_cast<DistributedImage<uint16_t>*>(p)) ||
								(dt == ImageUInt32 && dynamic_cast<DistributedImage<uint32_t>*>(p)) ||
								(dt == ImageUInt64 && dynamic_cast<DistributedImage<uint64_t>*>(p)) ||
								(dt == ImageFloat32 && dynamic_cast<DistributedImage<float32_t>*>(p)) ||
								(dt == ImageComplex32 && dynamic_cast<DistributedImage<complex32_t>*>(p)))
							{
								// Match
								if (doConversion)
									result.dimgval = p;
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

							if (dt == ImageUInt8)
							{
								result.img8 = new Image<uint8_t>();
								imgs[value] = result.img8;
							}
							else if (dt == ImageUInt16)
							{
								result.img16 = new Image<uint16_t>();
								imgs[value] = result.img16;
							}
							else if (dt == ImageUInt32)
							{
								result.img32 = new Image<uint32_t>();
								imgs[value] = result.img32;
							}
							else if (dt == ImageUInt64)
							{
								result.img64 = new Image<uint64_t>();
								imgs[value] = result.img64;
							}
							else if (dt == ImageFloat32)
							{
								result.imgf32 = new Image<float32_t>();
								imgs[value] = result.imgf32;
							}
							else if (dt == ImageComplex32)
							{
								result.imgc32 = new Image<complex32_t>();
								imgs[value] = result.imgc32;
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


							if (dt == ImageUInt8)
							{
								result.dimg8 = new DistributedImage<uint8_t>(value);
								distributedImgs[value] = result.dimg8;
							}
							else if (dt == ImageUInt16)
							{
								result.dimg16 = new DistributedImage<uint16_t>(value);
								distributedImgs[value] = result.dimg16;
							}
							else if (dt == ImageUInt32)
							{
								result.dimg32 = new DistributedImage<uint32_t>(value);
								distributedImgs[value] = result.dimg32;
							}
							else if (dt == ImageUInt64)
							{
								result.dimg64 = new DistributedImage<uint64_t>(value);
								distributedImgs[value] = result.dimg64;
							}
							else if (dt == ImageFloat32)
							{
								result.dimgf32 = new DistributedImage<float32_t>(value);
								distributedImgs[value] = result.dimgf32;
							}
							else if (dt == ImageComplex32)
							{
								result.dimgc32 = new DistributedImage<complex32_t>(value);
								distributedImgs[value] = result.dimgc32;
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
		int matchParameterTypes(const vector<CommandArgumentBase>& types, vector<string>& values, string& reason)
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
		Removes duplicate elements from the given list.
		*/
		template<typename T> static void removeDuplicates(vector<T>& v)
		{
			for (size_t n = 0; n < v.size(); n++)
			{
				for (size_t m = n + 1; m < v.size(); m++)
				{
					if (v[m] == v[n])
					{
						v.erase(v.begin() + m);
						m--;
					}
				}
			}
		}

		/**
		Finds some command of given priority from the given list, and returns count of items with given priority.
		*/
		size_t getByPriority(const vector<tuple<int, Command*> >& candidates, int allowedPriority, Command*& command)
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
		void executeCommand(const string& name, vector<string>& args)
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
			if(count <= 0)
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
					dist->runDistributedInternal(this, *distributor, convertedArgs);
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
		void parseLine(string& line)
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
		static string getDataTypeString(const ImageBase* p)
		{
			const Image<uint8_t>* p8 = dynamic_cast<const Image<uint8_t>*>(p);
			const Image<uint16_t>* p16 = dynamic_cast<const Image<uint16_t>*>(p);
			const Image<uint32_t>* p32 = dynamic_cast<const Image<uint32_t>*>(p);
			const Image<uint64_t>* p64 = dynamic_cast<const Image<uint64_t>*>(p);
			const Image<float32_t>* pf32 = dynamic_cast<const Image<float32_t>*>(p);
			const Image<complex32_t >* pc32 = dynamic_cast<const Image<complex32_t >*>(p);

			if (p8)
				return toString(ImageUInt8);
			if (p16)
				return toString(ImageUInt16);
			if (p32)
				return toString(ImageUInt32);
			if (p64)
				return toString(ImageUInt64);
			if (pf32)
				return toString(ImageFloat32);
			if (pc32)
				return toString(ImageComplex32);

			throw runtime_error("getDataTypeString not implemented for given data type.");
		}


		/**
		Gets a string representing data type of given image.
		*/
		static string getDataTypeString(const DistributedImageBase* p)
		{
			const DistributedImage<uint8_t>* p8 = dynamic_cast<const DistributedImage<uint8_t>*>(p);
			const DistributedImage<uint16_t>* p16 = dynamic_cast<const DistributedImage<uint16_t>*>(p);
			const DistributedImage<uint32_t>* p32 = dynamic_cast<const DistributedImage<uint32_t>*>(p);
			const DistributedImage<uint64_t>* p64 = dynamic_cast<const DistributedImage<uint64_t>*>(p);
			const DistributedImage<float32_t>* pf32 = dynamic_cast<const DistributedImage<float32_t>*>(p);
			const DistributedImage<complex32_t >* pc32 = dynamic_cast<const DistributedImage<complex32_t >*>(p);

			if (p8)
				return toString(ImageUInt8);
			if (p16)
				return toString(ImageUInt16);
			if (p32)
				return toString(ImageUInt32);
			if (p64)
				return toString(ImageUInt64);
			if (pf32)
				return toString(ImageFloat32);
			if (pc32)
				return toString(ImageComplex32);

			throw runtime_error("getDataTypeString not implemented for given data type.");
		}
		

	public:

		PISystem()
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
		}

		~PISystem()
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
		vector<string> getImageNames() const
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
		string getImageName(const ImageBase* img) const
		{
			for (auto const& item : imgs)
				if (item.second == img)
					return item.first;
			throw ITLException("Invalid image name query.");
		}

		/**
		Returns names of distributed images in the system.
		*/
		vector<string> getDistributedImageNames() const
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
		ImageBase* getImage(const string& name)
		{
			if (isDistributed())
				throw runtime_error("Images cannot be retrieved when distributed processing mode is activated.");

			return imgs[name];
		}

		/**
		Retrieve image having given name, do not throw exception but set last error if exception would be thrown.
		*/
		ImageBase* getImageNoThrow(const string& name)
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
				if (lastException == "")
					lastException = "Image not found.";
			}
			return 0;
		}

		/**
		Replace image with given image.
		Replace image by null pointer to remove it from the system.
		*/
		void replaceImage(const string& name, ImageBase* newImg)
		{
			if (imgs.find(name) != imgs.end())
			{
				if(imgs[name] != newImg)
					delete imgs[name];
				imgs.erase(name);
			}
			if (newImg)
				imgs[name] = newImg;
		}

		/**
		Get distributed image having given name.
		*/
		DistributedImageBase* getDistributedImage(const string& name)
		{
			return distributedImgs[name];
		}

		/**
		Replace distributed image with a new one.
		Set newImg to zero to remove the image with given name.
		*/
		void replaceDistributedImage(const string& name, DistributedImageBase* newImg, bool free = true)
		{
			if (distributedImgs.find(name) != distributedImgs.end())
			{
				if(distributedImgs[name] != newImg && free)
					delete distributedImgs[name];
				distributedImgs.erase(name);
			}
			if(newImg)
				distributedImgs[name] = newImg;
		}

		/**
		Gets error message for last exception that occured.
		*/
		const char* getLastErrorMessage() const
		{
			return lastException.c_str();
		}

		/**
		Gets line of code that caused last error.
		*/
		int getLastErrorLine() const
		{
			return lastExceptionLine;
		}

		/**
		Clear last error message.
		*/
		void clearLastError()
		{
			lastException = "";
			lastExceptionLine = 0;
		}

		/**
		Sets show commands flag.
		*/
		void showCommands(bool echo, bool timing)
		{
			showRunCommands = echo;
			showTiming = timing;
		}

		/**
		Gets a value indicating whethe distributed processing mode is active.
		*/
		bool isDistributed() const
		{
			return distributor != 0;
		}

		/**
		Enables or disables distributed processing.
		@param provider Provider name. Pass empty string to disable distributed processing.
		*/
		void distributedProcessing(string provider)
		{
			if (distributor)
			{
				delete distributor;
				distributor = 0;
				cout << "Disabled distributed computing mode." << endl;
			}

			if (provider.length() > 0)
			{
				toLower(provider);
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
		bool run(const string& commands)
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
		vector<string> getHelp(const string& cmdName) const
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
		string commandList(bool forUser) const
		{
			vector<tuple<string, bool> > names;
			for (size_t n = 0; n < commands.size(); n++)
			{
				names.push_back(make_tuple(commands[n]->toSimpleString(), commands[n]->canDistribute()));
			}
			removeDuplicates(names);
			sort(names.begin(), names.end());
			stringstream msg;
			if(forUser)
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
		string argumentToString(const CommandArgumentBase& argument, const ParamVariant& value, bool isDistributed) const
		{
			ArgumentDataType dt = argument.dataType();

			if (isDistributed)
			{
				switch (dt)
				{
				case ImageUInt8: dt = DImageUInt8; break;
				case ImageUInt16: dt = DImageUInt16; break;
				case ImageUInt32: dt = DImageUInt32; break;
				case ImageUInt64: dt = DImageUInt64; break;
				case ImageFloat32: dt = DImageFloat32; break;
				case ImageComplex32: dt = DImageComplex32; break;
				}
			}

			switch (dt)
			{
			case String: return get<string>(value);
			case Double: return itl2::toString(get<double>(value));
			case Int: return itl2::toString(get<coord_t>(value));
			case Size: return itl2::toString(get<size_t>(value));
			case NBType: return itl2::toString(get<NeighbourhoodType>(value));
			case BoundaryCond: return itl2::toString(get<BoundaryCondition>(value));
			case Bool: return itl2::toString(get<bool>(value));
			case Vect3d: return itl2::toString(get<Vec3d>(value));
			case Vect3c: return itl2::toString(get<Vec3c>(value));
			case ImageUInt8: return imageName(get<Image<uint8_t>* >(value));
			case ImageUInt16: return imageName(get<Image<uint16_t>* >(value));
			case ImageUInt32: return imageName(get<Image<uint32_t>* >(value));
			case ImageUInt64: return imageName(get<Image<uint64_t>* >(value));
			case ImageFloat32: return imageName(get<Image<float32_t>* >(value));
			case ImageComplex32: return imageName(get<Image<complex32_t>* >(value));
			case DImageUInt8: return distributedImageName(get<DistributedImage<uint8_t>* >(value));
			case DImageUInt16: return distributedImageName(get<DistributedImage<uint16_t>* >(value));
			case DImageFloat32: return distributedImageName(get<DistributedImage<float32_t>* >(value));
			case DImageComplex32: return distributedImageName(get<DistributedImage<complex32_t>* >(value));
			default: throw ITLException("Data type not configured.");
			}
		}

		
	};

	/**
	Convert command.
	(It is hard to find a place to define this command as it has to be in a header, and PISystem and Command must have been fully declared and this command is a template.)
	*/
	template<typename in_t> class ConvertCommand : public Command, public Distributable
	{
	public:
		ConvertCommand() : Command("convert", "Converts data type of input image.",
		{
			CommandArgument<Image<in_t> >(In, "input image", "Input image."),
			CommandArgument<string>(In, "output image", "Output image."),
			CommandArgument<string>(In, "data type", "Data type of the output image. Can be uint8, uint16, or float32.")
		})
		{
		}

		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const
		{
			Image<in_t>* in = pop<Image<in_t>*>(args);
			string outname = pop<string>(args);
			string dts = pop<string>(args);

			if (outname == system->getImageName(in))
				throw ITLException("Unable to convert image in-place. Please specify another target image name.");

			ImageDataType dt = fromString(dts);

			const vector<string>& names = system->getImageNames();
			auto it = find(names.begin(), names.end(), outname);
			ImageBase* pReplaceThis = 0;
			if (it != names.end())
			{
				ImageBase* p = system->getImage(outname);
				if (p->sizeEquals(in->dimensions()) && p->dataType() == dt)
				{
					pReplaceThis = p;
				}
			}

			if (dt == UInt8)
			{
				itl2::Image<uint8_t>* img = pReplaceThis != 0 ? (itl2::Image<uint8_t>*)pReplaceThis : new itl2::Image<uint8_t>();
				convert(*in, *img);
				system->replaceImage(outname, img);
			}
			else if (dt == UInt16)
			{
				itl2::Image<uint16_t>* img = pReplaceThis != 0 ? (itl2::Image<uint16_t>*)pReplaceThis : new itl2::Image<uint16_t>();
				convert(*in, *img);
				system->replaceImage(outname, img);
			}
			else if (dt == UInt32)
			{
				itl2::Image<uint32_t>* img = pReplaceThis != 0 ? (itl2::Image<uint32_t>*)pReplaceThis : new itl2::Image<uint32_t>();
				convert(*in, *img);
				system->replaceImage(outname, img);
			}
			else if (dt == UInt64)
			{
				itl2::Image<uint64_t>* img = pReplaceThis != 0 ? (itl2::Image<uint64_t>*)pReplaceThis : new itl2::Image<uint64_t>();
				convert(*in, *img);
				system->replaceImage(outname, img);
			}
			else if (dt == Float32)
			{
				itl2::Image<float32_t>* img = pReplaceThis != 0 ? (itl2::Image<float32_t>*)pReplaceThis : new itl2::Image<float32_t>();
				convert(*in, *img);
				system->replaceImage(outname, img);
			}
			else
				throw ITLException(string("Invalid data type: ") + dts);
		}

		virtual void runDistributedInternal(PISystem* system, Distributor& distributor, vector<ParamVariant>& args) const
		{
			DistributedImage<in_t>* in = get<DistributedImage<in_t>*>(args[0]);
			string outname = get<string>(args[1]);
			string dt = get<string>(args[2]);

			if (outname == in->varName())
				throw ITLException("Unable to convert image in-place. Please specify another target image name.");

			ParamVariant inImg;
			inImg.dimgval = in;
			ParamVariant outImg;

			vector<ParamVariant> dargs;
			dargs.push_back(inImg);

			if (dt == "uint8")
			{
				DistributedImage<uint8_t>* img = new DistributedImage<uint8_t>(outname, in->dimensions()[0], in->dimensions()[1], in->dimensions()[2]);
				system->replaceDistributedImage(outname, img);
				outImg.dimgval = img;
				dargs.push_back(outImg);

				// distribute in z, no overlap
				CopyCommand<in_t, uint8_t> cmd;
				distributor.distribute(&cmd, dargs, 2, Vec3c(0, 0, 0));
			}
			else if (dt == "uint16")
			{
				DistributedImage<uint16_t>* img = new DistributedImage<uint16_t>(outname, in->dimensions()[0], in->dimensions()[1], in->dimensions()[2]);
				system->replaceDistributedImage(outname, img);
				outImg.dimgval = img;
				dargs.push_back(outImg);

				// distribute in z, no overlap
				CopyCommand<in_t, uint16_t> cmd;
				distributor.distribute(&cmd, dargs, 2, Vec3c(0, 0, 0));
			}
			else if (dt == "uint32")
			{
				DistributedImage<uint32_t>* img = new DistributedImage<uint32_t>(outname, in->dimensions()[0], in->dimensions()[1], in->dimensions()[2]);
				system->replaceDistributedImage(outname, img);
				outImg.dimgval = img;
				dargs.push_back(outImg);

				// distribute in z, no overlap
				CopyCommand<in_t, uint32_t> cmd;
				distributor.distribute(&cmd, dargs, 2, Vec3c(0, 0, 0));
			}
			else if (dt == "uint64")
			{
				DistributedImage<uint64_t>* img = new DistributedImage<uint64_t>(outname, in->dimensions()[0], in->dimensions()[1], in->dimensions()[2]);
				system->replaceDistributedImage(outname, img);
				outImg.dimgval = img;
				dargs.push_back(outImg);

				// distribute in z, no overlap
				CopyCommand<in_t, uint64_t> cmd;
				distributor.distribute(&cmd, dargs, 2, Vec3c(0, 0, 0));
			}
			else if (dt == "float32")
			{
				DistributedImage<float32_t>* img = new DistributedImage<float32_t>(outname, in->dimensions()[0], in->dimensions()[1], in->dimensions()[2]);
				system->replaceDistributedImage(outname, img);
				outImg.dimgval = img;
				dargs.push_back(outImg);

				// distribute in z, no overlap
				CopyCommand<in_t, float32_t> cmd;
				distributor.distribute(&cmd, dargs, 2, Vec3c(0, 0, 0));
			}
			else
				throw ITLException(string("Invalid data type: ") + dt);

		}

		virtual void run(vector<ParamVariant>& args) const
		{
		}

		virtual void runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
		{
		}
	};

	template<typename pixel_t> class NewLikeCommand : public Command, public Distributable
	{

	public:
		NewLikeCommand() : Command("newlike", "Creates a new, empty image that has properties (dimensions, data type) similar to another image.",
			{
				CommandArgument<Image<pixel_t> >(In, "template image name", "Existing image where dimensions and data type will be copied from."),
				CommandArgument<string>(In, "image name", "Name of the new image in the system."),
				CommandArgument<string>(In, "data type", "Data type of the image. Can be uint8, uint16, uint32, uint64, or float32. Leave empty or set to Unknown to copy the value from the template image.", ""),
				CommandArgument<coord_t>(In, "width", "Width of the image. Set to zero to copy the value from the template image.", 0),
				CommandArgument<coord_t>(In, "height", "Height of the image. Set to zero to copy the value from the template image.", 0),
				CommandArgument<coord_t>(In, "depth", "Depth of the image. Set to zero to copy the value from the template image.", 0)
			})
		{
		}

		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const
		{
			Image<pixel_t>& templ = *pop<Image<pixel_t>* >(args);
			string name = pop<string>(args);
			string dts = pop<string>(args);
			coord_t w = pop<coord_t>(args);
			coord_t h = pop<coord_t>(args);
			coord_t d = pop<coord_t>(args);

			ImageDataType dt = fromString(dts);

			if (dt == Unknown)
				dt = templ.dataType();
			if (w <= 0)
				w = templ.width();
			if (h <= 0)
				h = templ.height();
			if (d <= 0)
				d = templ.depth();

			if (dt == UInt8)
				system->replaceImage(name, new itl2::Image<uint8_t>(w, h, d));
			else if (dt == UInt16)
				system->replaceImage(name, new itl2::Image<uint16_t>(w, h, d));
			else if (dt == UInt32)
				system->replaceImage(name, new itl2::Image<uint32_t>(w, h, d));
			else if (dt == UInt64)
				system->replaceImage(name, new itl2::Image<uint64_t>(w, h, d));
			else if (dt == Float32)
				system->replaceImage(name, new itl2::Image<float32_t>(w, h, d));
			else if (dt == Complex32)
				system->replaceImage(name, new itl2::Image<complex32_t>(w, h, d));
			else
				throw ParseException(string("Invalid data type: ") + dts);
		}

		virtual void run(vector<ParamVariant>& args) const
		{
		}

		virtual void runDistributedInternal(PISystem* system, Distributor& distributor, vector<ParamVariant>& args) const
		{
			DistributedImage<pixel_t>& templ = *pop<DistributedImage<pixel_t>* >(args);
			string name = pop<string>(args);
			string dts = pop<string>(args);
			coord_t w = pop<coord_t>(args);
			coord_t h = pop<coord_t>(args);
			coord_t d = pop<coord_t>(args);

			ImageDataType dt = fromString(dts);

			if (dt == Unknown)
				dt = templ.dataType();
			if (w <= 0)
				w = templ.width();
			if (h <= 0)
				h = templ.height();
			if (d <= 0)
				d = templ.depth();

			if (dt == UInt8)
				system->replaceDistributedImage(name, new DistributedImage<uint8_t>(name, w, h, d));
			else if (dt == UInt16)
				system->replaceDistributedImage(name, new DistributedImage<uint16_t>(name, w, h, d));
			else if (dt == UInt32)
				system->replaceDistributedImage(name, new DistributedImage<uint32_t>(name, w, h, d));
			else if (dt == UInt64)
				system->replaceDistributedImage(name, new DistributedImage<uint64_t>(name, w, h, d));
			else if (dt == Float32)
				system->replaceDistributedImage(name, new DistributedImage<float32_t>(name, w, h, d));
			else if (dt == Complex32)
				system->replaceDistributedImage(name, new DistributedImage<complex32_t>(name, w, h, d));
			else
				throw ParseException(string("Invalid data type: ") + dts);
		}

		virtual void runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
		{
		}
	};
}
