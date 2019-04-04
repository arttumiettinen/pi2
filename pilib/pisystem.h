#pragma once

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <sstream>
#include <iomanip>
#include <algorithm>

#include "command.h"
#include "parseexception.h"
#include "timer.h"
#include "distributor.h"
#include "distributable.h"
#include "distributedimage.h"
#include "slurmdistributor.h"
#include "localdistributor.h"
#include "pointprocesscommands.h"
#include "stringutils.h"
#include "conversions.h"

using namespace std;
using namespace itl2;

namespace pilib
{

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
		Parse line expected to contain function call
		funcname(param1, param2, param3, ...)
		*/
		static void parseFunctionCall(const string& line, string& name, vector<string>& args);

		/**
		Get all commands whose name is the given one.
		*/
		vector<Command*> getCommands(const string& name) const;

		/**
		Tests if the given value is valid image name.
		*/
		static bool isValidImageName(const string& value, string& reason);

		/**
		Converts image to variable name.
		*/
		string imageName(const ImageBase* img) const;

		/**
		Converts distributed image to variable name.
		*/
		string distributedImageName(const DistributedImageBase* img) const;

		/**
		Finds image with given name, checks that it has three pixels, and converts it to Vec3.
		*/
		bool getImageAsVec3(const string& name, Vec3d& v);

		/**
		Try to convert the given string to match the given argument.
		The conversion is done and, e.g., images are created if doConversion flag is set to true.
		The result of the conversion is assigned to the result variable.
		If conversion is not possible, reason string is assigned an explanation of the error.
		@return 0 if there is no match, 1 if the argument type and parameter type match, and 2 if they match after creation of new images.
		*/
		int tryConvert(string& value, const CommandArgumentBase& type, bool doConversion, ParamVariant& result, string& reason);

		/**
		Checks whether supplied string values can be converted to supplied types.
		Values are not changed but the reference is non-const as required by tryConvert.
		If 0 is returned, reason parameter is assigned an explanation why this match does not succeed.
		@return 0 if there is no match, 1 if there is match, and 2 if there is match after creation of new images.
		*/
		int matchParameterTypes(const vector<CommandArgumentBase>& types, vector<string>& values, string& reason);

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
		size_t getByPriority(const vector<tuple<int, Command*> >& candidates, int allowedPriority, Command*& command);

		/**
		Execute given command with given arguments.
		Searches the command based on name and arguments, and if found, executes it.
		If not found, throws ParseException.
		*/
		void executeCommand(const string& name, vector<string>& args);

		/**
		Parse one line (one command) of input code.
		*/
		void parseLine(string& line);


		/**
		Gets a string representing data type of given image.
		*/
		static string getDataTypeString(const ImageBase* p);


		/**
		Gets a string representing data type of given image.
		*/
		static string getDataTypeString(const DistributedImageBase* p);
		

	public:

		PISystem();

		~PISystem();

		/**
		Returns names of images in the system.
		*/
		vector<string> getImageNames() const;

		/**
		Returns name of the given image.
		*/
		string getImageName(const ImageBase* img) const;

		/**
		Returns names of distributed images in the system.
		*/
		vector<string> getDistributedImageNames() const;

		/**
		Retrieve image having given name.
		*/
		ImageBase* getImage(const string& name);

		/**
		Flush local changes to distributed image.
		Assumes that the distributed image has been retrieved using getImage.
		Deletes the local image.
		*/
		void flushIfDistributed(const string& imgName);

		/**
		Same than flushIfDistributed but sets lastException instead of throwing it.
		@return false if an error occurs.
		*/
		bool flushIfDistributedNoThrow(const string& name);

		/**
		Retrieve image having given name, do not throw exception but set last error and return 0 if exception would be thrown.
		*/
		ImageBase* getImageNoThrow(const string& name);

		/**
		Replace image with given image.
		Replace image by null pointer to remove it from the system.
		*/
		void replaceImage(const string& name, ImageBase* newImg);

		/**
		Get distributed image having given name.
		*/
		DistributedImageBase* getDistributedImage(const string& name);

		/**
		Replace distributed image with a new one.
		Set newImg to zero to remove the image with given name.
		*/
		void replaceDistributedImage(const string& name, DistributedImageBase* newImg, bool free = true);

		/**
		Gets error message for last exception that occured.
		*/
		const char* getLastErrorMessage() const;

		/**
		Gets line of code that caused last error.
		*/
		int getLastErrorLine() const;

		/**
		Clear last error message.
		*/
		void clearLastError();

		/**
		Sets show commands flag.
		*/
		void showCommands(bool echo, bool timing);

		/**
		Gets a value indicating whethe distributed processing mode is active.
		*/
		bool isDistributed() const;

		/**
		Enables or disables distributed processing.
		@param provider Provider name. Pass empty string to disable distributed processing.
		*/
		void distributedProcessing(string provider);

		/**
		Parses commands in the given string and runs them.
		*/
		bool run(const string& commands);

		/**
		Gets help topics for given command name.
		*/
		vector<string> getHelp(const string& cmdName) const;

		/**
		Get a list of commands and some basic information about each one, each command separated by newline.
		*/
		string commandList(bool forUser) const;

		/**
		Convert argument to string.
		@param isDistributed Set to true to convert arguments of type Image* to DistributedImage*
		*/
		string argumentToString(const CommandArgumentBase& argument, const ParamVariant& value, bool isDistributed) const;

		
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
			CommandArgument<Image<in_t> >(ParameterDirection::In, "input image", "Input image."),
			CommandArgument<string>(ParameterDirection::In, "output image", "Output image."),
			CommandArgument<string>(ParameterDirection::In, "data type", "Data type of the output image. Can be uint8, uint16, or float32.")
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

			ImageDataType dt = fromString<ImageDataType>(dts);

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

			if (dt == ImageDataType::UInt8)
			{
				itl2::Image<uint8_t>* img = pReplaceThis != 0 ? (itl2::Image<uint8_t>*)pReplaceThis : new itl2::Image<uint8_t>();
				convert(*in, *img);
				system->replaceImage(outname, img);
			}
			else if (dt == ImageDataType::UInt16)
			{
				itl2::Image<uint16_t>* img = pReplaceThis != 0 ? (itl2::Image<uint16_t>*)pReplaceThis : new itl2::Image<uint16_t>();
				convert(*in, *img);
				system->replaceImage(outname, img);
			}
			else if (dt == ImageDataType::UInt32)
			{
				itl2::Image<uint32_t>* img = pReplaceThis != 0 ? (itl2::Image<uint32_t>*)pReplaceThis : new itl2::Image<uint32_t>();
				convert(*in, *img);
				system->replaceImage(outname, img);
			}
			else if (dt == ImageDataType::UInt64)
			{
				itl2::Image<uint64_t>* img = pReplaceThis != 0 ? (itl2::Image<uint64_t>*)pReplaceThis : new itl2::Image<uint64_t>();
				convert(*in, *img);
				system->replaceImage(outname, img);
			}
			else if (dt == ImageDataType::Float32)
			{
				itl2::Image<float32_t>* img = pReplaceThis != 0 ? (itl2::Image<float32_t>*)pReplaceThis : new itl2::Image<float32_t>();
				convert(*in, *img);
				system->replaceImage(outname, img);
			}
			else
				throw ITLException(string("Invalid data type: ") + dts);
		}

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
		{
			DistributedImage<in_t>* in = get<DistributedImage<in_t>*>(args[0]);
			string outname = get<string>(args[1]);
			string dt = get<string>(args[2]);

			if (outname == in->varName())
				throw ITLException("Unable to convert image in-place. Please specify another target image name.");

			PISystem* system = distributor.getSystem();

			ParamVariant inImg;
			//inImg.dimgval = in;
			inImg = in;
			ParamVariant outImg;

			vector<ParamVariant> dargs;
			dargs.push_back(inImg);

			if (dt == "uint8")
			{
				DistributedImage<uint8_t>* img = new DistributedImage<uint8_t>(outname, in->dimensions()[0], in->dimensions()[1], in->dimensions()[2]);
				system->replaceDistributedImage(outname, img);
				//outImg.dimgval = img;
				outImg = img;
				dargs.push_back(outImg);

				// distribute in z, no overlap
				CopyCommand<in_t, uint8_t> cmd;
				distributor.distribute(&cmd, dargs, 2, Vec3c(0, 0, 0));
			}
			else if (dt == "uint16")
			{
				DistributedImage<uint16_t>* img = new DistributedImage<uint16_t>(outname, in->dimensions()[0], in->dimensions()[1], in->dimensions()[2]);
				system->replaceDistributedImage(outname, img);
				//outImg.dimgval = img;
				outImg = img;
				dargs.push_back(outImg);

				// distribute in z, no overlap
				CopyCommand<in_t, uint16_t> cmd;
				distributor.distribute(&cmd, dargs, 2, Vec3c(0, 0, 0));
			}
			else if (dt == "uint32")
			{
				DistributedImage<uint32_t>* img = new DistributedImage<uint32_t>(outname, in->dimensions()[0], in->dimensions()[1], in->dimensions()[2]);
				system->replaceDistributedImage(outname, img);
				//outImg.dimgval = img;
				outImg = img;
				dargs.push_back(outImg);

				// distribute in z, no overlap
				CopyCommand<in_t, uint32_t> cmd;
				distributor.distribute(&cmd, dargs, 2, Vec3c(0, 0, 0));
			}
			else if (dt == "uint64")
			{
				DistributedImage<uint64_t>* img = new DistributedImage<uint64_t>(outname, in->dimensions()[0], in->dimensions()[1], in->dimensions()[2]);
				system->replaceDistributedImage(outname, img);
				//outImg.dimgval = img;
				outImg = img;
				dargs.push_back(outImg);

				// distribute in z, no overlap
				CopyCommand<in_t, uint64_t> cmd;
				distributor.distribute(&cmd, dargs, 2, Vec3c(0, 0, 0));
			}
			else if (dt == "float32")
			{
				DistributedImage<float32_t>* img = new DistributedImage<float32_t>(outname, in->dimensions()[0], in->dimensions()[1], in->dimensions()[2]);
				system->replaceDistributedImage(outname, img);
				//outImg.dimgval = img;
				outImg = img;
				dargs.push_back(outImg);

				// distribute in z, no overlap
				CopyCommand<in_t, float32_t> cmd;
				distributor.distribute(&cmd, dargs, 2, Vec3c(0, 0, 0));
			}
			else
				throw ITLException(string("Invalid data type: ") + dt);

			return vector<string>();
		}

		virtual void run(vector<ParamVariant>& args) const
		{
		}
	};



	class NewLikeBase : public Command, public Distributable
	{
	protected:

		void createImage(const string& name, const string& dts, coord_t w, coord_t h, coord_t d,
							ImageDataType templDt, coord_t templw, coord_t templh, coord_t templd,
							PISystem* system) const
		{
			ImageDataType dt = fromString<ImageDataType>(dts);

			if (dt == ImageDataType::Unknown)
				dt = templDt;
			if (w <= 0)
				w = templw;
			if (h <= 0)
				h = templh;
			if (d <= 0)
				d = templd;

			if (dt == ImageDataType::UInt8)
				system->replaceImage(name, new itl2::Image<uint8_t>(w, h, d));
			else if (dt == ImageDataType::UInt16)
				system->replaceImage(name, new itl2::Image<uint16_t>(w, h, d));
			else if (dt == ImageDataType::UInt32)
				system->replaceImage(name, new itl2::Image<uint32_t>(w, h, d));
			else if (dt == ImageDataType::UInt64)
				system->replaceImage(name, new itl2::Image<uint64_t>(w, h, d));
			else if (dt == ImageDataType::Float32)
				system->replaceImage(name, new itl2::Image<float32_t>(w, h, d));
			else if (dt == ImageDataType::Complex32)
				system->replaceImage(name, new itl2::Image<complex32_t>(w, h, d));
			else
				throw ParseException(string("Invalid data type: ") + dts);
		}

		void createDistributedImage(const string& name, const string& dts, coord_t w, coord_t h, coord_t d,
			ImageDataType templDt, coord_t templw, coord_t templh, coord_t templd,
			PISystem* system) const
		{
			ImageDataType dt = fromString<ImageDataType>(dts);

			if (dt == ImageDataType::Unknown)
				dt = templDt;
			if (w <= 0)
				w = templw;
			if (h <= 0)
				h = templh;
			if (d <= 0)
				d = templd;

			if (dt == ImageDataType::UInt8)
				system->replaceDistributedImage(name, new DistributedImage<uint8_t>(name, w, h, d));
			else if (dt == ImageDataType::UInt16)
				system->replaceDistributedImage(name, new DistributedImage<uint16_t>(name, w, h, d));
			else if (dt == ImageDataType::UInt32)
				system->replaceDistributedImage(name, new DistributedImage<uint32_t>(name, w, h, d));
			else if (dt == ImageDataType::UInt64)
				system->replaceDistributedImage(name, new DistributedImage<uint64_t>(name, w, h, d));
			else if (dt == ImageDataType::Float32)
				system->replaceDistributedImage(name, new DistributedImage<float32_t>(name, w, h, d));
			else if (dt == ImageDataType::Complex32)
				system->replaceDistributedImage(name, new DistributedImage<complex32_t>(name, w, h, d));
			else
				throw ParseException(string("Invalid data type: ") + dts);
		}

	public:
		NewLikeBase(const string& name, const string& help, CommandArgumentBase templParam) : Command(name, help,
			{
				CommandArgument<string>(ParameterDirection::In, "image name", "Name of the new image in the system."),
				templParam,
				CommandArgument<string>(ParameterDirection::In, "data type", "Data type of the image. Can be uint8, uint16, uint32, uint64, float32, or complex32. Leave empty or set to Unknown to copy the value from the template image.", ""),
				CommandArgument<coord_t>(ParameterDirection::In, "width", "Width of the image. Set to zero to copy the value from the template image.", 0),
				CommandArgument<coord_t>(ParameterDirection::In, "height", "Height of the image. Set to zero to copy the value from the template image.", 0),
				CommandArgument<coord_t>(ParameterDirection::In, "depth", "Depth of the image. Set to zero to copy the value from the template image.", 0)
			})
		{
		}

		virtual void run(vector<ParamVariant>& args) const
		{
		}
	};


	template<typename pixel_t> class NewLikeCommand : public NewLikeBase
	{

	public:
		NewLikeCommand() : NewLikeBase("newlike",
										"Creates a new, empty image that has properties (dimensions, data type) similar to another image.",
										CommandArgument<Image<pixel_t> >(ParameterDirection::In, "template image", "Name of existing image where dimensions and data type will be copied from."))
		{
		}

		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const
		{
			string name = pop<string>(args);
			Image<pixel_t>& templ = *pop<Image<pixel_t>* >(args);
			string dts = pop<string>(args);
			coord_t w = pop<coord_t>(args);
			coord_t h = pop<coord_t>(args);
			coord_t d = pop<coord_t>(args);

			createImage(name, dts, w, h, d, templ.dataType(), templ.width(), templ.height(), templ.depth(), system);
		}

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
		{
			string name = pop<string>(args);
			DistributedImage<pixel_t>& templ = *pop<DistributedImage<pixel_t>* >(args);
			string dts = pop<string>(args);
			coord_t w = pop<coord_t>(args);
			coord_t h = pop<coord_t>(args);
			coord_t d = pop<coord_t>(args);

			PISystem* system = distributor.getSystem();

			createDistributedImage(name, dts, w, h, d, templ.dataType(), templ.width(), templ.height(), templ.depth(), system);

			return vector<string>();
		}
	};

	class NewLikeFileCommand : public NewLikeBase
	{

	public:
		NewLikeFileCommand() : NewLikeBase("newlikefile",
											"Creates a new, empty image that has properties (dimensions, data type) similar to another image that has been saved to disk.",
											CommandArgument<string>(ParameterDirection::In, "filename", "Name of existing image file where dimensions and data type will be copied from."))
		{
		}

		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const;

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const;
	};
}




