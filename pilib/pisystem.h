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
#include "stringutils.h"
#include "commandlist.h"
#include "pick.h"

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
		std::map<std::string, std::pair<ArgumentDataType, std::shared_ptr<ParamVariant>>> namedValues;

		/**
		Maps image variable names to distributed image objects.
		*/
		std::map<std::string, std::shared_ptr<DistributedImageBase> > distributedImgs;

		/**
		Temporary storage for distributed images that may be deleted during tryConvert calls.
		*/
		std::vector<std::shared_ptr<DistributedImageBase> > distributedImageStore;

		/**
		Last exception that has occured.
		*/
		std::string lastException;
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
		static void parseFunctionCall(const std::string& line, std::string& name, std::vector<std::string>& args);

		/**
		Finds some command of given priority from the given list, and returns count of items with given priority.
		*/
		static size_t getByPriority(const std::vector<std::tuple<int, Command*> >& candidates, int allowedPriority, Command*& command);

		/**
		Tests if the given value is valid image name.
		*/
		static bool isValidImageName(const std::string& value, std::string& reason);

		/**
		Finds image with given name, checks that it has three pixels, and converts it to Vec3.
		*/
		bool getImageAsVec3(const std::string& name, Vec3d& v);

		/**
		Try to convert the given string to match the given argument.
		The conversion is done and, e.g., images are created if doConversion flag is set to true.
		The result of the conversion is assigned to the result variable.
		If conversion is not possible, reason string is assigned an explanation of the error.
		@return 0 if there is no match, 1 if the argument type and parameter type match, and 2 if they match after creation of new images.
		*/
		int tryConvert(std::string& value, const CommandArgumentBase& type, bool doConversion, ParamVariant& result, std::shared_ptr<ParamVariant> resultPtr, std::string& reason);

		/**
		Checks whether supplied string values can be converted to supplied types.
		Values are not changed but the reference is non-const as required by tryConvert.
		If 0 is returned, reason parameter is assigned an explanation why this match does not succeed.
		@return 0 if there is no match, 1 if there is match, and 2 if there is match after creation of new images.
		*/
		int matchParameterTypes(const std::vector<CommandArgumentBase>& types, std::vector<std::string>& values, std::string& reason);

		/**
		Execute given command with given arguments.
		Searches the command based on name and arguments, and if found, executes it.
		If not found, throws ParseException.
		*/
		void executeCommand(const std::string& name, std::vector<std::string>& args);

		/**
		Parses one statement of input, e.g. read(img, abc);
		*/
		void parseStatement(std::string& statement);

		/**
		Parse one line of input code.
		Tests if the line is comment, if not, separates line to individual statements and passes them to parseStatement one by one.
		*/
		void parseLine(std::string& line);


		/**
		Gets a string representing data type of given image.
		*/
		static std::string getDataTypeString(const ImageBase* p);


		/**
		Gets a string representing data type of given image.
		*/
		static std::string getDataTypeString(const DistributedImageBase* p);
		

	public:

		PISystem();

		~PISystem();

		/**
		Converts image to variable name.
		*/
		std::string imageName(const ImageBase* img) const;

		/**
		Converts distributed image to variable name.
		*/
		//std::string distributedImageName(const DistributedImageBase* img) const;

		/**
		Returns names of images in the system.
		*/
		std::vector<std::string> getImageNames() const;

		/**
		Returns names of distributed images in the system.
		*/
		std::vector<std::string> getDistributedImageNames() const;

		/**
		Retrieve image having given name.
		*/
		ImageBase* getImage(const std::string& name);

		/**
		Gets smart pointer to given image. Used to store images during delayed distribution even if they are
		removed from the PISystem.
		*/
		std::shared_ptr<DistributedImageBase> getDistributedImagePointer(DistributedImageBase* img) const;


		/**
		Gets a value indicating whether the given distributed image is still accessible from this PI system.
		*/
		bool isDistributedImage(DistributedImageBase* img) const;

		/**
		Flush local changes to distributed image.
		Assumes that the distributed image has been retrieved using getImage.
		Deletes the local image.
		*/
		void flushIfDistributed(const std::string& imgName);

		/**
		Same than flushIfDistributed but sets lastException instead of throwing it.
		@return false if an error occurs.
		*/
		bool flushIfDistributedNoThrow(const std::string& name);

		/**
		Gets image information without reading distributed images to RAM.
		*/
		void getImageInfoNoThrow(const std::string& name, coord_t& width, coord_t& height, coord_t& depth, ImageDataType& dt);

		/**
		Retrieve image having given name, do not throw exception but set last error and return 0 if exception would be thrown.
		Causes reading of image data to RAM in the case of distributed processing.
		*/
		ImageBase* getImageNoThrow(const std::string& name);

		/**
		Replace image with given image.
		Replace image by null pointer to remove it from the system.
		*/
		void replaceImage(const std::string& name, std::shared_ptr<ParamVariant> img);

		void replaceNamedValue(const std::string& name, ArgumentDataType dt, std::shared_ptr<ParamVariant> newValue);

		/**
		Get distributed image having given name.
		*/
		DistributedImageBase* getDistributedImage(const std::string& name);

		/**
		Replace distributed image with a new one.
		Set newImg to zero to remove the image with given name.
		*/
		void replaceDistributedImage(const std::string& name, std::shared_ptr<DistributedImageBase> newImg);

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

		Distributor* getDistributor()
		{
			return distributor;
		}

		/**
		Enables or disables distributed processing.
		@param provider Provider name. Pass empty string to disable distributed processing.
		*/
		void distributedProcessing(std::string provider);

		/**
		Parses commands in the given string and runs them.
		*/
		bool run(const std::string& commands);


		
	};

	

	/**
	Functor that creates new image.
	*/
	template<typename pixel_t> struct CreateImage
	{
		static Image<pixel_t>* run(const Vec3c& dimensions, const std::string& imgName, PISystem* system)
		{
			//std::shared_ptr<Image<pixel_t>> img = std::make_shared<itl2::Image<pixel_t> >(dimensions);
			//system->replaceImage(imgName, img);

			ArgumentDataType argDataType = imageDataTypeToArgumentDataType(imageDataType<pixel_t>());
			std::shared_ptr<ParamVariant> p = std::make_shared<ParamVariant>(new itl2::Image<pixel_t>(dimensions));
			system->replaceNamedValue(imgName, argDataType, p);
			return std::get<Image<pixel_t>*>(*p.get());
		}
	};

	/**
	Functor that creates new distributed image.
	*/
	template<typename pixel_t> struct CreateEmptyDistributedImage
	{
		static void run(const std::string& imgName, const Vec3c& dimensions, PISystem* system)
		{
			std::shared_ptr<DistributedImageBase> ptr = std::make_shared<DistributedImage<pixel_t> >(*system->getDistributor(), imgName, dimensions);
			system->replaceDistributedImage(imgName, ptr);
		}
	};




	class NewLikeBase : public Command, public Distributable
	{
	protected:

		void createImage(const std::string& name, const std::string& dts, coord_t w, coord_t h, coord_t d,
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

			Vec3c dimensions(w, h, d);

			pick<CreateImage>(dt, dimensions, name, system);
		}

		void createDistributedImage(const std::string& name, const std::string& dts, coord_t w, coord_t h, coord_t d,
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

			Vec3c dimensions(w, h, d);

			pick<CreateEmptyDistributedImage>(dt, name, dimensions, system);
		}

		using Command::Command;

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
		}
	};


	inline std::string newlikeSeeAlso()
	{
		return "newlikefile, newimage, ensuresize";
	}

	template<typename pixel_t> class NewLikeCommand : public NewLikeBase
	{
	protected:
		friend class CommandList;

		NewLikeCommand() : NewLikeBase("newlike", "Creates a new, empty image that has properties (dimensions, data type) similar to another image.",
			{
				CommandArgument<std::string>(ParameterDirection::In, "image name", "Name of the new image in the system."),
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "template image", "Name of existing image where dimensions and data type will be copied from."),
				CommandArgument<std::string>(ParameterDirection::In, "data type", "Data type of the image. Can be " + listSupportedImageDataTypes() + ". Leave empty or set to Unknown to copy the value from the template image.", ""),
				CommandArgument<coord_t>(ParameterDirection::In, "width", "Width of the image. Set to zero to copy the value from the template image.", 0),
				CommandArgument<coord_t>(ParameterDirection::In, "height", "Height of the image. Set to zero to copy the value from the template image.", 0),
				CommandArgument<coord_t>(ParameterDirection::In, "depth", "Depth of the image. Set to zero to copy the value from the template image.", 0)
			},
			newlikeSeeAlso())
		{

		}

	public:
		virtual void runInternal(PISystem* system, std::vector<ParamVariant>& args) const override
		{
			std::string name = pop<std::string>(args);
			Image<pixel_t>& templ = *pop<Image<pixel_t>* >(args);
			std::string dts = pop<std::string>(args);
			coord_t w = pop<coord_t>(args);
			coord_t h = pop<coord_t>(args);
			coord_t d = pop<coord_t>(args);

			createImage(name, dts, w, h, d, templ.dataType(), templ.width(), templ.height(), templ.depth(), system);
		}

		virtual std::vector<std::string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const override
		{
			std::string name = pop<std::string>(args);
			DistributedImage<pixel_t>& templ = *pop<DistributedImage<pixel_t>* >(args);
			std::string dts = pop<std::string>(args);
			coord_t w = pop<coord_t>(args);
			coord_t h = pop<coord_t>(args);
			coord_t d = pop<coord_t>(args);

			PISystem* system = distributor.getSystem();

			createDistributedImage(name, dts, w, h, d, templ.dataType(), templ.width(), templ.height(), templ.depth(), system);

			return std::vector<std::string>();
		}
	};


	template<typename pixel_t> class NewLike2Command : public NewLikeBase
	{
	protected:
		friend class CommandList;

		NewLike2Command() : NewLikeBase("newlike", "Creates a new, empty image that has properties (dimensions, data type) similar to another image.",
			{
				CommandArgument<std::string>(ParameterDirection::In, "image name", "Name of the new image in the system."),
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "template image", "Name of existing image where dimensions and data type will be copied from."),
				CommandArgument<std::string>(ParameterDirection::In, "data type", "Data type of the image. Can be " + listSupportedImageDataTypes() + ". Leave empty or set to Unknown to copy the value from the template image.", ""),
				CommandArgument<Vec3c>(ParameterDirection::In, "dimensions", "Dimensions of the image. Set any component to zero to copy the value from the template image."), // NOTE: NewLikeCommand handles calls in form newlike(img1, template, uint32)
			},
			newlikeSeeAlso())
		{

		}

	public:
		virtual void runInternal(PISystem* system, std::vector<ParamVariant>& args) const override
		{
			std::string name = pop<std::string>(args);
			Image<pixel_t>& templ = *pop<Image<pixel_t>* >(args);
			std::string dts = pop<std::string>(args);
			Vec3c dim = pop<Vec3c>(args);

			createImage(name, dts, dim.x, dim.y, dim.z, templ.dataType(), templ.width(), templ.height(), templ.depth(), system);
		}

		virtual std::vector<std::string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const override
		{
			std::string name = pop<std::string>(args);
			DistributedImage<pixel_t>& templ = *pop<DistributedImage<pixel_t>* >(args);
			std::string dts = pop<std::string>(args);
			Vec3c dim = pop<Vec3c>(args);

			PISystem* system = distributor.getSystem();

			createDistributedImage(name, dts, dim.x, dim.y, dim.z, templ.dataType(), templ.width(), templ.height(), templ.depth(), system);

			return std::vector<std::string>();
		}
	};


	class NewLikeFileCommand : public NewLikeBase
	{
	protected:
		friend class CommandList;

		NewLikeFileCommand() : NewLikeBase("newlikefile", "Creates a new, empty image that has properties (dimensions, data type) similar to another image that has been saved to disk.",
			{
				CommandArgument<std::string>(ParameterDirection::In, "image name", "Name of the new image in the system."),
				CommandArgument<std::string>(ParameterDirection::In, "template filename", "Name of existing image file where dimensions and data type will be copied from."),
				CommandArgument<std::string>(ParameterDirection::In, "data type", "Data type of the image. Can be " + listSupportedImageDataTypes() + ". Leave empty or set to Unknown to copy the value from the template image.", ""),
				CommandArgument<coord_t>(ParameterDirection::In, "width", "Width of the image. Set to zero to copy the value from the template image.", 0),
				CommandArgument<coord_t>(ParameterDirection::In, "height", "Height of the image. Set to zero to copy the value from the template image.", 0),
				CommandArgument<coord_t>(ParameterDirection::In, "depth", "Depth of the image. Set to zero to copy the value from the template image.", 0)
			},
			newlikeSeeAlso())
		{
		}

	public:
		virtual void runInternal(PISystem* system, std::vector<ParamVariant>& args) const override;

		virtual std::vector<std::string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const override;
	};

	class NewLikeFile2Command : public NewLikeBase
	{
	protected:
		friend class CommandList;

		NewLikeFile2Command() : NewLikeBase("newlikefile", "Creates a new, empty image that has properties (dimensions, data type) similar to another image that has been saved to disk.",
			{
				CommandArgument<std::string>(ParameterDirection::In, "image name", "Name of the new image in the system."),
				CommandArgument<std::string>(ParameterDirection::In, "template filename", "Name of existing image file where dimensions and data type will be copied from."),
				CommandArgument<std::string>(ParameterDirection::In, "data type", "Data type of the image. Can be " + listSupportedImageDataTypes() + ". Leave empty or set to Unknown to copy the value from the template image.", ""),
				CommandArgument<Vec3c>(ParameterDirection::In, "dimensions", "Dimensions of the image. Set any component to zero to copy that from the template image.", Vec3c(0, 0, 0)),
			},
			newlikeSeeAlso())
		{
		}

	public:
		virtual void runInternal(PISystem* system, std::vector<ParamVariant>& args) const override;

		virtual std::vector<std::string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const override;
	};
}




