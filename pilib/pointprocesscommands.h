#pragma once

#include "command.h"
#include "commandsbase.h"
#include "distributable.h"
#include "pointprocess.h"
#include "math/mathutils.h"
#include "misc.h"
#include "autothreshold.h"
#include "distributedtempimage.h"
#include "projectioncommands.h"

#include <vector>

using namespace itl2;

namespace pilib
{

	/**
	Base class for point processes that process the image in-place.
	*/
	template<typename input_t> class InPlacePointProcess : public Command, public Distributable
	{
	protected:
		friend class CommandList;


		InPlacePointProcess(const string& name, const string& help, const std::vector<CommandArgumentBase>& extraArgs = {}, const std::string& seeAlso = "") :
			Command(name, help,
				concat({
					CommandArgument<Image<input_t> >(ParameterDirection::InOut, "image", "Image to process.")
					}, extraArgs),
				seeAlso
			)
		{
		}

	public:
		virtual void run(Image<input_t>& in, std::vector<ParamVariant>& args) const = 0;

		virtual void run(std::vector<ParamVariant>& args) const override
		{
			Image<input_t>& in = *pop<Image<input_t>* >(args);
			run(in, args);
		}

		using Distributable::runDistributed;

		virtual size_t getDistributionDirection2(const std::vector<ParamVariant>& args) const override
		{
			return 1;
		}

		virtual std::vector<string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const override
		{
			return distributor.distribute(this, args);
		}

		virtual JobType getJobType(const std::vector<ParamVariant>& args) const override
		{
			return JobType::Fast;
		}

		virtual bool canDelay(const std::vector<ParamVariant>& args) const override
		{
			return true;
		}
	};

	/**
	Base class for point process commands that need input and output image.
	*/
	template<typename input_t, typename output_t> class InputOutputPointProcess : public Command, public Distributable
	{
	protected:
		friend class CommandList;


		InputOutputPointProcess(const string& name, const string& help, const std::vector<CommandArgumentBase>& extraArgs = {}, const string& seeAlso = "") :
			Command(name, help,
				concat({
					CommandArgument<Image<input_t> >(ParameterDirection::In, "input image", "Input image."),
					CommandArgument<Image<output_t> >(ParameterDirection::Out, "output image", "Output image.")
					}, extraArgs),
				seeAlso
			)
		{
		}

	public:
		virtual void run(Image<input_t>& in, Image<output_t>& out, std::vector<ParamVariant>& args) const = 0;

		virtual void run(std::vector<ParamVariant>& args) const override
		{
			Image<input_t>& in = *pop<Image<input_t>* >(args);
			Image<output_t>& out = *pop<Image<output_t>* >(args);

			run(in, out, args);
		}

		using Distributable::runDistributed;

		virtual size_t getDistributionDirection2(const std::vector<ParamVariant>& args) const override
		{
			return 1;
		}

		virtual std::vector<string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const override
		{
			DistributedImage<input_t>& in = *std::get<DistributedImage<input_t>* >(args[0]);
			DistributedImage<output_t>& out = *std::get<DistributedImage<output_t>* >(args[1]);
			out.ensureSize(in);
			return distributor.distribute(this, args);
		}

		virtual JobType getJobType(const std::vector<ParamVariant>& args) const override
		{
			return JobType::Fast;
		}

		virtual bool canDelay(const std::vector<ParamVariant>& args) const override
		{
			return true;
		}
	};


	/**
	Base class for point processes that process image in-place and need a parameter image.
	*/
	template<typename input_t, typename param_t> class InputParamPointProcess : public InPlacePointProcess<input_t>
	{
	protected:
		friend class CommandList;

		InputParamPointProcess(const string& name, const string& help, const std::vector<CommandArgumentBase>& extraArgs = {}) :
			InPlacePointProcess<input_t>(name, help,
				concat({
					CommandArgument<Image<param_t> >(ParameterDirection::In, "parameter image", "Parameter image.")
					}, extraArgs)
				)
		{
		}

	public:
		virtual void run(Image<input_t>& in, std::vector<ParamVariant>& args) const override
		{
			Image<param_t>& out = *pop<Image<param_t>* >(args);
			run(in, out, args);
		}

		virtual void run(Image<input_t>& in, Image<param_t>& param, std::vector<ParamVariant>& args) const = 0;
	};





#define DEF_MATH_SINGLE_REAL2(classname, commandname, help, funcname) \
	template<typename pixel_t> class classname##Command : public InPlacePointProcess<pixel_t> \
	{ \
	protected: \
		friend class CommandList; \
	\
		classname##Command() : InPlacePointProcess<pixel_t>(#commandname, help, {}) {}\
	\
	public: \
		virtual void run(Image<pixel_t>& img, std::vector<ParamVariant>& args) const override \
		{ \
			itl2:: funcname (img); \
		} \
	}; \

#define DEF_MATH_SINGLE_REAL(classname, commandname, help) DEF_MATH_SINGLE_REAL2(classname, commandname, help, commandname)	

#define DEF_MATH_SINGLE_COMPLEX2(classname, commandname, help, funcname) \
	class classname##ComplexCommand : public InPlacePointProcess<complex32_t> \
	{ \
	protected: \
		friend class CommandList; \
	\
		classname##ComplexCommand() : InPlacePointProcess<complex32_t>(#commandname, help, {}) {} \
	\
	public: \
		virtual void run(Image<complex32_t>& img, std::vector<ParamVariant>& args) const override \
		{ \
			itl2:: funcname (img); \
		} \
	};

#define DEF_MATH_SINGLE_COMPLEX(classname, commandname, help) DEF_MATH_SINGLE_COMPLEX2(classname, commandname, help, commandname)

#define DEF_MATH_SINGLE_ALL2(classname, commandname, help, funcname) \
	DEF_MATH_SINGLE_REAL2(classname, commandname, help, funcname) \
	DEF_MATH_SINGLE_COMPLEX2(classname, commandname, help, funcname)

#define DEF_MATH_SINGLE_ALL(classname, commandname, help) \
	DEF_MATH_SINGLE_ALL2(classname, commandname, help, commandname)


	DEF_MATH_SINGLE_ALL2(SwapByteOrder, swapbyteorder, "Swaps byte order of each pixel value. Use this command to convert images read in wrong endianness to the correct one, or before saving an image if it should be saved in non-native byte order.", swapByteOrder)
	DEF_MATH_SINGLE_ALL(Negate, negate, "Negates pixel values.")
	DEF_MATH_SINGLE_ALL(Exponentiate, exponentiate, "Exponentates pixel values.")
	DEF_MATH_SINGLE_ALL(Abs, abs, "Calculates absolute value of pixel values.")
	DEF_MATH_SINGLE_ALL(Square, square, "Calculates square of pixel values.")
	DEF_MATH_SINGLE_ALL2(SquareRoot, squareroot, "Calculates square root of pixel values.", squareRoot)

	DEF_MATH_SINGLE_ALL(Log, log, "Calculates natural logarithm of pixel values.")
	DEF_MATH_SINGLE_ALL(Log10, log10, "Calculates base-10 logarithm of pixel values.")
	DEF_MATH_SINGLE_ALL(Sin, sin, "Calculates sine of pixel values.")
	DEF_MATH_SINGLE_ALL(Cos, cos, "Calculates cosine of pixel values.")
	DEF_MATH_SINGLE_ALL(Tan, tan, "Calculates tangent of pixel values.")
	DEF_MATH_SINGLE_ALL(Inv, inv, "Calculates inverse (1/x) of pixel values.")

	DEF_MATH_SINGLE_COMPLEX(Real, real, "Calculates real part of pixel values.")
	DEF_MATH_SINGLE_COMPLEX(Imag, imag, "Calculates imaginary part of pixel values.")
	DEF_MATH_SINGLE_COMPLEX(Conjugate, conjugate, "Calculates complex conjugate of pixel values.")
	DEF_MATH_SINGLE_COMPLEX(Normalize, normalize, "Makes norm of each pixel one.")
	DEF_MATH_SINGLE_COMPLEX(Arg, arg, "Calculates argument of pixel values.")
	DEF_MATH_SINGLE_COMPLEX2(NormSquared, normsquared, "Calculates squared norm of pixel values.", normSquared)

	DEF_MATH_SINGLE_REAL(Round, round, "Rounds pixel values.")
	DEF_MATH_SINGLE_REAL(Ceil, ceil, "Calculates ceiling of pixel values, i.e. the greatest integer less than or equal to the pixel value.")
	DEF_MATH_SINGLE_REAL(Floor, floor, "Calculates floor of pixel values, i.e. the least integer greater than or equal to the pixel value.")


	template<typename pixel_t> class SetConstantCommand : public InPlacePointProcess<pixel_t>
	{
	protected:
		friend class CommandList;

		SetConstantCommand() : InPlacePointProcess<pixel_t>("set", "Sets all pixels to the same value.",
			{
				CommandArgument<double>(ParameterDirection::In, "x", "Pixel value.")
			}) {}

	public:
		virtual void run(Image<pixel_t>& img, std::vector<ParamVariant>& args) const override
		{
			double param = pop<double>(args);
			itl2::setValue(img, param);
		}
	};


	template<typename pixel_t> class ThresholdRangeCommand : public InPlacePointProcess<pixel_t>
	{
	protected:
		friend class CommandList;

		ThresholdRangeCommand() :InPlacePointProcess<pixel_t>("thresholdrange", "Threshols a range from the image. Sets to one all pixels whose value is in ]min, max].",
			{
				CommandArgument<double>(ParameterDirection::In, "min", "Minimum of the threshold range."),
				CommandArgument<double>(ParameterDirection::In, "max", "Maximum of the threshold range.")
			})
		{

		}

	public:
		virtual void run(Image<pixel_t>& img, std::vector<ParamVariant>& args) const override
		{
			double min = pop<double>(args);
			double max = pop<double>(args);
			thresholdRange(img, Vec2d(min, max));
		}
	};

	template<typename pixel_t> class DoubleThresholdCommand : public InPlacePointProcess<pixel_t>
	{
	protected:
		friend class CommandList;

		DoubleThresholdCommand() : InPlacePointProcess<pixel_t>("doublethreshold", "Sets pixel to 0 if its value is less than the first threshold. Sets pixel to 1 if pixel its is larger than or equal to the first threshold and less than the second threshold (if any). Sets pixel to 2 if pixel its is larger than or equal to the second threshold.",
			{
				CommandArgument<double>(ParameterDirection::In, "first threshold", "Pixels whose value is larger than or equal to this threshold and less than the second threshold are set to 1."),
				CommandArgument<double>(ParameterDirection::In, "second threshold", "Pixels whose value is larger than or equal to this threshold are set to 2.")
			})
		{

		}

	public:
		virtual void run(Image<pixel_t>& img, std::vector<ParamVariant>& args) const override
		{
			double t1 = pop<double>(args);
			double t2 = pop<double>(args);
			std::vector<pixel_t> v = { pixelRound<pixel_t>(t1), pixelRound<pixel_t>(t2) };
			multiThreshold(img, v);
		}
	};

	template<typename pixel_t> class ThresholdPeriodicCommand : public InPlacePointProcess<pixel_t>
	{
	protected:
		friend class CommandList;

		ThresholdPeriodicCommand() : InPlacePointProcess<pixel_t>("thresholdperiodic", "Threshols a range from the image where pixel values are from periodic range (e.g. angle). Sets to one all pixels whose value is in ]min, max] mod period. If threshold range is larger than period, sets all pixels to 1. If threshold range start is larger than or equal to threshold range end, sets all pixels to 0.",
			{
				CommandArgument<double>(ParameterDirection::In, "period start", "Minimum of the period of pixel values."),
				CommandArgument<double>(ParameterDirection::In, "period end", "Maximum of the period of pixel values."),
				CommandArgument<double>(ParameterDirection::In, "min", "Minimum of the threshold range."),
				CommandArgument<double>(ParameterDirection::In, "max", "Maximum of the threshold range.")
			})
		{

		}

	public:
		virtual void run(Image<pixel_t>& img, std::vector<ParamVariant>& args) const override
		{
			double pmin = pop<double>(args);
			double pmax = pop<double>(args);
			double min = pop<double>(args);
			double max = pop<double>(args);
			thresholdPeriodic(img, Vec4d(pmin, pmax, min, max));
		}
	};





	template<typename pixel_t> class LinearMapCommand : public InPlacePointProcess<pixel_t>
	{
	protected:
		friend class CommandList;

		LinearMapCommand() : InPlacePointProcess<pixel_t>("linmap", "Maps pixel values linearly from one range to another. Maps position a in range [input min, input max] linearly to position in range [output min, output max]. E.g. if a = 0.5 and [input min, input max] = [0, 1] and [output min, output max] = [3, 4], result will be 3 + (0.5 - 0) / (1 - 0) * (4 - 3) = 3.5.",
			{
				CommandArgument<double>(ParameterDirection::In, "input min", "Minimum of the input range."),
				CommandArgument<double>(ParameterDirection::In, "input max", "Maximum of the input range."),
				CommandArgument<double>(ParameterDirection::In, "output min", "Minimum of the output range."),
				CommandArgument<double>(ParameterDirection::In, "output max", "Maximum of the output range.")
			})
		{

		}

	public:
		virtual void run(Image<pixel_t>& img, std::vector<ParamVariant>& args) const override
		{
			double imin = pop<double>(args);
			double imax = pop<double>(args);
			double omin = pop<double>(args);
			double omax = pop<double>(args);
			linearMap(img, Vec4d(imin, imax, omin, omax));
		}
	};

	template<typename pixel_t> class ReplaceCommand : public InPlacePointProcess<pixel_t>
	{
	protected:
		friend class CommandList;

		ReplaceCommand() :InPlacePointProcess<pixel_t>("replace", "Replaces pixel values a by value b.",
			{
				CommandArgument<double>(ParameterDirection::In, "a", "Value to be replaced by b."),
				CommandArgument<double>(ParameterDirection::In, "b", "Value that replaces a.")
			})
		{

		}

	public:
		virtual void run(Image<pixel_t>& img, std::vector<ParamVariant>& args) const override
		{
			double a = pop<double>(args);
			double b = pop<double>(args);
			replace(img, Vec2<pixel_t>(pixelRound<pixel_t>(a), pixelRound<pixel_t>(b)));
		}
	};















	template<typename input_t, typename output_t> class CopyCommand : public InputOutputPointProcess<input_t, output_t>
	{
	protected:
		friend class CommandList;

		CopyCommand() : InputOutputPointProcess<input_t, output_t>("copy", "Copies input image to output image. If pixel data types are not the same, performs conversion.", {})
		{
		}

	public:
		virtual void run(Image<input_t>& in, Image<output_t>& out, std::vector<ParamVariant>& args) const override
		{
			setValue(out, in);
		}
	};



	template<typename pixel_t> class SetEdgesCommand : public OneImageInPlaceCommand<pixel_t>
	{
	protected:
		friend class CommandList;

		SetEdgesCommand() :
			OneImageInPlaceCommand<pixel_t>("setedges", "Set edges of the image to specified value.",
				{
					CommandArgument<double>(ParameterDirection::In, "edge value", "The edges are set to this value.", 0),
					CommandArgument<coord_t>(ParameterDirection::In, "radius", "Pixels whose distance to the image edge is less than or equal to this value are set.", 1)
				}
				)
		{
		}

	public:
		virtual void run(Image<pixel_t>& img, std::vector<ParamVariant>& args) const override
		{
			double val = pop<double>(args);
			coord_t r = pop<coord_t>(args);

			setEdges(img, pixelRound<pixel_t>(val), r);
		}
	};


#define DEF_MATH_DUAL(classname, commandname, helpImageTopic, helpParamTopic, helpParamParam) \
	template<typename pixel_t, typename param_t> class classname##Command : public InputParamPointProcess<pixel_t, param_t> \
	{ \
	protected: \
		friend class CommandList; \
	\
		classname##Command() : InputParamPointProcess<pixel_t, param_t>(#commandname, helpImageTopic, \
		{ \
			CommandArgument<bool>(ParameterDirection::In, "allow broadcast", "Set to true to allow size of parameter image differ from size of input image. If there is a need to access pixel outside of parameter image, the nearest value inside the image is taken instead. If set to false, dimensions of input and parameter images must be equal. If set to true, the parameter image is always loaded in its entirety in distributed processing mode.", false)	\
		}) {} \
	public: \
		virtual void run(Image<pixel_t>& a, Image<param_t>& b, std::vector<ParamVariant>& args) const override \
		{ \
			bool allowBroadcast = pop<bool>(args);	\
			itl2:: commandname (a, b, allowBroadcast); \
		} \
		\
		virtual void getCorrespondingBlock(const std::vector<ParamVariant>& args, size_t argIndex, Vec3c& readStart, Vec3c& readSize, Vec3c& writeFilePos, Vec3c& writeImPos, Vec3c& writeSize) const override	\
		{	\
			bool allowBroadcast = std::get<bool>(args[2]);	\
			if(allowBroadcast && argIndex == 1)			\
			{	\
				DistributedImage<param_t>& img2 = *std::get<DistributedImage<param_t>*>(args[1]);	\
				readStart = Vec3c(0, 0, 0);		\
				readSize = img2.dimensions();	\
			}	\
		}	\
	};\
	template<typename pixel_t> class classname##ConstantCommand : public InPlacePointProcess<pixel_t> \
	{ \
	protected: \
		friend class CommandList; \
	\
		classname##ConstantCommand() : InPlacePointProcess<pixel_t>(#commandname, helpParamTopic, \
		{ \
			CommandArgument<double>(ParameterDirection::In, "x", helpParamParam) \
		}) {} \
	public: \
		virtual void run(Image<pixel_t>& img, std::vector<ParamVariant>& args) const override \
		{ \
			double param = pop<double>(args); \
			itl2:: commandname (img, param); \
		} \
	};

	DEF_MATH_DUAL(Add, add, "Adds two images. Output is placed to the first image. The operation is performed using saturation arithmetic.", "Adds and image and a constant. The operation is performed using saturation arithmetic.", "Constant to add to the image.")
	DEF_MATH_DUAL(Subtract, subtract, "Subtracts two images. Output is placed to the first image. The operation is performed using saturation arithmetic.", "Subtracts a constant from the image. The operation is performed using saturation arithmetic.", "Constant to subtract from the image.")
	DEF_MATH_DUAL(InvSubtract, invsubtract, "Subtracts two images in inverse order, and places output to the first image. In other words, first image = second image - first image. The operation is performed using saturation arithmetic.", "Subtracts the first image from a constant, and places the result to the first image. The operation is performed using saturation arithmetic.", "Constant from which the image is subtracted.")
	DEF_MATH_DUAL(Multiply, multiply, "Multiplies two images. Output is placed to the first image. The operation is performed using saturation arithmetic.", "Multiplies an image by a constant. The operation is performed using saturation arithmetic.", "Multiplier constant.")
	DEF_MATH_DUAL(Divide, divide, "Divides two images. Output is placed to the first image. The operation is performed using saturation arithmetic.", "Divides an image by a constant. The operation is performed using saturation arithmetic.", "Divider constant.")
	DEF_MATH_DUAL(Max, max, "Calculates pixelwise maximum of two images. Output is placed to the first image.", "Calculates maximum of an image and a constant.", "Constant value.")
	DEF_MATH_DUAL(Min, min, "Calculates pixelwise minimum of two images. Output is placed to the first image.", "Calculates minimum of an image and a constant.", "Constant value.")
	//DEF_MATH_DUAL(SetValue, setValue, "Copies pixel values from right image to left image.", "Sets all pixels to the same value.", "Pixel value.")
	DEF_MATH_DUAL(Threshold, threshold, "Thresholds left image, taking threshold values from right image. Sets pixel to 1 if pixel value > threshold and to 0 otherwise.", "Thresholds image. Sets pixel to 1 if pixel value > threshold and to 0 otherwise.", "Threshold value.")


	template<typename pixel_t> class SetCommand : public Command, public Distributable
	{
	protected:
		friend class CommandList;

		SetCommand() : Command("set", "Copies pixel values from the source image to the target image. Sets the size and data type of the target image to those of the source image. See also `copy` command.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::Out, "target image", "Image whose values are set."),
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "source image", "Image that is copied to the target image.")
			})
		{
		}

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& target = *pop<Image<pixel_t>* >(args);
			Image<pixel_t>& source = *pop<Image<pixel_t>* >(args);
			itl2::setValue(target, source);
		}

		virtual std::vector<string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const override
		{
			DistributedImage<pixel_t>& target = *std::get<DistributedImage<pixel_t>* >(args[0]);
			DistributedImage<pixel_t>& source = *std::get<DistributedImage<pixel_t>* >(args[1]);
			target.ensureSize(source.dimensions());

			// distribute in z, no overlap
			return distributor.distribute(this, args);
		}
	};

}

