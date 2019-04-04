#pragma once

#include "command.h"
#include "commandsbase.h"
#include "distributable.h"
#include "pointprocess.h"
#include "math/mathutils.h"
#include "misc.h"

#include <vector>

using namespace std;
using namespace itl2;

namespace pilib
{
	/**
	Base class for processes that can be distributed like point process.
	@param CMDBASE The parent class that is derived from Command.
	*/
	template<class CMDBASE>
	class PointProcessDistributable : public CMDBASE, public Distributable
	{
	public:

		PointProcessDistributable(const string& name, const string& help, const vector<CommandArgumentBase>& extraArgs) : CMDBASE(name, help, extraArgs)
		{
		}

		using Distributable::runDistributed;

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
		{
			// distribute in z, no overlap
			return distributor.distribute(this, args, 2, Vec3c(0, 0, 0));
		}

		virtual JobType getJobType() const
		{
			return JobType::Fast;
		}
	};



#define DEF_MATH_SINGLE_REAL2(classname, commandname, help, funcname) \
	template<typename pixel_t> class classname##Command : public PointProcessDistributable<OneImageInPlaceCommand<pixel_t> > \
	{ \
	public: \
		classname##Command() : PointProcessDistributable<OneImageInPlaceCommand<pixel_t> >(#commandname, help, {}) {}\
	\
		virtual void run(Image<pixel_t>& img, vector<ParamVariant>& args) const \
		{ \
			itl2:: funcname (img); \
		} \
	}; \

#define DEF_MATH_SINGLE_REAL(classname, commandname, help) DEF_MATH_SINGLE_REAL2(classname, commandname, help, commandname)	

#define DEF_MATH_SINGLE_COMPLEX2(classname, commandname, help, funcname) \
	class classname##ComplexCommand : public PointProcessDistributable<OneImageInPlaceCommand<complex32_t> > \
	{ \
	public: \
		classname##ComplexCommand() : PointProcessDistributable<OneImageInPlaceCommand<complex32_t> >(#commandname, help, {}) {}\
	\
		virtual void run(Image<complex32_t>& img, vector<ParamVariant>& args) const \
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
	DEF_MATH_SINGLE_REAL(Ceil, ceil, "Calculates floor of pixel values.")
	DEF_MATH_SINGLE_REAL(Floor, floor, "Calculates ceiling of pixel values.")


	template<typename pixel_t> class SetConstantCommand : public PointProcessDistributable<OneImageInPlaceCommand<pixel_t> >
	{
	public:
		SetConstantCommand() : PointProcessDistributable<OneImageInPlaceCommand<pixel_t> >("set", "Sets all pixels to the same value.",
			{
				CommandArgument<double>(ParameterDirection::In, "x", "Pixel value.")
			}) {}

		virtual void run(Image<pixel_t>& img, vector<ParamVariant>& args) const
		{
			double param = pop<double>(args);
			itl2::setValue(img, param);
		}
	};


	template<typename pixel_t> class ThresholdRangeCommand : public PointProcessDistributable<OneImageInPlaceCommand<pixel_t> >
	{
	public:
		ThresholdRangeCommand() :PointProcessDistributable<OneImageInPlaceCommand<pixel_t> >("thresholdrange", "Threshols a range from the image. Sets to one all pixels whose value is in ]min, max].",
			{
				CommandArgument<double>(ParameterDirection::In, "min", "Minimum of the threshold range."),
				CommandArgument<double>(ParameterDirection::In, "max", "Maximum of the threshold range.")
			})
		{

		}

		virtual void run(Image<pixel_t>& img, vector<ParamVariant>& args) const
		{
			double min = pop<double>(args);
			double max = pop<double>(args);
			thresholdRange(img, Vec2d(min, max));
		}
	};

	template<typename pixel_t> class DoubleThresholdCommand : public PointProcessDistributable<OneImageInPlaceCommand<pixel_t> >
	{
	public:
		DoubleThresholdCommand() : PointProcessDistributable<OneImageInPlaceCommand<pixel_t> >("doublethreshold", "Sets pixel to 0 if its value is less than the first threshold. Sets pixel to 1 if pixel its is larger than or equal to the first threshold and less than the second threshold (if any). Sets pixel to 2 if pixel its is larger than or equal to the second threshold.",
			{
				CommandArgument<double>(ParameterDirection::In, "first threshold", "Pixels whose value is larger than or equal to this threshold and less than the second threshold are set to 1."),
				CommandArgument<double>(ParameterDirection::In, "second threshold", "Pixels whose value is larger than or equal to this threshold are set to 2.")
			})
		{

		}

		virtual void run(Image<pixel_t>& img, vector<ParamVariant>& args) const
		{
			double t1 = pop<double>(args);
			double t2 = pop<double>(args);
			vector<pixel_t> v = { pixelRound<pixel_t>(t1), pixelRound<pixel_t>(t2) };
			multiThreshold(img, v);
		}
	};

	template<typename pixel_t> class ThresholdPeriodicCommand : public PointProcessDistributable<OneImageInPlaceCommand<pixel_t> >
	{
	public:
		ThresholdPeriodicCommand() : PointProcessDistributable<OneImageInPlaceCommand<pixel_t> >("thresholdperiodic", "Threshols a range from the image where pixel values are from periodic range (e.g. angle). Sets to one all pixels whose value is in ]min, max] mod period. If threshold range is larger than period, sets all pixels to 1. If threshold range start is larger than or equal to threshold range end, sets all pixels to 0.",
			{
				CommandArgument<double>(ParameterDirection::In, "period start", "Minimum of the period of pixel values."),
				CommandArgument<double>(ParameterDirection::In, "period end", "Maximum of the period of pixel values."),
				CommandArgument<double>(ParameterDirection::In, "min", "Minimum of the threshold range."),
				CommandArgument<double>(ParameterDirection::In, "max", "Maximum of the threshold range.")
			})
		{

		}

		virtual void run(Image<pixel_t>& img, vector<ParamVariant>& args) const
		{
			double pmin = pop<double>(args);
			double pmax = pop<double>(args);
			double min = pop<double>(args);
			double max = pop<double>(args);
			thresholdPeriodic(img, Vec4d(pmin, pmax, min, max));
		}
	};

	template<typename pixel_t> class LinearMapCommand : public PointProcessDistributable<OneImageInPlaceCommand<pixel_t> >
	{
	public:
		LinearMapCommand() : PointProcessDistributable<OneImageInPlaceCommand<pixel_t> >("linmap", "Maps pixel values linearly from one range to another. Maps position a in range [input min, input max] linearly to position in range [output min, output max]. E.g. if a = 0.5 and [input min, input max] = [0, 1] and [output min, output max] = [3, 4], result will be 3 + (0.5 - 0) / (1 - 0) * (4 - 3) = 3.5.",
			{
				CommandArgument<double>(ParameterDirection::In, "input min", "Minimum of the input range."),
				CommandArgument<double>(ParameterDirection::In, "input max", "Maximum of the input range."),
				CommandArgument<double>(ParameterDirection::In, "output min", "Minimum of the output range."),
				CommandArgument<double>(ParameterDirection::In, "output max", "Maximum of the output range.")
			})
		{

		}

		virtual void run(Image<pixel_t>& img, vector<ParamVariant>& args) const
		{
			double imin = pop<double>(args);
			double imax = pop<double>(args);
			double omin = pop<double>(args);
			double omax = pop<double>(args);
			linearMap(img, Vec4d(imin, imax, omin, omax));
		}
	};

	template<typename pixel_t> class ReplaceCommand : public PointProcessDistributable<OneImageInPlaceCommand<pixel_t> >
	{
	public:
		ReplaceCommand() :PointProcessDistributable<OneImageInPlaceCommand<pixel_t> >("replace", "Replaces pixel values a by value b.",
			{
				CommandArgument<double>(ParameterDirection::In, "a", "Value to be replaced by b."),
				CommandArgument<double>(ParameterDirection::In, "b", "Value that replaces a.")
			})
		{

		}

		virtual void run(Image<pixel_t>& img, vector<ParamVariant>& args) const
		{
			double a = pop<double>(args);
			double b = pop<double>(args);
			replace(img, Vec2<pixel_t>(pixelRound<pixel_t>(a), pixelRound<pixel_t>(b)));
		}
	};















	template<typename input_t, typename output_t> class CopyCommand : public PointProcessDistributable<TwoImageInputOutputCommand<input_t, output_t> >
	{
	public:
		CopyCommand() : PointProcessDistributable<TwoImageInputOutputCommand<input_t, output_t> >("copy", "Copies input image to output image.", {})
		{
		}

		virtual void run(Image<input_t>& in, Image<output_t>& out, vector<ParamVariant>& args) const
		{
			setValue(out, in);
		}
	};



	template<typename pixel_t> class SetEdgesCommand : public OneImageInPlaceCommand<pixel_t>
	{
	public:
		SetEdgesCommand() :
			OneImageInPlaceCommand<pixel_t>("setedges", "Set edges of the image to specified value.",
				{
					CommandArgument<double>(ParameterDirection::In, "edge value", "The edges are set to this value.", 0),
					CommandArgument<coord_t>(ParameterDirection::In, "radius", "Pixels whose distance to the image edge is less than or equal to this value are set.", 1)
				}
				)
		{
		}

		virtual void run(Image<pixel_t>& img, vector<ParamVariant>& args) const
		{
			double val = pop<double>(args);
			coord_t r = pop<coord_t>(args);

			setEdges(img, math::pixelRound<pixel_t>(val), r);
		}
	};


#define DEF_MATH_DUAL(classname, commandname, helpImageTopic, helpParamTopic, helpParamParam) \
	template<typename pixel_t, typename param_t> class classname##Command : public PointProcessDistributable<TwoImageInputParamCommand<pixel_t, param_t> > \
	{ \
	public: \
		classname##Command() : PointProcessDistributable<TwoImageInputParamCommand<pixel_t, param_t> >(#commandname, helpImageTopic, {}) {}\
	\
		virtual void run(Image<pixel_t>& a, Image<param_t>& b, vector<ParamVariant>& args) const \
		{ \
			itl2:: commandname (a, b); \
		} \
	};\
	template<typename pixel_t> class classname##ConstantCommand : public PointProcessDistributable<OneImageInPlaceCommand<pixel_t> > \
	{ \
	public: \
		classname##ConstantCommand() : PointProcessDistributable<OneImageInPlaceCommand<pixel_t> >(#commandname, helpParamTopic, \
		{ \
			CommandArgument<double>(ParameterDirection::In, "x", helpParamParam) \
		}) {}\
	\
		virtual void run(Image<pixel_t>& img, vector<ParamVariant>& args) const \
		{ \
			double param = pop<double>(args); \
			itl2:: commandname (img, param); \
		} \
	};

	DEF_MATH_DUAL(Add, add, "Adds two images. Output is placed to the first image.", "Adds image and constant", "Constant to add.")
	DEF_MATH_DUAL(Subtract, subtract, "Subtracts two images. Output is placed to the first image.", "Subtracts constant from image.", "Constant to subtract from the image.")
	DEF_MATH_DUAL(Multiply, multiply, "Multiplies two images. Output is placed to the first image. ", "Multiplies image by constant.", "Multiplier constant.")
	DEF_MATH_DUAL(Divide, divide, "Divides two images. Output is placed to the first image.", "Divides image by a constant.", "Divider constant.")
	DEF_MATH_DUAL(Max, max, "Calculates pixelwise maximum of two images. Output is placed to the first image.", "Calculates maximum of an image and a constant.", "Constant value.")
	DEF_MATH_DUAL(Min, min, "Calculates pixelwise minimum of two images. Output is placed to the first image.", "Calculates minimum of an image and a constant.", "Constant value.")
	//DEF_MATH_DUAL(SetValue, setValue, "Copies pixel values from right image to left image.", "Sets all pixels to the same value.", "Pixel value.")
	DEF_MATH_DUAL(Threshold, threshold, "Thresholds left image, taking threshold values from right image. Sets pixel to 1 if pixel value > threshold and to 0 otherwise.", "Thresholds image. Sets pixel to 1 if pixel value > threshold and to 0 otherwise.", "Threshold value.")


	template<typename pixel_t> class SetCommand : public Command, public Distributable
	{
	public:
		SetCommand() : Command("set", "Copies pixel values from right image to left image.",
			{
				CommandArgument<Image<pixel_t> >(ParameterDirection::Out, "target image", "Image whose values are set."),
				CommandArgument<Image<pixel_t> >(ParameterDirection::In, "source image", "Image that is copied to the target image.")
			})
		{
		}

		virtual void run(vector<ParamVariant>& args) const
		{
			Image<pixel_t>& target = *pop<Image<pixel_t>* >(args);
			Image<pixel_t>& source = *pop<Image<pixel_t>* >(args);
			itl2::setValue(target, source);
		}

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const
		{
			DistributedImage<pixel_t>& target = *get<DistributedImage<pixel_t>* >(args[0]);
			DistributedImage<pixel_t>& source = *get<DistributedImage<pixel_t>* >(args[1]);
			target.ensureSize(source.dimensions());

			// distribute in z, no overlap
			return distributor.distribute(this, args, 2, Vec3c(0, 0, 0));
		}
	};

	


}

