#pragma once

#include "command.h"
#include "commandsbase.h"
#include "distributable.h"

#include "eval.h"

namespace pilib
{
	class EvalBase : public Command, public Distributable
	{
	protected:
		EvalBase(const std::vector<CommandArgumentBase>& imageArgs) : Command("eval", "Evaluates a mathematical expression on each pixel of one or more images. The expression is given as a string.",
			concat({
				CommandArgument<string>(ParameterDirection::In, "expression", "The expression to evaluate on each pixel."),
			}, imageArgs)
		)
		{
		}
	public:

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

	template<typename output_t> class Eval1Command : public EvalBase
	{
	protected:
		friend class CommandList;

		Eval1Command() :
			EvalBase(
				{
					CommandArgument<Image<output_t> >(ParameterDirection::InOut, "image", "Image to process. The result will be assigned to this image. This image can be referenced by name x0 in the expression.")
				}
			)
		{
		}

	public:

		virtual void run(std::vector<ParamVariant>& args) const override
		{
			string expression = pop<string>(args);
			Image<output_t>& out = *pop<Image<output_t>* >(args);
			itl2::eval(expression, out, std::vector<ImageBase*>{&out});
		}

		virtual Vec3c getMargin(const std::vector<ParamVariant>& args) const
		{
			// Check sizes of argument images.
			return Vec3c(0, 0, 0);
		}
	};

	template<typename output_t, typename param1_t> class Eval2Command : public EvalBase
	{
	protected:
		friend class CommandList;

		Eval2Command() :
			EvalBase(
				{
					CommandArgument<Image<output_t> >(ParameterDirection::InOut, "image", "Image to process. The result will be assigned to this image. Pixel of this image can be referenced by name x0 in the expression."),
					CommandArgument<Image<param1_t> >(ParameterDirection::In, "argument image", "Argument image. Pixel of this image can be referenced by name x1 in the expression."),
				}
			)
		{
		}

	public:

		virtual void run(std::vector<ParamVariant>& args) const override
		{
			string expression = pop<string>(args);
			Image<output_t>& out = *pop<Image<output_t>* >(args);
			Image<param1_t>& param1 = *pop<Image<param1_t>* >(args);
			itl2::eval(expression, out, std::vector<ImageBase*>{&out, &param1});
		}

		virtual Vec3c getMargin(const std::vector<ParamVariant>& args) const
		{
			// Check sizes of argument images.
			DistributedImage<output_t>& out = *std::get<DistributedImage<output_t>* >(args[1]);
			DistributedImage<param1_t>& param1 = *std::get<DistributedImage<param1_t>* >(args[2]);
			out.checkSize(param1.dimensions());

			return Vec3c(0, 0, 0);
		}
	};


	template<typename output_t, typename param1_t, typename param2_t> class Eval3Command : public EvalBase
	{
	protected:
		friend class CommandList;

		Eval3Command() :
			EvalBase(
				{
					CommandArgument<Image<output_t> >(ParameterDirection::InOut, "image", "Image to process. The result will be assigned to this image. Pixel of this image can be referenced by name x0 in the expression."),
					CommandArgument<Image<param1_t> >(ParameterDirection::In, "argument image", "Argument image. Pixel of this image can be referenced by name x1 in the expression."),
					CommandArgument<Image<param2_t> >(ParameterDirection::In, "argument image 2", "Argument image. Pixel of this image can be referenced by name x2 in the expression."),
				}
				)
		{
		}

	public:

		virtual void run(std::vector<ParamVariant>& args) const override
		{
			string expression = pop<string>(args);
			Image<output_t>& out = *pop<Image<output_t>* >(args);
			Image<param1_t>& param1 = *pop<Image<param1_t>* >(args);
			Image<param2_t>& param2 = *pop<Image<param2_t>* >(args);
			itl2::eval(expression, out, std::vector<ImageBase*>{&out, &param1, &param2});
		}

		virtual Vec3c getMargin(const std::vector<ParamVariant>& args) const
		{
			// Check sizes of argument images.
			DistributedImage<output_t>& out = *std::get<DistributedImage<output_t>* >(args[1]);
			DistributedImage<param1_t>& param1 = *std::get<DistributedImage<param1_t>* >(args[2]);
			DistributedImage<param2_t>& param2 = *std::get<DistributedImage<param2_t>* >(args[3]);
			out.checkSize(param1.dimensions());
			out.checkSize(param2.dimensions());

			return Vec3c(0, 0, 0);
		}
	};

}
