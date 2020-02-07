#pragma once

#include "pisystem.h"
#include "conversions.h"
#include "pick.h"
#include "pointprocesscommands.h"

namespace pilib
{
	template<typename in_t> class ConvertCommandBase : public Command, public Distributable
	{
	private:
		template<typename out_t> struct CreateAndCopy
		{
			static void run(const string& outname, Image<in_t>* in, PISystem* system)
			{
				// Free memory if old output image is not in use anymore
				system->replaceImage(outname, nullptr);

				std::shared_ptr<ImageBase> ptr = std::make_shared<Image<out_t> >(in->dimensions());
				system->replaceImage(outname, ptr);

				convert<in_t, out_t>(*in, *(Image<out_t>*)ptr.get());
			}
		};

		template<typename out_t> struct CreateAndCopyDistributed
		{
			static void run(const string& outname, DistributedImage<in_t>* in, Distributor& distributor)
			{
				std::shared_ptr<DistributedImageBase> ptr = std::make_shared<DistributedImage<out_t> >(distributor, outname, in->dimensions());
				distributor.getSystem()->replaceDistributedImage(outname, ptr);

				auto& cmd = CommandList::get<CopyCommand<in_t, out_t> >();
				cmd.runDistributed(distributor, { in, (DistributedImage<out_t>*)ptr.get() });
			}
		};

    protected:
    
        // TODO: These wrappers are not really needed but I can't find any way to use template defined in the base class
        // as a template template argument for pick function, if the call is made from the derived class.
    
        void runCreateAndCopy(const string& dts, const string& outname, Image<in_t>* in, PISystem* system) const
        {
            pick<CreateAndCopy>(dts, outname, in, system);
        }
        
        void runCreateAndCopyDistributed(const string& dts, const string& outname, DistributedImage<in_t>* in, Distributor& distributor) const
        {
            pick<CreateAndCopyDistributed>(dts, outname, in, distributor);
        }

	public:
		using Command::Command;
	};
	

	/**
	Convert command.
	*/
	template<typename in_t> class ConvertCommand : public ConvertCommandBase<in_t>
	{
	protected:
		friend class CommandList;

		ConvertCommand() : ConvertCommandBase<in_t>("convert", "Converts data type of input image.",
			{
				CommandArgument<Image<in_t> >(ParameterDirection::In, "input image", "Input image."),
				CommandArgument<string>(ParameterDirection::In, "output image", "Output image."),
				CommandArgument<string>(ParameterDirection::In, "data type", "Data type of the output image. Can be " + listSupportedImageDataTypes() + ".")
			})
		{
		}

	public:
		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const override
		{
			Image<in_t>* in = pop<Image<in_t>*>(args);
			string outname = pop<string>(args);
			string dts = pop<string>(args);

            ConvertCommandBase<in_t>::runCreateAndCopy(dts, outname, in, system);

		}

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override
		{
			DistributedImage<in_t>* in = std::get<DistributedImage<in_t>*>(args[0]);
			string outname = std::get<string>(args[1]);
			string dts = std::get<string>(args[2]);

			ConvertCommandBase<in_t>::runCreateAndCopyDistributed(dts, outname, in, distributor);

			return vector<string>();
		}

		virtual void run(vector<ParamVariant>& args) const
		{
		}
	};

	/**
	Convert command.
	*/
	template<typename in_t> class ConvertInPlaceCommand : public ConvertCommandBase<in_t>
	{
	protected:
		friend class CommandList;

		ConvertInPlaceCommand() : ConvertCommandBase<in_t>("convert", "Converts data type of an image.",
			{
				CommandArgument<Image<in_t> >(ParameterDirection::InOut, "image", "Image whose data type is to be converted."),
				CommandArgument<string>(ParameterDirection::In, "data type", "Data type to convert to. Can be " + listSupportedImageDataTypes() + ".")
			})
		{
		}

	public:
		virtual void runInternal(PISystem* system, vector<ParamVariant>& args) const override
		{
			Image<in_t>* in = pop<Image<in_t>*>(args);
			string dts = pop<string>(args);

			ConvertCommandBase<in_t>::runCreateAndCopy(dts, system->imageName(in), in, system);
		}

		virtual vector<string> runDistributed(Distributor& distributor, vector<ParamVariant>& args) const override
		{
			DistributedImage<in_t>* in = std::get<DistributedImage<in_t>*>(args[0]);
			string dts = std::get<string>(args[1]);

			ConvertCommandBase<in_t>::runCreateAndCopyDistributed(dts, in->varName(), in, distributor);

			return vector<string>();
		}

		virtual void run(vector<ParamVariant>& args) const
		{
		}
	};

}
