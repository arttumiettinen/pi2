#pragma once

#include "command.h"
#include "distributable.h"
#include "io/imagedatatype.h"

namespace pilib
{

	namespace internals
	{
		/**
		Functor that creates a new image and stores it to images vector.
		*/
		template<typename pixel_t> struct AssignImage
		{
			static void run(ParamVariant& param, ImageBase* img)
			{
				param = (Image<pixel_t>*)img;
			}
		};
	}

	/**
	Identifies command to be trivially distributable, i.e. it runs the same way regardless of distributed processing state.
	Used for example for info and help commands.
	*/
	class TrivialDistributable : public Distributable, virtual public Command
	{
	protected:
		friend class CommandList;

		TrivialDistributable()
		{

		}

	public:
		using Distributable::runDistributed;

		virtual std::vector<std::string> runDistributed(Distributor& distributor, std::vector<ParamVariant>& args) const override
		{
			// Convert DistributedImage to Image.
			const std::vector<CommandArgumentBase>& argdefs = this->args();
			
			std::vector<ParamVariant> modifiedArgs, modifiedArgs2;

			std::vector<std::unique_ptr<ImageBase>> images;
			for (size_t n = 0; n < argdefs.size(); n++)
			{
				const CommandArgumentBase& argdef = argdefs[n];
				if (isImage(argdef.dataType()))
				{
					DistributedImageBase* dimage = getDistributedImage(args[n]);
					std::unique_ptr<ImageBase> nimage = dimage->toNormalImage();

					ParamVariant param;
					itl2::pick<internals::AssignImage>(dimage->dataType(), param, nimage.get());
					
					modifiedArgs.push_back(param);
					modifiedArgs2.push_back(param);

					images.push_back(std::move(nimage));
				}
				else
				{
					modifiedArgs.push_back(args[n]);
					modifiedArgs2.push_back(args[n]);
				}
			}

			// Run the command.
			runInternal(distributor.getSystem(), modifiedArgs2); // Note: this line usually clears argument list given as parameter.

			// Put data from Images to the DistributedImages.
			for (size_t n = 0; n < argdefs.size(); n++)
			{
				const CommandArgumentBase& argdef = argdefs[n];
				if (isImage(argdef.dataType()))
				{
					DistributedImageBase* dimage = getDistributedImage(args[n]);
					ImageBase* nimage = getImage(modifiedArgs[n]);
					dimage->setData(nimage);
				}
				else
				{
					args[n] = modifiedArgs[n];
				}
			}

			return std::vector<std::string>();
		}
	};
}
