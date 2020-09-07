#pragma once

#include "commandsbase.h"

#include "pathopening.h"

#include "pisystem.h"

namespace pilib
{

	template<typename pixel_t> class PathLengthCommand : public TwoImageInputOutputCommand<pixel_t, float32_t>
	{
	protected:
		friend class CommandList;

		PathLengthCommand() : TwoImageInputOutputCommand<pixel_t, float32_t>("pathlength", "Replaces value of each pixel by the length of the longest constrained path that goes through that pixel. Works only with binary input images. NOTE: This command creates temporary files to the current directory.",
			{
			},
			"")
		{

		}

	public:

		virtual void runInternal(PISystem* system, std::vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& in = *pop<Image<pixel_t>* >(args);
			Image<float32_t>& out = *pop<Image<float32_t>* >(args);

			std::string inName = system->imageName(&in);

			pathLength2Binary3dNormalOrChamferMemorySave(in, out, LengthType::Ones, std::string("TEMP-") + inName);
		}

		virtual void run(Image<pixel_t>& in, Image<float32_t>& out, std::vector<ParamVariant>& args) const override
		{
			
		}
	};

}
