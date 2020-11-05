#pragma once

#include "commandsbase.h"
#include "tomo/fbp.h"

namespace pilib
{

	class FBPPreprocessCommand : public TwoImageInputOutputCommand<float32_t>
	{
	protected:
		friend class CommandList;

		FBPPreprocessCommand() : TwoImageInputOutputCommand<float32_t>("fbppreprocess", "Performs preprocessing of transmission projection data for filtered backprojection. This command is experimental and may change in the near future.",
			{
				CommandArgument<std::string>(ParameterDirection::In, "reconstruction settings", "Settings for the reconstruction. If this string contains only a name of an existing file, the settings are read from that file. Otherwise, the string is treated as contents of the settings file.", "")
			},
			"fbp")
		{
		}

	public:
		virtual void run(Image<float32_t>& in, Image<float32_t>& out, std::vector<ParamVariant>& args) const override
		{
			std::string settings = pop<std::string>(args);

			if (fileExists(settings))
			{
				settings = readText(settings, true);
			}
			
			RecSettings sets = fromString<RecSettings>(settings);

			fbpPreprocess(in, out, sets);
		}
	};


	class FBPCommand : public TwoImageInputOutputCommand<float32_t>
	{
	protected:
		friend class CommandList;

		FBPCommand() : TwoImageInputOutputCommand<float32_t>("fbp", "Performs filtered backprojection of data for which fbppreprocess has been called. This command is experimental and may change in the near future.",
			{
				CommandArgument<std::string>(ParameterDirection::In, "reconstruction settings", "Settings for the reconstruction. If this string contains only a name of an existing file, the settings are read from that file. Otherwise, the string is treated as contents of the settings file.", ""),
				CommandArgument<bool>(ParameterDirection::In, "use GPU", "Set to true to allow processing on a GPU.", true)
			},
			"fbppreprocess")
		{
		}

	public:
		virtual void run(Image<float32_t>& in, Image<float32_t>& out, std::vector<ParamVariant>& args) const override
		{
			std::string settings = pop<std::string>(args);
			bool useGPU = pop<bool>(args);
			

			if (fileExists(settings))
			{
				settings = readText(settings, true);
			}

			RecSettings sets = fromString<RecSettings>(settings);

			if (useGPU)
			{
#if defined(USE_OPENCL)
				backprojectOpenCL(in, sets, out);
#else
				backproject(in, sets, out);
#endif
			}
			else
			{
				backproject(in, sets, out);
			}
		}
	};

	class CreateFBPFilterCommand : public Command
	{
	protected:
		friend class CommandList;

		CreateFBPFilterCommand() :
			Command("createfbpfilter", "Creates image containing values of filter function used in filtered backprojection. Can be used to visualize different filter functions.",
				{
					CommandArgument<Image<float32_t> >(ParameterDirection::Out, "image", "Image where the filter is to be placed."),
					CommandArgument<size_t>(ParameterDirection::In, "size", "Size of the filter. This corresponds to the padded image size in FBP.", 100),
					CommandArgument<std::string>(ParameterDirection::In, "filter type", "Type of the filter. Supported values are Ideal ramp, Ramp, Shepp-Logan, Cosine, Hamming, Hann, Blackman, Parze.", "Ramp"),
					CommandArgument<double>(ParameterDirection::In, "cut-off", "Filter cut-off frequency, 0 corresponds to DC and 1 corresponds to Nyquist.", 1.0f)
				}
			)
		{
		}

	public:
		virtual void run(std::vector<ParamVariant>& args) const override
		{
			Image<float32_t>& img = *pop<Image<float32_t>* >(args);
			size_t size = pop<size_t>(args);
			std::string s = pop<std::string>(args);
			double cutoff = pop<double>(args);

			FilterType t = fromString<FilterType>(s);

			if (size <= 0)
				size = img.width();

			img.ensureSize(size, 1, 1);

			createFilter(img, t, (float32_t)cutoff);
		}
	};
}
