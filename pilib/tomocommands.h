#pragma once

#include "commandsbase.h"
#include "tomo/fbp.h"

namespace pilib
{
	inline std::string fbpSeeAlso()
	{
		return "fbppreprocess, fbp, deadpixelremoval, createfbpfilter, defaultrecsettings";
	}

	class DefaultRecSettingsCommand : public OneImageInPlaceCommand<float32_t>
	{
	protected:
		friend class CommandList;

		DefaultRecSettingsCommand() : OneImageInPlaceCommand<float32_t>("defaultrecsettings", "Writes default reconstruction settings to image metadata.",
			{ },
			fbpSeeAlso())
		{

		}
	public:
		virtual void run(Image<float32_t>& in, std::vector<ParamVariant>& args) const override
		{
			RecSettings s;
			s.toMeta(in.metadata);
		}
	};

	class FBPPreprocessCommand : public TwoImageInputOutputCommand<float32_t>
	{
	protected:
		friend class CommandList;

		FBPPreprocessCommand() : TwoImageInputOutputCommand<float32_t>("fbppreprocess", "Performs preprocessing of transmission projection data for filtered backprojection. Reconstruction settings are read from the metadata of the input image. This command is experimental and may change in the near future.",
			{ },
			fbpSeeAlso())
		{
		}

	public:
		virtual void run(Image<float32_t>& in, Image<float32_t>& out, std::vector<ParamVariant>& args) const override
		{
			RecSettings sets = RecSettings::fromMeta(in.metadata);

			fbpPreprocess(in, out, sets);
		}
	};


	template<typename pixel_t> class DeadPixelRemovalCommand : public OneImageInPlaceCommand<pixel_t>
	{
	protected:
		friend class CommandList;

		DeadPixelRemovalCommand() : OneImageInPlaceCommand<pixel_t>("deadpixelremoval", "Removes dead pixels from projection dataset."
"Determines whether dead pixel removal algorithm should be applied to the projection images. "
"In the algorithm, each flat-field corrected projection $I$ is processed separately. Pixel at "
"position $x$ is classified as dead if its value $I(x)$ is $NaN$ or it satisfies "
"$|I(x) - m(x)| > M * std(|I - m|)$, "
"where "
"$m$ is a median filtering of $I$ with user-specified radius, "
"$M$ is a magnitude parameter, and "
"$std$ is standard deviation of whole image. "
"If a pixel is dead, its value is replaced by $m(x)$, otherwise the value is left unchanged."
			,
			{
				CommandArgument<size_t>(ParameterDirection::In, "radius", "Median filtering radius.", 1),
				CommandArgument<double>(ParameterDirection::In, "magnitude", "Magnitude parameter $M$.", 10.0)
			},
			fbpSeeAlso())
		{
		}

	public:
		virtual void run(Image<pixel_t>& img, std::vector<ParamVariant>& args) const override
		{
			size_t r = pop<size_t>(args);
			float32_t M = (float32_t)pop<double>(args);
			deadPixelRemoval(img, r, M);
		}
	};
	

	class FBPCommand : public TwoImageInputOutputCommand<float32_t>
	{
	protected:
		friend class CommandList;

		FBPCommand() : TwoImageInputOutputCommand<float32_t>("fbp", "Performs filtered backprojection of data for which fbppreprocess has been called. Reconstruction settings are read from the metadata of the input image. This command is experimental and may change in the near future.",
			{
				CommandArgument<bool>(ParameterDirection::In, "use GPU", "Set to true to allow processing on a GPU.", true)
			},
			fbpSeeAlso())
		{
		}

	public:
		virtual void run(Image<float32_t>& in, Image<float32_t>& out, std::vector<ParamVariant>& args) const override
		{
			bool useGPU = pop<bool>(args);
			
			RecSettings sets = RecSettings::fromMeta(in.metadata);

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
				},
				fbpSeeAlso()
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
