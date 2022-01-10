#pragma once

#include "inpaint.h"
#include "commandsbase.h"

namespace pilib
{
	template<typename pixel_t> class InpaintNearestCommand : public OneImageInPlaceCommand<pixel_t>
	{
	protected:
		friend class CommandList;

		InpaintNearestCommand() : OneImageInPlaceCommand<pixel_t>("inpaintn", "Replaces pixels that have specific flag value by the value of the nearest pixel that does not have the specific flag value.",
			{
				CommandArgument<double>(ParameterDirection::In, "flag", "Values of pixels that have this (flag) value are inpainted.", 0)
			},
			"inpaintn, inpaintg")
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const override
		{
			double value = pop<double>(args);
			
			inpaintNearest<pixel_t>(in, pixelRound<pixel_t>(value));
		}
	};

	template<typename pixel_t> class InpaintGarciaCommand : public OneImageInPlaceCommand<pixel_t>
	{
	protected:
		friend class CommandList;

		InpaintGarciaCommand() : OneImageInPlaceCommand<pixel_t>("inpaintg", "Replaces pixels that have specific flag value by a value interpolated from the neighbouring pixels that do not have the specific flag value..",
			{
				CommandArgument<double>(ParameterDirection::In, "flag", "Values of pixels that have this (flag) value are inpainted.", 0)
			},
			"inpaintn, inpaintg")
		{
		}

	public:
		virtual void run(Image<pixel_t>& in, vector<ParamVariant>& args) const override
		{
			double value = pop<double>(args);

			inpaintGarcia<pixel_t>(in, pixelRound<pixel_t>(value));
		}
	};
}