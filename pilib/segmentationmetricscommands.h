#pragma once

#include "command.h"
#include "commandsbase.h"
#include "segmentationmetrics.h"


namespace pilib
{

	inline std::string segmentationMetricsSeeAlso()
	{
		return "-, -";
	}


	inline std::string segmentationMetricsMethodsHelp()
	{
		return
			"**SegmentationMetrics**\n"
			"Calculates the segmentation metrics Accuracy, Sensitivity, Specificity, and the Dice Coefficent using confusion matrix."
			"\n"
			"\n";
		"\n";
	}


	template<typename pixel_t> class segmentationMetricsCommand : public OneImageCommand<pixel_t>
	{
	protected:
		friend class CommandList;

		segmentationMetricsCommand() : OneImageCommand<pixel_t>("segmentationMetrics",
			"Calculates Confusion matrix and then with it Acc, Se, Spe, Dice.\n\n" +
			segmentationMetricsMethodsHelp(),
			{
				CommandArgument<Image<pixel_t>>(ParameterDirection::In, "segImage", "Segmented Image to calculate the metrics for.")
			})
		{

		}

	public:
		virtual void run(const Image<pixel_t>& img, std::vector<ParamVariant>& args) const override
		{
			Image<pixel_t>& segImage = *pop<Image<pixel_t>* >(args);

			calculateSegmentationMetrics(img, segImage);

		}


	};


}