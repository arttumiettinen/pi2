
#include "histogramcommands.h"
#include "commandmacros.h"

namespace pilib
{
	void addHistogramCommands()
	{
		ADD_REAL2(HistogramCommand);
		ADD_REAL(WeightedHistogramCommand);
		ADD_REAL2(Histogram2Command);
		ADD_REAL2(WeightedHistogram2Command);
	}
}