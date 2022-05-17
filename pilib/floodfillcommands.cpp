
#include "floodfillcommands.h"
#include "commandmacros.h"

namespace pilib
{
	void addFloodfillCommands()
	{
		ADD_REAL(FloodFillBlockCommand);
		ADD_REAL(FloodFillCommand);
		ADD_REAL(GrowCommand);
		ADD_REAL2(GrowPriorityCommand);
		ADD_REAL(GrowLabelsCommand);
		ADD_REAL(DualThresholdCommand);
	}
}
