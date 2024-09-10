
#include "autothresholdcommands.h"
#include "commandmacros.h"

namespace pilib
{

	void addAutoThresholdCommands()
	{
		ADD_REAL(AutoThresholdValueCommand);
		ADD_REAL(AutoThresholdCommand);
		ADD_REAL(LocalThresholdCommand);
	}
}