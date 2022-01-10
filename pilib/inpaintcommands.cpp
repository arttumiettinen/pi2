
#include "inpaintcommands.h"
#include "commandmacros.h"

namespace pilib
{
	void addInpaintCommands()
	{
		ADD_ALL(InpaintNearestCommand);
		ADD_REAL(InpaintGarciaCommand);
	}
}