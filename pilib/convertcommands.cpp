
#include "convertcommands.h"
#include "commandmacros.h"

namespace pilib
{
	void addConvertCommands()
	{
		ADD_REAL(ConvertCommand);
		ADD_REAL(ConvertInPlaceCommand);
	}
}
