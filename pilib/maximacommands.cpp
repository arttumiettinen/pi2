
#include "maximacommands.h"
#include "commandmacros.h"

namespace pilib
{

	void addMaximaCommands()
	{
		ADD_REAL(LocalMaximaCommand);
		ADD_REAL(CleanMaximaCommand);
		ADD_REAL(DrawMaximaCommand);
	}
}