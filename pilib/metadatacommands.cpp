
#include "metadatacommands.h"
#include "commandmacros.h"

namespace pilib
{

	void addMetadataCommands()
	{
		ADD_ALL(SetMetadataCommand);
		ADD_ALL(GetMetadataCommand);
	}
}