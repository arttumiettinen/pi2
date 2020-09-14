
#include "metadatacommands.h"
#include "commandmacros.h"

namespace pilib
{

	void addMetadataCommands()
	{
		ADD_ALL(SetMetadataCommand);
		ADD_ALL(GetMetadataCommand);
		ADD_ALL(WriteMetadataCommand);
		ADD_ALL(ReadMetadataCommand);
		ADD_ALL(ClearMetadataCommand);
		ADD_ALL(ListMetadataCommand);
	}
}