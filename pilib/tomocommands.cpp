
#include "tomocommands.h"
#include "commandmacros.h"

namespace pilib
{
	void addTomoCommands()
	{
		CommandList::add<FBPPreprocessCommand>();
		CommandList::add<FBPCommand>();
		CommandList::add<CreateFBPFilterCommand>();
		CommandList::add<DefaultRecSettingsCommand>();
		ADD_REAL(DeadPixelRemovalCommand);
	}

}