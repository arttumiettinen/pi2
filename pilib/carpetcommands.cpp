
#include "carpetcommands.h"
#include "commandmacros.h"

namespace pilib
{

	void addCarpetCommands()
	{
		ADD_REAL(FindSurfaceCommand);
		ADD_REAL(FindSurface2Command);
		ADD_REAL(DrawHeightMapCommand);
		ADD_REAL(SetBeforeHeightMapCommand);
		ADD_REAL(SetAfterHeightMapCommand);
		ADD_REAL(ShiftZCommand);
	}

}