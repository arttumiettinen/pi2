
#include "dmapcommands.h"
#include "commandmacros.h"

namespace pilib
{
	void addDmapCommands()
	{
		ADD_REAL2(DistanceMap2Command);
		ADD_REAL2(DistanceMapCommand);
		ADD_REAL2(PrepareDistanceMapCommand);
		ADD_REAL(DistanceMapProcessDimensionCommand);
		ADD_REAL(SeededDistanceMapCommand);
		ADD_REAL(PolySeededDistanceMapCommand);
	}
}