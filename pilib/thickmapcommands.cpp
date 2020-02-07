
#include "thickmapcommands.h"
#include "commandmacros.h"

namespace pilib
{

	void addTmapCommands()
	{
		ADD_REAL(Danielsson2Command);

		ADD_REAL(RoundDistanceRidge2Command);
		
		//ADD_REAL(DrawSpheres2BlockCommand);
		ADD_REAL(DrawSpheres2ProcessDimensionCommand);
		ADD_REAL(DrawSpheres2Command);
		
		ADD_REAL2(FinalizeThicknessMapCommand);
		
		
		ADD_REAL2(ThicknessMapCommand);
	}

}
