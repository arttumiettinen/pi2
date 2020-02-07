
#include "structurecommands.h"
#include "commandmacros.h"

namespace pilib
{
	void addStructureCommands()
	{
		CommandList::add<CylindricalityCommand>();
		CommandList::add<CylinderOrientationCommand>();
		CommandList::add<PlateOrientationCommand>();
		ADD_REAL(MainOrientationColoringCommand);
		ADD_REAL(AxelssonColoringCommand);
		//new PlanarityCommand()
		ADD_REAL(SurfaceCurvatureCommand);
		ADD_REAL(MeanCurvatureCommand);
	}
}