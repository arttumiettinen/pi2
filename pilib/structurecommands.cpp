
#include "structurecommands.h"
#include "commandmacros.h"

namespace pilib
{
	void addStructureCommands()
	{
		CommandList::add<CylindricalityCommand>();
		CommandList::add<CylinderOrientationCommand>();
		CommandList::add<PlateOrientationCommand>();
		CommandList::add<OrientationDifferenceCommand>();
		ADD_REAL(MainOrientationColoringCommand);
		ADD_REAL(AxelssonColoringCommand);
		//new PlanarityCommand()
		ADD_REAL(SurfaceCurvatureCommand);
		ADD_REAL(MeanCurvatureCommand);
		ADD_REAL(SurfaceAreaCommand);
	}
}