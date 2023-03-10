
#include "othercommands.h"
#include "commandmacros.h"

namespace pilib
{
	void addOtherCommands()
	{
		ADD_REAL(BlockMatchCommand);
		ADD_REAL(BlockMatchNoDiskCommand);
		ADD_REAL(BlockMatchMultiCommand);
		CommandList::add<PointsToDeformedCommand>();
		CommandList::add<BlockMatchPartialLoadCommand>(); // Only one as this command determines data type itself
		CommandList::add<FilterDisplacementsCommand>(); // Only one no data type dependence
		ADD_REAL(PullbackCommand);
		ADD_REAL(PullbackNoDiskCommand);
		//ADD_REAL(StitchCommand)
		ADD_REAL(StitchVer2Command);
		ADD_REAL(StitchVer3Command);
		CommandList::add<DetermineWorldToLocalCommand>();

		ADD_REAL(NormalizeZCommand);

		ADD_REAL(CannyCommand);
		ADD_REAL(CannyPart1Command);
		ADD_REAL(CannyPart2Command);

		ADD_REAL(NoiseCommand);

		ADD_ALL(MontageCommand);
	}
}