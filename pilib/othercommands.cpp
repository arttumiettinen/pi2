
#include "othercommands.h"
#include "commandmacros.h"

namespace pilib
{
	void addOtherCommands(vector<Command*>& commands)
	{
		commands.insert(commands.end(),
			{
			ADD_REAL(BlockMatchCommand),
			new BlockMatchPartialLoadCommand(), // Only one, as this command determines data type itself
			new FilterDisplacementsCommand(), // Only one, no data type dependence
			ADD_REAL(PullbackCommand),
			//ADD_REAL(StitchCommand)
			ADD_REAL(StitchVer2Command),
			new DetermineWorldToLocalCommand(),

			ADD_REAL(DistanceMapCommand),
			ADD_REAL(HistogramCommand),

			ADD_REAL(MinAllPixelsCommand),
			ADD_REAL(MaxAllPixelsCommand),
			ADD_REAL(SumAllPixelsCommand),
			ADD_REAL(MeanAllPixelsCommand),
			ADD_REAL(MinProjectCommand),
			ADD_REAL(MaxProjectCommand),
			ADD_REAL(SumProjectCommand),
			ADD_REAL(MeanProjectCommand),

			ADD_REAL(FloodFillCommand),

			ADD_REAL(NormalizeZCommand),
			ADD_REAL(RegionRemovalCommand),

			ADD_REAL(RampCommand)
			}
		);
	}
}