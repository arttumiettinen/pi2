
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
			ADD_REAL(PullbackNoDiskCommand),
			//ADD_REAL(StitchCommand)
			ADD_REAL(StitchVer2Command),
			new DetermineWorldToLocalCommand(),

			ADD_REAL(DistanceMapCommand),
			ADD_REAL(HistogramCommand),

			ADD_REAL(FloodFillCommand),
			ADD_REAL(LabelCommand),

			ADD_REAL(NormalizeZCommand),

			ADD_REAL(CannyCommand),
			ADD_REAL(CannyPart1Command),
			ADD_REAL(CannyPart2Command),

			new CylindricalityCommand(),
			//new PlanarityCommand()

			ADD_REAL(DualThresholdCommand),
			ADD_REAL(GrowCommand),

			ADD_REAL(NoiseCommand)
			}
		);
	}
}