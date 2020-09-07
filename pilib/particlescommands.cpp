
#include "particlescommands.h"
#include "commandmacros.h"

namespace pilib
{
	void addParticlesCommands()
	{
		ADD_REAL(RegionRemovalCommand);
		ADD_REAL(AnalyzeParticlesBlockCommand);
		ADD_REAL(PrepareAnalyzeParticlesCommand);
		ADD_REAL(AnalyzeParticlesCommand);
		ADD_REAL(LabelCommand);
		ADD_REAL(AnalyzeLabelsCommand);
		CommandList::add<HeadersCommand>();
		CommandList::add<ListAnalyzersCommand>();
		ADD_REAL(FillParticlesCommand);
		ADD_REAL(DrawEllipsoidsCommand);
		ADD_REAL(GreedyColoringCommand);

		ADD_REAL(CSACommand);
	}
}