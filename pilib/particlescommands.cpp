
#include "particlescommands.h"
#include "commandmacros.h"

namespace pilib
{
	void addParticlesCommands(vector<Command*>& commands)
	{
		commands.insert(commands.end(),
			{
			ADD_REAL(RegionRemovalCommand),
			ADD_REAL(AnalyzeParticlesBlockCommand),
			ADD_REAL(AnalyzeParticlesCommand),
			new HeadersCommand(),
			new ListAnalyzersCommand(),
			ADD_REAL(FillParticlesCommand)
			}
		);
	}
}