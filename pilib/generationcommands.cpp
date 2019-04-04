
#include "generationcommands.h"
#include "commandmacros.h"

namespace pilib
{
	void addGenerationCommands(vector<Command*>& commands)
	{
		commands.insert(commands.end(),
			{
			ADD_REAL(RampCommand),
			ADD_REAL(SphereCommand),
			ADD_REAL(BoxCommand),
			ADD_REAL(LineCommand)
			}
		);
	}
}