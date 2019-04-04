
#include "projectioncommands.h"
#include "commandmacros.h"

namespace pilib
{
	void addProjectionCommands(vector<Command*>& commands)
	{
		commands.insert(commands.end(),
			{
			
			ADD_REAL(MinAllPixelsCommand),
			ADD_REAL(MaxAllPixelsCommand),
			ADD_REAL(SumAllPixelsCommand),
			ADD_REAL(MeanAllPixelsCommand),
			//ADD_REAL(VarianceAllPixelsCommand),
			ADD_REAL(StdDevAllPixelsCommand),
			ADD_REAL(MinProjectCommand),
			ADD_REAL(MaxProjectCommand),
			ADD_REAL(SumProjectCommand),
			ADD_REAL(MeanProjectCommand)
			//ADD_REAL(VarianceProjectCommand),
			//ADD_REAL(StdDevProjectCommand)
			}
		);
	}
}