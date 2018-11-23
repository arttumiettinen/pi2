
#include "thinandskeletoncommands.h"
#include "commandmacros.h"

namespace pilib
{
	void addThinAndSkeletonCommands(vector<Command*>& commands)
	{
		commands.insert(commands.end(),
			{
			ADD_REAL(HybridThinCommand),
			ADD_REAL(HybridSkeletonCommand),
			ADD_REAL(LineThinCommand),
			ADD_REAL(LineSkeletonCommand),
			ADD_REAL(ClassifySkeletonCommand),
			ADD_REAL(TraceLineSkeletonCommand),
			ADD_REAL(TraceLineSkeletonBlockCommand),
			new CleanSkeletonCommand()
			}
		);
	}
}
