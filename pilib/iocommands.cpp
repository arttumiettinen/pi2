
#include "iocommands.h"
#include "commandmacros.h"

namespace pilib
{
	void addIOCommands(vector<Command*>& commands)
	{
		commands.insert(commands.end(),
			{
			ADD_ALL(NopSingleImageCommand),
			new FileInfoCommand(),
			new RawInfoCommand(),
			new SequenceInfoCommand(),
			ADD_ALL(WriteRawCommand),
			ADD_ALL(WriteRawBlockCommand),
			ADD_ALL(WriteSequenceCommand),
			ADD_ALL(WriteSequenceBlockCommand),
			}
		);
	}
}