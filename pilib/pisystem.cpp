
#include "pisystem.h"
#include "commandmacros.h"

namespace pilib
{
	void addSpecialSystemCommands(vector<Command*>& commands)
	{
		commands.insert(commands.end(),
			{
			ADD_REAL(ConvertCommand),
			ADD_ALL(NewLikeCommand)
			});
	}
}
