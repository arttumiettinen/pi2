
#include "projectioncommands.h"
#include "commandmacros.h"

namespace pilib
{
	void addProjectionCommands()
	{
		ADD_REAL(MinAllPixelsCommand);
		ADD_REAL(MaxAllPixelsCommand);
		ADD_REAL(SumAllPixelsCommand);
		ADD_REAL(SquareSumAllPixelsCommand);
		ADD_REAL(MeanAllPixelsCommand);
		ADD_REAL(MaskedMeanAllPixelsCommand);
		//ADD_REAL(VarianceAllPixelsCommand);
		ADD_REAL(StdDevAllPixelsCommand);
		ADD_REAL(MinProjectCommand);
		ADD_REAL(MaxProjectCommand);
		ADD_REAL(SumProjectCommand);
		//ADD_REAL(SquareSumProjectCommand);
		ADD_REAL(MeanProjectCommand)
		ADD_REAL2(MinProject2ImageCommand);
		ADD_REAL2(MaxProject2ImageCommand);
		//ADD_REAL(VarianceProjectCommand);
		//ADD_REAL(StdDevProjectCommand)
	}
}