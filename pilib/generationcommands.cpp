
#include "generationcommands.h"
#include "commandmacros.h"

namespace pilib
{
	void addGenerationCommands()
	{
		ADD_REAL(RampCommand);
		ADD_REAL(Ramp3Command);
		ADD_REAL(SetPixelCommand);
		ADD_REAL(SetPixelsCommand);
		ADD_REAL(GetPixelsCommand);
		ADD_REAL(GetPixelsToTempFileCommand);
		ADD_REAL(SphereCommand);
		ADD_REAL(EllipsoidCommand);
		ADD_REAL(BoxCommand);
		ADD_REAL(GenericBoxCommand);
		ADD_REAL(LineCommand);
		ADD_REAL(CapsuleCommand);
		ADD_REAL(DrawGraphCommand);
		ADD_REAL(DrawGraph2Command);
	}
}