
#include "transformcommands.h"
#include "commandmacros.h"

namespace pilib
{
	void addTransformCommands()
	{
		ADD_ALL(Rotate90CWCommand);
		ADD_ALL(Rotate90CCWCommand);
		ADD_ALL(FlipCommand);
		ADD_ALL(RotateCommand);
		ADD_ALL(Rotate2Command);
		ADD_ALL(ResliceCommand);
		ADD_REAL2(BinCommand);
		ADD_REAL(MaskedBinCommand);
		ADD_ALL(CropCommand);
		ADD_ALL(ScaleCommand);
		ADD_REAL(ScaleLabelsCommand);
		ADD_ALL(GenericTransformCommand);
		ADD_ALL(TranslateCommand);
		ADD_ALL(Copy2Command);
		ADD_REAL(CartesianToCylindricalCommand);
	}

}
