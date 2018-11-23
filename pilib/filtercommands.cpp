
#include "filtercommands.h"
#include "commandmacros.h"

namespace pilib
{
	void addFilterCommands(vector<Command*>& commands)
	{
		commands.insert(commands.end(),
			{
			new BandpassFilterCommand(),
			new FFTCommand(),
			new InverseFFTCommand(),
			ADD_REAL(MinFilterCommand),
			ADD_REAL(MaxFilterCommand),
			ADD_REAL(VarianceFilterCommand),
			ADD_REAL(StddevFilterCommand),
			ADD_REAL(VaWeFilterCommand),
			ADD_REAL(OpeningFilterCommand),
			ADD_REAL(ClosingFilterCommand),
			ADD_REAL(OpeningFilter2ParamCommand),
			ADD_REAL(ClosingFilter2ParamCommand),
			ADD_REAL(BilateralFilterCommand),
			ADD_REAL(GaussianFilterCommand),
			ADD_REAL(HighpassFilterCommand),
			ADD_SIGNED(DerivativeCommand),
			ADD_REAL(SatoLineFilterCommand),
			ADD_REAL(FrangiLineFilterCommand)
			}
		);
	}
}