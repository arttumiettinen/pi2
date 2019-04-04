
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
			//ADD_REAL(OpeningFilter2ParamCommand), These are not intuitive as they destroy the input image.
			//ADD_REAL(ClosingFilter2ParamCommand),
			ADD_REAL(BilateralFilterCommand),
			ADD_REAL(GaussianFilterCommand),
			ADD_REAL(HighpassFilterCommand),
			
			new DerivativeCommand<uint8_t, float32_t>(),
			new DerivativeCommand<uint16_t, float32_t>(),
			new DerivativeCommand<uint32_t, float32_t>(),
			new DerivativeCommand<uint64_t, float32_t>(),
			new DerivativeCommand<float32_t, float32_t>(),

			new GradientMagnitudeCommand<uint8_t, float32_t>(),
			new GradientMagnitudeCommand<uint16_t, float32_t>(),
			new GradientMagnitudeCommand<uint32_t, float32_t>(),
			new GradientMagnitudeCommand<uint64_t, float32_t>(),
			new GradientMagnitudeCommand<float32_t, float32_t>(),

			ADD_REAL(SatoLineFilterCommand),
			ADD_REAL(FrangiLineFilterCommand),
			ADD_REAL2(MorphoRecCommand)
			}
		);
	}
}