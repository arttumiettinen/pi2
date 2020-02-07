
#include "filtercommands.h"
#include "commandmacros.h"

namespace pilib
{
	void addFilterCommands()
	{
		CommandList::add<BandpassFilterCommand>();
		CommandList::add<FFTCommand>();
		CommandList::add<InverseFFTCommand>();
		ADD_REAL(MinFilterCommand);
		ADD_REAL(MaxFilterCommand);
		ADD_REAL(MedianFilterCommand);
		ADD_REAL(VarianceFilterCommand);
		ADD_REAL(StddevFilterCommand);
		ADD_REAL(VaWeFilterCommand);
		ADD_REAL(OpeningFilterCommand);
		ADD_REAL(ClosingFilterCommand);
		//ADD_REAL(OpeningFilter2ParamCommand); These are not intuitive as they destroy the input image.
		//ADD_REAL(ClosingFilter2ParamCommand);
		ADD_REAL(BilateralFilterCommand);
		ADD_REAL(ApproxBilateralFilterCommand);
		ADD_REAL(GaussianFilterCommand);
		ADD_REAL(HighpassFilterCommand);
			
		CommandList::add<DerivativeCommand<uint8_t, float32_t> >();
		CommandList::add<DerivativeCommand<uint16_t, float32_t> >();
		CommandList::add<DerivativeCommand<uint32_t, float32_t> >();
		CommandList::add<DerivativeCommand<uint64_t, float32_t> >();
		CommandList::add<DerivativeCommand<int8_t, float32_t> >();
		CommandList::add<DerivativeCommand<int16_t, float32_t> >();
		CommandList::add<DerivativeCommand<int32_t, float32_t> >();
		CommandList::add<DerivativeCommand<int64_t, float32_t> >();
		CommandList::add<DerivativeCommand<float32_t, float32_t> >();

		CommandList::add<GradientCommand<uint8_t, float32_t> >();
		CommandList::add<GradientCommand<uint16_t, float32_t> >();
		CommandList::add<GradientCommand<uint32_t, float32_t> >();
		CommandList::add<GradientCommand<uint64_t, float32_t> >();
		CommandList::add<GradientCommand<int8_t, float32_t> >();
		CommandList::add<GradientCommand<int16_t, float32_t> >();
		CommandList::add<GradientCommand<int32_t, float32_t> >();
		CommandList::add<GradientCommand<int64_t, float32_t> >();
		CommandList::add<GradientCommand<float32_t, float32_t> >();

		CommandList::add<GradientMagnitudeCommand<uint8_t, float32_t> >();
		CommandList::add<GradientMagnitudeCommand<uint16_t, float32_t> >();
		CommandList::add<GradientMagnitudeCommand<uint32_t, float32_t> >();
		CommandList::add<GradientMagnitudeCommand<uint64_t, float32_t> >();
		CommandList::add<GradientMagnitudeCommand<int8_t, float32_t> >();
		CommandList::add<GradientMagnitudeCommand<int16_t, float32_t> >();
		CommandList::add<GradientMagnitudeCommand<int32_t, float32_t> >();
		CommandList::add<GradientMagnitudeCommand<int64_t, float32_t> >();
		CommandList::add<GradientMagnitudeCommand<float32_t, float32_t> >();

		ADD_REAL(SatoLineFilterCommand);
		ADD_REAL(FrangiLineFilterCommand);
		ADD_REAL2(MorphoRecCommand);
	}
}