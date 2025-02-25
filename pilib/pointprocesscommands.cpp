
#include "pointprocesscommands.h"
#include "commandmacros.h"

namespace pilib
{
	void addPointProcessCommands()
	{
		ADD_ALL(SwapByteOrderCommand);
		ADD_ALL(NegateCommand);
		ADD_ALL(ExponentiateCommand);
		ADD_ALL(SquareCommand);
		ADD_ALL(SquareRootCommand);
		ADD_ALL(AbsCommand);
		ADD_ALL(LogCommand);
		ADD_ALL(Log10Command);
		ADD_ALL(SinCommand);
		ADD_ALL(CosCommand);
		ADD_ALL(TanCommand);
		ADD_ALL(InvCommand);
		ADD_REAL(RoundCommand);
		ADD_REAL(CeilCommand);
		ADD_REAL(FloorCommand);

		ADD_REAL(ReplaceCommand);


		CommandList::add<ConjugateComplexCommand>();
		CommandList::add<NormalizeComplexCommand>();
		CommandList::add<RealComplexCommand>();
		CommandList::add<ImagComplexCommand>();
		CommandList::add<ArgComplexCommand>();
		CommandList::add<NormSquaredComplexCommand>();

		ADD_ALL2(AddCommand);
		ADD_ALL2(SubtractCommand);
		ADD_ALL2(InvSubtractCommand);
		ADD_ALL2(DivideCommand);
		ADD_ALL2(MultiplyCommand);
		ADD_ALL(SetCommand);
		ADD_REAL2(ThresholdCommand);
		ADD_REAL2(MaxCommand);
		ADD_REAL2(MinCommand);

		ADD_REAL(AddConstantCommand);
		ADD_REAL(SubtractConstantCommand);
		ADD_REAL(InvSubtractConstantCommand);
		ADD_REAL(DivideConstantCommand);
		ADD_REAL(MultiplyConstantCommand);
		ADD_REAL(SetConstantCommand);
		ADD_REAL(ThresholdConstantCommand);
		ADD_REAL(MaxConstantCommand);
		ADD_REAL(MinConstantCommand);

		ADD_REAL(ThresholdRangeCommand);
		ADD_REAL(ThresholdPeriodicCommand);
		ADD_REAL(DoubleThresholdCommand);
		ADD_REAL(LinearMapCommand);

		ADD_ALL2(CopyCommand);
		ADD_REAL(SetEdgesCommand);

		ADD_REAL(ClampCommand);
		CommandList::add<EnsureFiniteCommand<float32_t> >();
	}
}