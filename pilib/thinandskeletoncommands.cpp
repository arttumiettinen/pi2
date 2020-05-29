
#include "thinandskeletoncommands.h"
#include "commandmacros.h"

namespace pilib
{
	void addThinAndSkeletonCommands()
	{
		//ADD_REAL(HybridThinCommand);
		//ADD_REAL(HybridSkeletonCommand);
		//ADD_REAL(ThinCommand);
		ADD_REAL(LineThinCommand);
		ADD_REAL(LineSkeletonCommand);
		ADD_REAL(SurfaceThinCommand);
		ADD_REAL(SurfaceSkeletonCommand);
		ADD_REAL(ClassifySkeletonCommand);
		ADD_REAL(ClassifyForTracingCommand);
		ADD_REAL(TraceLineSkeletonCommand);
		ADD_REAL(TraceLineSkeleton2Command);
		ADD_REAL(TraceLineSkeletonBlockCommand);
		ADD_REAL(TraceLineSkeletonBlock2Command);
		CommandList::add<CleanSkeletonCommand>();
		CommandList::add<PruneSkeletonCommand>();
		CommandList::add<GetPointsAndLinesCommand>();
		CommandList::add<WriteVtkCommand>();
		CommandList::add<WriteVtk2Command>();
		CommandList::add<WriteVtk3Command>();
		ADD_REAL(RemoveEdgesCommand);
		ADD_REAL(FillSkeletonCommand);
	}
}
