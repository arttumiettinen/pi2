
#include "iocommands.h"
#include "commandmacros.h"

namespace pilib
{
	void addIOCommands()
	{
		ADD_ALL(NopSingleImageCommand);
		CommandList::add<FileInfoCommand>();
		CommandList::add<IsImageFileCommand>();
		CommandList::add<ShowFileInfoCommand>();
		CommandList::add<ShowRawInfoCommand>();
		CommandList::add<ShowSequenceInfoCommand>();
		ADD_ALL(WriteTiffCommand);
		ADD_ALL(WriteNRRDCommand);
		ADD_ALL(WriteRawCommand);
		ADD_ALL(WriteLZ4Command);
		ADD_ALL(WriteNN5Command);
		//ADD_ALL(WriteRawBlockCommand);
		ADD_ALL(WriteRawBlock2Command);
		CommandList::add<WriteRGBRawCommand>();
		ADD_ALL(WriteSequenceCommand);
		//ADD_ALL(WriteSequenceBlockCommand);
		ADD_ALL(WriteSequenceBlock2Command);
		//ADD_ALL(WriteNN5BlockCommand);
		ADD_ALL(WriteNN5Block2Command);
		CommandList::add<EndConcurrentWriteCommand>();
	}
}