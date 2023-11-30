
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
		ADD_ALL(WritePngCommand);
		ADD_ALL(WriteNRRDCommand);
		ADD_ALL(WriteRawCommand);
		ADD_ALL(WriteLZ4Command);
		ADD_ALL(WriteNN5Command);
		ADD_ALL(WriteRawBlockCommand);
		CommandList::add<WriteRGBRawCommand>();
		ADD_ALL(WriteSequenceCommand);
		ADD_ALL(WriteSequenceBlockCommand);
		ADD_ALL(WriteNN5BlockCommand);
		CommandList::add<EndConcurrentWriteCommand>();
	}
}