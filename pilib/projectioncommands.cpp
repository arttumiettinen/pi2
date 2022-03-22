
#include "projectioncommands.h"
#include "commandmacros.h"

namespace pilib
{
	

/*
This macro adds commands CMD<src_type, target_type, intermediate_type> such that target_type is equal to or wider than source type.
*/
#define ADD_REAL2_WIDER_TYPE(cmd) \
	CommandList::add<cmd<uint8_t, uint8_t> >(); \
	CommandList::add<cmd<uint8_t, uint16_t> >(); \
	CommandList::add<cmd<uint8_t, uint32_t> >(); \
	CommandList::add<cmd<uint8_t, uint64_t> >(); \
	CommandList::add<cmd<uint8_t, float32_t> >(); \
	CommandList::add<cmd<uint16_t, uint16_t> >(); \
	CommandList::add<cmd<uint16_t, uint32_t> >(); \
	CommandList::add<cmd<uint16_t, uint64_t> >(); \
	CommandList::add<cmd<uint16_t, float32_t> >(); \
	CommandList::add<cmd<uint32_t, uint32_t> >(); \
	CommandList::add<cmd<uint32_t, uint64_t> >(); \
	CommandList::add<cmd<uint32_t, float32_t> >(); \
	CommandList::add<cmd<uint64_t, uint64_t> >(); \
	CommandList::add<cmd<uint64_t, float32_t> >(); \
	CommandList::add<cmd<float32_t, float32_t> >(); \
	\
	CommandList::add<cmd<int8_t, int8_t> >(); \
	CommandList::add<cmd<int8_t, int16_t> >(); \
	CommandList::add<cmd<int8_t, int32_t> >(); \
	CommandList::add<cmd<int8_t, int64_t> >(); \
	CommandList::add<cmd<int8_t, float32_t> >(); \
	CommandList::add<cmd<int16_t, int16_t> >(); \
	CommandList::add<cmd<int16_t, int32_t> >(); \
	CommandList::add<cmd<int16_t, int64_t> >(); \
	CommandList::add<cmd<int16_t, float32_t> >(); \
	CommandList::add<cmd<int32_t, int32_t> >(); \
	CommandList::add<cmd<int32_t, int64_t> >(); \
	CommandList::add<cmd<int32_t, float32_t> >(); \
	CommandList::add<cmd<int64_t, int64_t> >(); \
	CommandList::add<cmd<int64_t, float32_t> >();


	void addProjectionCommands()
	{
		ADD_REAL2_WIDER_TYPE(MinAllPixelsCommand);
		ADD_REAL2_WIDER_TYPE(MaxAllPixelsCommand);
		ADD_REAL2_WIDER_TYPE(SumAllPixelsCommand);
		ADD_REAL2_WIDER_TYPE(SquareSumAllPixelsCommand);

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