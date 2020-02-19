
#include "projectioncommands.h"
#include "commandmacros.h"

namespace pilib
{
	/**
	This is used to select suitable distributed accumulator type for sum projections.
	The difference to itl2::sum_intermediate_type is that we must use float32_t instead of double
	in the final accumulation for floating-point images as pi2 does not currently support double/float64_t pixels.
	*/
	template<class pixel_t> struct sum_distributed_intermediate_type {
		using type = typename std::conditional <
			std::is_floating_point_v<pixel_t>,
			float32_t,  // floating-point pixels -> float32_t accumulator
			typename std::conditional<std::is_signed_v<pixel_t>,
				int64_t, // signed integer pixels -> int64 accumulator
				uint64_t // unsigned integer pixels -> uint64 accumulator
			>::type
		>::type;
	};

/*
This macro adds commands CMD<src_type, target_type, intermediate_type> such that target_type is equal to or wider than source type.
*/
#define ADD_REAL2_WIDER_TYPE_SUM(cmd) \
	CommandList::add<cmd<uint8_t, uint8_t, sum_distributed_intermediate_type<uint8_t>::type> >(); \
	CommandList::add<cmd<uint8_t, uint16_t, sum_distributed_intermediate_type<uint8_t>::type> >(); \
	CommandList::add<cmd<uint8_t, uint32_t, sum_distributed_intermediate_type<uint8_t>::type> >(); \
	CommandList::add<cmd<uint8_t, uint64_t, sum_distributed_intermediate_type<uint8_t>::type> >(); \
	CommandList::add<cmd<uint8_t, float32_t, sum_distributed_intermediate_type<uint8_t>::type> >(); \
	CommandList::add<cmd<uint16_t, uint16_t, sum_distributed_intermediate_type<uint16_t>::type> >(); \
	CommandList::add<cmd<uint16_t, uint32_t, sum_distributed_intermediate_type<uint16_t>::type> >(); \
	CommandList::add<cmd<uint16_t, uint64_t, sum_distributed_intermediate_type<uint16_t>::type> >(); \
	CommandList::add<cmd<uint16_t, float32_t, sum_distributed_intermediate_type<uint16_t>::type> >(); \
	CommandList::add<cmd<uint32_t, uint32_t, sum_distributed_intermediate_type<uint32_t>::type> >(); \
	CommandList::add<cmd<uint32_t, uint64_t, sum_distributed_intermediate_type<uint32_t>::type> >(); \
	CommandList::add<cmd<uint32_t, float32_t, sum_distributed_intermediate_type<uint32_t>::type> >(); \
	CommandList::add<cmd<uint64_t, uint64_t, sum_distributed_intermediate_type<uint64_t>::type> >(); \
	CommandList::add<cmd<uint64_t, float32_t, sum_distributed_intermediate_type<uint64_t>::type> >(); \
	CommandList::add<cmd<float32_t, float32_t, sum_distributed_intermediate_type<float32_t>::type> >(); \
	\
	CommandList::add<cmd<int8_t, int8_t, sum_distributed_intermediate_type<int8_t>::type> >(); \
	CommandList::add<cmd<int8_t, int16_t, sum_distributed_intermediate_type<int8_t>::type> >(); \
	CommandList::add<cmd<int8_t, int32_t, sum_distributed_intermediate_type<int8_t>::type> >(); \
	CommandList::add<cmd<int8_t, int64_t, sum_distributed_intermediate_type<int8_t>::type> >(); \
	CommandList::add<cmd<int8_t, float32_t, sum_distributed_intermediate_type<int8_t>::type> >(); \
	CommandList::add<cmd<int16_t, int16_t, sum_distributed_intermediate_type<int16_t>::type> >(); \
	CommandList::add<cmd<int16_t, int32_t, sum_distributed_intermediate_type<int16_t>::type> >(); \
	CommandList::add<cmd<int16_t, int64_t, sum_distributed_intermediate_type<int16_t>::type> >(); \
	CommandList::add<cmd<int16_t, float32_t, sum_distributed_intermediate_type<int16_t>::type> >(); \
	CommandList::add<cmd<int32_t, int32_t, sum_distributed_intermediate_type<int32_t>::type> >(); \
	CommandList::add<cmd<int32_t, int64_t, sum_distributed_intermediate_type<int32_t>::type> >(); \
	CommandList::add<cmd<int32_t, float32_t, sum_distributed_intermediate_type<int32_t>::type> >(); \
	CommandList::add<cmd<int64_t, int64_t, sum_distributed_intermediate_type<int64_t>::type> >(); \
	CommandList::add<cmd<int64_t, float32_t, sum_distributed_intermediate_type<int64_t>::type> >();


	void addProjectionCommands()
	{
		ADD_REAL(MinAllPixelsCommand);
		ADD_REAL(MaxAllPixelsCommand);
		ADD_REAL2_WIDER_TYPE_SUM(SumAllPixelsCommand);
		ADD_REAL2_WIDER_TYPE_SUM(SquareSumAllPixelsCommand);
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