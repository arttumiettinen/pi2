
#include "evalcommands.h"
#include "commandmacros.h"

namespace pilib
{
	void addEvalCommands()
	{
		ADD_REAL(Eval1Command);
		ADD_REAL2(Eval2Command);

		// NOTE: We need to do something to avoid the data type explosion below!

		CommandList::add<Eval3Command<uint8_t, uint8_t, uint8_t> >();
		CommandList::add<Eval3Command<uint8_t, uint8_t, uint16_t> >();
		CommandList::add<Eval3Command<uint8_t, uint8_t, uint32_t> >();
		CommandList::add<Eval3Command<uint8_t, uint8_t, uint64_t> >();
		CommandList::add<Eval3Command<uint8_t, uint8_t, float32_t> >();
		CommandList::add<Eval3Command<uint8_t, uint16_t, uint8_t> >();
		CommandList::add<Eval3Command<uint8_t, uint16_t, uint16_t> >();
		CommandList::add<Eval3Command<uint8_t, uint16_t, uint32_t> >();
		CommandList::add<Eval3Command<uint8_t, uint16_t, uint64_t> >();
		CommandList::add<Eval3Command<uint8_t, uint16_t, float32_t> >();
		CommandList::add<Eval3Command<uint8_t, uint32_t, uint8_t> >();
		CommandList::add<Eval3Command<uint8_t, uint32_t, uint16_t> >();
		CommandList::add<Eval3Command<uint8_t, uint32_t, uint32_t> >();
		CommandList::add<Eval3Command<uint8_t, uint32_t, uint64_t> >();
		CommandList::add<Eval3Command<uint8_t, uint32_t, float32_t> >();
		CommandList::add<Eval3Command<uint8_t, uint64_t, uint8_t> >();
		CommandList::add<Eval3Command<uint8_t, uint64_t, uint16_t> >();
		CommandList::add<Eval3Command<uint8_t, uint64_t, uint32_t> >();
		CommandList::add<Eval3Command<uint8_t, uint64_t, uint64_t> >();
		CommandList::add<Eval3Command<uint8_t, uint64_t, float32_t> >();
		CommandList::add<Eval3Command<uint8_t, float32_t, uint8_t> >();
		CommandList::add<Eval3Command<uint8_t, float32_t, uint16_t> >();
		CommandList::add<Eval3Command<uint8_t, float32_t, uint32_t> >();
		CommandList::add<Eval3Command<uint8_t, float32_t, uint64_t> >();
		CommandList::add<Eval3Command<uint8_t, float32_t, float32_t> >();


		CommandList::add<Eval3Command<uint16_t, uint8_t, uint8_t> >();
		CommandList::add<Eval3Command<uint16_t, uint8_t, uint16_t> >();
		CommandList::add<Eval3Command<uint16_t, uint8_t, uint32_t> >();
		CommandList::add<Eval3Command<uint16_t, uint8_t, uint64_t> >();
		CommandList::add<Eval3Command<uint16_t, uint8_t, float32_t> >();
		CommandList::add<Eval3Command<uint16_t, uint16_t, uint8_t> >();
		CommandList::add<Eval3Command<uint16_t, uint16_t, uint16_t> >();
		CommandList::add<Eval3Command<uint16_t, uint16_t, uint32_t> >();
		CommandList::add<Eval3Command<uint16_t, uint16_t, uint64_t> >();
		CommandList::add<Eval3Command<uint16_t, uint16_t, float32_t> >();
		CommandList::add<Eval3Command<uint16_t, uint32_t, uint8_t> >();
		CommandList::add<Eval3Command<uint16_t, uint32_t, uint16_t> >();
		CommandList::add<Eval3Command<uint16_t, uint32_t, uint32_t> >();
		CommandList::add<Eval3Command<uint16_t, uint32_t, uint64_t> >();
		CommandList::add<Eval3Command<uint16_t, uint32_t, float32_t> >();
		CommandList::add<Eval3Command<uint16_t, uint64_t, uint8_t> >();
		CommandList::add<Eval3Command<uint16_t, uint64_t, uint16_t> >();
		CommandList::add<Eval3Command<uint16_t, uint64_t, uint32_t> >();
		CommandList::add<Eval3Command<uint16_t, uint64_t, uint64_t> >();
		CommandList::add<Eval3Command<uint16_t, uint64_t, float32_t> >();
		CommandList::add<Eval3Command<uint16_t, float32_t, uint8_t> >();
		CommandList::add<Eval3Command<uint16_t, float32_t, uint16_t> >();
		CommandList::add<Eval3Command<uint16_t, float32_t, uint32_t> >();
		CommandList::add<Eval3Command<uint16_t, float32_t, uint64_t> >();
		CommandList::add<Eval3Command<uint16_t, float32_t, float32_t> >();


		CommandList::add<Eval3Command<uint32_t, uint8_t, uint8_t> >();
		CommandList::add<Eval3Command<uint32_t, uint8_t, uint16_t> >();
		CommandList::add<Eval3Command<uint32_t, uint8_t, uint32_t> >();
		CommandList::add<Eval3Command<uint32_t, uint8_t, uint64_t> >();
		CommandList::add<Eval3Command<uint32_t, uint8_t, float32_t> >();
		CommandList::add<Eval3Command<uint32_t, uint16_t, uint8_t> >();
		CommandList::add<Eval3Command<uint32_t, uint16_t, uint16_t> >();
		CommandList::add<Eval3Command<uint32_t, uint16_t, uint32_t> >();
		CommandList::add<Eval3Command<uint32_t, uint16_t, uint64_t> >();
		CommandList::add<Eval3Command<uint32_t, uint16_t, float32_t> >();
		CommandList::add<Eval3Command<uint32_t, uint32_t, uint8_t> >();
		CommandList::add<Eval3Command<uint32_t, uint32_t, uint16_t> >();
		CommandList::add<Eval3Command<uint32_t, uint32_t, uint32_t> >();
		CommandList::add<Eval3Command<uint32_t, uint32_t, uint64_t> >();
		CommandList::add<Eval3Command<uint32_t, uint32_t, float32_t> >();
		CommandList::add<Eval3Command<uint32_t, uint64_t, uint8_t> >();
		CommandList::add<Eval3Command<uint32_t, uint64_t, uint16_t> >();
		CommandList::add<Eval3Command<uint32_t, uint64_t, uint32_t> >();
		CommandList::add<Eval3Command<uint32_t, uint64_t, uint64_t> >();
		CommandList::add<Eval3Command<uint32_t, uint64_t, float32_t> >();
		CommandList::add<Eval3Command<uint32_t, float32_t, uint8_t> >();
		CommandList::add<Eval3Command<uint32_t, float32_t, uint16_t> >();
		CommandList::add<Eval3Command<uint32_t, float32_t, uint32_t> >();
		CommandList::add<Eval3Command<uint32_t, float32_t, uint64_t> >();
		CommandList::add<Eval3Command<uint32_t, float32_t, float32_t> >();



		CommandList::add<Eval3Command<uint64_t, uint8_t, uint8_t> >();
		CommandList::add<Eval3Command<uint64_t, uint8_t, uint16_t> >();
		CommandList::add<Eval3Command<uint64_t, uint8_t, uint32_t> >();
		CommandList::add<Eval3Command<uint64_t, uint8_t, uint64_t> >();
		CommandList::add<Eval3Command<uint64_t, uint8_t, float32_t> >();
		CommandList::add<Eval3Command<uint64_t, uint16_t, uint8_t> >();
		CommandList::add<Eval3Command<uint64_t, uint16_t, uint16_t> >();
		CommandList::add<Eval3Command<uint64_t, uint16_t, uint32_t> >();
		CommandList::add<Eval3Command<uint64_t, uint16_t, uint64_t> >();
		CommandList::add<Eval3Command<uint64_t, uint16_t, float32_t> >();
		CommandList::add<Eval3Command<uint64_t, uint32_t, uint8_t> >();
		CommandList::add<Eval3Command<uint64_t, uint32_t, uint16_t> >();
		CommandList::add<Eval3Command<uint64_t, uint32_t, uint32_t> >();
		CommandList::add<Eval3Command<uint64_t, uint32_t, uint64_t> >();
		CommandList::add<Eval3Command<uint64_t, uint32_t, float32_t> >();
		CommandList::add<Eval3Command<uint64_t, uint64_t, uint8_t> >();
		CommandList::add<Eval3Command<uint64_t, uint64_t, uint16_t> >();
		CommandList::add<Eval3Command<uint64_t, uint64_t, uint32_t> >();
		CommandList::add<Eval3Command<uint64_t, uint64_t, uint64_t> >();
		CommandList::add<Eval3Command<uint64_t, uint64_t, float32_t> >();
		CommandList::add<Eval3Command<uint64_t, float32_t, uint8_t> >();
		CommandList::add<Eval3Command<uint64_t, float32_t, uint16_t> >();
		CommandList::add<Eval3Command<uint64_t, float32_t, uint32_t> >();
		CommandList::add<Eval3Command<uint64_t, float32_t, uint64_t> >();
		CommandList::add<Eval3Command<uint64_t, float32_t, float32_t> >();



		CommandList::add<Eval3Command<float32_t, uint8_t, uint8_t> >();
		CommandList::add<Eval3Command<float32_t, uint8_t, uint16_t> >();
		CommandList::add<Eval3Command<float32_t, uint8_t, uint32_t> >();
		CommandList::add<Eval3Command<float32_t, uint8_t, uint64_t> >();
		CommandList::add<Eval3Command<float32_t, uint8_t, float32_t> >();
		CommandList::add<Eval3Command<float32_t, uint16_t, uint8_t> >();
		CommandList::add<Eval3Command<float32_t, uint16_t, uint16_t> >();
		CommandList::add<Eval3Command<float32_t, uint16_t, uint32_t> >();
		CommandList::add<Eval3Command<float32_t, uint16_t, uint64_t> >();
		CommandList::add<Eval3Command<float32_t, uint16_t, float32_t> >();
		CommandList::add<Eval3Command<float32_t, uint32_t, uint8_t> >();
		CommandList::add<Eval3Command<float32_t, uint32_t, uint16_t> >();
		CommandList::add<Eval3Command<float32_t, uint32_t, uint32_t> >();
		CommandList::add<Eval3Command<float32_t, uint32_t, uint64_t> >();
		CommandList::add<Eval3Command<float32_t, uint32_t, float32_t> >();
		CommandList::add<Eval3Command<float32_t, uint64_t, uint8_t> >();
		CommandList::add<Eval3Command<float32_t, uint64_t, uint16_t> >();
		CommandList::add<Eval3Command<float32_t, uint64_t, uint32_t> >();
		CommandList::add<Eval3Command<float32_t, uint64_t, uint64_t> >();
		CommandList::add<Eval3Command<float32_t, uint64_t, float32_t> >();
		CommandList::add<Eval3Command<float32_t, float32_t, uint8_t> >();
		CommandList::add<Eval3Command<float32_t, float32_t, uint16_t> >();
		CommandList::add<Eval3Command<float32_t, float32_t, uint32_t> >();
		CommandList::add<Eval3Command<float32_t, float32_t, uint64_t> >();
		CommandList::add<Eval3Command<float32_t, float32_t, float32_t> >();
	}
}