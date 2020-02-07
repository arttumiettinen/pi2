#pragma once

#include "commandlist.h"

//#define ADD_SIGNED(cmd) \
//	CommandList::add<cmd<float32_t>()

#define ADD_REAL(cmd) \
	CommandList::add<cmd<uint8_t> >(); \
	CommandList::add<cmd<uint16_t> >(); \
	CommandList::add<cmd<uint32_t> >(); \
	CommandList::add<cmd<uint64_t> >(); \
	CommandList::add<cmd<int8_t> >(); \
	CommandList::add<cmd<int16_t> >(); \
	CommandList::add<cmd<int32_t> >(); \
	CommandList::add<cmd<int64_t> >(); \
	CommandList::add<cmd<float32_t> >();

#define ADD_COMPLEX(cmd) \
	CommandList::add<cmd<complex32_t> >();

#define ADD_ALL(cmd) \
	ADD_REAL(cmd) \
	ADD_COMPLEX(cmd)


// TODO: Should I add signed int versions here?
// If we support all possible argument combinations, there will be a bit too many commands!
// I will add cross-type versions only for uints and floats, and signed ints are going to get only the non-cross type versions for now.
#define ADD_REAL2(cmd) \
	CommandList::add<cmd<uint8_t, uint8_t> >(); \
	CommandList::add<cmd<uint8_t, uint16_t> >(); \
	CommandList::add<cmd<uint8_t, uint32_t> >(); \
	CommandList::add<cmd<uint8_t, uint64_t> >(); \
	CommandList::add<cmd<uint8_t, float32_t> >(); \
	CommandList::add<cmd<uint16_t, uint8_t> >(); \
	CommandList::add<cmd<uint16_t, uint16_t> >(); \
	CommandList::add<cmd<uint16_t, uint32_t> >(); \
	CommandList::add<cmd<uint16_t, uint64_t> >(); \
	CommandList::add<cmd<uint16_t, float32_t> >(); \
	CommandList::add<cmd<uint32_t, uint8_t> >(); \
	CommandList::add<cmd<uint32_t, uint16_t> >(); \
	CommandList::add<cmd<uint32_t, uint32_t> >(); \
	CommandList::add<cmd<uint32_t, uint64_t> >(); \
	CommandList::add<cmd<uint32_t, float32_t> >(); \
	CommandList::add<cmd<uint64_t, uint8_t> >(); \
	CommandList::add<cmd<uint64_t, uint16_t> >(); \
	CommandList::add<cmd<uint64_t, uint32_t> >(); \
	CommandList::add<cmd<uint64_t, uint64_t> >(); \
	CommandList::add<cmd<uint64_t, float32_t> >(); \
	\
	CommandList::add<cmd<int8_t, int8_t> >(); \
	CommandList::add<cmd<int16_t, int16_t> >(); \
	CommandList::add<cmd<int32_t, int32_t> >(); \
	CommandList::add<cmd<int64_t, int64_t> >(); \
	\
	CommandList::add<cmd<float32_t, uint8_t> >(); \
	CommandList::add<cmd<float32_t, uint16_t> >(); \
	CommandList::add<cmd<float32_t, uint32_t> >(); \
	CommandList::add<cmd<float32_t, uint64_t> >(); \
	CommandList::add<cmd<float32_t, float32_t> >();

#define ADD_COMPLEX2(cmd) \
	CommandList::add<cmd<complex32_t, complex32_t> >();

#define ADD_ALL2(cmd) \
	ADD_REAL2(cmd) \
	ADD_COMPLEX2(cmd)
	
	
	
	
