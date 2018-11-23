#pragma once

#define ADD_SIGNED(cmd) \
	new cmd<float32_t>()

#define ADD_REAL(cmd) \
	new cmd<uint8_t>(), \
	new cmd<uint16_t>(), \
	new cmd<uint32_t>(), \
	new cmd<uint64_t>(), \
	new cmd<float32_t>()

#define ADD_COMPLEX(cmd) \
	new cmd<complex32_t>()

#define ADD_ALL(cmd) \
	ADD_REAL(cmd), \
	ADD_COMPLEX(cmd)



#define ADD_REAL2(cmd) \
	new cmd<uint8_t, uint8_t>(), \
	new cmd<uint8_t, uint16_t>(), \
	new cmd<uint8_t, uint32_t>(), \
	new cmd<uint8_t, uint64_t>(), \
	new cmd<uint8_t, float32_t>(), \
	new cmd<uint16_t, uint8_t>(), \
	new cmd<uint16_t, uint16_t>(), \
	new cmd<uint16_t, uint32_t>(), \
	new cmd<uint16_t, uint64_t>(), \
	new cmd<uint16_t, float32_t>(), \
	new cmd<uint32_t, uint8_t>(), \
	new cmd<uint32_t, uint16_t>(), \
	new cmd<uint32_t, uint32_t>(), \
	new cmd<uint32_t, uint64_t>(), \
	new cmd<uint32_t, float32_t>(), \
	new cmd<uint64_t, uint8_t>(), \
	new cmd<uint64_t, uint16_t>(), \
	new cmd<uint64_t, uint32_t>(), \
	new cmd<uint64_t, uint64_t>(), \
	new cmd<uint64_t, float32_t>(), \
	new cmd<float32_t, uint8_t>(), \
	new cmd<float32_t, uint16_t>(), \
	new cmd<float32_t, uint32_t>(), \
	new cmd<float32_t, uint64_t>(), \
	new cmd<float32_t, float32_t>()

#define ADD_COMPLEX2(cmd) \
	new cmd<complex32_t, complex32_t>()

#define ADD_ALL2(cmd) \
	ADD_REAL2(cmd), \
	ADD_COMPLEX2(cmd)
	
	
	
	
