#pragma once

#include <stdlib.h>
#include <cstddef>
#include <stdint.h>
#include <complex>
#include <utility>

using namespace std::literals;
using namespace std::complex_literals;
using namespace std::literals::complex_literals;

namespace itl2
{

    /**
    Type for 32bit float images.
    */
    typedef float float32_t;

	/**
	Complex number consisting of two 32-bit floats.
	*/
	typedef std::complex<float32_t> complex32_t;
	
	/**
    Defines type for image coordinates.
    This type is signed integer.
    */
    typedef int64_t coord_t;

}
