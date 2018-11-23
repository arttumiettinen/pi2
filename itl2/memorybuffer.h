#pragma once

#include "buffer.h"

#include "fftw3.h"

namespace itl2
{

	/**
	Memory-resident buffer, compatible with fftw.
	*/
	template<typename pixel_t> class MemoryBuffer : public Buffer<pixel_t>
	{
	private:

		/**
		The memory buffer.
		*/
		pixel_t* pBuffer;

	public:

		/**
		Constructor
		@param size Count of T:s to allocate.
		*/
		MemoryBuffer(size_t size)
		{
			pBuffer = (pixel_t*)fftwf_malloc(size * sizeof(pixel_t));
			if (!pBuffer)
				throw ITLException("Out of memory.");
		}

		virtual ~MemoryBuffer()
		{
			fftwf_free(pBuffer);
		}

		virtual pixel_t* getBufferPointer()
		{
			return pBuffer;
		}

		virtual void prefetch(size_t start, size_t end) const
		{
			// Do nothing, there's nothing to prefetch.
		}
	};

}
