#pragma once

namespace itl2
{
	/**
	Base class for data storage buffers.
	*/
	template<typename pixel_t> class Buffer
	{
	public:
		/**
		* Virtual destructor
		*/
		virtual ~Buffer()
		{
		}

		/**
		Get pointer to the contents of the buffer.
		*/
		virtual pixel_t* getBufferPointer() = 0;

		/*
		Prefetch a range of buffer.
		Start and end are given as pixel indices relative to buffer start.
		*/
		virtual void prefetch(size_t start, size_t end) const = 0;
	};
}
