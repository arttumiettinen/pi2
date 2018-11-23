#pragma once

#include <string>

#include "buffer.h"
#include "itlexception.h"
#include "io/fileutils.h"

#if defined(__linux__)

	#include <sys/mman.h>
	#include <sys/types.h>
	#include <sys/stat.h>
	#include <fcntl.h>
	#include <unistd.h>

#elif defined(_WIN32)

	#define WIN32_LEAN_AND_MEAN
	#define NOMINMAX
	#include <Windows.h>

#else

	#error diskmappedbuffer.h not configured for this platform.

#endif


namespace itl2
{


#if defined(__linux__)

	/**
	Disk-mapped buffer.
	*/
	template<typename T> class DiskMappedBuffer : public Buffer<T>
	{
	private:

		/**
		Pointer to mapped memory.
		*/
		T* pBuffer;

		/**
		Pointer to mapped memory, passed to the user.
		*/
		T* pUserBuffer;

		/*
		Size of mapped region.
		*/
		size_t mappedSize;

		/**
		Dummy variable that will cause memory leak if the handles are not closed properly.
		*/
		int* pDummy;

		/**
		Closes all handles.
		*/
		void close()
		{
			munmap(pBuffer, mappedSize);

			delete pDummy;
			pDummy = 0;

			pUserBuffer = 0;
		}

	public:

		/**
		Constructor
		@param size Size to map, pass zero to map the whole file.
		@param filename Name of file to map to.
		*/
		DiskMappedBuffer(size_t size, const string& filename, size_t bytesToSkip = 0) :
			pDummy(0),
			mappedSize(size * sizeof(T))
		{
			createFoldersFor(filename);

			// Open file (if the file is created, user has read and write permissions, group has read permissions.
			int fd = open(filename.c_str(), O_RDWR | O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP);
			if (fd == -1)
				throw ITLException(std::string("Unable to open file ") + filename);

			// Resize
			if (ftruncate(fd, mappedSize) == -1)
			{
				::close(fd);
				throw ITLException(std::string("Unable to set size of ") + filename);
			}

			// Create memory map
			pBuffer = (T*)mmap(0, mappedSize, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
			if (pBuffer == MAP_FAILED)
			{
				::close(fd);
				throw ITLException(std::string("Unable to create view for file ") + filename);
			}

			::close(fd);

			pUserBuffer = (T*)((uint8_t*)(pBuffer)+bytesToSkip);

			pDummy = new int();
		}

		/**
		Destructor
		*/
		virtual ~DiskMappedBuffer()
		{
			close();
		}

		virtual T* getBufferPointer()
		{
			return pUserBuffer;
		}

		virtual void prefetch(size_t start, size_t end) const
		{
			madvise(pUserBuffer + start, end - start, MADV_WILLNEED);
		}
	};

#elif defined(_WIN32)

	/**
	Disk-mapped buffer.
	*/
	template<typename T> class DiskMappedBuffer : public Buffer<T>
	{
	private:

		/**
		Pointer to mapped memory.
		*/
		T* pBuffer;

		/**
		Pointer to mapped memory, passed to the user.
		*/
		T* pUserBuffer;

		/**
		The required handles.
		*/
		HANDLE fileHandle;
		HANDLE fileMappingHandle;

		/**
		Dummy variable that will cause memory leak if the handles are not closed properly.
		*/
		int* pDummy;

		/**
		Closes all handles.
		*/
		void close()
		{
			if(fileMappingHandle != INVALID_HANDLE_VALUE)
			{
				CloseHandle(fileMappingHandle);
				fileMappingHandle = INVALID_HANDLE_VALUE;
			}

			if(fileHandle != INVALID_HANDLE_VALUE)
			{
				CloseHandle(fileHandle);
				fileHandle = INVALID_HANDLE_VALUE;
			}

			delete pDummy;
			pDummy = 0;

			pUserBuffer = 0;
		}

	public:

		/**
		Constructor
		Resizes input file to specified size and maps it to memory.
		@param size Size to map, pass zero to map the whole file.
		@param filename Name of file to map to.
		*/
		DiskMappedBuffer(size_t size, const std::string& filename, size_t bytesToSkip = 0) :
		  fileHandle(INVALID_HANDLE_VALUE),
		  fileMappingHandle(INVALID_HANDLE_VALUE),
		  pDummy(0)
		{
			try
			{
				// TODO: Add GetLastError messages to exception text.

				createFoldersFor(filename);

				fileHandle = CreateFileA(filename.c_str(),					// filename
										GENERIC_READ | GENERIC_WRITE,		// access mode
										FILE_SHARE_READ | FILE_SHARE_WRITE,	// share mode
										NULL,								// security attributes
										OPEN_ALWAYS,						// creation/opening mode
										FILE_ATTRIBUTE_NORMAL,				// attributes
										NULL);								// template file

				if(fileHandle == INVALID_HANDLE_VALUE)
					throw ITLException("Unable to open file " + filename + " for disk mapping.");

				LARGE_INTEGER datasize;
				datasize.QuadPart = size * sizeof(T);

				// Resize file to requested size
				SetFilePointer(fileHandle, datasize.LowPart, &datasize.HighPart, FILE_BEGIN);
				SetEndOfFile(fileHandle);
				SetFilePointer(fileHandle, 0, 0, FILE_BEGIN);

				// Create mapping
				fileMappingHandle = CreateFileMappingA(fileHandle,			// file handle
													  NULL,					// security attributes
													  PAGE_READWRITE,		// mapping attributes
													  datasize.HighPart,	// maximum size h-dword
													  datasize.LowPart,		// maximum size l-dword
													  NULL);				// name

				if(fileMappingHandle == NULL)	// NULL, not INVALID_HANDLE_VALUE
					throw ITLException("Unable to open file mapping for " + filename);

				pBuffer = (T*)MapViewOfFile(fileMappingHandle,						// file mapping
											FILE_MAP_ALL_ACCESS,					// access
											0,										// offset h-dword
											0,										// offset l-dword
											0);										// size of view, 0-whole mapping

				if(pBuffer == NULL)
					throw ITLException("Unable to create view for file " + filename);

                pUserBuffer = (T*)((uint8_t*)(pBuffer) + bytesToSkip);

				pDummy = new int();
			}
			catch(...)
			{
				close();
				throw;
			}

		}

		/**
		Destructor
		*/
		virtual ~DiskMappedBuffer()
		{
			close();
		}

		virtual T* getBufferPointer()
		{
			return pUserBuffer;
		}

		virtual void prefetch(size_t start, size_t end) const
		{
			// TODO: This does not work in Windows 7 (no support => runtime error)
			// Enable when no Win7 support is required anymore.
			//WIN32_MEMORY_RANGE_ENTRY entry;
			//entry.VirtualAddress = pUserBuffer + start;
			//entry.NumberOfBytes = end - start;
			//PrefetchVirtualMemory(GetCurrentProcess(), 1, &entry, 0);
		}
	};

#else
	#error Disk mapped buffer not configured for this platform.
#endif
}
