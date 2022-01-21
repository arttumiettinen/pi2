#pragma once

#include <string>

#include "buffer.h"
#include "itlexception.h"
#include "io/fileutils.h"
#include "utilities.h"

#if defined(__linux__) || defined(__APPLE__)

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

#if defined(__linux__) || defined(__APPLE__)

    /**
	Reads errno and converts it to error message.
	*/
    std::string getErrNoMessage();

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

			// TODO: Set modification and access times.
		}

		/**
		Reads file size from file descriptor.
		*/
		size_t getFileSize(int fd)
		{
			struct stat s;
			if (fstat(fd, &s) == -1)
				throw ITLException(std::string("Unable to read file size. ") + " (" + getErrNoMessage() + ")");
			return s.st_size;
		}

	public:

		DiskMappedBuffer(const DiskMappedBuffer&) = delete;
		DiskMappedBuffer& operator=(DiskMappedBuffer const&) = delete;

		/**
		Constructor
		If readOnly is false, resizes input file to specified size and maps it to memory.
		Otherwise just maps the file to memory.
		@param size Size to map, pass zero to map the whole file.
		@param filename Name of file to map to.
		@param readOnly Specify true to map a read-only view. Writing to read-only view leads to undefined behaviour.
		*/
		DiskMappedBuffer(size_t size, const std::string& filename, bool readOnly, size_t bytesToSkip = 0) :
			pDummy(0),
			mappedSize(size * sizeof(T))
		{
    		int fd;
		    if(!readOnly)
		    {
		
    			createFoldersFor(filename);

			    // Open file (if the file is created, user has read and write permissions, group has read permissions.
			    fd = open(filename.c_str(), O_RDWR | O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP);
			    
			    if (fd == -1)
				    throw ITLException(std::string("Unable to open file for reading and writing: ") + filename + " (" + getErrNoMessage() + ")");
			
				if (mappedSize <= 0)
					mappedSize = getFileSize(fd);

			    // Resize to correct size
    		    if (ftruncate(fd, mappedSize) == -1)
			    {
			    	::close(fd);
			    	throw ITLException(std::string("Unable to set size of ") + filename + " (" + getErrNoMessage() + ")");
			    }

			    // Create memory map
			    pBuffer = (T*)mmap(0, mappedSize, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
			}
			else
			{
			    // Open read-only
			    fd = open(filename.c_str(), O_RDONLY, S_IRUSR | S_IWUSR | S_IRGRP);
			    
			    if (fd == -1)
				    throw ITLException(std::string("Unable to open file for reading: ") + filename + " (" + getErrNoMessage() + ")");
				
				if (mappedSize <= 0)
					mappedSize = getFileSize(fd);

				// Create read-only mapping
    			pBuffer = (T*)mmap(0, mappedSize, PROT_READ, MAP_SHARED, fd, 0);
    		}
			
			if (pBuffer == MAP_FAILED)
			{
				::close(fd);
				throw ITLException(std::string("Unable to create view for file ") + filename + ". " + getErrNoMessage() + " Size = " + itl2::toString(mappedSize) + ".");
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

		virtual T* getBufferPointer() override
		{
			return pUserBuffer;
		}

		virtual void prefetch(size_t start, size_t end) const override
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
		Indicates if the file mapping is read-only.
		*/
		bool readOnly;

		/**
		Closes all handles.
		*/
		void close()
		{
			pUserBuffer = 0;

			if (pBuffer)
			{
				UnmapViewOfFile(pBuffer);
				pBuffer = 0;
			}

			if(fileMappingHandle != INVALID_HANDLE_VALUE)
			{
				CloseHandle(fileMappingHandle);
				fileMappingHandle = INVALID_HANDLE_VALUE;
			}

			if(fileHandle != INVALID_HANDLE_VALUE)
			{
				FILETIME now;
				SYSTEMTIME st;
				GetSystemTime(&st);
				SystemTimeToFileTime(&st, &now);

				// Set access time
				SetFileTime(fileHandle, NULL, &now, NULL);
				if (!readOnly)
				{
					// Set modification time
					SetFileTime(fileHandle, NULL, NULL, &now);
				}

				CloseHandle(fileHandle);
				fileHandle = INVALID_HANDLE_VALUE;
			}

			delete pDummy;
			pDummy = 0;
		}

	public:

		DiskMappedBuffer(const DiskMappedBuffer&) = delete;
		DiskMappedBuffer& operator=(DiskMappedBuffer const&) = delete;

		/**
		Constructor
		If readOnly is false, resizes input file to specified size and maps it to memory.
		Otherwise just maps the file to memory.
		@param size Size to map, pass zero to map the whole file.
		@param filename Name of file to map to.
		@param readOnly Specify true to map a read-only view. Writing to read-only view leads to undefined behaviour.
		*/
		DiskMappedBuffer(size_t size, const std::string& filename, bool readOnly, size_t bytesToSkip = 0) :
			pBuffer(0),
			pUserBuffer(0),
			fileHandle(INVALID_HANDLE_VALUE),
			fileMappingHandle(INVALID_HANDLE_VALUE),
			pDummy(0),
			readOnly(readOnly)
		{
			try
			{
				// TODO: Add GetLastError messages to exception text.

				if (!readOnly)
				{
					createFoldersFor(filename);

					// Copy string to wstring.
					//std::wstring wfilename(filename.length(), L' '); // Make room for characters
					//std::copy(filename.begin(), filename.end(), wfilename.begin());

					fileHandle = CreateFileA(filename.c_str(),					// filename
						GENERIC_READ | GENERIC_WRITE | FILE_WRITE_ATTRIBUTES,	// access mode
						FILE_SHARE_READ | FILE_SHARE_WRITE,	// share mode
						NULL,								// security attributes
						OPEN_ALWAYS,						// creation/opening mode
						FILE_ATTRIBUTE_NORMAL,				// attributes
						NULL);								// template file

					if (fileHandle == INVALID_HANDLE_VALUE)
						throw ITLException("Unable to open file " + filename + " for disk mapping.");

					//// Get file size
					//LARGE_INTEGER filesize;
					//if(!GetFileSizeEx(fileHandle, &filesize))
					//	throw ITLException("Unable to get file size " + filename);

					//// Make sure file size is not zero.
					//if(filesize.QuadPart == 0)
					//{
					//	T buffer = T();
					//	DWORD written;
					//	WriteFile(fileHandle, &buffer, sizeof(T), &written, NULL);
					//	SetFilePointer(fileHandle, 0, 0, FILE_BEGIN);
					//}

					LARGE_INTEGER datasize;
					datasize.QuadPart = size * sizeof(T);

					// Resize file to requested size
					SetFilePointer(fileHandle, datasize.LowPart, &datasize.HighPart, FILE_BEGIN);
					SetEndOfFile(fileHandle);
					SetFilePointer(fileHandle, 0, 0, FILE_BEGIN);

					// Create mapping
					fileMappingHandle = CreateFileMappingA(fileHandle,			// file handle
						NULL,					// security attributes
						PAGE_READWRITE,			// mapping attributes
						datasize.HighPart,		// maximum size h-dword
						datasize.LowPart,		// maximum size l-dword
						NULL);					// name

					if (fileMappingHandle == NULL)	// NULL, not INVALID_HANDLE_VALUE
						throw ITLException("Unable to open file mapping for " + filename + ". Is there enough free disk space available?");

					pBuffer = (T*)MapViewOfFile(fileMappingHandle,						// file mapping
						FILE_MAP_ALL_ACCESS,					// access
						0,										// offset h-dword
						0,										// offset l-dword
						0);										// size of view, 0-whole mapping
				}
				else
				{
					fileHandle = CreateFileA(filename.c_str(),					// filename
						GENERIC_READ,						// access mode
						FILE_SHARE_READ | FILE_SHARE_WRITE,	// share mode
						NULL,								// security attributes
						OPEN_ALWAYS,						// creation/opening mode
						FILE_ATTRIBUTE_NORMAL,				// attributes
						NULL);								// template file

					if (fileHandle == INVALID_HANDLE_VALUE)
						throw ITLException("Unable to open file " + filename + " for read-only disk mapping.");

					LARGE_INTEGER datasize;
					datasize.QuadPart = size * sizeof(T);

					// Create mapping
					fileMappingHandle = CreateFileMappingA(fileHandle,			// file handle
						NULL,					// security attributes
						PAGE_READONLY,			// mapping attributes
						datasize.HighPart,		// maximum size h-dword
						datasize.LowPart,		// maximum size l-dword
						NULL);					// name

					if (fileMappingHandle == NULL)	// NULL, not INVALID_HANDLE_VALUE
						throw ITLException("Unable to open read-only file mapping for " + filename + ".");

					pBuffer = (T*)MapViewOfFile(fileMappingHandle,						// file mapping
						FILE_MAP_READ,							// access
						0,										// offset h-dword
						0,										// offset l-dword
						0);										// size of view, 0-whole mapping
				}

				if(pBuffer == NULL)
					throw ITLException(string("Unable to create view for file ") + filename);

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

		virtual T* getBufferPointer() override
		{
			return pUserBuffer;
		}

		virtual void prefetch(size_t start, size_t end) const override
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
