
#include "diskmappedbuffer.h"

#if defined(__linux__) || defined(__APPLE__)

#include <errno.h>

#endif


namespace itl2
{

#if defined(__linux__) || defined(__APPLE__)

    std::string getErrNoMessage()
	{
		std::string msg;
		switch (errno)
		{
		case EACCES: msg = "A file descriptor refers to a non - regular file. Or a file	mapping was requested, but fd is not open for reading.Or MAP_SHARED was requested and PROT_WRITE is set, but fd is not open in read / write(O_RDWR) mode.Or PROT_WRITE is set, but the file is append - only."; break;
		case EAGAIN: msg = "The file has been locked, or too much memory has been locked (see setrlimit(2))."; break;
		case EBADF: msg = "fd is not a valid file descriptor(and MAP_ANONYMOUS was not set)."; break;
		case EEXIST: msg = "MAP_FIXED_NOREPLACE was specified in flags, and the range covered by addr and length clashes with an existing mapping."; break;
		case EINVAL: msg = "We don't like addr, length, or offset (e.g., they are too large, or not aligned on a page boundary). Maybe length was 0 or flags contained none of MAP_PRIVATE, MAP_SHARED or MAP_SHARED_VALIDATE."; break;
		case ENFILE: msg = "The system-wide limit on the total number of open files has been reached."; break;
		case ENODEV: msg = "The underlying filesystem of the specified file does not support memory mapping."; break;
		case ENOMEM: msg = "No memory is available. The process's maximum number of mappings would have been exceeded. This error can also occur for munmap(), when unmapping a region in the middle of an existing mapping, since this results in two smaller mappings on either side of the region being unmapped. The process's RLIMIT_DATA limit, described in getrlimit(2), would have been exceeded."; break;
		case EOVERFLOW: msg = "On 32 - bit architecture together with the large file extension (i.e., using 64 - bit off_t) : the number of pages used for length plus number of pages used for offset would overflow unsigned long(32 bits)."; break;
		case EPERM: msg = "The prot argument asks for PROT_EXEC but the mapped area belongs to a file on a filesystem that was mounted no-exec. The operation was prevented by a file seal; see fcntl(2)."; break;
		case ETXTBSY: msg = "MAP_DENYWRITE was set but the object specified by fd is open for writing."; break;
		default: msg = std::string("Unknown error: ") + itl2::toString(errno); break;
		}
		return msg;
	}
	
#endif

}
