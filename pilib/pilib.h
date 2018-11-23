#pragma once

#if defined(__linux__)

#define PILIB_API extern

#elif defined(_WIN32)

// The following ifdef block is the standard way of creating macros which make exporting 
// from a DLL simpler. All files within this DLL are compiled with the PILIB_EXPORTS
// symbol defined on the command line. This symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see 
// PILIB_API functions as being imported from a DLL, whereas this DLL sees symbols
// defined with this macro as being exported.
#ifdef PILIB_EXPORTS
#define PILIB_API __declspec(dllexport)
#else
#define PILIB_API __declspec(dllimport)
#endif

#else
#error pilib.h not configured for this platform.
#endif

#include "datatypes.h"

using namespace itl2;

extern "C"
{
	/**
	Create the pi library object.
	@return Handle that must be passed to all other methods.
	*/
	PILIB_API void* createPI();

	/**
	Close the pi library object.
	*/
	PILIB_API void destroyPI(void* pi);

	/**
	Run commands.
	@param pi Handle returned by create method.
	@param commands String containig one or more commands that should be executed.
	@return True if command execution was successful; false otherwise.
	*/
	PILIB_API uint8_t run(void* pi, const char* commands);

	/**
	Get error message identifying the last error that has occured.
	The returned pointer is valid until the next call run or clearLastError.
	*/
	PILIB_API const char* lastErrorMessage(void* pi);

	/**
	Get line of script code that caused the last error.
	*/
	PILIB_API int32_t lastErrorLine(void* pi);

	/**
	Clear last error.
	*/
	PILIB_API void clearLastError(void* pi);

	/**
	Gets list commands available in the system.
	The returned pointer is valid until the library is unloaded.
	*/
	PILIB_API const char* commandList(void* pi);

	/**
	Gets help text for given command.
	The returned pointer is valid until the next call of help.
	*/
	PILIB_API const char* help(void* pi, const char* commandName);

	/**
	Gets pointer to data storing the given image.
	Stores the size of the image into the values pointed by the three last arguments.
	If an error occurs, returns zero and sets width, height and depth to zero, and sets dataType to Unknown (zero).
	@param dataType The system sets this int to 1 to signify uint8 image, 2 for uint16 image, 3 for float32 image and 4 for complex32 image.
	*/
	PILIB_API void* getImage(void* pi, const char* imgName, int64_t* width, int64_t* height, int64_t* depth, int32_t* dataType);
}

