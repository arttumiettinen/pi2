#pragma once

#if defined(__linux__) || defined(__APPLE__)

#define PILIB_API __attribute__ ((visibility("default"))) extern

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
	@param pi Pi object created using createPI() function.
	*/
	PILIB_API void destroyPI(void* pi);

	/**
	Run commands.
	@param pi Pi object created using createPI() function.
	@param commands String containig one or more commands that should be executed.
	@return True if command execution was successful; false otherwise.
	*/
	PILIB_API uint8_t run(void* pi, const char* commands);

	/**
	Get error message identifying the last error that has occured.
	The returned pointer is valid until the next call run or clearLastError.
	@param pi Pi object created using createPI() function.
	*/
	PILIB_API const char* lastErrorMessage(void* pi);

	/**
	Get line of script code that caused the last error.
	@param pi Pi object created using createPI() function.
	*/
	PILIB_API int32_t lastErrorLine(void* pi);

	/**
	Clear last error.
	@param pi Pi object created using createPI() function.
	*/
	PILIB_API void clearLastError(void* pi);

	/**
	Gets list commands available in the system.
	The returned pointer is valid until the library is unloaded.
	@param pi Pi object created using createPI() function.
	*/
	PILIB_API const char* commandList(void* pi);

	/**
	Gets help text for given command.
	The returned pointer is valid until the next call of help.
	@param pi Pi object created using createPI() function.
	@param commandName Name of command whose help is to be returned.
	@param Help text.
	*/
	PILIB_API const char* help(void* pi, const char* commandName);

	/**
	Gets information of an image given its name.
	Does not read distributed images into RAM, and does not return pointer to image data.
	If the image does not exists, sets width, height, and depth to zero and data type to Unknown (zero).
	@param pi Pi object created using createPI() function.
	@param imgName Name of image.
	@param width, height, depth Dimensions of the image will be filled into these values.
	@param dataType Data type of the image will be filled into this value.
	*/
	PILIB_API void getImageInfo(void* pi, const char* imgName, int64_t* width, int64_t* height, int64_t* depth, int32_t* dataType);

	/**
	Gets pointer to data storing the given image.
	In distributed mode the image is read to RAM.
	Stores the size of the image into the values pointed by the three last arguments.
	If an error occurs, returns zero and sets width, height and depth to zero, and sets dataType to Unknown (zero).
	@param pi Pi object created using createPI() function.
	@param imgName Name of image.
	@param width, height, depth Dimensions of the image will be filled into these values.
	@param dataType The system sets this to a value describing pixel data type. See ImageDataType.
	*/
	PILIB_API void* getImage(void* pi, const char* imgName, int64_t* width, int64_t* height, int64_t* depth, int32_t* dataType);

	/**
	Gets value of a string object.
	The returned pointer is valid until the object is destroyed.
	@param pi Pi object created using createPI() function.
	@param name Name of object whose contents are to be returned.
	@param Pointer to the value or zero if no object with given name exists, or if the object is not of correct type.
	*/
	PILIB_API const char* getString(void* pi, const char* name);

	/**
	Gets value of a integer object.
	The returned pointer is valid until the object is destroyed.
	@param pi Pi object created using createPI() function.
	@param name Name of object whose contents are to be returned.
	@param The value, or minimum possible value if an object with the given name does not exist, or if the object is of wrong type.
	*/
	PILIB_API const int64_t getInt(void* pi, const char* name);

	/**
	Gets value of a real number object.
	The returned pointer is valid until the object is destroyed.
	@param pi Pi object created using createPI() function.
	@param name Name of object whose contents are to be returned.
	@param The value, or minimum possible value if an object with the given name does not exist, or if the object is of wrong type.
	*/
	PILIB_API const float32_t getReal(void* pi, const char* name);

	/**
	Gets value of a boolean object as 0 or 1.
	The returned pointer is valid until the object is destroyed.
	@param pi Pi object created using createPI() function.
	@param name Name of object whose contents are to be returned.
	@param The value (0 for false or 1 for true), or maximum possible value if an object with the given name does not exist, or if the object is of wrong type.
	*/
	PILIB_API const uint8_t getBool(void* pi, const char* name);

	/**
	In distributed computing mode, flushes changes made to image data through pointers returned by getData to files.
	Does nothing in normal mode.
	@param pi Pi object created using createPI() function.
	*/
	PILIB_API uint8_t finishUpdate(void* pi, const char* imgName);
}

