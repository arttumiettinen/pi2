#pragma once

#include <string>
#include <fstream>
#include <cstring>

namespace itl2
{
	/**
	Returns size of given file.
	*/
	std::ifstream::pos_type fileSize(const std::string& filename);

	/**
	Gets error message corresponding to the current errno variable.
	*/
	inline const std::string getStreamErrorMessage()
	{
#if defined(_WIN32)
		char buf[1024];
		strerror_s(&buf[0], 1024, errno);
		return std::string(&buf[0]);
#else
		return std::string(strerror(errno));
#endif
	}

	/**
	Sets size of file to given value.
	*/
	void setFileSize(const std::string& filename, size_t size);

	/**
	Tests if a file exists.
	*/
	bool fileExists(const std::string& filename);

	/**
	Deletes a file.
	*/
	void deleteFile(const std::string& filename);

	/**
	Copies a file.
	The destination file is overwritten.
	*/
	void copyFile(const std::string& sourceName, const std::string& destinationName);

	/**
	Moves a file.
	The destination file is overwritten.
	*/
	void moveFile(const std::string& sourceName, const std::string& destinationName);

	/**
	Creates folders in the path of a given file if they don't exist.
	Does not create the file, only directories.
	*/
	void createFoldersFor(const std::string& filename);

}
