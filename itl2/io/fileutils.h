#pragma once

#include <string>
#include <fstream>

namespace itl2
{
	/**
	Returns size of given file.
	*/
	std::ifstream::pos_type fileSize(const std::string& filename);

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