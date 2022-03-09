#pragma once

#include <string>
#include <fstream>

#include "math/vec3.h"
#include "filesystem.h"

namespace itl2
{
	/**
	Returns size of given file in bytes.
	*/
	std::ifstream::pos_type fileSize(const std::string& filename);

	/**
	Gets error message corresponding to the current errno variable.
	*/
	const std::string getStreamErrorMessage();

	/**
	Sets size of file to given value.
	*/
	void setFileSize(const std::string& filename, size_t size);

	/**
	Tests if a file exists.
	*/
	bool fileExists(const std::string& filename);
	bool fileExists(const fs::path& filename);

	/**
	Deletes a file.
	*/
	void deleteFile(const std::string& filename);

	/**
	Copies a file.
	The destination file is overwritten.
	*/
	void copyFile(const std::string& sourceName, const std::string& destinationName, bool showProgressInfo);

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

	/**
	Adds image dimensions and ".raw" suffix to file name.
	*/
	std::string concatDimensions(const std::string& baseName, const Vec3c& dimensions);


	/**
	Removes .raw image dimensions from file name.
	*/
	std::string getPrefix(std::string filename);

}
