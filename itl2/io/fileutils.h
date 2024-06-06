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

	/**
	Adds image dimensions and ".raw" suffix to file name.
	*/
	std::string concatDimensions(const std::string& baseName, const Vec3c& dimensions);


	/**
	Removes .raw image dimensions from file name.
	*/
	std::string getPrefix(std::string filename);


    /**
	Builds naturally sorted list of files that match the given template.
	The file name part of the template may contain wildcards and @ as discussed in matches(...) function above,
	but directory must not contain wildcards.
	If the template is a directory (no file name specified), all files in the directory are listed.
	*/
	std::vector<std::string> buildFileList(const std::string& templ);

	/**
	Works as buildFileList but removed non-image files from the list.
	*/
	std::vector<std::string> buildFilteredFileList(const std::string& templ);

    /**
	Tests if given str matches template.
	In the template
	* matches to sequence of any characters
	? matches to any single character
	@ matches to one or more numerical digits
	*/
	bool matches(const std::string& str, const std::string& templ);

    /**
	Separates directory and filename parts of a sequence template.
	*/
	void separatePathAndFileTemplate(const std::string& templ, fs::path& dir, std::string& fileTemplate);

	/**
	Finds out the host name of the current system.
	*/
	std::string getHostname();
}
