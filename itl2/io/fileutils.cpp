
#include "fileutils.h"
#include "itlexception.h"
#include "utilities.h"
#include "stringutils.h"
#include "io/alphanum.h"

#include <cstring>

#include "filesystem.h"

#if defined(__linux__) || defined(__APPLE__)

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#elif defined(_WIN32)

#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <Windows.h>

#else

#error fileutils.cpp not configured for this platform.

#endif

using namespace std;

namespace itl2
{
	const std::string getStreamErrorMessage()
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
	Returns size of given file in bytes.
	*/
	std::ifstream::pos_type fileSize(const string& filename)
	{
		ifstream in(filename.c_str(), ifstream::ate | ifstream::binary);
		return in.tellg();
	}

	/**
	Sets size of file to given value.
	*/
	void setFileSize(const std::string& filename, size_t size)
	{
#if defined(__linux__) || defined(__APPLE__)

		// Open file (if the file is created, user has read and write permissions, group has read permissions.
		int fd = open(filename.c_str(), O_RDWR | O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP);
		if (fd == -1)
			throw ITLException(std::string("Unable to open file ") + filename);

		// Resize
		if (ftruncate(fd, size) == -1)
		{
			::close(fd);
			throw ITLException(std::string("Unable to set size of ") + filename);
		}

		::close(fd);

#elif defined(_WIN32)

		HANDLE fileHandle = CreateFileA(filename.c_str(),					// filename
			GENERIC_READ | GENERIC_WRITE,		// access mode
			FILE_SHARE_READ | FILE_SHARE_WRITE,	// share mode
			NULL,								// security attributes
			OPEN_ALWAYS,						// creation/opening mode
			FILE_ATTRIBUTE_NORMAL,				// attributes
			NULL);								// template file

		if (fileHandle == INVALID_HANDLE_VALUE)
			throw ITLException("Unable to open file " + filename);

		LARGE_INTEGER datasize;
		datasize.QuadPart = size;

		// Resize file to requested size
		SetFilePointer(fileHandle, datasize.LowPart, &datasize.HighPart, FILE_BEGIN);
		SetEndOfFile(fileHandle);
		//SetFilePointer(fileHandle, 0, 0, FILE_BEGIN);
		CloseHandle(fileHandle);
#else

#error fileutils.cpp not configured for this platform.

#endif
	}

	/**
	Tests if a file exists.
	*/
	bool fileExists(const std::string& filename)
	{
		//std::ifstream infile(filename);
		//return (bool)infile;
		return fs::exists(filename);
	}

	bool fileExists(const fs::path& filename)
	{
		return fileExists(filename.string());
	}

	/**
	Deletes a file.
	*/
	void deleteFile(const std::string& filename)
	{
		remove(filename.c_str());
	}

#if defined(_WIN32)
	DWORD progressRoutine(
		LARGE_INTEGER TotalFileSize,
		LARGE_INTEGER TotalBytesTransferred,
		LARGE_INTEGER StreamSize,
		LARGE_INTEGER StreamBytesTransferred,
		DWORD dwStreamNumber,
		DWORD dwCallbackReason,
		HANDLE hSourceFile,
		HANDLE hDestinationFile,
		LPVOID lpData
	)
	{
		//showProgress(TotalBytesTransferred.QuadPart, TotalFileSize.QuadPart + 1);
		cout << bytesToString((double)TotalBytesTransferred.QuadPart) << " / " << bytesToString((double)TotalFileSize.QuadPart) << "\r";
		return PROGRESS_CONTINUE;
	};
#endif

	/**
	Copies a file.
	The destination file is overwritten.
	*/
	void copyFile(const std::string& sourceName, const std::string& destinationName, bool showProgressInfo)
	{
		fs::path p1(sourceName);
		fs::path p2(destinationName);
		std::error_code code; // This is used to choose the nothrow overload of fs::equivalent
		if (fs::equivalent(p1, p2, code))
			return;

		if (fileExists(sourceName))
			deleteFile(destinationName);

#if defined(__linux__) || defined(__APPLE__)

        if(system(0) == 0)
            throw ITLException("Command processor is not available.");

		// TODO: This does not show progress bar
		if(system((string("cp \"") + sourceName + string("\" \"") + destinationName + string("\"")).c_str()) != 0)
            throw ITLException(string("Unable to copy ") + sourceName + " to " + destinationName);

#elif defined(_WIN32)

		BOOL cancel = FALSE;
		LPPROGRESS_ROUTINE progress = nullptr;
		if (showProgressInfo)
			progress = &progressRoutine;
		
		if (CopyFileExA(sourceName.c_str(), destinationName.c_str(), progress, NULL, &cancel, 0) == 0)
			throw ITLException(string("Unable to copy ") + sourceName + " to " + destinationName);

		if(showProgressInfo)
			cout << "                          \r" << flush;
#else

#error fileutils.cpp not configured for this platform.

#endif
	}

	/**
	Creates folders in the path of a given file if they don't exist.
	Does not create the file, only directories.
	*/
	void createFoldersFor(const std::string& filename)
	{
		fs::path folder(filename);
		folder = folder.remove_filename();
		if (folder != filename && folder != "")
			fs::create_directories(folder);
	}

	/**
	Moves a file.
	The destination file is overwritten.
	*/
	void moveFile(const std::string& sourceName, const std::string& destinationName)
	{
		fs::path p1(sourceName);
		fs::path p2(destinationName);
		std::error_code code; // This is used to choose the nothrow overload of fs::equivalent
		if (fs::equivalent(p1, p2, code))
			return;

		createFoldersFor(destinationName);

#if defined(__linux__) || defined(__APPLE__)

        if(system(0) == 0)
            throw ITLException("Command processor is not available.");

		if(system((string("mv \"") + sourceName + string("\" \"") + destinationName + string("\"")).c_str()) != 0)
            throw ITLException(string("Unable to move ") + sourceName + " to " + destinationName);

#elif defined(_WIN32)
		if (MoveFileWithProgressA(sourceName.c_str(), destinationName.c_str(), &progressRoutine, NULL, MOVEFILE_COPY_ALLOWED | MOVEFILE_REPLACE_EXISTING) == 0)
			throw ITLException(string("Unable to move ") + sourceName + " to " + destinationName);
#else

#error fileutils.cpp not configured for this platform.

#endif
	}

	std::string concatDimensions(const std::string& baseName, const Vec3c& dimensions)
	{
		std::stringstream suffix;
		suffix << "_" << dimensions.x << "x" << dimensions.y << "x" << dimensions.z << ".raw";

		if (endsWithIgnoreCase(baseName, suffix.str()))
			return baseName;

		std::stringstream name;
		name << baseName << suffix.str();
		return name.str();
	}

	std::string getPrefix(std::string filename)
	{
		size_t pos = filename.find_last_of('_');
		if (pos != std::string::npos)
			filename.erase(filename.begin() + pos, filename.end());
		return filename;
	}


	bool matches(const string& str, const string& templ)
	{
		size_t strPos = 0;
		size_t templPos = 0;
		while (true)
		{
			if (strPos >= str.length())
				return templPos == templ.length();

			if (templPos >= templ.length())
				return false;

			if (templ[templPos] == '?')
			{
				// Any character is ok, no need to test.
				templPos++;
				strPos++;
			}
			else if (templ[templPos] == '*')
			{
				if (templPos >= templ.length() - 1)
					return true;

				for (size_t n = strPos; n < str.length(); n++)
					if (matches(str.substr(n), templ.substr(templPos + 1)))
						return true;

				return false;
			}
			else if (templ[templPos] == '@')
			{
				if (!isdigit(str[strPos]))
					return false;

				while(strPos < str.length() && isdigit(str[strPos]))
					strPos++;

				templPos++;
			}
			else
			{
				// Characters must be equal.
				if (str[strPos] != templ[templPos])
					return false;
				strPos++;
				templPos++;
			}
		}
	}

    /**
	Separates directory and filename parts of a sequence template.
	*/
	void separatePathAndFileTemplate(const string& templ, fs::path& dir, string& fileTemplate)
	{
		fs::path p(templ);

		if (fs::is_directory(p))
		{
			dir = p;
			fileTemplate = "";
		}
		else
		{
			dir = p.parent_path();
			fileTemplate = p.filename().string();
			if (fileTemplate == ".")
				fileTemplate = "";
		}
	}


    vector<string> buildFileList(const string& templ)
	{
		// Separate directory and file name template
		string fileTemplate;
		fs::path dir;
		separatePathAndFileTemplate(templ, dir, fileTemplate);

		if (dir == "")
			dir = ".";

		if(fileTemplate == "")
			fileTemplate = "*";
		
		//cout << "Directory: " << dir << endl;
		//cout << "Template: " << fileTemplate << endl;

		// Get those files in directory that match the template
		vector<string> filenames;

		if (fs::is_directory(dir)) // Note: This is required in Linux, or otherwise we get an exception for non-existing directories.
		{
			for (auto & p : fs::directory_iterator(dir))
			{
				if (p.is_regular_file())
				{
					string filename = p.path().filename().string();
					if (matches(filename, fileTemplate))
						filenames.push_back(p.path().string());
				}
			}

			// Sort to natural order
			sort(filenames.begin(), filenames.end(), doj::alphanum_less<std::string>());
		}

		//for (size_t n = 0; n < filenames.size(); n++)
		//	cout << filenames[n] << endl;
		return filenames;
	}

	vector<string> buildFilteredFileList(const string& templ)
	{
		vector<string> results = buildFileList(templ);

		// Remove non-images
		//coord_t w, h;
		//ImageDataType dt;
		for (size_t n = 0; n < results.size(); n++)
		{
			// This is too slow!
			//if (!getInfo2D(result[n], w, h, dt))
			fs::path p(results[n]);
			fs::path ext = p.extension();
			string exts = ext.string();
			toLower(exts);

			// TODO: Add other formats here.
			if(ext != ".tif" && ext != ".tiff" && ext != ".png")
			{
				results.erase(results.begin() + n);
				n--;
			}
		}

		return results;
	}

}
