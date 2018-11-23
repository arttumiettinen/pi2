#pragma once

#include "image.h"
#include "io/imagedatatype.h"
#include "io/itlpng.h"
#include "utilities.h"

#include <cmath>
#include <experimental/filesystem>

namespace fs = std::experimental::filesystem;

namespace itl2
{

	namespace sequence
	{
		namespace internals
		{
			/**
			Builds naturally sorted list of files that match the given template.
			The file name part of the template may contain wildcards and @ as discussed in matches(...) function above,
			but directory must not contain wildcards.
			If the template is a directory (no file name specified), all files in the directory are listed.
			*/
			vector<string> buildFileList(const string& templ);

			/**
			Separates directory and filename parts of a sequence template.
			*/
			void separatePathAndFileTemplate(const string& templ, fs::path& dir, string& fileTemplate);

			/**
			Returns metadata of a 2D image.
			*/
			inline bool getInfo2D(const string& filename, coord_t& width, coord_t& height, ImageDataType& dataType)
			{
				// TODO: Add other formats here.
				return png::getInfo(filename, width, height, dataType);
			}

			/**
			Reads a 2D image.
			*/
			template<typename pixel_t> void read2D(Image<pixel_t>& image, const string& filename, coord_t z)
			{
				// TODO: Add other formats here.
				png::read(image, filename, z);
			}

			/**
			Writes a 2D image.
			*/
			template<typename pixel_t> void write2D(const Image<pixel_t>& image, const string& filename, coord_t z)
			{
				// TODO: Add other formats here.
				png::write(image, filename, z);
			}
			
		}

		/**
		Get general information of an image sequence that is assumed to form a 3D image.
		@param filename Directory name and possibly file name template including wildcards *, ?, @ that identifies files that belong to the sequence. (@ corresponds to one or more numerical digits)
		@param width, height, depth Dimensions of the image
		@param dataType Pixel data type of the image.
		@return True if the sequence seems to be a valid image with supported pixel data type.
		*/
		inline bool getInfo(const string& filename, coord_t& width, coord_t& height, coord_t& depth, ImageDataType& dataType)
		{
			width = 0;
			height = 0;
			depth = 0;
			dataType = Unknown;

			vector<string> files = internals::buildFileList(filename);

			if (files.size() <= 0)
				return false;

			depth = files.size();
			return internals::getInfo2D(files[0], width, height, dataType);
		}

		inline bool getInfo(const string& filename, math::Vec3c& dims, ImageDataType& dataType)
		{
			return getInfo(filename, dims.x, dims.y, dims.z, dataType);
		}

		/**
		Tests if the given template corresponds to a valid image sequence.
		*/
		inline bool isSequence(const string& filename)
		{
			math::Vec3c dims;
			ImageDataType dt;
			return getInfo(filename, dims, dt);
		}

		/**
		Read an image sequence to a 3D image.
		@param img Image where the data will be placed. The size of the image is set automatically. NOTE: Data type of the image must be correct - use getInfo to determine that.
		@param filename Directory name and possibly file name template including wildcards *, ?, @ that identifies files that belong to the sequence. (@ corresponds to one or more numerical digits)
		@param firstSlice Index of first slice to be read.
		@param lastSlice Index of last slice to be read. Defaults to all slices. The value is clamped to actual count of slices.
		*/
		template<typename pixel_t> void read(Image<pixel_t>& img, const string& filename, size_t firstSlice = 0, size_t lastSlice = numeric_limits<size_t>::max())
		{
			vector<string> files = internals::buildFileList(filename);

			if (files.size() <= 0)
				throw ITLException("Sequence contains no images.");

			math::clamp<size_t>(firstSlice, 0, files.size() - 1);
			math::clamp<size_t>(lastSlice, 0, files.size() - 1);

			if (lastSlice < firstSlice)
				throw ITLException("No slices to be read, last slice index is less than first slice index.");

			coord_t depth = (size_t)(lastSlice - firstSlice + 1);
			coord_t width, height;
			ImageDataType dataType;

			if (!internals::getInfo2D(files[0], width, height, dataType))
				throw ITLException("Unable to get basic information of the image sequence. Possibly image type cannot be opened.");

			// Check that image data type matches
			if (dataType != imageDataType<pixel_t>())
				throw ITLException("The data types of the image and the sequence do not match.");

			img.ensureSize(width, height, depth);

			size_t counter = 0;
			#pragma omp parallel for if(!omp_in_parallel())
			for (coord_t z = (coord_t)firstSlice; z <= (coord_t)lastSlice; z++)
			{
				internals::read2D(img, files[z], z - firstSlice);
				showThreadProgress(counter, img.depth() + 1);
			}

		}

		/**
		Reads part of an image sequence to given image.
		Supports multiple computers accessing the same image at the same time only if the file system where the sequence is saved supports file locking correctly among all the computers accessing the data.
		NOTE: Does not support out-of-bounds start position.
		@param img Image where the data is placed. The size of the image defines the size of the block that is read.
		@param filename Directory name and possibly file name template including wildcards *, ?, @ that identifies files that belong to the sequence. (@ corresponds to one or more numerical digits)
		@param start Start location of the read. The size of the image defines the size of the block that is read.
		*/
		template<typename pixel_t> void readBlock(Image<pixel_t>& img, const std::string& filename, const math::Vec3c& start, bool showProgressInfo = false)
		{
			vector<string> files = internals::buildFileList(filename);

			if (files.size() <= 0)
				throw ITLException("Sequence contains no images.");

			math::Vec3c dimensions(0, 0, files.size());
			ImageDataType dataType;

			if (!internals::getInfo2D(files[0], dimensions.x, dimensions.y, dataType))
				throw ITLException("Unable to get basic information of the image sequence. Possibly image type cannot be opened.");

			// Check that image data type matches
			if (dataType != imageDataType<pixel_t>())
				throw ITLException("The data types of the image and the sequence do not match.");

			if (start.x < 0 || start.y < 0 || start.z < 0 || start.x >= dimensions.x || start.y >= dimensions.y || start.z >= dimensions.z)
				throw ITLException("Out of bounds start position in sequence::readBlock.");

			math::Vec3c cStart = start;
			clamp(cStart, math::Vec3c(0, 0, 0), dimensions);
			math::Vec3c cEnd = start + img.dimensions();
			clamp(cEnd, math::Vec3c(0, 0, 0), dimensions);

			size_t counter = 0;
			#pragma omp parallel for if(!omp_in_parallel())
			for (coord_t z = cStart.z; z < cEnd.z; z++)
			{
				// Read whole file
				Image<pixel_t> slice(dimensions.x, dimensions.y, 1);
				{
					//FileLock lock(files[z]);
					internals::read2D(slice, files[z], 0);
				}

				// Pick the part of the file we need
				for (coord_t y = cStart.y; y < cEnd.y; y++)
				{
					for (coord_t x = cStart.x; x < cEnd.x; x++)
					{
						pixel_t pix = slice(x, y, 0);
						img(x - cStart.x, y - cStart.y, z - cStart.z) = pix;
					}
				}

				showThreadProgress(counter, img.depth() + 1, showProgressInfo);
			}

		}

		namespace internals
		{
			inline void prepareTemplate(const string& filename, size_t sliceCount, fs::path& dir, string& fileTempl, int& fieldWidth, size_t& atPos)
			{
				internals::separatePathAndFileTemplate(filename, dir, fileTempl);

				if (fileTempl.find('@') >= fileTempl.size())
				{
					dir /= fileTempl;
					fileTempl = "";
				}

				if (fileTempl == "")
					fileTempl = "@.png";

				// Find out width of field.
				fieldWidth = 0;
				atPos = fileTempl.find("@(");
				if (atPos < fileTempl.size())
				{
					// We found @(, take all digits until )
					size_t endPos = fileTempl.find(')', atPos);
					string countStr = fileTempl.substr(atPos + 2, endPos - atPos - 2);
					if (countStr == "-")
					{
						fieldWidth = -1;
					}
					else
					{
						fieldWidth = fromString<int>(countStr);
					}

					// Remove count specification
					fileTempl.erase(atPos, endPos - atPos + 1);
				}
				else
				{
					atPos = fileTempl.find('@');
					fileTempl.erase(atPos, 1);
				}

				//if (atPos > fileTempl.size())
				//	throw ITLException("File name template does not contain letter @ that signifies the position where the slice number should be inserted.");

				// If field width is negative, set it automatically.
				if (fieldWidth < 0)
				{
					fieldWidth = (int)std::log10(sliceCount + 1) + 1;
				}
			}

			inline string insertSliceNumber(const string& fileTempl, int fieldWidth, size_t atPos, size_t z)
			{
				string zStr = toString(z);
				while (zStr.size() < fieldWidth)
					zStr = "0" + zStr;

				string filename = fileTempl;
				filename.insert(atPos, zStr);

				return filename;
			}
		}

		/**
		Write an image to multiple files, one file per z-slice.
		@param img Image to write.
		@param filename Filename template. Character @ is replaced by the z-coordinate of the slice.
		Sequence @(n) is replaced by the z-coordinate of the slice and leading zeroes so that total width of the numeric field is n.
		Specify @(-) to set field width automatically to a suitable value.
		Empty file name defaults to template "@.png" that corresponds to file names like 0.png, 1.png, 2.png, etc.
		If the file name template does not contain letter @, it is assumed to correspond to a folder and /@.png is added to the end of the template.
		@param outputFirstSlice Index that the slice corresponding to first output slice will get in the output file names. Typically 0 or 1 unless partial sequence is being saved.
		@param firstSlice Index of first slice to write.
		@param lastSlice Index of last slice to write. This value is clamped to the size of the image.
		*/
		template<typename pixel_t> void write(const Image<pixel_t>& img, const string& filename, coord_t outputFirstSlice = 0, size_t firstSlice = 0, size_t lastSlice = numeric_limits<size_t>::max())
		{
			math::clamp<size_t>(firstSlice, 0, img.depth() - 1);
			math::clamp<size_t>(lastSlice, 0, img.depth() - 1);

			if (lastSlice < firstSlice)
				throw ITLException("No slices to be written, last slice index is less than first slice index.");

			fs::path dir;
			string fileTempl;
			int fieldWidth;
			size_t atPos;
			internals::prepareTemplate(filename, lastSlice - firstSlice, dir, fileTempl, fieldWidth, atPos);

			// Write all slices
			fs::create_directories(dir);
			size_t counter = 0;
			#pragma omp parallel for if(!omp_in_parallel())
			for (coord_t z = firstSlice; z <= (coord_t)lastSlice; z++)
			{
				string filename = internals::insertSliceNumber(fileTempl, fieldWidth, atPos, z - firstSlice + outputFirstSlice);

				internals::write2D(img, (dir / filename).string(), z);

				showThreadProgress(counter, img.depth());
			}

		}

		/**
		Writes a block of image to specified location in a raw file.
		The output files are not truncated if they exist.
		If some output file does not exist, it is created.
		Part of image extending beyond [0, fileDimensions[ is not written.
		@param img Image to write.
		@param filename Filename template. Character @ is replaced by the z-coordinate of the slice.
		Sequence @(n) is replaced by the z-coordinate of the slice and leading zeroes so that total width of the numeric field is n.
		Specify @(-) to set field width automatically to a suitable value.
		Empty file name defaults to template "@.png" that corresponds to file names like 0.png, 1.png, 2.png, etc.
		If the file name template does not contain letter @, it is assumed to correspond to a folder and /@.png is added to the end of the template.
		@param filePosition Position in the file to write to.
		@param fileDimension Total dimensions of the output file.
		@param imagePosition Position in the image where the block to be written starts.
		@param imageDimensions Dimensions of the block of the source image to write.
		@param showProgressInfo Set to true to show a progress bar.
		*/
		template<typename pixel_t> void writeBlock(const Image<pixel_t>& img, const std::string& filename, const math::Vec3c& filePosition, const math::Vec3c& fileDimensions,
			const math::Vec3c& imagePosition, const math::Vec3c& imageDimensions,
			bool showProgressInfo = false)
		{
			constexpr size_t READ_WRITE_TRIALS = 5;

			fs::path dir;
			string fileTempl;
			int fieldWidth;
			size_t atPos;
			internals::prepareTemplate(filename, fileDimensions.z, dir, fileTempl, fieldWidth, atPos);

			math::Vec3c cStart = filePosition;
			clamp(cStart, math::Vec3c(0, 0, 0), fileDimensions);
			math::Vec3c cEnd = filePosition + imageDimensions;
			clamp(cEnd, math::Vec3c(0, 0, 0), fileDimensions);

			if (!img.isInImage(imagePosition))
				throw ITLException("Block start position must be inside the image.");
			if (!img.isInImage(imagePosition + imageDimensions - math::Vec3c(1, 1, 1)))
				throw ITLException("Block end position must be inside the image.");

			const pixel_t* pBuffer = img.getData();

			fs::create_directories(dir);
			#pragma omp parallel for if(!omp_in_parallel())
			for (coord_t z = cStart.z; z < cEnd.z; z++)
			{
				string filename = internals::insertSliceNumber(fileTempl, fieldWidth, atPos, z/* - cStart.z + imagePosition.z*/);
				filename = (dir / filename).string();

				Image<pixel_t> slice(fileDimensions.x, fileDimensions.y, 1);

				{
					//FileLock lock(filename);

					// Read existing data
					if (fs::exists(filename))
					{
						for (size_t n = 0; n < READ_WRITE_TRIALS; n++)
						{
							try
							{
								coord_t fw, fh;
								ImageDataType fd;
								if (internals::getInfo2D(filename, fw, fh, fd))
								{
									if (fw == fileDimensions.x && fh == fileDimensions.y)
									{
										internals::read2D(slice, filename, 0);
										break;
									}
								}
							}
							catch (ITLException&)
							{
								if (n >= READ_WRITE_TRIALS - 1)
									throw;
								sleep(100);
							}
						}
					}

					// Copy block
					for (coord_t y = cStart.y; y < cEnd.y; y++)
					{
						for (coord_t x = cStart.x; x < cEnd.x; x++)
						{
							slice(x, y, 0) = img(x - cStart.x + imagePosition.x, y - cStart.y + imagePosition.y, z - cStart.z + imagePosition.z);
						}
					}

					// Write slice back to disk
					for (size_t n = 0; n < READ_WRITE_TRIALS; n++)
					{
						try
						{
							internals::write2D(slice, filename, 0);
							break;
						}
						catch (ITLException&)
						{
							if (n >= READ_WRITE_TRIALS - 1)
								throw;
							sleep(100);
						}
					}
				}

				if (showProgressInfo)
					showProgress(z - cStart.z, cEnd.z - cStart.z);
			}
		}

		/**
		Writes an image to a specified location in a .raw file.
		The output file is not truncated if it exists.
		If the output file does not exist, it is created.
		Part of image extending beyond [0, fileDimensions[ is not written.
		@param img Image to write.
		@param filename Name of file to write.
		@param filePosition Position in the file to write to.
		@param fileDimension Total dimensions of the output file.
		@param showProgressInfo Set to true to show a progress bar.
		*/
		template<typename pixel_t> void writeBlock(const Image<pixel_t>& img, const std::string& filename, const math::Vec3c& filePosition, const math::Vec3c& fileDimensions, bool showProgressInfo = false)
		{
			writeBlock(img, filename, filePosition, fileDimensions, math::Vec3c(0, 0, 0), img.dimensions(), showProgressInfo);
		}

		/**
		Move or copy image sequence from one location to another.
		*/
		inline void moveOrCopySequence(const std::string& inputName, const std::string& outputName, bool move)
		{
			// TODO: Check that input and output sequences have the same file type.

			// Get list of input files
			vector<string> inputFiles = internals::buildFileList(inputName);

			if (inputFiles.size() <= 0)
				throw ITLException("Sequence contains no images.");

			// Prepare output file name template
			fs::path dir;
			string fileTempl;
			int fieldWidth;
			size_t atPos;
			internals::prepareTemplate(outputName, inputFiles.size(), dir, fileTempl, fieldWidth, atPos);

			fs::create_directories(dir);
			#pragma omp parallel for if(!omp_in_parallel())
			for (coord_t z = 0; z < (coord_t)inputFiles.size(); z++)
			{
				string infile = inputFiles[z];

				string outfile = internals::insertSliceNumber(fileTempl, fieldWidth, atPos, z);
				outfile = (dir / outfile).string();

				if (move)
					fs::rename(infile, outfile);
				else
					fs::copy(infile, outfile);
			}
		}

		inline void moveSequence(const std::string& inputName, const std::string& outputName)
		{
			moveOrCopySequence(inputName, outputName, true);
		}

		inline void copySequence(const std::string& inputName, const std::string& outputName)
		{
			moveOrCopySequence(inputName, outputName, false);
		}

		namespace tests
		{
			void match();
			void sequence();
			void readWriteBlock();
		}
	}

}
