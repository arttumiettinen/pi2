#pragma once

#include "image.h"
#include "io/imagedatatype.h"
#include "io/itlpng.h"
#include "io/itltiff.h"
#include "utilities.h"
#include "ompatomic.h"

#include <cmath>
#include "filesystem.h"

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
			std::vector<std::string> buildFileList(const std::string& templ);

			/**
			Works as buildFileList but removed non-image files from the list.
			*/
			std::vector<std::string> buildFilteredFileList(const std::string& templ);

			/**
			Separates directory and filename parts of a sequence template.
			*/
			void separatePathAndFileTemplate(const std::string& templ, fs::path& dir, std::string& fileTemplate);

			/**
			Returns metadata of a 2D image.
			*/
			bool getInfo2D(const std::string& filename, coord_t& width, coord_t& height, ImageDataType& dataType, std::string& reason);

			/**
			Reads a 2D image.
			*/
			template<typename pixel_t> void read2D(Image<pixel_t>& image, const std::string& filename, coord_t z)
			{
				// TODO: Add other formats here.

				Vec3c dimensions;
				coord_t width, height;
				ImageDataType dataType;
				std::string reason;

				if (png::getInfo(filename, width, height, dataType, reason))
					png::read(image, filename, z);
				else if (tiff::getInfo(filename, dimensions, dataType, reason))
					tiff::read2D(image, filename, z);
				else
					throw ITLException(std::string("Unsupported sequence input file format (") + filename + ").");
			}

			/**
			Writes a 2D image.
			*/
			template<typename pixel_t> void write2D(const Image<pixel_t>& image, const std::string& filename, coord_t z)
			{
				// TODO: Add other formats here.

				if (endsWithIgnoreCase(filename, ".png"))
					png::write(image, filename, z);
				else if (endsWithIgnoreCase(filename, ".tif") || endsWithIgnoreCase(filename, ".tiff"))
					tiff::write2D(image, filename, z);
				else
					throw ITLException(std::string("Unsupported sequence output file format (") + filename + ").");
			}
			

		}

		/**
		Get general information of an image sequence that is assumed to form a 3D image.
		@param filename Directory name and possibly file name template including wildcards *, ?, @ that identifies files that belong to the sequence. (@ corresponds to one or more numerical digits)
		@param width, height, depth Dimensions of the image
		@param dataType Pixel data type of the image.
		@return True if the sequence seems to be a valid image with supported pixel data type.
		*/
		bool getInfo(const std::string& filename, coord_t& width, coord_t& height, coord_t& depth, ImageDataType& dataType, std::string& reason);

		inline bool getInfo(const std::string& filename, Vec3c& dims, ImageDataType& dataType, std::string& reason)
		{
			return getInfo(filename, dims.x, dims.y, dims.z, dataType, reason);
		}

		/**
		Tests if the given template corresponds to a valid image sequence.
		*/
		inline bool isSequence(const std::string& filename)
		{
			Vec3c dims;
			ImageDataType dt;
			std::string reason;
			return getInfo(filename, dims, dt, reason);
		}

		/**
		Read an image sequence to a 3D image.
		@param img Image where the data will be placed. The size of the image is set automatically. NOTE: Data type of the image must be correct - use getInfo to determine that.
		@param filename Directory name and possibly file name template including wildcards *, ?, @ that identifies files that belong to the sequence. (@ corresponds to one or more numerical digits)
		@param firstSlice Index of first slice to be read.
		@param lastSlice Index of last slice to be read. Defaults to all slices. The value is clamped to actual count of slices.
		*/
		template<typename pixel_t> void read(Image<pixel_t>& img, const std::string& filename, size_t firstSlice = 0, size_t lastSlice = std::numeric_limits<size_t>::max())
		{
			std::vector<std::string> files = internals::buildFilteredFileList(filename);

			if (files.size() <= 0)
				throw ITLException("Sequence contains no images.");

			clamp<size_t>(firstSlice, 0, files.size() - 1);
			clamp<size_t>(lastSlice, 0, files.size() - 1);

			if (lastSlice < firstSlice)
				throw ITLException("No slices to be read, last slice index is less than first slice index.");

			coord_t depth = (size_t)(lastSlice - firstSlice + 1);
			coord_t width, height;
			ImageDataType dataType;
			std::string reason;

			if (!internals::getInfo2D(files[0], width, height, dataType, reason))
				throw ITLException(reason);
				//throw ITLException("Unable to read dimensions of the image sequence.");

			// Check that image data type matches
			if (dataType != imageDataType<pixel_t>())
				throw ITLException(std::string("Expected data type is ") + toString(imageDataType<pixel_t>()) + " but the sequence contains data of type " + toString(dataType) + ".");

			img.ensureSize(width, height, depth);

			std::string errorMessage;
			OmpAtomic<bool> broken = false;

			size_t counter = 0;
			#pragma omp parallel for if(!omp_in_parallel())
			for (coord_t z = (coord_t)firstSlice; z <= (coord_t)lastSlice; z++)
			{
				if (!broken)
				{
					try
					{
						internals::read2D(img, files[z], z - firstSlice);
					}
					catch (const ITLException& ex)
					{
						broken = true;
						
						#pragma omp critical
						{
							errorMessage = ex.message();
						}
					}
					catch (const std::exception& ex)
					{
						broken = true;

						#pragma omp critical
						{
							errorMessage = ex.what();
						}
					}
					catch (...)
					{
						broken = true;

						#pragma omp critical
						{
							errorMessage = "Unknown error.";
						}
					}
				}

				showThreadProgress(counter, img.depth() + 1);
			}

			if (broken)
				throw ITLException(errorMessage);
		}

		/**
		Reads part of an image sequence to given image.
		Supports multiple computers accessing the same image at the same time only if the file system where the sequence is saved supports file locking correctly among all the computers accessing the data.
		NOTE: Does not support out-of-bounds start position.
		@param img Image where the data is placed. The size of the image defines the size of the block that is read.
		@param filename Directory name and possibly file name template including wildcards *, ?, @ that identifies files that belong to the sequence. (@ corresponds to one or more numerical digits)
		@param start Start location of the read. The size of the image defines the size of the block that is read.
		*/
		template<typename pixel_t> void readBlock(Image<pixel_t>& img, const std::string& filename, const Vec3c& start, bool showProgressInfo = false)
		{
			std::vector<std::string> files = internals::buildFilteredFileList(filename);

			if (files.size() <= 0)
				throw ITLException("Sequence contains no images.");

			Vec3c dimensions(0, 0, files.size());
			ImageDataType dataType;
			std::string reason;
			if (!internals::getInfo2D(files[0], dimensions.x, dimensions.y, dataType, reason))
				throw ITLException(std::string("Unable to read sequence information. ") + reason);

			// Check that image data type matches
			if (dataType != imageDataType<pixel_t>())
				throw ITLException("The data types of the image and the sequence do not match.");

			if (start.x < 0 || start.y < 0 || start.z < 0 || start.x >= dimensions.x || start.y >= dimensions.y || start.z >= dimensions.z)
				throw ITLException("Out of bounds start position in sequence::readBlock.");

			Vec3c cStart = start;
			clamp(cStart, Vec3c(0, 0, 0), dimensions);
			Vec3c cEnd = start + img.dimensions();
			clamp(cEnd, Vec3c(0, 0, 0), dimensions);

			std::string errorMessage;
			OmpAtomic<bool> broken = false;
			
			size_t counter = 0;
			#pragma omp parallel for if(!omp_in_parallel())
			for (coord_t z = cStart.z; z < cEnd.z; z++)
			{
				if (!broken)
				{
					// Read whole file
					Image<pixel_t> slice(dimensions.x, dimensions.y, 1);
					try
					{
						internals::read2D(slice, files[z], 0);

						// Pick the part of the file we need
						for (coord_t y = cStart.y; y < cEnd.y; y++)
						{
							for (coord_t x = cStart.x; x < cEnd.x; x++)
							{
								pixel_t pix = slice(x, y, 0);
								img(x - cStart.x, y - cStart.y, z - cStart.z) = pix;
							}
						}
					}
					catch (const ITLException& ex)
					{
						broken = true;

						#pragma omp critical
						{
							errorMessage = ex.message();
						}
					}
				}

				showThreadProgress(counter, img.depth() + 1, showProgressInfo);
			}

			if (broken)
				throw ITLException(errorMessage);
		}

		namespace internals
		{
			inline void prepareTemplate(const std::string& filename, size_t sliceCount, fs::path& dir, std::string& fileTempl, int& fieldWidth, size_t& atPos)
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
					std::string countStr = fileTempl.substr(atPos + 2, endPos - atPos - 2);
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

				// No extension given => select .png
				if (!endsWithIgnoreCase(fileTempl, ".tif") && !endsWithIgnoreCase(fileTempl, ".tiff") && !endsWithIgnoreCase(fileTempl, ".png"))
					fileTempl += ".png";
			}

			inline std::string insertSliceNumber(const std::string& fileTempl, int fieldWidth, size_t atPos, size_t z)
			{
				std::string zStr = toString(z);
				while (zStr.size() < fieldWidth)
					zStr = "0" + zStr;

				std::string filename = fileTempl;
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
		template<typename pixel_t> void write(const Image<pixel_t>& img, const std::string& filename, coord_t outputFirstSlice = 0, size_t firstSlice = 0, size_t lastSlice = std::numeric_limits<size_t>::max())
		{
			clamp<size_t>(firstSlice, 0, img.depth() - 1);
			clamp<size_t>(lastSlice, 0, img.depth() - 1);

			if (lastSlice < firstSlice)
				throw ITLException("No slices to be written, last slice index is less than first slice index.");

			fs::path dir;
			std::string fileTempl;
			int fieldWidth;
			size_t atPos;
			internals::prepareTemplate(filename, lastSlice - firstSlice, dir, fileTempl, fieldWidth, atPos);

			// Write all slices
			fs::create_directories(dir);
			size_t counter = 0;
			std::string errorMessage;
			OmpAtomic<bool> broken = false;

			#pragma omp parallel for if(!omp_in_parallel())
			for (coord_t z = firstSlice; z <= (coord_t)lastSlice; z++)
			{
				if (!broken)
				{
					try
					{
						std::string filename = internals::insertSliceNumber(fileTempl, fieldWidth, atPos, z - firstSlice + outputFirstSlice);

						internals::write2D(img, (dir / filename).string(), z);

					}
					catch (ITLException& ex)
					{
						broken = true;
						#pragma omp critical
						{
							errorMessage = ex.message();
						}
					}
				}

				showThreadProgress(counter, lastSlice - firstSlice + 1);
			}

			if (broken)
				throw ITLException(errorMessage);
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
		template<typename pixel_t> void writeBlock(const Image<pixel_t>& img, const std::string& filename, const Vec3c& filePosition, const Vec3c& fileDimensions,
			const Vec3c& imagePosition, const Vec3c& imageDimensions,
			bool showProgressInfo = false)
		{
			constexpr size_t READ_WRITE_TRIALS = 12;
			constexpr int INITIAL_WAIT_TIME = 1;

			fs::path dir;
			std::string fileTempl;
			int fieldWidth;
			size_t atPos;
			internals::prepareTemplate(filename, fileDimensions.z, dir, fileTempl, fieldWidth, atPos);

			Vec3c cStart = filePosition;
			clamp(cStart, Vec3c(0, 0, 0), fileDimensions);
			Vec3c cEnd = filePosition + imageDimensions;
			clamp(cEnd, Vec3c(0, 0, 0), fileDimensions);

			if (!img.isInImage(imagePosition))
				throw ITLException("Block start position must be inside the image.");
			if (!img.isInImage(imagePosition + imageDimensions - Vec3c(1, 1, 1)))
				throw ITLException("Block end position must be inside the image.");

			fs::create_directories(dir);

			std::string errorMessage;
			OmpAtomic<bool> broken = false;

			size_t counter = 0;
			#pragma omp parallel for if(!omp_in_parallel())
			for (coord_t z = cStart.z; z < cEnd.z; z++)
			{
				if (!broken)
				{
					try
					{
						std::string filename = internals::insertSliceNumber(fileTempl, fieldWidth, atPos, z/* - cStart.z + imagePosition.z*/);
						filename = (dir / filename).string();


						if (cStart.x == 0 && cStart.y == 0 && cEnd.x == img.width() && cEnd.y == img.height() &&
							imagePosition.x == 0 && imagePosition.y == 0)
						{
							// We are writing full slices. Reading phase can be skipped (see below).
							internals::write2D(img, filename, z - cStart.z + imagePosition.z);
						}
						else
						{
							// We are writing only a part of each slice.

							//FileLock lock(filename); File locking does not really work on network shares so... :(

							Image<pixel_t> slice(fileDimensions.x, fileDimensions.y, 1);

							// Read existing data
							if (fs::exists(filename))
							{
								int waitTime = INITIAL_WAIT_TIME;
								for (size_t n = 0; n < READ_WRITE_TRIALS; n++)
								{
									try
									{
										coord_t fw, fh;
										ImageDataType fd;
										std::string reason;
										if (internals::getInfo2D(filename, fw, fh, fd, reason))
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
										sleep(waitTime);
										waitTime *= 2;
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
							int waitTime = INITIAL_WAIT_TIME;
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
									sleep(waitTime);
									waitTime *= 2;
								}
							}
						}

					}
					catch (const ITLException& ex)
					{
						broken = true;
						
						#pragma omp critical
						{
							errorMessage = ex.message();
						}
					}
				}

				showThreadProgress(counter, cEnd.z - cStart.z, showProgressInfo);
			}

			if (broken)
				throw ITLException(errorMessage);
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
		template<typename pixel_t> void writeBlock(const Image<pixel_t>& img, const std::string& filename, const Vec3c& filePosition, const Vec3c& fileDimensions, bool showProgressInfo = false)
		{
			writeBlock(img, filename, filePosition, fileDimensions, Vec3c(0, 0, 0), img.dimensions(), showProgressInfo);
		}

		/**
		Move or copy image sequence from one location to another.
		*/
		inline void moveOrCopySequence(const std::string& inputName, const std::string& outputName, bool move)
		{
			// TODO: Check that input and output sequences have the same file type.

			// Get list of input files
			std::vector<std::string> inputFiles = internals::buildFileList(inputName);

			if (inputFiles.size() <= 0)
				throw ITLException("Sequence contains no images.");

			// Prepare output file name template
			fs::path dir;
			std::string fileTempl;
			int fieldWidth;
			size_t atPos;
			internals::prepareTemplate(outputName, inputFiles.size(), dir, fileTempl, fieldWidth, atPos);

			fs::create_directories(dir);
			#pragma omp parallel for if(!omp_in_parallel())
			for (coord_t z = 0; z < (coord_t)inputFiles.size(); z++)
			{
				std::string infile = inputFiles[z];

				std::string outfile = internals::insertSliceNumber(fileTempl, fieldWidth, atPos, z);
				outfile = (dir / outfile).string();

				if (move)
					fs::rename(infile, outfile);
				else
					fs::copy(infile, outfile, fs::copy_options::overwrite_existing);
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
			void readWriteBlockOptimization();
			void fileFormats();
		}
	}

}
