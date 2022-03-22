#pragma once

#include <string>

#include "datatypes.h"
#include "itlexception.h"
#include "image.h"
#include "utilities.h"
#include "timer.h"
#include "math/mathutils.h"
#include "io/imagedatatype.h"
#include "math/vec3.h"
#include "progress.h"

namespace itl2
{
	/**
	Methods for manipulating raw files.
	*/
	namespace raw
	{
		namespace internals
		{
			/**
			Tries to find pixel data type based on file size and dimensions.
			*/
			inline ImageDataType estimateDataType(size_t fileSize, const Vec3c& dimensions, size_t& pixelSizeBytes)
			{
				size_t pixelCount = (size_t)dimensions.x * (size_t)dimensions.y * (size_t)dimensions.z;

				pixelSizeBytes = fileSize / pixelCount;

				if (fileSize == pixelCount * sizeof(uint8_t))
					return ImageDataType::UInt8;

				if (fileSize == pixelCount * sizeof(uint16_t))
					return ImageDataType::UInt16;

				// uint32 and float32 cannot be differentiated based on file size only.
				// Prefer float32 as that is more common.

				if (fileSize == pixelCount * sizeof(float32_t))
					return ImageDataType::Float32;

				if (fileSize == pixelCount * sizeof(uint64_t))
					return ImageDataType::UInt64;

				// complex32 and uint64 cannot be differentiated based on file size only.
				// Prefer uint64 as that is more common.

				//if (fileSize == pixelCount * sizeof(complex32_t))
				//	return ImageDataType::Complex32;

				return ImageDataType::Unknown;
			}


			/**
			Tries to find pixel data type based on file name, size and dimensions.
			*/
			inline ImageDataType estimateDataType(const std::string& filename, const Vec3c& dimensions, size_t& pixelSizeBytes)
			{
				// First try to find data type from file name.
				ImageDataType dt = ImageDataType::Unknown;
				if (containsIgnoreCase(filename, toString(ImageDataType::UInt8)))
					dt = ImageDataType::UInt8;
				else if (containsIgnoreCase(filename, toString(ImageDataType::UInt16)))
					dt = ImageDataType::UInt16;
				else if (containsIgnoreCase(filename, toString(ImageDataType::UInt32)))
					dt = ImageDataType::UInt32;
				else if (containsIgnoreCase(filename, toString(ImageDataType::UInt64)))
					dt = ImageDataType::UInt64;
				else if (containsIgnoreCase(filename, toString(ImageDataType::Int8)))
					dt = ImageDataType::Int8;
				else if (containsIgnoreCase(filename, toString(ImageDataType::Int16)))
					dt = ImageDataType::Int16;
				else if (containsIgnoreCase(filename, toString(ImageDataType::Int32)))
					dt = ImageDataType::Int32;
				else if (containsIgnoreCase(filename, toString(ImageDataType::Int64)))
					dt = ImageDataType::Int64;
				else if (containsIgnoreCase(filename, toString(ImageDataType::Float32)))
					dt = ImageDataType::Float32;
				else if (containsIgnoreCase(filename, toString(ImageDataType::Complex32)))
					dt = ImageDataType::Complex32;

				if (dt != ImageDataType::Unknown)
				{
					pixelSizeBytes = pixelSize(dt);
					return dt;
				}

				// No data type information found from the file name, so try to estimate from file size.
				size_t fileSize = (size_t)itl2::fileSize(filename);
				return estimateDataType(fileSize, dimensions, pixelSizeBytes);
			}

			/**
			Adds .raw image dimensions to file name.
			*/
			template<typename pixel_t> std::string concatDimensions(const std::string& baseName, const Image<pixel_t>& img)
			{
				return itl2::concatDimensions(baseName, img.dimensions());
			}

			/**
			If a file with given name exists, does nothing. If it does not exist, tests if the name matches any files with template "{filename}*.raw", and if
			it matches a single file, sets filename to name of that file. If it matches none, throws an ITLException.
			If it matches multiple files, tries again with template "{filename}_@x@x@.raw", where @ is one or more numerical digits.
			If unique match is still not found, throws and ITLException.
			*/
			inline void expandRawFilename(::std::string& filename)
			{
				if (!fs::exists(filename))
				{
					::std::vector<::std::string> candidates = buildFileList(filename + "*.raw");

					if (candidates.size() == 0)
						throw ITLException(::std::string("No file found matching to template ") + filename + "*.raw");

					if (candidates.size() == 1)
						filename = candidates[0];
					else
					{
						// Try with more restricting templates.
						
						// First, try with the template as-is.
						candidates = buildFileList(filename + "_@x@x@.raw");

						if (candidates.size() == 1)
							filename = candidates[0];
						else
						{
							// None or multiple candidates.

							// Try removing "_" from the end of the template, as the user might have added that already.
							::std::string reducedFilename = filename;
							if (endsWith(reducedFilename, "_"))
							{
								reducedFilename = reducedFilename.substr(0, reducedFilename.length() - 1);
								candidates = buildFileList(reducedFilename + "_@x@x@.raw");
							}

							if (candidates.size() == 0)
								throw ITLException(string("Multiple .raw files match template ") + filename + "*.raw but no files match template " + reducedFilename + "_@x@x@.raw with one or more underscores.");

							if (candidates.size() == 1)
								filename = candidates[0];
							else
								throw ITLException(string("Multiple .raw files match the templates ") + filename + "*.raw and " + reducedFilename + "_@x@x@.raw with one or more underscores.");
						}
					}
				}
			}

			/**
			Parses dimensions from .raw file name.
			Supports names like
			image_name_100x200x300.raw
			and
			image_name_100x200.raw
			*/
			inline bool parseDimensions(const ::std::string& filename, Vec3c& dimensions)
			{
				::std::string::size_type temppos = filename.find_last_of('.');
				if (temppos != ::std::string::npos)
					temppos = temppos - 1;
				::std::string::size_type startpos = filename.find_last_not_of("0123456789x", temppos);
				if (startpos == ::std::string::npos)
					return false;

				::std::string::size_type endpos = filename.find('x', startpos);
				if (endpos == ::std::string::npos)
					return false;

				::std::string ws = filename.substr(startpos + 1, endpos - startpos - 1);
				size_t w = itl2::fromString<size_t>(ws);

				::std::string::size_type startpos2 = endpos;
				::std::string::size_type endpos2 = filename.find('x', startpos2 + 1);
				if (endpos2 == ::std::string::npos)
				{
					// 2D image
					endpos2 = filename.find('.', startpos);
					if (endpos2 == ::std::string::npos)
						return false;


					::std::string hs = filename.substr(startpos2 + 1, endpos2 - startpos2 - 1);
					size_t h = itl2::fromString<size_t>(hs);
					dimensions = Vec3c((coord_t)w, (coord_t)h, 1);

					return true;
				}

				::std::string hs = filename.substr(startpos2 + 1, endpos2 - startpos2 - 1);
				size_t h = itl2::fromString<size_t>(hs);

				::std::string::size_type startpos3 = endpos2;
				::std::string::size_type endpos3 = filename.find('.', startpos3 + 1);
				if (endpos3 == ::std::string::npos)
					return false;

				::std::string ds = filename.substr(startpos3 + 1, endpos3 - startpos3 - 1);
				size_t d = itl2::fromString<size_t>(ds);

				dimensions = Vec3c((coord_t)w, (coord_t)h, (coord_t)d);
				return true;
			}
		}

		/**
		Gets dimension vector and data type from a file name like
		image_name_100x200x300.raw or
		image_name_100x200.raw
		@return True if the file name could be parsed and the file exists.
		*/
		inline bool getInfo(std::string filename, Vec3c& dimensions, ImageDataType& dataType, size_t& pixelSizeBytes, std::string& reason)
		{
			try
			{
				internals::expandRawFilename(filename);
			}
			catch (ITLException e)
			{
				reason = e.message();
				return false;
			}

			if (!internals::parseDimensions(filename, dimensions))
			{
				reason = "Unable to find dimensions from file name in [name]_[width]x[height]x[depth].raw format.";
				return false;
			}

			if (dimensions.min() <= 0)
			{
				reason = string("Invalid image dimensions: ") + toString(dimensions);
				return false;
			}

			if (!fs::exists(filename))
			{
				reason = "Input file does not exist.";
				return false;
			}

			dataType = internals::estimateDataType(filename, dimensions, pixelSizeBytes);
			//return dataType != ImageDataType::Unknown;
			return true;
		}

		inline bool getInfo(std::string filename, Vec3c& dimensions, ImageDataType& dataType, std::string& reason)
		{
			size_t pixelSizeBytes;
			return getInfo(filename, dimensions, dataType, pixelSizeBytes, reason);
		}

        /**
		Writes any trivially copyable value to stream.
		*/
		template<typename T>
		typename std::enable_if<std::is_trivially_copyable_v<T> >::type writePixel(std::ofstream& out, const T& item)
		{
			out.write((char*)&item, sizeof(T));
		}

		/**
		Reads any trivially copyable value from stream.
		*/
		template<typename T>
		typename std::enable_if<std::is_trivially_copyable_v<T> >::type readPixel(std::ifstream& in, T& item)
		{
			in.read((char*)&item, sizeof(T));
		}


		/**
		Reads a raw file to the given image.
		If the pixel data type is trivially copyable, the pixels are read directly from the file to the image memory buffer.
		In this case the pixels are assumed to be stored in native byte order.
		If the pixel data type is NOT trivially copyable, suitable function readPixel(ifstream& in, pixel_t& target)
		must be passed as parameter, and the function must read one pixel from the stream and assign its value to target.
		@param img Image where the data is placed. The size of the image must be correct.
		@param filename The name of the file to read.
		@param bytesToSkip Skip this many bytes from the beginning of the file.
		@param readPixel Pixel reading function. Relevant only for non-trivially copyable pixel data types.
		*/
		template<typename pixel_t, typename ReadPixel = decltype(raw::readPixel<pixel_t>)> void readNoParse(Image<pixel_t>& img, const std::string& filename, size_t bytesToSkip = 0, ReadPixel readPixel = raw::readPixel<pixel_t>)
		{
			std::ifstream in(filename.c_str(), std::ios_base::in | std::ios_base::binary);

			if(!in)
            {
                throw ITLException(std::string("Unable to open ") + filename + std::string(", ") + getStreamErrorMessage());
            }

			in.seekg(bytesToSkip, std::ios::beg);

			if constexpr (std::is_trivially_copyable_v<pixel_t>)
			{
				// Load directly to the buffer
				size_t READ_SIZE = 200 * 1024 * 1024;
				size_t read_start = 0;
				char* pBuffer = (char*)img.getData();
				while (read_start < img.pixelCount() * sizeof(pixel_t))
				{
					size_t read_size = READ_SIZE;
					if (read_start + read_size > img.pixelCount() * sizeof(pixel_t))
						read_size = img.pixelCount() * sizeof(pixel_t) - read_start;

					in.read(pBuffer + read_start, read_size);

					if (in.bad())
					{
						throw ITLException(std::string("Read failed, file size might be incorrect: ") + filename);
					}

					read_start += read_size;
				}
			}
			else
			{
				for (coord_t n = 0; n < img.pixelCount(); n++)
				{
					readPixel(in, img(n));
				}
			}
		}
		
		/**
		Reads a part of a .raw file to the given image.
		NOTE: Does not support out of bounds start position.
		@param img Image where the data is placed. The size of the image defines the size of the block that is read.
		@param filename The name of the file to read.
		@param fileDimensions Dimensions of the whole file.
		@param fileStartPos Start location of the read. The size of the image defines the size of the block that is read.
		*/
		template<typename pixel_t> void readBlockNoParse(Image<pixel_t>& img, const std::string& filename, const Vec3c& fileDimensions, const Vec3c& fileStart, bool showProgressInfo = false, size_t bytesToSkip = 0)
		{
			if (fileStart.x < 0 || fileStart.y < 0 || fileStart.z < 0 || fileStart.x >= fileDimensions.x || fileStart.y >= fileDimensions.y || fileStart.z >= fileDimensions.z)
				throw ITLException("Out of bounds start position in raw::readBlock.");

			Vec3c cStart = fileStart;
			clamp(cStart, Vec3c(0, 0, 0), fileDimensions);
			Vec3c cEnd = fileStart + img.dimensions();
			clamp(cEnd, Vec3c(0, 0, 0), fileDimensions);

			if (cStart == Vec3c(0, 0, 0) && cEnd == fileDimensions && fileDimensions == img.dimensions())
			{
				// Reading whole file, use the whole file reading function.
				raw::readNoParse(img, filename, bytesToSkip);
				return;
			}


			std::ifstream in(filename.c_str(), std::ios_base::in | std::ios_base::binary);

			if (!in)
			{
				throw ITLException(std::string("Unable to open ") + filename + std::string(", ") + getStreamErrorMessage());
			}

			in.seekg(bytesToSkip, std::ios::beg);

			pixel_t* pBuffer = img.getData();

			if (cStart.x == 0 && cEnd.x == fileDimensions.x)
			{
				// Reading whole scan lines.
				// We can read one slice per one read call.
				for (coord_t z = cStart.z; z < cEnd.z; z++)
				{
					size_t filePos = (z * fileDimensions.x * fileDimensions.y + cStart.y * fileDimensions.x + cStart.x) * sizeof(pixel_t);
					in.seekg(filePos);

					size_t imgPos = img.getLinearIndex(0, 0, z - cStart.z);
					in.read((char*)&pBuffer[imgPos], (cEnd.x - cStart.x) * (cEnd.y - cStart.y) * sizeof(pixel_t));

					if (in.bad())
					{
						showProgress(cEnd.z - cStart.z, cEnd.z - cStart.z, showProgressInfo);
						throw ITLException(std::string("Failed to read (slice at a time) block of ") + filename + std::string(", ") + getStreamErrorMessage());
					}

					showProgress(z - cStart.z, cEnd.z - cStart.z, showProgressInfo);
				}
			}
			else
			{
				// Reading partial scan lines.
				for (coord_t z = cStart.z; z < cEnd.z; z++)
				{
					for (coord_t y = cStart.y; y < cEnd.y; y++)
					{
						size_t filePos = (z * fileDimensions.x * fileDimensions.y + y * fileDimensions.x + cStart.x) * sizeof(pixel_t);
						in.seekg(filePos);

						size_t imgPos = img.getLinearIndex(0, y - cStart.y, z - cStart.z);
						in.read((char*)&pBuffer[imgPos], (cEnd.x - cStart.x) * sizeof(pixel_t));

						if (in.bad())
						{
							showProgress(cEnd.z - cStart.z, cEnd.z - cStart.z, showProgressInfo);
							throw ITLException(std::string("Failed to read block of ") + filename + std::string(", ") + getStreamErrorMessage());
						}
					}

					
					showProgress(z - cStart.z, cEnd.z - cStart.z, showProgressInfo);
				}
			}

		}

		template<typename pixel_t> void getInfoAndCheck(const std::string& filename, Vec3c& dimensions)
		{
			ImageDataType dataType;
			size_t pixelSizeBytes;

			std::string reason;
			if (!getInfo(filename, dimensions, dataType, pixelSizeBytes, reason))
				throw ITLException("Unable to read " + filename + ". " + reason);

			// Only check pixel size if the pixels are trivially copyable.
			if constexpr (std::is_trivially_copyable_v<pixel_t>)
			{
				// Check that pixel size is correct.
				if (pixelSizeBytes != sizeof(pixel_t))
					throw ITLException(std::string("Expected pixel size is ") + toString(sizeof(pixel_t)) + " bytes but the .raw file contains pixels of size " + toString(pixelSizeBytes) + " bytes.");
			}

			// NOTE: This does not work for e.g. int32_t images as getInfo assumes they are float32_t images.
			//if (dataType != ImageDataType::Unknown)
			//{
			//	if (dataType != imageDataType<pixel_t>())
			//		throw ITLException(string("Expected data type is ") + toString(imageDataType<pixel_t>()) + " but the .raw file contains data of type " + toString(dataType) + ".");
			//}
			//else
			//{
			//	// Unknown data types can be read if the pixel size is correct.
			//	if(pixelSizeBytes != sizeof(pixel_t))
			//		throw ITLException(string("Expected pixel size is ") + toString(sizeof(pixel_t)) + " but the .raw file contains pixels of size " + toString(pixelSizeBytes) + ".");
			//}
		}
		
		/**
		Reads a .raw file to the given image, initializes the image to correct size read from file name.
		The file name must be in format image_name_100x200x300.raw or image_name_100x200.raw.
		@param img Image where the data is placed. The size of the image will be set based on the .raw file name.
		@param filename The name of the file to read.
		@param bytesToSkip Skip this many bytes from the beginning of the file.
		@param readPixel Function that reads one pixel. Relevant only for non-trivially copyable pixel data types.
		*/
		template<typename pixel_t, typename ReadPixel = decltype(raw::readPixel<pixel_t>)> void read(Image<pixel_t>& img, std::string filename, size_t bytesToSkip = 0, ReadPixel readPixel = raw::readPixel<pixel_t>)
		{
			Vec3c dimensions;
			getInfoAndCheck<pixel_t>(filename, dimensions);

			img.ensureSize(dimensions);
			internals::expandRawFilename(filename);
			readNoParse<pixel_t, ReadPixel>(img, filename, bytesToSkip, readPixel);
		}

		/**
		Reads a part of a .raw file to the given image.
		Reads image size from .raw file name.
		NOTE: Does not support out of bounds start position.
		@param img Image where the data is placed. The size of the image defines the size of the block that is read.
		@param filename The name of the file to read.
		@param fileStartPos Start location of the read. The size of the image defines the size of the block that is read.
		*/
		template<typename pixel_t> void readBlock(Image<pixel_t>& img, std::string filename, const Vec3c& fileStart, bool showProgressInfo = false)
		{
			Vec3c dimensions;
			getInfoAndCheck<pixel_t>(filename, dimensions);

			internals::expandRawFilename(filename);
			readBlockNoParse(img, filename, dimensions, fileStart, showProgressInfo);
		}


		/**
		Writes a block of an image to the specified location in a raw file.
		The output file is not truncated if it exists.
		If the output file does not exist, it is created.
		Part of image extending beyond [0, fileDimensions[ is not written.
		@param img Image to write.
		@param filename Name of file to write.
		@param filePosition Position in the file to write to.
		@param fileDimension Total dimensions of the entire output file.
		@param imagePosition Position in the image where the block to be written starts.
		@param blockDimensions Dimensions of the block of the source image to write.
		@param showProgressInfo Set to true to show a progress bar.
		*/
		template<typename pixel_t> void writeBlock(const Image<pixel_t>& img, const std::string& filename,
			const Vec3c& filePosition, const Vec3c& fileDimensions,
			const Vec3c& imagePosition,
			const Vec3c& blockDimensions,
			bool showProgressInfo = false)
		{
			Vec3c fileStartPos = filePosition;
			clamp(fileStartPos, Vec3c(0, 0, 0), fileDimensions);
			Vec3c fileEndPos = filePosition + blockDimensions;
			clamp(fileEndPos, Vec3c(0, 0, 0), fileDimensions);

			if (!img.isInImage(imagePosition))
				throw ITLException("Block start position must be inside the image.");
			if (!img.isInImage(imagePosition + blockDimensions - Vec3c(1, 1, 1)))
				throw ITLException("Block end position must be inside the image.");

			createFoldersFor(filename);

			// Create file if it does not exist, otherwise set file size to the correct value.
			setFileSize(filename, fileDimensions.x * fileDimensions.y * fileDimensions.z * sizeof(pixel_t));

			std::ofstream out(filename.c_str(), std::ios_base::in | std::ios_base::out | std::ios_base::binary);

			if(!out)
				throw ITLException(std::string("Unable to open ") + filename + std::string(", ") + getStreamErrorMessage());

			const pixel_t* pBuffer = img.getData();
			{
				ProgressIndicator prog(fileEndPos.z - fileStartPos.z, showProgressInfo);

				for (coord_t z = fileStartPos.z; z < fileEndPos.z; z++)
				{
					if (fileStartPos.x == 0 && fileEndPos.x == fileDimensions.x &&
						fileDimensions.x == blockDimensions.x &&
						img.width() == fileDimensions.x)
					{
						// Writing whole scanlines.
						// Write all scanlines in region [fileStartPos.y, fileEndPos.y[ at once in order to increase write speed.
						size_t filePos = (z * fileDimensions.x * fileDimensions.y + fileStartPos.y * fileDimensions.x + fileStartPos.x) * sizeof(pixel_t);
						out.seekp(filePos);

						if (!out)
							throw ITLException(std::string("Seek failed (fast write) for file ") + filename + std::string(", ") + getStreamErrorMessage());

						size_t imgPos = img.getLinearIndex(imagePosition.x, imagePosition.y, z - fileStartPos.z + imagePosition.z);
						out.write((char*)&pBuffer[imgPos], (fileEndPos.x - fileStartPos.x) * (fileEndPos.y - fileStartPos.y) * sizeof(pixel_t));

						if (!out)
							throw ITLException(std::string("Unable to write (fast) to ") + filename + std::string(", ") + getStreamErrorMessage());
					}
					else
					{
						// Writing partial scanlines.

						for (coord_t y = fileStartPos.y; y < fileEndPos.y; y++)
						{
							size_t filePos = (z * fileDimensions.x * fileDimensions.y + y * fileDimensions.x + fileStartPos.x) * sizeof(pixel_t);
							out.seekp(filePos);

							if (!out)
								throw ITLException(std::string("Seek failed for file ") + filename + std::string(", ") + getStreamErrorMessage());

							size_t imgPos = img.getLinearIndex(imagePosition.x, y - fileStartPos.y + imagePosition.y, z - fileStartPos.z + imagePosition.z);
							out.write((char*)&pBuffer[imgPos], (fileEndPos.x - fileStartPos.x) * sizeof(pixel_t));

							if (!out)
								throw ITLException(std::string("Unable to write to ") + filename + std::string(", ") + getStreamErrorMessage());
						}
					}

					prog.step();
				}
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
		template<typename pixel_t> void writeBlock(const Image<pixel_t>& img, const std::string& filename, const Vec3c& filePosition, const Vec3c& fileDimensions, bool showProgressInfo = false)
		{
			writeBlock(img, filename, filePosition, fileDimensions, Vec3c(0, 0, 0), img.dimensions(), showProgressInfo);
		}

		/**
		Writes an image to a .raw file.
		If the pixel data type is trivially copyable, the pixels are written directly to the file from the image memory buffer.
		In this case the pixels are written in native byte order.
		If the pixel data type is not trivially copyable, suitable overload of itl2::raw::write(ofstream& out, pixel_type& p)
		is called, and that overload should write the pixel to the file.
		Therefore, by providing the overload of itl2::raw::write to each non-trivially copyable pixel type any type of pixel data
		can be written to the .raw file.
		@param img Image to write.
		@param filename Name of file to write.
		@param truncate Set to false to add to existing file.
		*/
		template<typename pixel_t, typename WritePixel = decltype(raw::writePixel<pixel_t>)> void write(const Image<pixel_t>& img, const std::string& filename, bool truncate = true, WritePixel writePixel = raw::writePixel<pixel_t>)
		{
			createFoldersFor(filename);

			std::ios::openmode mode;
			if(truncate)
				mode = std::ios_base::out | std::ios_base::trunc | std::ios_base::binary;
			else
				mode = std::ios_base::out | std::ios_base::app | std::ios_base::binary;

			std::ofstream out(filename.c_str(), mode);

			// If the file does not exist, create it.
			if(!out.is_open())
				out.open(filename.c_str(), std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);

			if(!out)
                throw ITLException(std::string("Unable to write to ") + filename + std::string(", ") + getStreamErrorMessage());


			if constexpr (std::is_trivially_copyable_v<pixel_t>)
			{

				// Write directly from the native buffer
				// Writing everything at once does not work for big files.

				size_t WRITE_SIZE = 200 * 1024 * 1024;
				size_t write_start = 0;
				const char* pBuffer = (const char*)img.getData();
				while (write_start < img.pixelCount() * sizeof(pixel_t))
				{
					size_t write_size = WRITE_SIZE;
					if (write_start + write_size > img.pixelCount() * sizeof(pixel_t))
						write_size = img.pixelCount() * sizeof(pixel_t) - write_start;

					out.write(pBuffer + write_start, write_size);
					write_start += write_size;
				}
			}
			else
			{
				for (coord_t n = 0; n < img.pixelCount(); n++)
				{
					writePixel(out, img(n));
				}
			}
		}

		/**
		Writes an image to a .raw file.
		Concatenates image dimensions to file name automatically.
		@param img Image to write.
		@param filename Template of file name. The full file name will be [template]_[width]x[height]x[depth].raw.
		@param truncate Set to false to add to existing file.
		@return The name of the file that was written.
		*/
		template<typename pixel_t, typename WritePixel = decltype(raw::writePixel<pixel_t>)> std::string writed(const Image<pixel_t>& img, const std::string& filename, bool truncate = true, WritePixel writePixel = raw::writePixel<pixel_t>)
		{
			std::string cf = internals::concatDimensions(filename, img);
			write(img, cf, truncate, writePixel);
			return cf;
		}



		/**
		Writes an RGB image to a .raw file.
		@param r, g, b Red, green and blue color component images.
		@param filename Name of file to write.
		@param truncate Set to false to add to existing file.
		*/
		void write(const Image<uint8_t>& r, const Image<uint8_t>& g, const Image<uint8_t>& b, const std::string& filename, bool truncate = true);

		/**
		Writes an RGB image to .raw file.
		Concatenates image dimensions to file name automatically.
		@param r, g, b Red, green and blue color component images.
		@param filename Template of file name. The full file name will be [template]_[width]x[height]x[depth].raw.
		@param truncate Set to false to add to existing file.
		@return The name of the file that was written.
		*/
		std::string writed(const Image<uint8_t>& r, const Image<uint8_t>& g, const Image<uint8_t>& b, const std::string& filename, bool truncate = true);


		namespace tests
		{
			void parseDimensions();
			void raw();
			void writeBlock();
			void writeBlockFast();
			void expandFilename();
		}
	}

}
