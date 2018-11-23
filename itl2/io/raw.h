#pragma once

#include <string>

#include "datatypes.h"
#include "itlexception.h"
#include "image.h"
#include "utilities.h"
#include "timer.h"
#include "math/mathutils.h"
#include "io/imagedatatype.h"
#include "io/sequence.h"
#include "io/fileutils.h"


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
			inline ImageDataType estimateDataType(const std::string& filename, const math::Vec3c& dimensions)
			{
				size_t fileSize = (size_t)itl2::fileSize(filename);
				size_t pixelCount = (size_t)dimensions.x * (size_t)dimensions.y * (size_t)dimensions.z;
				if (fileSize == pixelCount * sizeof(uint8_t))
					return UInt8;
				if (fileSize == pixelCount * sizeof(uint16_t))
					return UInt16;
				if (fileSize == pixelCount * sizeof(float32_t))
					return Float32;
				if (fileSize == pixelCount * sizeof(complex32_t))
					return Complex32;
				return Unknown;
			}

			/**
			Gets dimension vector and data type from a file name like
			image_name_100x200x300.raw or
			image_name_100x200.raw
			*/
			inline bool parseDimensions(const std::string& filename, math::Vec3c& dimensions, ImageDataType& dataType)
			{
				std::string::size_type startpos = filename.find_last_of('_');
				if (startpos == std::string::npos)
					return false;

				std::string::size_type endpos = filename.find('x', startpos);
				if (endpos == std::string::npos)
					return false;

				std::string ws = filename.substr(startpos + 1, endpos - startpos - 1);
				size_t w = itl2::fromString<size_t>(ws);

				std::string::size_type startpos2 = endpos;
				std::string::size_type endpos2 = filename.find('x', startpos2 + 1);
				if (endpos2 == std::string::npos)
				{
					// 2D image
					endpos2 = filename.find('.', startpos);
					if (endpos2 == std::string::npos)
						return false;


					std::string hs = filename.substr(startpos2 + 1, endpos2 - startpos2 - 1);
					size_t h = itl2::fromString<size_t>(hs);
					dimensions = math::Vec3c((coord_t)w, (coord_t)h, 1);

					dataType = estimateDataType(filename, dimensions);
					return true;
				}

				std::string hs = filename.substr(startpos2 + 1, endpos2 - startpos2 - 1);
				size_t h = itl2::fromString<size_t>(hs);

				std::string::size_type startpos3 = endpos2;
				std::string::size_type endpos3 = filename.find('.', startpos3 + 1);
				if (endpos3 == std::string::npos)
					return false;

				std::string ds = filename.substr(startpos3 + 1, endpos3 - startpos3 - 1);
				size_t d = itl2::fromString<size_t>(ds);

				dimensions = math::Vec3c((coord_t)w, (coord_t)h, (coord_t)d);
				dataType = estimateDataType(filename, dimensions);
				return true;
			}

		}

		/**
		Reads part of a .raw file to given image.
		NOTE: Does not support out of bounds start position.
		@param img Image where the data is placed. The size of the image defines the size of the block that is read.
		@param filename The name of the file to read.
		@param dimensions Dimensions of the whole file.
		@param start Start location of the read. The size of the image defines the size of the block that is read.
		*/
		template<typename pixel_t> void readBlock(Image<pixel_t>& img, const std::string& filename, const math::Vec3c& dimensions, const math::Vec3c& start, bool showProgressInfo = false)
		{
			ifstream in(filename.c_str(), ios_base::in | ios_base::binary);

			//if (!in.good())
			if(!in)
			{
				throw ITLException(std::string("File not found: ") + filename);
			}

			if (start.x < 0 || start.y < 0 || start.z < 0 || start.x >= dimensions.x || start.y >= dimensions.y || start.z >= dimensions.z)
				throw ITLException("Out of bounds start position in raw::readBlock.");

			math::Vec3c cStart = start;
			clamp(cStart, math::Vec3c(0, 0, 0), dimensions);
			math::Vec3c cEnd = start + img.dimensions();
			clamp(cEnd, math::Vec3c(0, 0, 0), dimensions);

			pixel_t* pBuffer = img.getData();
			for (coord_t z = cStart.z; z < cEnd.z; z++)
			{
				for (coord_t y = cStart.y; y < cEnd.y; y++)
				{
					size_t filePos = (z * dimensions.x * dimensions.y + y * dimensions.x + cStart.x) * sizeof(pixel_t);
					in.seekg(filePos);

					size_t imgPos = img.getLinearIndex(0, y - cStart.y, z - cStart.z);
					in.read((char*)&pBuffer[imgPos], (cEnd.x - cStart.x) * sizeof(pixel_t));

					if (in.bad())
					{
						throw ITLException(std::string("Read raw block failed, size of file is incorrect: ") + filename);
					}
				}

				if (showProgressInfo)
					showProgress(z - cStart.z, cEnd.z - cStart.z);
			}
		}
		
		

		/**
		Reads raw file to given image.
		@param img Image where the data is placed. The size of the image must be correct.
		@param filename The name of the file to read.
		@param bytesToSkip Skip this many bytes from the beginning of the file.
		*/
		template<typename pixel_t> void read(Image<pixel_t>& img, const std::string& filename, size_t bytesToSkip = 0)
		{
			ifstream in(filename.c_str(), ios_base::in | ios_base::binary);

            //if(!in.good())
			if(!in)
            {
                throw ITLException(std::string("File not found: ") + filename);
            }

			in.seekg(bytesToSkip, ios::beg);

			// Load directly to the buffer
			size_t READ_SIZE = 200 * 1024 * 1024;
			size_t read_start = 0;
			char* pBuffer = (char*)img.getData();
			while(read_start < img.pixelCount() * sizeof(pixel_t))
			{
				size_t read_size = READ_SIZE;
				if(read_start + read_size > img.pixelCount() * sizeof(pixel_t))
					read_size = img.pixelCount() * sizeof(pixel_t) - read_start;

				in.read(pBuffer + read_start, read_size);

				if(in.bad())
				{
					throw ITLException(std::string("Read raw failed, size of file is incorrect: ") + filename);
				}

				read_start += read_size;
			}
			/*
			// Swap byte order
			if(fileIsBigEndian)
			{
				pixel_t* pData = img.getData();
				size_t count = img.pixelCount();
				#pragma omp parallel for
				for(index_t n = 0; n < count; n++)
				{
					pData[n] = toLittleEndian(pData[n]);
				}
            }
            */
		}

		namespace internals
		{
			/**
			Add dimensions to file name.
			*/
			inline std::string concatDimensions(const std::string& baseName, const math::Vec3c& dimensions)
			{
				std::stringstream name;
				name << baseName << "_" << dimensions.x << "x" << dimensions.y << "x" << dimensions.z << ".raw";
				return name.str();
			}

			/**
			* Add image dimensions to file name.
			*/
			template<typename pixel_t> std::string concatDimensions(const std::string& baseName, const Image<pixel_t>& img)
			{
				return concatDimensions(baseName, img.dimensions());
			}

		}
		
		/**
		Reads raw file to given image, initializes the image to correct size read from file name.
		The file name must be in format image_name_100x200x300.raw or image_name_100x200.raw.
		@param img Image where the data is placed. The size of the image will be set based on the .raw file name.
		@param filename The name of the file to read.
		@param bytesToSkip Skip this many bytes from the beginning of the file.
		*/
		template<typename pixel_t> void readd(Image<pixel_t>& img, std::string filename, size_t bytesToSkip = 0)
		{
			math::Vec3c dimensions;
			ImageDataType dt;

			if (!fileExists(filename))
			{
				vector<string> candidates = sequence::internals::buildFileList(filename + "*.raw");

				if (candidates.size() == 0)
					throw ITLException(string("No file found matching to template ") + filename);

				if (candidates.size() == 1)
					filename = candidates[0];
				else
					throw ITLException(string("Multiple .raw files match the template ") + filename);
			}

			if (!internals::parseDimensions(filename, dimensions, dt))
				throw ITLException(std::string("No dimensions found from file name ") + filename);

			img.ensureSize(dimensions);
			read(img, filename);
		}


		/**
		Writes a block of image to specified location in a raw file.
		The output file is not truncated if it exists.
		If the output file does not exist, it is created.
		Part of image extending beyond [0, fileDimensions[ is not written.
		@param img Image to write.
		@param filename Name of file to write.
		@param filePosition Position in the file to write to.
		@param fileDimension Total dimensions of the output file.
		@param imagePosition Position in the image where the block to be written starts.
		@param imageDimensions Dimensions of the block of the source image to write.
		@param showProgressInfo Set to true to show a progress bar.
		*/
		template<typename pixel_t> void writeBlock(const Image<pixel_t>& img, const std::string& filename, const math::Vec3c& filePosition, const math::Vec3c& fileDimensions,
			const math::Vec3c& imagePosition, const math::Vec3c& imageDimensions,
			bool showProgressInfo = false)
			//bool waitOtherWriters = true)
		{
			createFoldersFor(filename);

			// Create file if it does not exist, otherwise set file size to correct value.
			setFileSize(filename, fileDimensions.x * fileDimensions.y * fileDimensions.z * sizeof(pixel_t));

			ofstream out(filename.c_str(), ios_base::in | ios_base::out | ios_base::binary);

			//if (!out.good())
			if(!out)
				throw ITLException(std::string("Write raw block: unable to write to ") + filename);

			math::Vec3c cStart = filePosition;
			clamp(cStart, math::Vec3c(0, 0, 0), fileDimensions);
			math::Vec3c cEnd = filePosition + imageDimensions;
				//img.dimensions();
			clamp(cEnd, math::Vec3c(0, 0, 0), fileDimensions);

			if (!img.isInImage(imagePosition))
				throw ITLException("Block start position must be inside the image.");
			if (!img.isInImage(imagePosition + imageDimensions - math::Vec3c(1, 1, 1)))
				throw ITLException("Block end position must be inside the image.");

			const pixel_t* pBuffer = img.getData();
			for (coord_t z = cStart.z; z < cEnd.z; z++)
			{
				for (coord_t y = cStart.y; y < cEnd.y; y++)
				{
					size_t filePos = (z * fileDimensions.x * fileDimensions.y + y * fileDimensions.x + cStart.x) * sizeof(pixel_t);
					out.seekp(filePos);

					size_t imgPos = img.getLinearIndex(imagePosition.x, y - cStart.y + imagePosition.y, z - cStart.z + imagePosition.z);
					out.write((char*)&pBuffer[imgPos], (cEnd.x - cStart.x) * sizeof(pixel_t));

					//if (!out.good())
					if(!out)
						throw ITLException(std::string("Write raw block: unable to write to ") + filename);
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
		Write image to .raw file.
		The pixels are written in native byte order.
		@param img Image to write.
		@param filename Name of file to write.
		@param truncate Set to false to add to existing file.
		*/
		template<typename pixel_t> void write(const Image<pixel_t>& img, const std::string& filename, bool truncate = true)
		{
			createFoldersFor(filename);

			ios::openmode mode;
			if(truncate)
				mode = ios_base::out | ios_base::trunc | ios_base::binary;
			else
				mode = ios_base::out | ios_base::app | ios_base::binary;

			ofstream out(filename.c_str(), mode);

			// If the file does not exist, create it.
			if(!out.is_open())
				out.open(filename.c_str(), ios_base::out | ios_base::trunc | ios_base::binary);

            //if(!out.good())
			if(!out)
                throw ITLException(std::string("Write raw: unable to write to ") + filename);


			// Write directly from the native buffer
			// Writing everything at once does not work for big files.

			size_t WRITE_SIZE = 200 * 1024 * 1024;
			size_t write_start = 0;
			const char* pBuffer = (const char*)img.getData();
			while(write_start < img.pixelCount() * sizeof(pixel_t))
			{
				size_t write_size = WRITE_SIZE;
				if(write_start + write_size > img.pixelCount() * sizeof(pixel_t))
					write_size = img.pixelCount() * sizeof(pixel_t) - write_start;

				out.write(pBuffer + write_start, write_size);
				write_start += write_size;
			}
		}

		/**
		 Write image to .raw file.
		 Concatenates image dimensions to file name automatically.
		 @param img Image to write.
		 @param filename Template of file name. The full file name will be [template]_[width]x[height]x[depth].raw.
		 @param truncate Set to false to add to existing file.
		 @return The name of the file that was written.
		 */
		template<typename pixel_t> std::string writed(const Image<pixel_t>& img, const std::string& filename, bool truncate = true)
		{
			std::string cf = internals::concatDimensions(filename, img);
			write(img, cf, truncate);
			return cf;
		}

		/*
		 * Write image to .raw file.
		 * @param r, g, b The red, green and blue component images.
		 */
		//template<typename pixel_t> void writeRGB(const Image<pixel_t>& r, const Image<pixel_t>& g, const Image<pixel_t>& b, const std::string& filename)
		//{
		//	r.checkSize(g);
		//	r.checkSize(b);

		//	ofstream out(filename.c_str(), ios_base::out | ios_base::trunc | ios_base::binary);

		//	for(coord_t n = 0; n < r.pixelCount(); n++)
		//	{
		//		out.write((const char*)&r(n), sizeof(pixel_t));
		//		out.write((const char*)&g(n), sizeof(pixel_t));
		//		out.write((const char*)&b(n), sizeof(pixel_t));
		//	}
		//}


		namespace tests
		{
			void raw();
			void readWriteBlock();
		}
	}

}
