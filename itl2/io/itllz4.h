#pragma once

#include <iostream>

#include "image.h"
#include "byteorder.h"
#include "transform.h"

#include "lz4/lz4frame.h"

namespace itl2
{
	namespace lz4
	{

		namespace internals
		{
			/**
			Writes the given value into the stream in little-endian format.
			*/
			template<typename T> void writeSafe(std::ofstream& out, T val)
			{
				if (isBigEndian())
					swapByteOrder(val);

				out.write((char*)&val, sizeof(T));
			}

			/**
			Reads a value from the given file and swaps endianness if the computer is big endian.
			*/
			template<typename T> T readSafe(std::ifstream& in)
			{
				T result;
				in.read((char*)&result, sizeof(T));
				if (isBigEndian())
					swapByteOrder(result);
				return result;
			}

			static const size_t LZ4_CHUNK_SIZE = 64 * 1024;

			static const LZ4F_preferences_t lz4Prefs = {
				{ LZ4F_max256KB,
					LZ4F_blockLinked,
					LZ4F_noContentChecksum,
					LZ4F_frame,
					0, // unknown content size
					0, // no dictID
					LZ4F_noBlockChecksum },
				0,   // compression level; 0 == default
				0,   // autoflush
				0,   // favor decompression speed
				{ 0, 0, 0 },  // reserved, must be set to 0
			};

			/**
			Get LZ4Raw info from file stream.
			*/
			bool getInfo(std::ifstream& in, Vec3c& dimensions, ImageDataType& dataType, string& reason);

			/**
			Decompress LZ4 compressed bytes from given stream to the dst buffer.
			Dst buffer size is dstSizeBytes.
			File name is used only for error messages.
			*/
			void decompress(std::ifstream& in, uint8_t* dst, size_t dstSizeBytes, const string& filenameForErrorMessages);

		}

		/**
		Get information of .lz4raw image file.
		@param dimensions Dimensions of the image
		@param dataType Pixel data type of the image.
		@return True if the file seems to be an existing, valid .lz4raw file with supported pixel data type.
		*/
		bool getInfo(const std::string& filename, Vec3c& dimensions, ImageDataType& dataType, string& reason);

		/**
		Tests if the given path points to a .lz4raw file.
		*/
		bool isFile(const std::string& filename);
		
		
		/**
		Reads an .lz4raw image from disk.
		*/
		template<typename pixel_t> void read(Image<pixel_t>& target, const std::string& filename)
		{
std::cout << "DEBUG: LZ4 read file " << filename << std::endl;

			// This is required in some Linux systems to differentiate files from directories.
			if (!fs::is_regular_file(filename))
				throw ITLException("Not a file.");

			std::ifstream in(filename.c_str(), std::ios_base::in | std::ios_base::binary);

			if (!in)
				throw ITLException(std::string("Unable to open file. ") + getStreamErrorMessage());

			Vec3c dimensions;
			ImageDataType fileDT;
			std::string reason;
			if (!internals::getInfo(in, dimensions, fileDT, reason))
				throw ITLException(reason);

			if (fileDT != target.dataType())
				throw ITLException(std::string("Image data type is ") + toString(target.dataType()) + std::string(" but the file contains data of type ") + toString(fileDT));

			if (!in)
				throw ITLException(std::string("Unable to read from ") + filename);

			target.ensureSize(dimensions);

			internals::decompress(in, (uint8_t*)target.getData(), target.pixelCount() * target.pixelSize(), filename);
		}

		

		/**
		Reads a part of a .lz4raw file to the given image.
		NOTE: Does not support out of bounds start position.
		@param img Image where the data is placed. The size of the image defines the size of the block that is read.
		@param filename The name of the file to read.
		@param filePos Start location of the read. The size of the image defines the size of the block that is read.
		@param temp Temporary image. If this image has the same dimensions than the file, no temporary memory allocations for image data are made. At output, this image will contain the entire decompressed file.
		*/
		template<typename pixel_t> void readBlock(Image<pixel_t>& img, std::string filename, const Vec3c& filePos, Image<pixel_t>& temp)
		{
std::cout << "DEBUG: LZ4 read block from file " << filename << ", filePos = " << filePos << ", size = " << img.dimensions() << std::endl;

			Vec3c fileDimensions;
			ImageDataType fileDT;
			std::string reason;
			if (!getInfo(filename, fileDimensions, fileDT, reason))
				throw ITLException(reason);

			if (filePos == Vec3c(0, 0, 0) && img.dimensions() == fileDimensions)
			{
				// Read the entire file.
				read(img, filename);
				return;
			}

			// Read and crop, as the compressed data must be decompressed entirely before getting access to the required block.
			read(temp, filename);
			crop(temp, img, filePos);
		}

        /**
		Reads a part of a .lz4raw file to the given image.
		NOTE: Does not support out of bounds start position.
		@param img Image where the data is placed. The size of the image defines the size of the block that is read.
		@param filename The name of the file to read.
		@param filePos Start location of the read. The size of the image defines the size of the block that is read.
		*/
		template<typename pixel_t> void readBlock(Image<pixel_t>& img, std::string filename, const Vec3c& filePos)
		{
			Image<pixel_t> temp;
			readBlock(img, filename, filePos, temp);
		}

		/**
		Writes an image to an .lz4raw file.
		@param img Image to write.
		@param filename Name of file to write.
		*/
		template<typename pixel_t> void write(const Image<pixel_t>& source, const std::string& filename)
		{
std::cout << "DEBUG: LZ4 writing file " << filename << std::endl;

			createFoldersFor(filename);

			std::ofstream out(filename, std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);

			if (!out)
				throw ITLException(std::string("Unable to write to ") + filename + std::string(", ") + getStreamErrorMessage());

			
			internals::writeSafe(out, source.width());
			internals::writeSafe(out, source.height());
			internals::writeSafe(out, source.depth());
			internals::writeSafe(out, (int32_t)source.dataType());

			LZ4F_cctx* ctx;
			LZ4F_errorCode_t err = LZ4F_createCompressionContext(&ctx, LZ4F_VERSION);
			if (LZ4F_isError(err))
				throw ITLException(string("Unable to create LZ4 compression context: ") + LZ4F_getErrorName(err));
			std::unique_ptr<LZ4F_cctx, decltype(LZ4F_freeCompressionContext)*> pCtx(ctx, LZ4F_freeCompressionContext);
			
			size_t outputCapacity = LZ4F_compressBound(internals::LZ4_CHUNK_SIZE, &internals::lz4Prefs);
			std::unique_ptr<uint8_t[]> pDest = std::make_unique<uint8_t[]>(outputCapacity);

			// Frame header
			{
				size_t headerSize = LZ4F_compressBegin(ctx, pDest.get(), outputCapacity, &internals::lz4Prefs);
				if (LZ4F_isError(headerSize))
					throw ITLException(string("Unable to init LZ4 compression: ") + LZ4F_getErrorName(headerSize));

				out.write((char*)pDest.get(), headerSize);
			}

			// Compress data
			size_t pos = 0;
			size_t imgBytes = source.pixelCount() * source.pixelSize();
			while(pos < imgBytes)
			{
				size_t readSize = internals::LZ4_CHUNK_SIZE;
				if (pos + readSize > imgBytes)
					readSize = imgBytes - pos;
				uint8_t* sourceAddr = (uint8_t*)source.getData() + pos;
				pos += readSize;

				size_t compressedSize = LZ4F_compressUpdate(ctx,
					pDest.get(), outputCapacity,
					sourceAddr, readSize,
					NULL);
				if (LZ4F_isError(compressedSize))
					throw ITLException(string("Unable to perform LZ4 compression: ") + LZ4F_getErrorName(compressedSize));
				
				out.write((char*)pDest.get(), compressedSize);
			}

			// Finalize compression
			{   
				size_t const compressedSize = LZ4F_compressEnd(ctx, pDest.get(), outputCapacity, NULL);
				if (LZ4F_isError(compressedSize))
					throw ITLException(string("Unable to finalize LZ4 compression: ") + LZ4F_getErrorName(compressedSize));

				out.write((char*)pDest.get(), compressedSize);
			}

		}


		namespace internals
		{
			template<typename pixel_t> void writeBlockToFullFile(const Image<pixel_t>& img, const std::string& filename,
				const Vec3c& imagePosition,
				const Vec3c& blockDimensions)
			{
				createFoldersFor(filename);

				std::ofstream out(filename, std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);

				if (!out)
					throw ITLException(std::string("Unable to write to ") + filename + std::string(", ") + getStreamErrorMessage());

				internals::writeSafe(out, blockDimensions.x);
				internals::writeSafe(out, blockDimensions.y);
				internals::writeSafe(out, blockDimensions.z);
				internals::writeSafe(out, (int32_t)img.dataType());

				LZ4F_cctx* ctx;
				LZ4F_errorCode_t err = LZ4F_createCompressionContext(&ctx, LZ4F_VERSION);
				if (LZ4F_isError(err))
					throw ITLException(string("Unable to create LZ4 compression context: ") + LZ4F_getErrorName(err));
				std::unique_ptr<LZ4F_cctx, decltype(LZ4F_freeCompressionContext)*> pCtx(ctx, LZ4F_freeCompressionContext);

				size_t lzChunkSize = internals::LZ4_CHUNK_SIZE;
				if (lzChunkSize < (size_t)blockDimensions.x)
					lzChunkSize = (size_t)blockDimensions.x;

				size_t outputCapacity = LZ4F_compressBound(lzChunkSize, &internals::lz4Prefs);
				std::unique_ptr<uint8_t[]> pDest = std::make_unique<uint8_t[]>(outputCapacity);

				// Frame header
				{
					size_t headerSize = LZ4F_compressBegin(ctx, pDest.get(), outputCapacity, &internals::lz4Prefs);
					if (LZ4F_isError(headerSize))
						throw ITLException(string("Unable to init LZ4 compression: ") + LZ4F_getErrorName(headerSize));

					out.write((char*)pDest.get(), headerSize);
				}

				// Compress data one x-directional scan line at time.
				for (coord_t z = imagePosition.z; z < imagePosition.z + blockDimensions.z; z++)
				{
					for (coord_t y = imagePosition.y; y < imagePosition.y + blockDimensions.y; y++)
					{
						size_t readSize = blockDimensions.x * img.pixelSize();
						uint8_t* sourceAddr = (uint8_t*)&img(imagePosition.x, y, z);

						size_t compressedSize = LZ4F_compressUpdate(ctx,
							pDest.get(), outputCapacity,
							sourceAddr, readSize,
							NULL);
						if (LZ4F_isError(compressedSize))
							throw ITLException(string("Unable to perform LZ4 compression: ") + LZ4F_getErrorName(compressedSize));

						out.write((char*)pDest.get(), compressedSize);
					}
				}

				// Finalize compression
				{
					size_t const compressedSize = LZ4F_compressEnd(ctx, pDest.get(), outputCapacity, NULL);
					if (LZ4F_isError(compressedSize))
						throw ITLException(string("Unable to finalize LZ4 compression: ") + LZ4F_getErrorName(compressedSize));

					out.write((char*)pDest.get(), compressedSize);
				}
			}
		}

		/**
		Writes block of an image into a .lz4raw file.
		Replaces existing file.
		If the output file exists, writing a block that does not cover it entirely, triggers reading of the entire file to a temporary image.
		@param img Image to write.
		@param filename Name of file to write.
		@param filePosition Position in the file where the first pixel should be written.
		@param fileDimensions Dimensions of the entire output file. If zero, 
		*/
		template<typename pixel_t> void writeBlock(const Image<pixel_t>& img, const std::string& filename,
			const Vec3c& filePosition, const Vec3c& fileDimensions,
			const Vec3c& imagePosition,
			const Vec3c& blockDimensions)
		{
std::cout << "DEBUG: LZ4 write block to file " << filename << ", filePos = " << filePosition << ", size = " << blockDimensions << std::endl;


			if (!img.isInImage(imagePosition))
				throw ITLException("Block start position must be inside the image.");
			if (!img.isInImage(imagePosition + blockDimensions - Vec3c(1, 1, 1)))
				throw ITLException("Block end position must be inside the image.");

			
			if (filePosition == Vec3c(0, 0, 0) && (fileDimensions == blockDimensions || fileDimensions == Vec3c(0, 0, 0)))
			{
				// We are overwriting entire output file if it exists.
				internals::writeBlockToFullFile(img, filename, imagePosition, blockDimensions);
			}
			else
			{
				// We are writing a block of the image to a block of the output file.
				// As the file is compressed, we must read the entire file, modify the data, and write the entire file again.
				Image<pixel_t> temp;
				if (isFile(filename))
				{
					// We have existing data, so read it.
					lz4::read(temp, filename);
					if (fileDimensions != Vec3c(0, 0, 0))
						temp.ensureSize(fileDimensions);
				}
				else
				{
					// No existing data; just ensure that temp image size is correct.
					if (fileDimensions != Vec3c(0, 0, 0))
						temp.ensureSize(fileDimensions);
					else
						temp.ensureSize(filePosition + blockDimensions);
				}
				copyValues(temp, img, filePosition, imagePosition, blockDimensions);
				lz4::write(temp, filename);
			}

			
		}


		/**
		Reads an lz4raw image from the disk.
		Appends .lz4raw suffix to the file name if the file does not exist and does not have that suffix.
		*/
		template<typename pixel_t> void readd(Image<pixel_t>& target, const std::string& filename)
		{
			if (fs::exists(filename))
				read(target, filename);
			else if (fs::exists(filename + ".lz4raw"))
				read(target, filename + ".lz4raw");
			else
				throw ITLException(string("File not found: ") + filename);
		}

		/**
		Writes an image to an .lz4raw file.
		@param img Image to write.
		@param filename Template of file name. The full file name will be [template].lz4, if the template does not contain the .lz4 suffix.
		*/
		template<typename pixel_t> void writed(const Image<pixel_t>& img, const std::string& filename)
		{
			if (endsWithIgnoreCase(filename, ".lz4raw"))
				write(img, filename);
			else
				write(img, filename + ".lz4raw");
		}

		namespace tests
		{
			void lz4io();
			void lz4blockIo();
		}
	}

	
}
