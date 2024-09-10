#pragma once

#include <string>

#include "io/sequence.h"
#include "io/raw.h"
#include "io/vol.h"
#include "io/nrrd.h"
#include "io/pcr.h"
#include "io/nn5.h"
#include "io/zarr.h"

namespace itl2
{
	namespace io
	{
		namespace internals
		{
			inline std::string combineReasons(
				const std::string& rawReason,
				const std::string& tiffReason,
				const std::string& sequenceReason,
				const std::string& volReason,
				const std::string& nrrdReason,
				const std::string& pcrReason,
				const std::string& nn5Reason,
				const std::string& zarrReason,
				const std::string& lz4Reason,
				const std::string& dicomReason)
			{
				return std::string() +
					"raw: " + rawReason + "\n" +
					"tiff: " + tiffReason + "\n" +
					"sequence: " + sequenceReason + "\n" +
					"vol: " + volReason + "\n" +
					"nrrd: " + nrrdReason + "\n" +
					"pcr: " + pcrReason + "\n" +
					"nn5: " + nn5Reason + "\n" +
					"zarr: " + zarrReason + "\n" +
					"lz4raw: " + lz4Reason + "\n" +
					"DICOM: " + dicomReason;
			}

			inline std::string combineReasonsNoDICOM(
				const std::string& rawReason,
				const std::string& tiffReason,
				const std::string& sequenceReason,
				const std::string& volReason,
				const std::string& nrrdReason,
				const std::string& pcrReason,
				const std::string& nn5Reason,
				const std::string& zarrReason,
				const std::string& lz4Reason)
			{
				return std::string() +
					"raw: " + rawReason + "\n" +
					"tiff: " + tiffReason + "\n" +
					"sequence: " + sequenceReason + "\n" +
					"vol: " + volReason + "\n" +
					"nrrd: " + nrrdReason + "\n" +
					"pcr: " + pcrReason + "\n" +
					"nn5: " + nn5Reason + "\n" +
					"zarr: " + zarrReason + "\n" +
					"lz4raw: " + lz4Reason;
			}
		}

		/**
		Finds data format of given file and reads it to the given image.
		The data type of the target image must be correct but its size is set automatically.
		*/
		template<typename pixel_t> void read(Image<pixel_t>& img, const std::string& filename)
		{
			Vec3c dimensions;
			ImageDataType dt;
			
			std::string volReason, tiffReason, nrrdReason, sequenceReason, rawReason, pcrReason, nn5Reason, zarrReason, lz4Reason, dicomReason;
			if (vol::getInfo(filename, dimensions, dt, volReason))
			{
				vol::read(img, filename);
			}
			else if (tiff::getInfo(filename, dimensions, dt, tiffReason))
			{
				tiff::read(img, filename);
			}
			else if (nrrd::getInfo(filename, dimensions, dt, nrrdReason))
			{
				nrrd::read(img, filename);
			}
			else if (sequence::getInfo(filename, dimensions, dt, sequenceReason))
			{
				sequence::read(img, filename);
			}
			else if (pcr::getInfo(filename, dimensions, dt, pcrReason))
			{
				pcr::read(img, filename);
			}
			else if (raw::getInfo(filename, dimensions, dt, rawReason))
			{
				raw::read(img, filename);
			}
			else if (nn5::getInfo(filename, dimensions, dt, nn5Reason))
			{
				nn5::read(img, filename);
			}
			else if (zarr::getInfo(filename, dimensions, dt, zarrReason))
			{
				zarr::read(img, filename);
			}
			else if (lz4::getInfo(filename, dimensions, dt, lz4Reason))
			{
				lz4::read(img, filename);
			}
			else if (dicom::getInfo(filename, dimensions, dt, dicomReason))
			{
				dicom::read(img, filename);
			}
			else
			{
				throw ITLException(std::string("Unsupported file type, file not found, or cannot be read: ") + filename + "\n" +
					internals::combineReasons(rawReason, tiffReason, sequenceReason, volReason, nrrdReason, pcrReason, nn5Reason, zarrReason, lz4Reason, dicomReason));
			}
		}

		/**
		Parses the given file and finds the dimensions and pixel data type of image stored in the file.
		@return True if the parsing succeeded, false if the file cannot be loaded.
		*/
		bool getInfo(const std::string& filename, Vec3c& dimensions, ImageDataType& dataType, std::string& reason);

		/**
		Finds data format of given file and reads part of it to the given image.
		The data type of the target image must be correct but its size is set automatically.
		The size of the block to read is determined by the size of the image.
		@param blockStart The start position of the block.
		*/
		template<typename pixel_t> void readBlock(Image<pixel_t>& img, const std::string& filename, const Vec3c& blockStart)
		{
			Vec3c dimensions;
			ImageDataType dt;

			std::string volReason, tiffReason, nrrdReason, sequenceReason, rawReason, pcrReason, nn5Reason, zarrReason, lz4Reason;
			if (vol::getInfo(filename, dimensions, dt, volReason))
			{
				vol::readBlock(img, filename, blockStart);
			}
			else if (tiff::getInfo(filename, dimensions, dt, tiffReason))
			{
				tiff::readBlock(img, filename, blockStart);
			}
			else if (nrrd::getInfo(filename, dimensions, dt, nrrdReason))
			{
				throw ITLException(std::string("The file ") + filename + " is identified to be an NRRD file, but reading blocks of NRRD files is not supported at the moment.");
			}
			else if (sequence::getInfo(filename, dimensions, dt, sequenceReason))
			{
				sequence::readBlock(img, filename, blockStart);
			}
			else if (pcr::getInfo(filename, dimensions, dt, pcrReason))
			{
				pcr::readBlock(img, filename, blockStart);
			}
			else if (raw::getInfo(filename, dimensions, dt, rawReason))
			{
				raw::readBlock(img, filename, blockStart);
			}
			else if (nn5::getInfo(filename, dimensions, dt, nn5Reason))
			{
				nn5::readBlock(img, filename, blockStart);
			}
			else if (zarr::getInfo(filename, dimensions, dt, zarrReason))
			{
				zarr::readBlock(img, filename, blockStart);
			}
			else if (lz4::getInfo(filename, dimensions, dt, lz4Reason)) //TODO: is this supposed to be lz4Reason?
			{
				lz4::readBlock(img, filename, blockStart);
			}
			// TODO: DICOM read block
			else
			{
				throw ITLException(std::string("Unsupported file type, file not found, or cannot be read: ") + filename + "\n" +
					internals::combineReasonsNoDICOM(rawReason, tiffReason, sequenceReason, volReason, nrrdReason, pcrReason, nn5Reason, zarrReason, lz4Reason));
			}
		}

		namespace tests
		{
			void readWrite();
			void badnn5();
		}
	}
}
