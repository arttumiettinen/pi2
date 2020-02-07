#pragma once

#include <string>

#include "math/vec3.h"
#include "io/imagedatatype.h"
#include "stringutils.h"
#include "image.h"
#include "io/raw.h"
#include "io/io.h"


using itl2::Vec3c;
using itl2::toString;
using itl2::imageDataType;
using itl2::ImageDataType;
using itl2::coord_t;
using itl2::ImageBase;
using itl2::ITLException;
using itl2::endsWith;
using itl2::Image;

namespace pilib
{
	class Distributor;

	/**
	Base class for distributed images.
	*/
	class DistributedImageBase
	{
	private:

		/**
		Dimensions of the image.
		*/
		Vec3c dims;

		/**
		Name of the image in the system, i.e., name of variable holding this image.
		*/
		std::string name;

		/**
		Unique name of this image, i.e., any other image will not have the same unique name in a PISystem session.
		*/
		std::string uniqName;

		/**
		Filename where the image data should be read.
		*/
		std::string readSource;

		/**
		Filename where modified image data should be saved.
		*/
		std::string writeTarget;

		/**
		Names of temporary files where image data is stored temporarily.
		These files are system-generated, so they will be deleted when the image is deleted.
		*/
		std::string tempFilename1, tempFilename2;

		/**
		Flag that indicates that this image is new (should not be read from temp file)
		*/
		bool isNewImage;

		/**
		Pixel data type
		*/
		ImageDataType pixelDataType;

		/**
		Generate filename for temporary storage.
		If source file is raw, the temp file is raw.
		Otherwise, the temp file is folder containing an image sequence.
		*/
		void createTempFilenames();

		/**
		Changes the location where the image is read from.
		@param check Set to false if calling from a constructor to avoid pure virtual call.
		*/
		void setReadSource(const std::string& filename, bool check);

		/**
		The distributor that owns this image.
		*/
		Distributor* distributor;

	protected:
		/**
		Calls distributor.flush()
		*/
		void flush() const;

	public:

        /**
        Returns true if the image has been saved to disk (it is not just created
        unprocessed new image)
        */
		bool isSavedToDisk() const
		{
			return !isNewImage;
		}

		/**
		Creates distributed image that points to given file.
		*/
		DistributedImageBase(Distributor& distributor, const std::string& name, const Vec3c& dimensions, ImageDataType dataType, const std::string& sourceFilename);

		/**
		Creates distributed image using temporary file as storage.
		*/
		DistributedImageBase(Distributor& distributor, const std::string& name, const Vec3c& dimensions, ImageDataType dataType) :
			DistributedImageBase(distributor, name, dimensions, dataType, "")
		{

		}

		virtual ~DistributedImageBase();

		/**
		Returns the dimensionality of the image.
		*/
		size_t dimensionality() const
		{
			if (dims.x <= 1 && dims.y <= 1 && dims.z <= 1)
				return 0;
			if (dims.x > 1 && dims.y <= 1 && dims.z <= 1)
				return 1;
			if (dims.x > 1 && dims.y > 1 && dims.z <= 1)
				return 2;

			return 3;
		}

		/**
		Gets dimensions of this image.
		*/
		Vec3c dimensions() const
		{
			return dims;
		}

		coord_t width() const
		{
			return dims.x;
		}

		coord_t height() const
		{
			return dims.y;
		}

		coord_t depth() const
		{
			return dims.z;
		}

		/**
		Gets count of pixels in the image.
		*/
		coord_t pixelCount() const
		{
			return dims.x * dims.y * dims.z;
		}

		/**
		Gets data type of this image.
		*/
		ImageDataType dataType() const
		{
			return pixelDataType;
		}

		/**
		Gets pixel size in bytes.
		*/
		virtual size_t pixelSize() const = 0;

		/**
		Gets piece of pi2 code to read a block of this image.
		@param dataNeeded Set to true if the image is used as input data. In that case the image data should be read from disk.
		*/
		std::string emitReadBlock(const Vec3c& filePos, const Vec3c& blockSize, bool dataNeeded) const;

		/**
		Gets piece of pi2 code to write a block of this image.
		*/
		std::string emitWriteBlock(const Vec3c& filePos, const Vec3c& imagePos, const Vec3c& blockSize);

		/**
		Call when all blocks of this image have been written.
		*/
		void writeComplete();

		/**
		Gets the file path where the image data should be read.
		*/
		std::string currentReadSource() const
		{
			return readSource;
		}

		/**
		Gets the file path where the changed image will be saved.
		*/
		std::string currentWriteTarget() const
		{
			return writeTarget;
		}

		/**
		Enures that image writing process does not override old data (that might be needed by another processes if doing distribution with overlapping blocks).
		Swaps between two write targets.
		*/
		void newWriteTarget();

		/**
		Gets a value indicating whether the current file read location is a temporary file.
		*/
		bool isSavedToTemp() const
		{
			return currentReadSource() == tempFilename1 || currentReadSource() == tempFilename2;
		}

		/**
		Changes the location where the image is read from.
		*/
		void setReadSource(const std::string& filename)
		{
			setReadSource(filename, true);
		}
		
		/**
		Changes the file or location where the image is saved to.
		*/
		void setWriteTarget(const std::string& filename)
		{
		    writeTarget = filename;   
		}

		/**
		Gets name of variable that stores this image.
		*/
		const std::string& varName() const
		{
			return name;
		}

		/**
		Gets unique name of this image.
		*/
		const std::string& uniqueName() const
		{
			return uniqName;
		}

		/**
		Creates new Image<pixel_t> object and reads the distributed image data to that image.
		*/
		virtual std::shared_ptr<ImageBase> toNormalImage() const = 0;

		/**
		Set data of this distributed image from given normal image.
		*/
		virtual void setData(const ImageBase* pImage) = 0;

		/**
		Makes sure that the size of this image equals to what is given as an argument.
		*/
		void ensureSize(const Vec3c& newDimensions);
		
		void ensureSize(coord_t w, coord_t h = 1, coord_t d = 1)
		{
			ensureSize(Vec3c(w, h, d));
		}

		void ensureSize(const DistributedImageBase& p)
		{
			ensureSize(p.dimensions());
		}

		/**
		Tests whether the size of this image equals the size of the given image.
		*/
		inline bool sizeEquals(const Vec3c& dims) const
		{
			return dims.equals(dimensions());
		}

		/**
		Tests whether the size of this image and the given image are equal and throws an
		exception if they are not.
		*/
		void checkSize(const Vec3c& r) const
		{
			if (!sizeEquals(r))
				throw ITLException("Image sizes are not equal.");
		}

		/**
		Tests whether the size of this image and the given image are equal and throws an
		exception if they are not.
		*/
		void checkSize(const ImageBase& r) const
		{
			checkSize(r.dimensions());
		}
		
		/**
		Tests if output file is .raw.
		*/
		bool isOutputRaw() const
		{
		    return endsWith(currentWriteTarget(), ".raw"); // The file does not need to exist, and write target is always .raw or sequence.
			//Vec3c dims;
			//ImageDataType dt;
			//return raw::getInfo(currentWriteTarget(), dims, dt);
		}

        /**
        Tests if input file is raw.
        */
		bool isRaw() const
		{
			//return endsWith(currentReadSource(), ".raw");
			Vec3c dims;
			ImageDataType dt;
			std::string reason;
			return itl2::raw::getInfo(currentReadSource(), dims, dt, reason);
		}

        /**
        Tests if input file is sequence.
        */
		bool isSequence() const
		{
			//return !isRaw();
			Vec3c dims;
			ImageDataType dt;
			std::string reason;
			return itl2::sequence::getInfo(currentReadSource(), dims, dt, reason);
		}


	};

	template<typename pixel_t> class DistributedImage : public DistributedImageBase
	{
	public:
		/**
		Creates distributed image whose source points to the given file.
		*/
		DistributedImage(Distributor& distributor, const std::string& name, coord_t width, coord_t height, coord_t depth, const std::string& filename) :
			DistributedImage(distributor, name, Vec3c(width, height, depth), filename)
		{
		}

		DistributedImage(Distributor& distributor, const std::string& name, const Vec3c& dimensions, const std::string& filename) :
			DistributedImageBase(distributor, name, dimensions, imageDataType<pixel_t>(), filename)
		{
		}

		/**
		Creates distributed image without source file.
		*/
		DistributedImage(Distributor& distributor, const std::string& name, coord_t width = 1, coord_t height = 1, coord_t depth = 1) :
			DistributedImage(distributor, name, Vec3c(width, height, depth))
		{
		}

		/**
		Creates distributed image without source file.
		*/
		DistributedImage(Distributor& distributor, const std::string& name, const Vec3c& dimensions) :
			DistributedImageBase(distributor, name, dimensions, imageDataType<pixel_t>())
		{
		}

		virtual size_t pixelSize() const override
		{
			return sizeof(pixel_t);
		}

		/**
		Converts this distributed image to non-distributed image.
		*/
		virtual std::shared_ptr<ImageBase> toNormalImage() const override
		{
			DistributedImageBase::flush();
			std::shared_ptr<itl2::Image<pixel_t>> pNormalImg = std::make_shared<itl2::Image<pixel_t> >(dimensions());
			readTo(*pNormalImg);
			return pNormalImg;
		}

		/**
		Gets value of the first pixel in this image.
		Causes the image to be read to RAM.
		*/
		pixel_t getValue() const
		{
			Image<pixel_t> img(dimensions());
			readTo(img);
			return img(0);
		}

		virtual void setData(const ImageBase* pImage) override
		{
			const Image<pixel_t>* pi = dynamic_cast<const Image<pixel_t>*>(pImage);
			if (!pi)
				throw ITLException("The data type of the normal image is not the same than the data type of the distributed image.");

			setData(*pi);
		}

		/**
		Reads the data of this distributed image to the given normal image.
		*/
		void readTo(itl2::Image<pixel_t>& img) const
		{
			// Flush so that there are no pending writes to this image.
			DistributedImageBase::flush();

			img.ensureSize(dimensions());

			if (isSavedToDisk())
			{
				itl2::io::read(img, currentReadSource());
			}
		}

		/**
		Copies data of given normal image to this distributed image.
		*/
		void setData(const itl2::Image<pixel_t>& img)
		{
			// Flush so that there are no pending accesses to this image.
			DistributedImageBase::flush();

			ensureSize(img.dimensions());

			if (isOutputRaw())
			{
				itl2::raw::write(img, currentWriteTarget());
			}
			else
			{
				itl2::sequence::write(img, currentWriteTarget());
			}

			writeComplete();
		}


		/**
		Throws exception if the given image is the same than this image.
		*/
		template<typename pixel2_t> void mustNotBe(const DistributedImage<pixel2_t>& other) const
		{
			if ((void*)this == (void*)&other)
				throw ITLException("This operation cannot be executed if the parameter images are the same.");
		}
	};
}
