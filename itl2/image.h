#pragma once

#include <omp.h>

#include "buildsettings.h"
#include "datatypes.h"
#include "math/mathutils.h"
#include "math/vec2.h"
#include "math/vec3.h"
#include "test.h"
#include "memorybuffer.h"
#include "diskmappedbuffer.h"
#include "io/imagedatatype.h"

namespace itl2
{

	/**
	Base class for all images.
	*/
	class ImageBase
	{
	protected:

		math::Vec3c dims;

	public:
		virtual ~ImageBase()
		{

		}

		/**
		Return raw pointer to the image data.
		*/
		virtual void* getRawData() = 0;

		/**
		Get data type of this image.
		*/
		virtual ImageDataType dataType() const = 0;

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
		Returns dimensions of the image.
		*/
		math::Vec3c dimensions() const
		{
			return dims;
		}

		/**
		Returns width of the image.
		*/
		coord_t width() const
		{
			return dims.x;
		}

		/**
		Returns height of the image.
		*/
		coord_t height() const
		{
			return dims.y;
		}

		/**
		Returns depth of the image.
		*/
		coord_t depth() const
		{
			return dims.z;
		}

		/**
		Gets dimension (width, height or depth) of the image in some direction.
		*/
		coord_t dimension(size_t dim) const
		{
			return dims[dim];
		}

		/**
		Gets count of pixels in the image.
		*/
		coord_t pixelCount() const
		{
			return width() * height() * depth();
		}

		/**
		Removes singleton dimensions.
		*/
		void squeeze()
		{
			if (width() <= 1)
			{
				dims.x = dims.y;
				dims.y = dims.z;
				dims.z = 1;
			}

			if (dims.y <= 1)
			{
				dims.y = dims.z;
				dims.z = 1;
			}
		}

		/**
		Tests whether the size of this image equals the size of the given image.
		*/
		inline bool sizeEquals(const math::Vec3c& dims) const
		{
			return dims.equals(dimensions());
		}
	};

	/*
	Converts linear index to coordinates.
	*/
	inline math::Vec3c indexToCoords(coord_t ind, const math::Vec3c& size)
	{
		coord_t z = ind / (size[0] * size[1]);
		coord_t y = (ind - z * size[0] * size[1]) / size[0];
		coord_t x = ind - z * size[0] * size[1] - y * size[0];

		return math::Vec3c(x, y, z);
	}

	/**
	0-, 1-, 2- or 3-dimensional image.
	*/
	template<typename pixel_t> class Image : public ImageBase
	{
	private:

		/**
		Pointer to storage for image data.
		This points to a buffer managed by *pBufferObject and is used to access image data.
		*/
		pixel_t* pData;
		const pixel_t* pDataConst;

		/**
		Pointer to buffer object.
		*/
		Buffer<pixel_t>* pBufferObject;

		/*
		File mapped to this image.
		*/
		string mapFile;

		/**
		Used in constructors to allocate memory.
		Does not set pixel values.
		*/
		void initBuffer(coord_t width, coord_t height, coord_t depth)
		{
			dims.x = math::max<coord_t>(1, width);
			dims.y = math::max<coord_t>(1, height);
			dims.z = math::max<coord_t>(1, depth);

			if (mapFile.length() <= 0)
			{
				// Create memory buffer
				pBufferObject = new MemoryBuffer<pixel_t>(pixelCount());
			}
			else
			{
				// Create buffer mapped to file
				pBufferObject = new DiskMappedBuffer<pixel_t>(pixelCount(), mapFile, 0);
			}

			pData = pBufferObject->getBufferPointer();
			pDataConst = pData;
		}

		/**
		Used in constructors to allocate memory.
		Sets pixel values to the given initial value.
		*/
		void initBuffer(coord_t width, coord_t height, coord_t depth, const pixel_t initialValue)
		{
			initBuffer(width, height, depth);

			// Zero memory
#pragma omp parallel for if(pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
			for (coord_t n = 0; n < pixelCount(); n++)
			{
				new (&pData[n]) pixel_t();
				pData[n] = initialValue;
			}
		}


		/**
		Disable copy constructor
		*/
		Image(const Image<pixel_t>& right);

		/**
		Disable assignment operator.
		*/
		const Image<pixel_t>& operator=(const Image<pixel_t>& rhs);

	public:

		/**
		Default constructor, creates memory-resident image.
		*/
		explicit Image()
		{
			initBuffer(1, 1, 1, pixel_t());
		}

		/**
		Constructor, creates memory-resident image.
		*/
		Image(coord_t width, coord_t height = 0, coord_t depth = 0, const pixel_t val = pixel_t())
		{
			initBuffer(width, height, depth, val);
		}

		/**
		Constructor, creates memory-resident image.
		*/
		Image(const math::Vec3c& dimensions, const pixel_t val = pixel_t())
		{
			initBuffer(dimensions.x, dimensions.y, dimensions.z, val);
		}

		/**
		Constructor, creates disk-mapped image.
		*/
		Image(const string& filename, coord_t width = 0, coord_t height = 0, coord_t depth = 0)
		{
			this->mapFile = filename;
			initBuffer(width, height, depth);
		}

		/**
		Constructor, creates image that points to a z-range in another image.
		@param startZ z-coordinate of the first slice to include in the view.
		@param endZ z-coordinate of the last slice to include in the view.
		*/
		Image(Image<pixel_t>& source, coord_t startZ, coord_t endZ)
		{
			pBufferObject = 0;
			init(source, startZ, endZ);
		}

		/**
		Constructor, creates image that points to a z-range in another image.
		@param startZ z-coordinate of the first slice to include in the view.
		@param endZ z-coordinate of the last slice to include in the view.
		*/
		Image(const Image<pixel_t>& source, coord_t startZ, coord_t endZ)
		{
			pBufferObject = 0;
			init(source, startZ, endZ);
		}

		/**
		Destructor
		*/
		virtual ~Image()
		{
			deleteData();
		}

		/**
		Re-init the image to specified size and value.
		*/
		void init(coord_t width, coord_t height = 0, coord_t depth = 0, const pixel_t val = pixel_t())
		{
			deleteData();
			initBuffer(width, height, depth, val);
		}

		void init(const math::Vec3c& dims)
		{
			init(dims.x, dims.y, dims.z);
		}

		/**
		Re-init the image to point to specified buffer file.
		*/
		void init(const string& filename, coord_t width, coord_t height = 0, coord_t depth = 0)
		{
			deleteData();
			this->mapFile = filename;
			initBuffer(width, height, depth);
		}

		void init(Image<pixel_t>& source, coord_t startZ, coord_t endZ)
		{
			if (startZ < 0 || endZ < 0 || startZ >= source.depth() || endZ >= source.depth() || startZ > endZ)
				throw ITLException("Invalid z range.");

			deleteData();

			dims.x = source.width();
			dims.y = source.height();
			dims.z = endZ - startZ + 1;
			pBufferObject = 0;
			pData = &source(0, 0, startZ);
			pDataConst = pData;
		}

		void init(const Image<pixel_t>& source, coord_t startZ, coord_t endZ)
		{
			if (startZ < 0 || endZ < 0 || startZ >= source.depth() || endZ >= source.depth() || startZ > endZ)
				throw ITLException("Invalid z range.");

			deleteData();

			dims.x = source.width();
			dims.y = source.height();
			dims.z = endZ - startZ + 1;
			pBufferObject = 0;
			pData = 0;
			pDataConst = &source(0, 0, startZ);
		}

		/**
		Delete image data resident in memory.
		Use this to free large images before normal destruction (stack walk) takes place.
		*/
		void deleteData()
		{
			if (pBufferObject)
			{
#pragma omp parallel for if(pixelCount() > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
				for (coord_t n = 0; n < pixelCount(); n++)
					pData[n].~pixel_t();

				pData = 0;
				pDataConst = 0;
				delete pBufferObject;
				pBufferObject = 0;
				// Don't reset mapFile here so that init methods can re-init to same file.
			}
		}

		/*
		Lets the image know that data in cube [start, end] will be processed soon.
		*/
		//void prefetch(math::Vec3c start, math::Vec3c end) const
		//{
		//	clamp(start, math::Vec3c(0, 0, 0), dimensions() - math::Vec3c(1, 1, 1));
		//	clamp(end, math::Vec3c(0, 0, 0), dimensions() - math::Vec3c(1, 1, 1));

		//	for (coord_t z = start.z; z <= end.z; z++)
		//	{
		//		for (coord_t y = start.y; y <= end.y; y++)
		//		{
		//			pBufferObject->prefetch(getLinearIndex(start.x, y, z), getLinearIndex(end.x, y, z));
		//		}
		//	}
		//}

		

		virtual void* getRawData()
		{
			return getData();
		}

		virtual ImageDataType dataType() const
		{
			return imageDataType<pixel_t>();
		}

		/**
		Gets pointer to the pixel data.
		*/
		pixel_t* getData()
		{
			return pData;
		}

		const pixel_t* getData() const
		{
			return pDataConst;
		}

		/**
		Get index of pixel at (x, y, z) in the array returned by getData() method.
		*/
		size_t getLinearIndex(coord_t x, coord_t y = 0, coord_t z = 0) const
		{
			size_t ind = (size_t)z * (size_t)width() * (size_t)height() + (size_t)y * (size_t)width() + (size_t)x;
#if _DEBUG
			//if (x < 0 || y < 0 || z < 0 || x >= w || y >= h || z >= d)
			if(ind >= (size_t)pixelCount())
				throw runtime_error("Bounds check failure.");
#endif
			return ind;
		}

		/*
		Gets index of pixel at p in the array returned by the getData() method.
		*/
		size_t getLinearIndex(const math::Vec3c& p) const
		{
			return getLinearIndex(p.x, p.y, p.z);
		}

		/**
		Converts linear index to image coordinates.
		*/
		math::Vec3c getCoords(size_t linearIndex) const
		{
			return indexToCoords(linearIndex, dimensions());

			// NOTE: This is a dimensionality-independent version (for use if Vec3 is replaced by dimensionality-independent vector)
			//math::Vec3c tempSpace;
			//math::Vec3c result;
			//math::Vec3c dims = dimensions();

			//size_t currVal = 1;
			//for (size_t n = 0; n < dims.size(); n++)
			//{
			//	tempSpace[n] = currVal;
			//	currVal *= dims[n];
			//}


			//size_t mn1 = linearIndex;
			//for (size_t i = 0; i < dims.size(); i++)
			//{
			//	size_t n = dims.size() - 1 - i;

			//	size_t mn;
			//	if (mn1 > 0)
			//		mn = mn1 % tempSpace[n];
			//	else
			//		mn = 0;
			//	size_t xn = (mn1 - mn) / tempSpace[n];
			//	mn1 = mn;

			//	result[n] = (coord_t)xn;
			//}

			//return result;
		}

		/**
		Get a pixel at specified location.
		No bounds checking is performed.
		*/
		const pixel_t& operator()(coord_t x, coord_t y = 0, coord_t z = 0) const
		{
			return pDataConst[getLinearIndex(x, y, z)];
		}

		/**
		Get a reference to pixel at the specified location.
		No bounds checking is performed.
		*/
		pixel_t& operator()(coord_t x, coord_t y = 0, coord_t z = 0)
		{
			return pData[getLinearIndex(x, y, z)];
		}

		/**
		Get a pixel at specified location.
		No bounds checking is performed.
		*/
		const pixel_t& operator()(const math::Vec3c& p) const
		{
			return operator()(p.x, p.y, p.z);
		}

		/**
		Get a reference to pixel at the specified location.
		No bounds checking is performed.
		*/
		pixel_t& operator()(const math::Vec3c& p)
		{
			return operator()(p.x, p.y, p.z);
		}
		
		/**
		Get a pixel at specified location.
		No bounds checking is performed.
		*/
		const pixel_t& operator()(const math::Vec3sc& p) const
		{
			return operator()(p.x, p.y, p.z);
		}

		/**
		Get a reference to pixel at the specified location.
		No bounds checking is performed.
		*/
		pixel_t& operator()(const math::Vec3sc& p)
		{
			return operator()(p.x, p.y, p.z);
		}

		/**
		Calculates distance of pos to the nearest edge of the image.
		The distance is defined as
		min(min_i(|pos_i|), min_i(|dimensions_i-1-pos_i|))
		@param pos Position
		@param dimensions Dimensions of the image.
		*/
		inline coord_t edgeDistance(const math::Vec3c& pos) const
		{
			math::Vec3c dims = dimensions();

			coord_t dist = ::abs(pos[0]);

			for (size_t n = 0; n < pos.size(); n++)
			{
				coord_t cd = ::abs(pos[n]);
				if (cd < dist)
					dist = cd;

				cd = ::abs(dims[n] - 1 - pos[n]);
				if (cd < dist)
					dist = cd;
			}

			return dist;
		}

		inline coord_t edgeDistance(const math::Vec3sc& pos) const
		{
            return edgeDistance(math::Vec3c(pos));
		}
		
		/**
		Gets a value indicating whether the given position is on the edge of the image.
		*/
		inline bool isOnEdge(const math::Vec3c& pos) const
		{
			return edgeDistance(pos) <= 0;
		}

		/**
		Gets a value indicating whether the given position is on the edge of the image.
		*/
		inline bool isOnEdge(const math::Vec3sc& pos) const
		{
			return edgeDistance(pos) <= 0;
		}

		/**
		Test whether the given coordinates are inside this image.
		*/
		template<typename val_t> bool isInImage(val_t x, val_t y = 0, val_t z = 0) const
		{
			return x >= 0 && y >= 0 && z >= 0 && (coord_t)x < width() && (coord_t)y < height() && (coord_t)z < depth();
		}

		/**
		Test whether the given position is inside this image.
		*/
		template<typename val_t> bool isInImage(const math::Vec3<val_t>& pos) const
		{
			return isInImage(pos.x, pos.y, pos.z);
		}

		/**
		Test whether the size of this image is the same than size of the given image.
		*/
		template<typename pixel2_t> bool sizeEquals(const Image<pixel2_t>& r) const
		{
			return width() == r.width() && height() == r.height() && depth() == r.depth();
		}

		/**
		Tests whether the size of this image and the given image are equal and throws an
		exception if they are not.
		*/
		template<typename pixel2_t> void checkSize(const Image<pixel2_t>& r) const
		{
			if(!sizeEquals(r))
				throw ITLException("Image sizes are not equal.");
		}

		/**
		Makes sure that the size of this image equals given dimensions.
		Initializes image again if required.
		*/
		inline void ensureSize(const math::Vec3c& dims)
		{
			if (!ImageBase::sizeEquals(dims))
				init(dims);
		}

		/**
		Makes sure that the size of this image equals given dimensions.
		Initializes image again if required.
		*/
		inline void ensureSize(coord_t w, coord_t h = 1, coord_t d = 1)
		{
			ensureSize(math::Vec3c(w, h, d));
		}

		/**
		Makes sure that the size of this image equals the size of the given image.
		Initializes image again if required.
		*/
		template<typename pixel2_t> void ensureSize(const Image<pixel2_t>& r)
		{
			ensureSize(r.dimensions());
		}

		/**
		Throws exception if the given image is the same than this image.
		*/
		template<typename pixel2_t> void mustNotBe(const Image<pixel2_t>& other) const
		{
			if ((void*)this == (void*)&other)
				throw ITLException("This operation cannot be executed if the parameter images are the same.");
		}
	};

	namespace tests
	{
		void image();
		void buffers();
	}
}
