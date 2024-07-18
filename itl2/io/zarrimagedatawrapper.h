#pragma once

#include <string>
#include <sstream> // Ensure you include this for stringstream

#include "json.h"
#include "utilities.h"

namespace itl2
{
	using std::cout, std::endl;
	namespace zarr::internals
	{
		template<typename pixel_t>
		class ImageDataWrapper
		{
			Vec3c transposeOrder;
			const Vec3c datasetChunkShape;


		 public:
			Image<pixel_t>& img;

			//img is virtual and its shape is corresponding to the shape in pi2
			//physical coords are the transposed data written on disk by zarr
			//datasetShape and datasetChunkShape are equal to zarr.json file (virtual)
			ImageDataWrapper(Image<pixel_t>& img, const Vec3c& datasetChunkShape)
				: img(img), transposeOrder(0, 1, 2),  datasetChunkShape(datasetChunkShape)
			{
			}

			Vec3c virtualToPhysicalCoords(const Vec3c& p) const{
				return _transpose(p, transposeOrder);
			}

			Vec3c physicalToVirtualCoords(const Vec3c& p) const{
				return _transpose(p, _inverseOrder(transposeOrder));
			}

			void transpose(const Vec3c& order)
			{
				this->transposeOrder = order;
			}

			Vec3c dims() const
			{
				return virtualToPhysicalCoords(img.dimensions());
			}

			Vec3c physicalChunkShape() const
			{
				return virtualToPhysicalCoords(datasetChunkShape);
			}

			/**
			   Get a reference to a pixel at the specified location.
			   No bounds checking is performed.
		   */
			pixel_t& operator()(const Vec3c& p)
			{
				Vec3c tp = physicalToVirtualCoords(p);
				cout << "Accessing pixel at " << p << " transposed to " << tp << " value: " << img(tp) << endl;
				return img(tp); // Assuming img provides operator() to access pixel data
			}

			/**
			Get a reference to a pixel at the specified location.
			No bounds checking is performed.
			*/
			pixel_t& operator()(coord_t x, coord_t y, coord_t z)
			{
				return operator()(Vec3c(x, y, z));
			}
		};
	}
}