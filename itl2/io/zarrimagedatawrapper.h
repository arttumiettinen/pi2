#pragma once

#include <string>
#include <sstream> // Ensure you include this for stringstream

#include "json.h"
#include "utilities.h"

namespace itl2
{
	using std::cout, std::endl;
	namespace zarr::internals{
		template<typename pixel_t>
		class ImageDataWrapper
		{
			Vec3c transposeOrder;

			Vec3c unTransposedCoords(Vec3c p) const
			{
				Vec3c tp(0, 0, 0);
				tp[transposeOrder.x] = p.x;
				tp[transposeOrder.y] = p.y;
				tp[transposeOrder.z] = p.z;
				return tp;
			}

			Vec3c transposedCoords(Vec3c p) const
			{
				Vec3c tp(0, 0, 0);
				tp.x = p[transposeOrder.x];
				tp.y = p[transposeOrder.y];
				tp.z = p[transposeOrder.z];
				return tp;
			}

		 public:
			Image<pixel_t>& img;

			// Initialize img using the initializer list
			ImageDataWrapper(Image<pixel_t>& img)
				: img(img), transposeOrder(0, 1, 2)
			{
			}

			void transpose(const Vec3c& order)
			{
				this->transposeOrder = order;
			}

			Vec3c dims() const
			{
				return transposedCoords(img.dimensions());
			}

			/**
			   Get a reference to a pixel at the specified location.
			   No bounds checking is performed.
		   */
			pixel_t& operator()(const Vec3c& p)
			{
				Vec3c tp = unTransposedCoords(p);
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
	}}