#pragma once

#include "image.h"
#include "boundarycondition.h"
#include "math/vec3.h"
#include "progress.h"

#include <vector>

namespace itl2
{
	namespace internals
	{
		/**
		Calculates minimum or maximum filtering (defined by function op) of array of values using van Herk algorithm.
		See van Herk - A fast algorithm for local minimum and maximum filters on rectangular and octagonal kernels.
		*/
		template<typename pixel_t, class Operation> void lineOp(std::vector<pixel_t>& img, coord_t r, BoundaryCondition bc, Operation op, pixel_t padValue, std::vector<pixel_t>& g, std::vector<pixel_t>& h)
		{
			// This function is implemented with img as separate array. It is very easy to change that to pointer + stride if that is
			// good for performance.

			while (img.size() % (2 * r + 1) != 0)
				img.push_back(padValue);

			while (g.size() < img.size())
				g.push_back(0);

			while (h.size() < img.size())
				h.push_back(0);

			coord_t count = (coord_t)img.size();

			// Build g and h arrays
			// TODO: This could be made by incrementing pointers
			for (coord_t x = 0; x < count; x += 2 * r + 1)
			{
				g[x] = img[x];
				coord_t xx = x + 1;
				coord_t xxmax = std::min(count, x + 2 * r + 1);
				while (xx < xxmax)
				{
					g[xx] = op(g[xx - 1], img[xx]);
					xx++;
				}

				xx--;
				h[xx] = img[xx];
				xx--;
				while (xx >= x)
				{
					h[xx] = op(h[xx + 1], img[xx]);
					xx--;
				}
			}

			// Store output
			if (bc == BoundaryCondition::Nearest)
			{
				// The code below is optimization of this loop:
				//for (coord_t x = 0; x < count; x++)
				//{
				//	pixel_t res;
				//	if (x + r < count && x - r >= 0)
				//		res = op(g[x + r], h[x - r]);
				//	else if (x + r < count)
				//		res = g[x + r];
				//	else if (x - r >= 0)
				//		res = op(h[x - r], img(count - 1, y, z));
				//	else
				//		res = op(g[g.size() - 1], h[0]);
				//	img(x, y, z) = res;
				//}

				// Left boundary
				for (coord_t x = 0; x < std::min<coord_t>(r, count); x++)
				{
					pixel_t res;
					if (x + r < count && x - r >= 0)
						res = op(g[x + r], h[x - r]);
					else if (x + r < count)
						res = g[x + r];
					else if (x - r >= 0)
						res = op(h[x - r], img[count - 1]);
					else
						res = op(g[count - 1], h[0]);

					img[x] = res;
				}

				// Center
				for (coord_t x = r; x < count - r; x++)
					img[x] = op(g[x + r], h[x - r]);

				// Right boundary
				for (coord_t x = std::max<coord_t>(0, count - r); x < count; x++)
				{
					pixel_t res;
					if (x + r < count && x - r >= 0)
						res = op(g[x + r], h[x - r]);
					else if (x + r < count)
						res = g[x + r];
					else if (x - r >= 0)
						res = op(h[x - r], img[count - 1]);
					else
						res = op(g[count - 1], h[0]);

					img[x] = res;
				}
			}
			else if (bc == BoundaryCondition::Zero)
			{
				// The code below is optimization of this loop:
				//for (coord_t x = 0; x < count; x++)
				//{
				//	pixel_t res;
				//	if (x + r < count && x - r >= 0)
				//		res = op(g[x + r], h[x - r]);
				//	else if (x + r < count)
				//		res = op(g[x + r], 0);
				//	else if (x - r >= 0)
				//		res = op(h[x - r], 0);
				//	else
				//		res = op(0, op(g[g.size() - 1], h[0]));
				//	img(x, y, z) = res;
				//}

				// Left boundary
				for (coord_t x = 0; x < std::min<coord_t>(r, count); x++)
				{
					pixel_t res;
					if (x + r < count && x - r >= 0)
						res = op(g[x + r], h[x - r]);
					else if (x + r < count)
						res = op(g[x + r], 0);
					else if (x - r >= 0)
						res = op(h[x - r], 0);
					else
						res = op(0, op(g[count - 1], h[0]));

					img[x] = res;
				}

				// Center
				for (coord_t x = r; x < count - r; x++)
					img[x] = op(g[x + r], h[x - r]);

				// Right boundary
				for (coord_t x = std::max<coord_t>(0, count - r); x < count; x++)
				{
					pixel_t res;
					if (x + r < count && x - r >= 0)
						res = op(g[x + r], h[x - r]);
					else if (x + r < count)
						res = op(g[x + r], 0);
					else if (x - r >= 0)
						res = op(h[x - r], 0);
					else
						res = op(0, op(g[count - 1], h[0]));

					img[x] = res;
				}
			}
			else
			{
				throw ITLException("Unsupported boundary condition.");
			}
		}

		/**
		Calculates minimum or maximum filtering (defined by function op) of image with periodic line structuring element, using van Herk algorithm.
		See van Herk - A fast algorithm for local minimum and maximum filters on rectangular and octagonal kernels
		and Jones - Periodic lines Definition, cascades, and application to granulometries
		@param r Half the length of the structuring element.
		@param step Vector that gives the direction of the periodic line. If the step is anything else than single pixel step (including diagonal steps), the periodic line may not be continuous (as expected, see Jones' paper).
		*/
		template<typename pixel_t, typename Operation> void lineOp(Image<pixel_t>& img, coord_t r, const Vec3c& step, pixel_t padValue, BoundaryCondition bc, Operation op)
		{
			// Find all start points and process each, making sure each of them is processed only once.
			std::vector<Vec3c> startPoints;
			for (coord_t z = 0; z < img.depth(); z++)
			{
				for (coord_t y = 0; y < img.height(); y++)
				{
					coord_t xmin, xmax;
					if ((step.z > 0 && z < step.z) ||
						(step.z < 0 && z >= img.depth() - abs(step.z)))
					{
						// Process all x and y
						xmin = 0;
						xmax = img.width();
					}
					else
					{
						if ((step.y > 0 && y < step.y) ||
							(step.y < 0 && y >= img.height() - abs(step.y)))
						{
							// Process all x
							xmin = 0;
							xmax = img.width();
						}
						else
						{
							if (step.x > 0)
							{
								// Process [0, step.x[
								xmin = 0;
								xmax = step.x;
							}
							else if (step.x < 0)
							{
								// Process [width - abs(step.x), width[
								xmin = img.width() - abs(step.x);
								xmax = img.width();
							}
							else
							{
								// Process nothing
								xmin = 0;
								xmax = 0;
							}
						}
					}


					for (coord_t x = xmin; x < xmax; x++)
					{
						startPoints.push_back(Vec3c(x, y, z));
					}
				}
			}

			ProgressIndicator progress(startPoints.size());
			#pragma omp parallel if(!omp_in_parallel() && img.pixelCount() > PARALLELIZATION_THRESHOLD)
			{
				std::vector<pixel_t> g(img.width() + img.height() + img.depth(), 0);
				std::vector<pixel_t> h(g.size(), 0);
				std::vector<pixel_t> row;

				#pragma omp for
				for (coord_t m = 0; m < (coord_t)startPoints.size(); m++)
				{
					Vec3c start = startPoints[m];

					// Gather pixels. This should be good for cache coherence (compared to start pointer + stride) in further
					// processing, but no measurements have been made in this particular case.
					row.clear();
					Vec3c p = start;
					while (img.isInImage(p))
					{
						row.push_back(img(p));
						//img(start)++;
						p += step;
					}

					// Process
					lineOp<pixel_t, Operation>(row, r, bc, op, padValue, g, h);

					// Scatter pixels back
					size_t n = 0;
					p = start;
					while (img.isInImage(p))
					{
						img(p) = row[n];
						p += step;
						n++;
					}

					progress.step();
				}
			}

		}
	}

	/**
	Calculates maximum filtering of image with periodic line structuring element.
	@param r Half the length of the structuring element.
	@param step Vector that gives the direction of the periodic line.
	*/
	template<typename pixel_t> void lineMax(Image<pixel_t>& img, coord_t r, const Vec3c& step, BoundaryCondition bc)
	{
		pixel_t padValue;
		if (bc == BoundaryCondition::Zero)
			padValue = 0;
		else
			padValue = std::numeric_limits<pixel_t>::lowest();
		internals::lineOp<pixel_t, const pixel_t&(const pixel_t&, const pixel_t&)>(img, r, step, padValue, bc, std::max<pixel_t>);
	}

	/**
	Calculates minimum filtering of image with periodic line structuring element.
	@param r Half the length of the structuring element.
	@param step Vector that gives the direction of the periodic line.
	*/
	template<typename pixel_t> void lineMin(Image<pixel_t>& img, coord_t r, const Vec3c& step, BoundaryCondition bc)
	{
		pixel_t padValue;
		if (bc == BoundaryCondition::Zero)
			padValue = 0;
		else
			padValue = std::numeric_limits<pixel_t>::max();
		internals::lineOp<pixel_t, const pixel_t&(const pixel_t&, const pixel_t&)>(img, r, step, padValue, bc, std::min<pixel_t>);
	}



	namespace internals
	{
		
		/**
		Represents spherical structuring element approximately decomposed into set of lines with different orientations.
		Stores direction vectors and corresponding line radii.
		TODO: This could be made in a clearer way without separate Ns array.
		*/
		struct DecomposedSphere
		{
			/**
			Direction vectors.
			*/
			::std::vector<Vec3c> dirs;

			/**
			Indices to dirs array giving end of each direction group forming a set of symmetrical lines (that all have the same length).
			*/
			::std::vector<coord_t> Ns;

			/**
			Radius value for each direction group.
			*/
			::std::vector<coord_t> rls;
		};

		/**
		Calculates approximate maximum or minimum filter with spherical structuring element.
		The structuring element is set of lines and it is defined by the given DecomposedSphere structure.
		*/
		template<typename pixel_t, typename Operation> void sphereOpApprox(Image<pixel_t>& img, const DecomposedSphere& directions, BoundaryCondition bc, Operation op, pixel_t padValue)
		{
			coord_t j = 0;
			ProgressIndicator progress(directions.dirs.size());
			for (size_t i = 0; i < directions.Ns.size(); i++)
			{
				coord_t rl = directions.rls[i];
				coord_t N = directions.Ns[i];

				while (j < N)
				{
					if (rl > 0)
						lineOp<pixel_t, Operation>(img, rl, directions.dirs[j], padValue, bc, op);
					j++;
					progress.step();
				}
			}
		}

		/**
		Calculates maximum filtering with spherical structuring element by decomposing the structuring element to a set of periodic lines.
		*/
		template<typename pixel_t> void maxFilterSphereApprox(Image<pixel_t>& img, const DecomposedSphere& directions, BoundaryCondition bc)
		{
			pixel_t padValue;
			if (bc == BoundaryCondition::Zero)
				padValue = 0;
			else
				padValue = std::numeric_limits<pixel_t>::lowest();

			sphereOpApprox<pixel_t, const pixel_t&(const pixel_t&, const pixel_t&)>(img, directions, bc, std::max<pixel_t>, padValue);
		}

		/**
		Calculates minimum filtering with spherical structuring element by decomposing the structuring element to a set of periodic lines.
		*/
		template<typename pixel_t> void minFilterSphereApprox(Image<pixel_t>& img, const DecomposedSphere& directions, BoundaryCondition bc)
		{
			pixel_t padValue;
			if (bc == BoundaryCondition::Zero)
				padValue = 0;
			else
				padValue = std::numeric_limits<pixel_t>::max();

			sphereOpApprox<pixel_t, const pixel_t&(const pixel_t&, const pixel_t&)>(img, directions, bc, std::min<pixel_t>, padValue);
		}

		/**
		Finds optimal decomposition of sphere of radius r from cache, or if the cache does not contain it, calculates the decomposition.
		*/
		DecomposedSphere optimizeStructuringElementCached(coord_t r);
	}


	/**
	Calculates approximation of maximum filtering by sphere.
	Decomposes spherical structuring element approximately to a set of periodic lines.
	It is not advisable to use this approximation for radiuses less than 5.
	The filtering will give wrong results where distance from image edge is less than r.
	Consider enlarging the image by r to all directions before processing.
	Enlarging in the z-direction is especially important for 2D images.
	*/
	template<typename pixel_t> void maxFilterSphereApprox(Image<pixel_t>& img, coord_t r, BoundaryCondition bc)
	{
		internals::DecomposedSphere elem = internals::optimizeStructuringElementCached(r);

		std::cout << "Using " << elem.dirs.size() << " linear elements to approximate a sphere." << std::endl;

		internals::maxFilterSphereApprox(img, elem, bc);
	}

	/**
	Calculates approximation of minimum filtering by sphere.
	Decomposes spherical structuring element approximately to a set of periodic lines.
	It is not advisable to use this approximation for radiuses less than 5.
	The filtering will give wrong results where distance from image edge is less than r.
	Consider enlarging the image by r to all directions before processing.
	Enlarging in the z-direction is especially important for 2D images.
	*/
	template<typename pixel_t> void minFilterSphereApprox(Image<pixel_t>& img, coord_t r, BoundaryCondition bc)
	{
		internals::DecomposedSphere elem = internals::optimizeStructuringElementCached(r);

		std::cout << "Using " << elem.dirs.size() << " linear elements to approximate a sphere." << std::endl;

		internals::minFilterSphereApprox(img, elem, bc);
	}


	namespace tests
	{
		void lineMax();
		void lineMin();
		void sphereMaxSpeed();
		void sphereMaxApprox();
		void sphereMinApprox();
	}
}
