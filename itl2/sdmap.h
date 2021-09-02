#pragma once

#include "image.h"
#include "connectivity.h"
#include "iteration.h"

#include <vector>
#include <queue>
#include <iostream>

namespace itl2
{

	namespace internals
	{
		/**
		Storage class for priority queue in seeded dmap algorithm.
		*/
		template<typename Tlabel> class DMapSeed
		{
		private:
			Vec3sc pos;
			Tlabel targetLabel;
			float32_t dist;


		public:

			/**
			Constructor
			*/
			DMapSeed(const Vec3sc& p, Tlabel label, float32_t dist)
			{
				this->pos = p;
				this->dist = dist;
				this->targetLabel = label;
			}

			/**
			Gets position.
			*/
			const Vec3sc& position() const
			{
				return pos;
			}

			/**
			Compares weights.
			*/
			bool operator < (const DMapSeed& right) const
			{
				return dist > right.dist;
			}

			/**
			Gets label.
			*/
			const Tlabel region() const
			{
				return targetLabel;
			}

			const float32_t distance() const
			{
				return dist;
			}

		};
	}

	/**
	Calculates seeded distance map.
	@param seeds Seed image containing the set where the distance is zero. The set is marked with nonzero values, i.e. the distance map will propagate to pixels that have zero value in this image. This image is not modified.
	@param geometry Image containing the geometry. The distance transform will only proceed to pixels whose color in this image matches the color of the seed point (in this image). This image is not modified.
	@param distance Will contain distance to the nearest seed region.
	*/
	template<typename Tseed, typename Tregion> void seededDistanceMap(const Image<Tseed>& seeds, const Image<Tregion>& geometry, Image<float32_t>& distance, Connectivity connectivity = Connectivity::AllNeighbours)
	{
		seeds.checkSize(geometry);
		
		distance.ensureSize(seeds);

		std::priority_queue<internals::DMapSeed<Tregion> > points;

		Vec3sc currentRelativePosition;

		// Add neighbours of seed points to the priority queue
		{
			forAllPixels(seeds, [&](coord_t x, coord_t y, coord_t z)
			{
				Vec3sc p((int32_t)x, (int32_t)y, (int32_t)z);
				if (seeds(p) != 0) // if seed != 0
				{
					// The color of the seed region
					Tregion region = geometry(p);

					/*
					// This is ok for connectivity == Nearest
					for(size_t n = 0; n < p.size(); n++)
					{

						if(p[n] > 0)
						{
							vector<coord_t> np = p;
							np[n]--;

							Tregion lbl = geometry.getPixel(np);
							if(lbl == region && seeds.getPixel(np) == 0)
							{
								points.push(internal::DMapCompare<Tregion>(np, region, 1));
							}
						}

						if((size_t)p[n] < geometry.getDimension(n) - 1)
						{
							vector<coord_t> np = p;
							np[n]++;

							Tregion lbl = geometry.getPixel(np);
							if(lbl == region && seeds.getPixel(np) == 0)
							{
								points.push(internal::DMapCompare<Tregion>(np, region, 1));
							}
						}
					}
					*/

					if (connectivity == Connectivity::NearestNeighbours)
					{
						// Only nearest neighbours

						for (size_t n = 0; n < seeds.dimensionality(); n++)
						{
							if (p[n] > 0)
							{
								Vec3sc np = p;
								np[n]--;

								Tregion lbl = geometry(np);
								if (lbl == region && seeds(np) == 0)
								{
									points.push(internals::DMapSeed<Tregion>(np, region, (p - np).norm<float32_t>()));
								}
							}

							if (p[n] < seeds.dimension(n) - 1)
							{
								Vec3sc np = p;
								np[n]++;

								Tregion lbl = geometry(np);
								if (lbl == region && seeds(np) == 0)
								{
									points.push(internals::DMapSeed<Tregion>(np, region, (p - np).norm<float32_t>()));
								}
							}
						}
					}
					else
					{
						// All neighbours
						for (size_t n = 0; n < seeds.dimensionality(); n++)
							currentRelativePosition[n] = -1;

						do
						{
							Vec3sc np = p + currentRelativePosition;
							if(seeds.isInImage(np))
							{
								Tregion lbl = geometry(np);
								if (lbl == region && seeds(np) == 0)
								{
									points.push(internals::DMapSeed<Tregion>(np, region, (p - np).norm<float32_t>()));
								}
							}

							// Proceed in first dimension
							currentRelativePosition[0]++;

							// If first dimension is out of bounds, proceed one pixel in second dimension.
							// Cascade updates to upper dimensions, if required.
							for (size_t n = 0; n < currentRelativePosition.size() - 1; n++)
							{
								if (currentRelativePosition[n] > 1)
								{
									currentRelativePosition[n] = -1;
									currentRelativePosition[n + 1]++;
								}
								else
								{
									break;
								}
							}
						} while (currentRelativePosition[currentRelativePosition.size() - 1] <= 1);

					}


					distance(p) = 0;
				}
				else
				{
					distance(p) = std::numeric_limits<float32_t>::infinity();
				}
			});
		}

		/*
		// Plot start points
		while(!points.empty())
		{
			const internal::DMapCompare<Tregion>& obj = points.top();
			const vector<coord_t> p = obj.Position();
			Tregion region = obj.Region();
			points.pop();

			distance.setPixel(p, 1);
		}
		return;
		*/

		long filledCount = 0;
		size_t dispStep = 30000;
		size_t k = 0;
		//long round = 0;
		//long savecounter = 0;

		// Grow from the point p to all directions if they are not filled yet.
		while (!points.empty())
		{
			const internals::DMapSeed<Tregion>& obj = points.top();
			const Vec3sc p = obj.position();
			Tregion region = obj.region();
			float32_t currDistance = obj.distance();
			points.pop();

			if (currDistance < distance(p))
			{
				distance(p) = currDistance;

				/*
				// This is ok for connectivity == Nearest

				uint16_t newDistance = currDistance + 1;

				// Insert neighbours into the priority queue.
				for(size_t n = 0; n < p.size(); n++)
				{
					if(p[n] > 0)
					{
						vector<coord_t> np = p;
						np[n]--;

						Tregion lbl = geometry.getPixel(np);
						if(lbl == region && newDistance < distance.getPixel(np))
						{
							points.push(internal::DMapCompare<Tregion>(np, region, newDistance));
						}
					}

					if((size_t)p[n] < geometry.getDimension(n) - 1)
					{
						vector<coord_t> np = p;
						np[n]++;

						Tregion lbl = geometry.getPixel(np);
						if(lbl == region && newDistance < distance.getPixel(np))
						{
							points.push(internal::DMapCompare<Tregion>(np, region, newDistance));
						}
					}
				}
				*/

				if (connectivity == Connectivity::NearestNeighbours)
				{
					// Only nearest neighbours

					for (size_t n = 0; n < p.size(); n++)
					{
						if (p[n] > 0)
						{
							Vec3sc np = p;
							np[n]--;

							Tregion lbl = geometry(np);
							float32_t newDistance = currDistance + (p - np).norm<float32_t>();
							if (lbl == region && newDistance < distance(np))
							{
								points.push(internals::DMapSeed<Tregion>(np, region, newDistance));
							}
						}

						if (p[n] < seeds.dimension(n) - 1)
						{
							Vec3sc np = p;
							np[n]++;

							Tregion lbl = geometry(np);
							float32_t newDistance = currDistance + (p - np).norm();
							if (lbl == region && newDistance < distance(np))
							{
								points.push(internals::DMapSeed<Tregion>(np, region, newDistance));
							}
						}
					}
				}
				else
				{
					// All neighbours
					currentRelativePosition = Vec3sc(0, 0, 0);
					for (size_t n = 0; n < seeds.dimensionality(); n++)
						currentRelativePosition[n] = -1;
					

					do
					{
						Vec3sc np = p + currentRelativePosition;
						if (seeds.isInImage(np))
						{
							Tregion lbl = geometry(np);
							float32_t newDistance = currDistance + (p - np).norm();
							if (lbl == region && newDistance < distance(np))
							{
								points.push(internals::DMapSeed<Tregion>(np, region, newDistance));
							}
						}

						// Proceed in first dimension
						currentRelativePosition[0]++;

						// If first dimension is out of bounds, proceed one pixel in second dimension.
						// Cascade updates to upper dimensions, if required.
						for (size_t n = 0; n < currentRelativePosition.size() - 1; n++)
						{
							if (currentRelativePosition[n] > 1)
							{
								currentRelativePosition[n] = -1;
								currentRelativePosition[n + 1]++;
							}
							else
							{
								break;
							}
						}
					} while (currentRelativePosition[currentRelativePosition.size() - 1] <= 1);

				}
			}

			
			k++;
			if (k > dispStep)
			{
				k = 0;
				std::cout << "Queue size: " << points.size() << "                       \r" << std::flush;

				// This saves a movie of progress.
				//round++;
				//if(round > 10)
				//{
				//	round = 0;
				//	savecounter++;

				//	stringstream name;
				//	name << "movie/labels" << savecounter;
				//	raw::write(labels, raw::concatDimensions(name.str(), labels.dimensions()));

				//}
			}
		}


		//savecounter++;
		//stringstream name;
		//name << "movie/labels" << savecounter;
		//raw::write(labels, raw::concatDimensions(name.str(), labels.dimensions()));

	}

	namespace tests
	{
		void seededDMap();
	}
}
