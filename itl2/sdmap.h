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

		/**
		Merge two priority queues.
		*/
		template<typename T> void merge_pq(std::priority_queue<T>& dest, std::priority_queue<T>& src) {
			if (dest.size() < src.size()) {
				std::swap(dest, src);
			}
			while (!src.empty()) {
				dest.push(src.top());
				src.pop();
			}
		}
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
		#pragma omp parallel if(seeds.pixelCount() > PARALLELIZATION_THRESHOLD)
		{
			std::priority_queue<internals::DMapSeed<Tregion> > pointsLocal;
			
			#pragma omp for
			for (coord_t z = 0; z < seeds.depth(); z++)
			{
				for (coord_t y = 0; y < seeds.height(); y++)
				{
					for (coord_t x = 0; x < seeds.width(); x++)
					{
						Vec3sc p((int32_t)x, (int32_t)y, (int32_t)z);
						if (seeds(p) != 0)
						{
							// The color of the seed region
							Tregion region = geometry(p);

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
											pointsLocal.push(internals::DMapSeed<Tregion>(np, region, (p - np).norm<float32_t>()));
										}
									}

									if (p[n] < seeds.dimension(n) - 1)
									{
										Vec3sc np = p;
										np[n]++;

										Tregion lbl = geometry(np);
										if (lbl == region && seeds(np) == 0)
										{
											pointsLocal.push(internals::DMapSeed<Tregion>(np, region, (p - np).norm<float32_t>()));
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
									if (seeds.isInImage(np))
									{
										Tregion lbl = geometry(np);
										if (lbl == region && seeds(np) == 0)
										{
											pointsLocal.push(internals::DMapSeed<Tregion>(np, region, (p - np).norm<float32_t>()));
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
					}
				}
			}

			#pragma omp critical(sdmap_init)
			{
				internals::merge_pq(points, pointsLocal);
			}
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
