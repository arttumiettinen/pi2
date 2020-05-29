#pragma once

#include "image.h"
#include "floodfill.h"
#include "pointprocess.h"
#include "neighbourhood.h"
#include "math/vec2.h"
#include "math/vec3.h"
#include "utilities.h"
#include "network.h"
#include "aabox.h"
#include "math/matrix3x3.h"
#include "math/vectoroperations.h"
#include "interpolation.h"
#include "projections.h"
#include "indexforest.h"
#include "lineskeleton.h"
#include "misc.h"

#include <omp.h>
#include <vector>
#include <tuple>

namespace itl2
{

	namespace internals
	{
		/**
		Skeleton point types.
		*/
		constexpr int UNCLASSIFIED = 1;
		constexpr int ENDPOINT = 2;
		constexpr int CURVE = 3;
		constexpr int BRANCHING = 4;
		constexpr int JUNCTION = 5;
		constexpr int INTERNAL = 6;
		constexpr int EDGE = 7;
		constexpr int UNKNOWN = 8;
		

		/**
		Classifies all pixels in the image. Overwrites the image.
		*/
		template<typename pixel_t, bool test(Image<pixel_t>& nb), int value>
		void classify(Image<pixel_t>& img, bool showProgress)
		{
			size_t totalProcessed = 0;
			#pragma omp parallel if(!omp_in_parallel())
			{

				Image<pixel_t> nb(Vec3c(3, 3, 3));

				#pragma omp for
				for (coord_t z = 0; z < img.depth(); z++)
				{
					for (coord_t y = 0; y < img.height(); y++)
					{
						for (coord_t x = 0; x < img.width(); x++)
						{
							if (img(x, y, z) == UNCLASSIFIED) // Background pixels do not need to be classified.
							{
								getNeighbourhood(img, Vec3c(x, y, z), Vec3c(1, 1, 1), nb, BoundaryCondition::Zero);

								if(test(nb))
									img(x, y, z) = (pixel_t)value;
							}
						}
					}

					showThreadProgress(totalProcessed, img.depth(), showProgress);
				}
			}
		}

		/**
		Resets neighbourhood so that nonzero pixels become 1.
		*/
		template<typename pixel_t> void reset(Image<pixel_t>& nb)
		{
			threshold(nb, (pixel_t)0);
		}

		/**
		Counts object components in the given neighbourhood.
		*/
		template<typename pixel_t> coord_t countObjectComponents(Image<pixel_t>& nb)
		{
			coord_t count = 0;

			for (coord_t z = 0; z < nb.depth(); z++)
			{
				for (coord_t y = 0; y < nb.height(); y++)
				{
					for (coord_t x = 0; x < nb.width(); x++)
					{
						if (nb(x, y, z) != (pixel_t)0)
						{
							pixel_t fillColor = (pixel_t)(count + 2);
							floodfill(nb, Vec3c(x, y, z), fillColor, fillColor, Connectivity::AllNeighbours);
							count++;
						}
					}
				}
			}

			return count;
		}

		/**
		Endpoint voxel:
		An object voxel p with exactly one object neighbor is a curve end point voxel.
		*/
		template<typename pixel_t> bool testCurveEnd(Image<pixel_t>& nb)
		{
			nb(1, 1, 1) = 0;
			reset(nb);

			coord_t numberOfNeighbors = 0;
			for (coord_t n = 0; n < nb.pixelCount(); n++)
			{
				if (nb(n) != (pixel_t)0)
					numberOfNeighbors++;
			}

			if (numberOfNeighbors == 1)
				return true;

			return false;
		}

		/**
		Curve voxel:
		An object voxel p with at most two disjoint object neighbors is a curve voxel.
		*/
		template<typename pixel_t> bool testCurve(Image<pixel_t>& nb)
		{
			nb(1, 1, 1) = 0;
			reset(nb);

			coord_t numberOfNeighbors = 0;
			for (coord_t n = 0; n < nb.pixelCount(); n++)
			{
				if (nb(n) != (pixel_t)0)
					numberOfNeighbors++;
			}

			if (numberOfNeighbors == 1)
				return true;

			if (numberOfNeighbors == 2)
			{
				if (countObjectComponents(nb) == 2)
					return true;
			}

			return false;
		}

		/**
		Branching voxel:
		An object voxel p with a neighboring curve voxel is a branching voxel.
		*/
		template<typename pixel_t> bool testBranching(Image<pixel_t>& nb)
		{
			nb(1, 1, 1) = 0;
			
			for (coord_t n = 0; n < nb.pixelCount(); n++)
			{
				if (nb(n) == CURVE)
					return true;
			}
			

			return false;
		}


		/**
		Counts 6-connected background components in the given neighbourhood.
		*/
		template<typename pixel_t> coord_t count6BackgroundComponents(Image<pixel_t>& nb)
		{
			coord_t count = 0;

			for (coord_t z = 0; z < nb.depth(); z++)
			{
				for (coord_t y = 0; y < nb.height(); y++)
				{
					for (coord_t x = 0; x < nb.width(); x++)
					{
						if (nb(x, y, z) == 0)
						{
							pixel_t fillColor = (pixel_t)(count + 2);
							floodfill(nb, Vec3c(x, y, z), fillColor, fillColor, Connectivity::NearestNeighbours);
							count++;
						}
					}
				}
			}

			return count;
		}

		/**
		Junction voxel:
		An object voxel p with at least three 6-connected components of background voxels
		in N(p), or placed in a 2x2x2 configuration of object voxels is a junction voxel.
		*/
		template<typename pixel_t> bool testJunction(Image<pixel_t>& nb)
		{
			if (nb(0) == UNCLASSIFIED && nb(1) == UNCLASSIFIED && nb(3) == UNCLASSIFIED && nb(4) == UNCLASSIFIED && nb(9) == UNCLASSIFIED && nb(10) == UNCLASSIFIED && nb(12) == UNCLASSIFIED)
				return true;
			if (nb(1) == UNCLASSIFIED && nb(2) == UNCLASSIFIED && nb(4) == UNCLASSIFIED && nb(5) == UNCLASSIFIED && nb(10) == UNCLASSIFIED && nb(11) == UNCLASSIFIED && nb(14) == UNCLASSIFIED)
				return true;
			if (nb(9) == UNCLASSIFIED && nb(10) == UNCLASSIFIED && nb(12) == UNCLASSIFIED && nb(18) == UNCLASSIFIED && nb(19) == UNCLASSIFIED && nb(21) == UNCLASSIFIED && nb(22) == UNCLASSIFIED)
				return true;
			if (nb(10) == UNCLASSIFIED && nb(11) == UNCLASSIFIED && nb(14) == UNCLASSIFIED && nb(19) == UNCLASSIFIED && nb(20) == UNCLASSIFIED && nb(22) == UNCLASSIFIED && nb(23) == UNCLASSIFIED)
				return true;
			if (nb(3) == UNCLASSIFIED && nb(4) == UNCLASSIFIED && nb(6) == UNCLASSIFIED && nb(7) == UNCLASSIFIED && nb(12) == UNCLASSIFIED && nb(15) == UNCLASSIFIED && nb(16) == UNCLASSIFIED)
				return true;
			if (nb(4) == UNCLASSIFIED && nb(5) == UNCLASSIFIED && nb(7) == UNCLASSIFIED && nb(8) == UNCLASSIFIED && nb(14) == UNCLASSIFIED && nb(16) == UNCLASSIFIED && nb(17) == UNCLASSIFIED)
				return true;
			if (nb(12) == UNCLASSIFIED && nb(15) == UNCLASSIFIED && nb(16) == UNCLASSIFIED && nb(21) == UNCLASSIFIED && nb(22) == UNCLASSIFIED && nb(24) == UNCLASSIFIED && nb(25) == UNCLASSIFIED)
				return true;
			if (nb(14) == UNCLASSIFIED && nb(16) == UNCLASSIFIED && nb(17) == UNCLASSIFIED && nb(22) == UNCLASSIFIED && nb(23) == UNCLASSIFIED && nb(25) == UNCLASSIFIED && nb(26) == UNCLASSIFIED)
				return true;

			if (count6BackgroundComponents(nb) >= 3)
				return true;

			return false;
		}

		/**
		Helper function for c*_p calculation in testInternal.
		*/
		template<typename pixel_t> int cHelper(Image<pixel_t>& nb, coord_t x, coord_t y, coord_t z, int count)
		{
			x += 1;
			y += 1;
			z += 1;
			if (nb(x, y, z) == 0)
			{
				pixel_t fillColor = (pixel_t)(count + 2);
				floodfill(nb, Vec3c(x, y, z), fillColor, fillColor, Connectivity::NearestNeighbours);
				return count + 1;
			}

			return count;
		}

		/**
		Internal voxel:
		An object voxel p such that c*_p != 1 is an internal voxel.

		The number of 6-connected components of background voxels having p as
		face-neighbor and computed in N*(p) is denoted by by c*_p.

		The set N*(p) includes only the 6 face- and the 12 edge-neighbors of p.
		*/
		template<typename pixel_t> bool testInternal(Image<pixel_t>& nb)
		{
			nb(0) = UNCLASSIFIED;
			nb(2) = UNCLASSIFIED;
			nb(6) = UNCLASSIFIED;
			nb(8) = UNCLASSIFIED;
			nb(18) = UNCLASSIFIED; 
			nb(20) = UNCLASSIFIED; 
			nb(24) = UNCLASSIFIED; 
			nb(26) =  UNCLASSIFIED;

			int count = 0;
			count = cHelper(nb, -1, 0, 0, count);
			count = cHelper(nb, 1, 0, 0, count);
			count = cHelper(nb, 0, -1, 0, count);
			count = cHelper(nb, 0, 1, 0, count);
			count = cHelper(nb, 0, 0, -1, count);
			count = cHelper(nb, 0, 0, 1, count);

			return count != 1;
		}

		/**
		Used to test for unclassified voxel.

		Edge voxel:
		Any other object voxel p is an edge voxel.
		*/
		template<typename pixel_t> bool testAlwaysTrue(Image<pixel_t>& nb)
		{
			// No need to test anything!
			return true;
		}


		/**
		Set non-background edge voxels to given value.
		See also setEdges function.
		*/
		template<typename pixel_t> void handleEdges(Image<pixel_t>& img, pixel_t edgeValue = internals::BRANCHING)
		{
			for (size_t skip = 0; skip < img.dimensionality(); skip++)
			{
				Vec3c reducedsize = img.dimensions();
				reducedsize[skip] = 1;

				for (coord_t z = 0; z < reducedsize.z; z++)
				{
					for (coord_t y = 0; y < reducedsize.y; y++)
					{
						for (coord_t x = 0; x < reducedsize.x; x++)
						{
							Vec3c coords(x, y, z);

							coords[skip] = 0;
							
							if (img(coords) != 0)
								img(coords) = edgeValue;

							coords[skip] = img.dimension(skip) - 1;
							if (img(coords) != 0)
								img(coords) = edgeValue;
						}
					}
				}
			}
		}
	}

	/**
	Classifies line skeleton points to end points, normal points, intersection points, and edge points.
	Background is assumed to have value zero.
	See Arcelli - From 3D Discrete Surface Skeletons to Curve Skeletons.
	Call classifySkeleton(img, false, false) corresponds to Arcelli classification, except points at image edge are given
	UNKNOWN label.
	@param img Image containing the skeleton.
	@param curveEnds Set to true to classify curve endpoints separately from curve points.
	@param curveAndIntersectionOnly Set to true to skip classes JUNCTION, INTERNAL, and EDGE, and to set all those to BRANCHING.
	*/
	template<typename pixel_t> void classifySkeleton(Image<pixel_t>& img, bool curveEnds, bool curveAndIntersectionOnly, bool showProgress, pixel_t edgeValue = internals::BRANCHING)
	{
		threshold(img, (pixel_t)0);

		if(curveEnds)
			internals::classify<pixel_t, internals::testCurveEnd<pixel_t>, internals::ENDPOINT>(img, showProgress);

		internals::classify<pixel_t, internals::testCurve<pixel_t>, internals::CURVE>(img, showProgress);

		if (!curveAndIntersectionOnly)
		{
			// Proceed with classification normally
			internals::classify<pixel_t, internals::testBranching<pixel_t>, internals::BRANCHING>(img, showProgress);
			internals::classify<pixel_t, internals::testJunction<pixel_t>, internals::JUNCTION>(img, showProgress);
			internals::classify<pixel_t, internals::testInternal<pixel_t>, internals::INTERNAL>(img, showProgress);
			internals::classify<pixel_t, internals::testAlwaysTrue<pixel_t>, internals::EDGE>(img, showProgress);
		}
		else
		{
			// Set all unclassified points to BRANCHING
			internals::classify<pixel_t, internals::testAlwaysTrue<pixel_t>, internals::BRANCHING>(img, showProgress);
		}

		internals::handleEdges(img, edgeValue);
	}


	/**
	Smooths a line using anchored convolution with Gaussian. Does not change start and end point of the line.
	For the algorithm, see Suhadolnik - An anchored discrete convolution algorithm for measuring length in digital images.
	@param points The list of points that will be smoothed.
	@param sigma Standard deviation of Gaussian.
	@param delta The distance between original location of a point and new location won't be more than this value.
	*/
	void smoothLine(std::vector<Vec3f>& points, double sigma, double delta = 0.5);


	namespace internals
	{
		class BranchStartInfo
		{
		public:
			Vec3sc start;
			coord_t startVertexIndex = 0;
		};

		/**
		Checks if any of the points in the 'points' list is at the edge of an image whose dimensions are 'imageDimensions'.
		*/
		inline bool isOnEdge(const std::vector<Vec3sc>& points, const Vec3c& imageDimensions)
		{
			for (const Vec3sc& p : points)
			{
				if (itl2::isOnEdge(p, imageDimensions))
					return true;
			}

			return false;
		}

		template<typename pixel_t> std::vector<BranchStartInfo> addIntersection(Image<pixel_t>& img, const Vec3sc& pos, Network& network)
		{
			//if (pos == Vec3sc(161, 67, 0))
			//{
			//	std::cout << "Here" << std::endl;
			//}

			// Find all points in the intersection region
			std::vector<Vec3sc> filledPoints;
			floodfill<pixel_t>(img, Vec3c(pos), 0, 0, Connectivity::AllNeighbours, nullptr, &filledPoints);

			if (filledPoints.size() <= 0)
				throw ITLException("Invalid position passed to intersection tracer. There is no intersection at that location.");

			// Calculate average position of the intersection
			Vec3f center(0, 0, 0);
			for (const Vec3sc& p : filledPoints)
				center += Vec3f(p);
			center /= (float32_t)filledPoints.size();

			// Insert the new vertex to the network
			size_t vertexIndex = network.vertices.size();
			network.vertices.push_back(center);

			// If the vertex is on edge or any of the branch starts is on edge, insert point list to the network, too.
			bool onEdge = isOnEdge(filledPoints, img.dimensions());

			// Find all branches that start in point neighbouring the intersection region.
			// Add them to processing list.
			std::vector<BranchStartInfo> branchStarts;
			for (const Vec3sc& p : filledPoints)
			{
				for (int32_t z = -1; z <= 1; z++)
				{
					for (int32_t y = -1; y <= 1; y++)
					{
						for (int32_t x = -1; x <= 1; x++)
						{
							Vec3sc pp = p + Vec3sc(x, y, z);
							if (img.isInImage(pp) && img(pp) != 0)
							{
								// pp is a start of branch
								BranchStartInfo si;
								si.start = pp;
								si.startVertexIndex = vertexIndex;
								branchStarts.push_back(si);

								// If branch start is on edge, the intersection must be classified as incomplete
								// so that it participates in edge-ending branch resolving process.
								if (img.isOnEdge(pp))
									onEdge = true;
							}
						}
					}
				}
			}

			if (onEdge)
				network.incompleteVertices.push_back(IncompleteVertex(vertexIndex, filledPoints));

			return branchStarts;
		}


		/**
		Finds the next point on a curve starting at given start location.
		*/
		template<typename pixel_t> Vec3sc findNextOnCurve(Image<pixel_t>& img, Vec3sc start)
		{
			
			// Find if there is unique point where we could continue
			for (int32_t z = -1; z <= 1; z++)
			{
				for (int32_t y = -1; y <= 1; y++)
				{
					for (int32_t x = -1; x <= 1; x++)
					{
						Vec3sc pp = start + Vec3sc(x, y, z);
						if (img.isInImage(pp) && img(pp) == CURVE)
						{
							return pp;
						}
					}
				}
			}

			// No next point found
			return start;
		}

		/**
		Trace single branch of the skeleton and return the other end point of the branch, and erase it from the image.
		@param count The count of pixels traversed is stored in this variable on output.
		*/
		template<typename pixel_t> std::tuple<std::vector<Vec3sc>, Vec3sc> traceBranch(Image<pixel_t>& img, Vec3sc start, size_t& count)
		{
			std::vector<Vec3sc> points;
			
			count = 0;
			while (true)
			{
				points.push_back(start);
				img(start) = 0;
				count++;

				//if (start == Vec3sc(116, 130, 24))
				//{
				//	std::cout << "Here" << std::endl;
				//}

				Vec3sc next = findNextOnCurve(img, start);
				if(next == start)
					return make_tuple(points, next);
				start = next;
			}
		}

		/**
		Estimates length of line that goes through pixels whose locations are given.
		Uses anchored convolution approach from
		Suhadolnik - An anchored discrete convolution algorithm for measuring length in digital images
		@param points Positions of pixels through which the line goes. At output, contains smoothed (and anchored) points that can be used, e.g., to calculate tangent of the line using finite differences.
		@param pStraightLength Pointer to float where distance between smoothed and anchored end points of the line will be stored.
		@param sigma Amount of smoothing. Set to zero to calculate point-to-point polygon length.
		@return Length estimate.
		*/
		float32_t lineLength(std::vector<Vec3f>& points, /*float32_t* pStraightLength, */double sigma = 1, double maxDisplacement = 1);

		/**
		Extracts 2D slice from 3D image.
		@param img Image where the pixel data is extracted.
		@param pos Position of the center point of the slice.
		@param dir Direction of slice normal.
		@param slice Slice pixels will be set to this image. The size of the image must be the required size of the slice.
		*/
		template<typename pixel_t> void getSlice(const Image<pixel_t>& img, const Vec3d& pos, Vec3d dir, Image<pixel_t>& slice)
		{
			const Interpolator<pixel_t, pixel_t, double>& interpolator = LinearInterpolator<pixel_t, pixel_t, double, double>(BoundaryCondition::Zero);

			// Create orthogonal base by choosing two more directions
			dir = dir.normalized();
			Vec3d up(0, 0, 1);
			Vec3d right = dir.cross(up);
			if (right.norm() < 0.0001)
			{
				up = Vec3d(1, 0, 0);
				right = dir.cross(up);
			}
			right = right.normalized();
			up = right.cross(dir);

			// Build rotation matrix
			Matrix3x3d rot(up.x, right.x, dir.x,
								up.y, right.y, dir.y,
								up.z, right.z, dir.z);

			if (!NumberUtils<double>::equals(rot.det(), 1.0, 100 * NumberUtils<double>::tolerance()))
				throw ITLException("Invalid rotation matrix.");

			Vec2d sliceRadius = Vec2d((double)slice.width() - 1, (double)slice.height() - 1) / 2.0;

			for (coord_t yi = 0; yi < (coord_t)slice.height(); yi++)
			{
				for (coord_t xi = 0; xi < (coord_t)slice.width(); xi++)
				{
					Vec3d p(xi - (double)sliceRadius.x, (double)yi - (double)sliceRadius.y, 0);
					Vec3d x = rot * p + pos;
					Vec3c xc = round(x);
					
					if (img.isInImage(xc))
					{
						pixel_t pixel = interpolator(img, x);
						slice(xi, yi) = pixel;
					}
					else
					{
						slice(xi, yi) = 0;
					}
				}
			}
		}

		/**
		Estimates properties of skeleton edge from the points that belong to the edge.
		@param origShift Position of image origin in pOriginal image.
		@param areaReplacement Value that will be used to fill EdgeMeasurements::area if pOriginal is nullptr.
		*/
		template<typename orig_t> EdgeMeasurements measureEdge(const Image<orig_t>* pOriginal, const Vec3d& origShift, std::vector<Vec3f> edgePoints, const std::vector<Vec3sc>& edgePointsNoIntersections, float32_t areaReplacement, bool storeAllEdgePoints, double smoothingSigma, double maxDisplacement)
		{
			EdgeMeasurements result;

			result.pointCount = (float32_t)edgePoints.size();
			result.length = lineLength(edgePoints, smoothingSigma, maxDisplacement);

			// Default value for area used if area measurements fail or original image is not available.
			result.area = areaReplacement;

			if (edgePointsNoIntersections.size() > 0)
			{
				if (!storeAllEdgePoints)
					result.edgePoints.push_back(edgePointsNoIntersections[0]);
				else
					result.edgePoints.insert(result.edgePoints.end(), edgePointsNoIntersections.begin(), edgePointsNoIntersections.end());
			}

			// Calculate area by extracting a slice through each point on the path
			// If there are less than three points, we cannot measure area. (two points would be possible but a special case)
			if (edgePointsNoIntersections.size() >= 2 && pOriginal)
			{
				//bool write = false;
				//if (edgePointsNoIntersections[0] == Vec3sc(157, 189, 53))
				//{
				//	write = true;
				//	std::cout << "STOP" << std::endl;
				//}

				// Calculate smoothed version for determination of tangent vectors
				std::vector<Vec3f> smoothed;
				for (const Vec3sc& v : edgePointsNoIntersections)
					smoothed.push_back(Vec3f(v));
				lineLength(smoothed, smoothingSigma, maxDisplacement);

				if (smoothed.size() < 3)
				{
					// If we have only two points, duplicate the second one so that
					// when calculating tangent (we have n=1), the tangent becomes
					// smoothed[2] - smoothed[0] == second point - first point.
					smoothed.push_back(*smoothed.rbegin());
				}

				Image<orig_t> slice(75, 75, 0); // TODO: Hardcoded slice size
				std::vector<float32_t> areaSamples;
				areaSamples.reserve(smoothed.size());
				for (size_t n = 1; n < smoothed.size() - 1; n++)
				{
					Vec3d center = Vec3d(smoothed[n]);
					//Vec3d tangent = Vec3d(edgePointsNoIntersections[n + 1] - edgePointsNoIntersections[n - 1]);
					Vec3d tangent = Vec3d(smoothed[n + 1] - smoothed[n - 1]);

					try
					{
						getSlice(*pOriginal, center + origShift, tangent, slice);
					}
					catch (ITLException)
					{
						// The slice cannot be calculated as tangent vector is too close to zero.
						continue;
					}

					//if(write)
					//	raw::writed(slice, "./slice");

					orig_t M = max(slice);
					threshold(slice, M / 2);

					//raw::writed(slice, "./skeleton/tracing/slice_th");

					std::vector<Vec3sc> filledPoints;
					floodfill<orig_t>(slice, slice.dimensions() / 2, (orig_t)0, (orig_t)0, Connectivity::AllNeighbours, nullptr, &filledPoints);

					//raw::writed(slice, "./skeleton/tracing/slice_fill");

					if (filledPoints.size() > 0)
					{
						areaSamples.push_back((float32_t)filledPoints.size());
					}
					//else
					//{
					//	raw::writed(slice, "./slice_filled");
					//}
				}

				// NOTE: mode has been changed to mean as we may have sample lists like [7, 6, 8, 9, 10, 1, 1], and then mode will be 1.
				// TODO: Maybe outlier-discarding mean would be even better?
				//if(areaSamples.size() > 0)
				//	result.area = (float32_t)mode(areaSamples);
				result.area = mean(areaSamples);

				//if (abs(mode(areaSamples) - mean(areaSamples)) > 10)
				//{
				//	for (float32_t s : areaSamples)
				//		std::cout << s << " ";
				//	std::cout << std::endl;
				//	std::cout << "mean == " << mean(areaSamples) << std::endl;
				//	std::cout << "mode == " << mode(areaSamples) << std::endl;
				//}

				//if (std::isnan(result.area) || result.area > 1000)
				//{
				//	raw::writed(slice, "./slice_filled");
				//	std::cout << "Bad area" << std::endl;
				//}
			}

			return result;
		}

		inline void buildEdgePointList(const std::vector<Vec3sc>& edgePoints, coord_t startVertexIndex, coord_t endVertexIndex, Network& net, std::vector<Vec3f>& points)
		{
			points.clear();
			points.reserve(edgePoints.size() + 2);
			if (startVertexIndex >= 0)
				points.push_back(net.vertices[startVertexIndex]);

			//points.insert(points.end(), edgePoints.begin(), edgePoints.end());
			for(auto& e : edgePoints)
			    points.push_back(Vec3f(e));

			if (endVertexIndex >= 0)
				points.push_back(net.vertices[endVertexIndex]);
		}

		/**
		Add edge to network.
		@param edgePoints List of points on edge, not including intersection points.
		@param startVertexIndex, endVertexIndex Indices of the start and the end vertices.
		@param areaReplacement This value is used as edge area if pOriginal is nullptr.
		*/
		template<typename orig_t> void addEdge(const std::vector<Vec3sc>& edgePoints, coord_t startVertexIndex, coord_t endVertexIndex, Network& net, const Image<orig_t>* pOriginal, const Vec3d& origShift, bool storeAllEdgePoints, double smoothingSigma, double maxDisplacement, float32_t areaReplacement = std::numeric_limits<float32_t>::signaling_NaN())
		{
			//Vec3sc stopPoint(138, 1, 94);
			//if (edgePoints.front() == stopPoint || edgePoints.back() == stopPoint)
			////if(isNeighbour(edgePoints.front(), stopPoint) || isNeighbour(edgePoints.back(), stopPoint))
			//{
			//	std::cout << "stop point" << std::endl;
			//}

			//if (startVertexIndex == 1 && endVertexIndex == 3)
			//{
			//	std::cout << "stop point" << std::endl;
			//}

			std::vector<Vec3f> currEdge;
			buildEdgePointList(edgePoints, startVertexIndex, endVertexIndex, net, currEdge);

			EdgeMeasurements props = internals::measureEdge(pOriginal, origShift, currEdge, edgePoints, areaReplacement, storeAllEdgePoints, smoothingSigma, maxDisplacement);

			if (startVertexIndex >= 0 && endVertexIndex >= 0 && !net.isIncompleteVertex(startVertexIndex) && !net.isIncompleteVertex(endVertexIndex))
			{

				//if ((props.pointCount > 4 && std::isnan(props.area)) || props.area > 1000)
				//{
				//	std::cout << "Bad area" << std::endl;
				//}


				Edge e(startVertexIndex, endVertexIndex, props);
				net.edges.push_back(e);
			}
			else
			{
				net.incompleteEdges.push_back(IncompleteEdge(startVertexIndex, endVertexIndex, edgePoints, props));
			}
		}

		/**
		Trace branches in the given seeds list and continue recursively to all branches that are connected to them.
		Add everything to the network.

		Single pixel long branch connecting to single intersection is missed.
		@param img Classified skeleton image.
		@param orig Original (non-skeletonized) image that is used for branch area and shape measurements.
		*/
		template<typename pixel_t, typename orig_t> void trace(Image<pixel_t>& img, const Image<orig_t>* pOriginal, const Vec3d& origShift, bool storeAllEdgePoints, double smoothingSigma, double maxDisplacement, const std::vector<BranchStartInfo>& seeds, Network& net)
		{
			// Add branches to the tracing list
			std::deque<internals::BranchStartInfo> tracingList;
			for (const BranchStartInfo& f : seeds)
				tracingList.push_front(f);

			// Process branches until the tracing list is empty
			while (!tracingList.empty())
			{
				BranchStartInfo info = tracingList.front();
				tracingList.pop_front();

				//if (info.start == Vec3sc(87, 135, 25))
				//{
				//	std::cout << "Here" << std::endl;
				//}


				// Only trace if the branch has not been traced yet. This is done to remove duplicate points in branches.
				// Duplicates may occur (at least) if edge ends in the edge of the image
				// and contacts multiple intersection regions. This is possible at least for edges of length 1:
				// ------------------------ image edge -----------------
				// i - i - x - i - i
				//   /          \
				// i              i
				// Points marked with i are all intersection points, point marked with x has two neighbours so
				// it is not an intersection point, but it is a branch start point for both of the intersections.
				// This way it ends up in incomplete edge list two times if we allow trace to start if the start
				// pixel is already set to zero.
				if (img(info.start) != 0)
				{
					size_t count = 0;
					std::tuple<std::vector<Vec3sc>, Vec3sc> result = traceBranch(img, info.start, count);
					const std::vector<Vec3sc>& edgePoints = std::get<0>(result);
					Vec3sc end = std::get<1>(result);

					bool found = false;
					// Find known edges where the end point is a start point (to resolve loops)
					for (size_t n = 0; n < tracingList.size(); n++)
					{
						if (tracingList[n].start == end)
						{
							found = true;

							// Found start point corresponding to our end point
							addEdge(edgePoints, info.startVertexIndex, tracingList[n].startVertexIndex, net, pOriginal, origShift, storeAllEdgePoints, smoothingSigma, maxDisplacement);

							tracingList.erase(tracingList.begin() + n);
							n--;
						}
					}

					
					// Find untraced intersections around the end point (there might be many of them)
					for (int32_t z = -1; z <= 1; z++)
					{
						for (int32_t y = -1; y <= 1; y++)
						{
							for (int32_t x = -1; x <= 1; x++)
							{
								Vec3sc pp = end + Vec3sc(x, y, z);
								//if (img.isInImage(pp) && (img(pp) == BRANCHING || img(pp) == internals::ENDPOINT))
								if (img.isInImage(pp) && img(pp) != 0 && img(pp) != internals::CURVE)
								{
									found = true;

									std::vector<BranchStartInfo> newBranches = addIntersection(img, pp, net);

									// The current branch is already erased so it is not going to be in the newBranches list.
									for (BranchStartInfo& f : newBranches)
									{
										tracingList.push_front(f);
									}

									// Measure edge only if both end point vertices are complete.
									// Otherwise add the edge to incompleteEdges list (with start and end vertex indices) and process it later.
									size_t endVertexIndex = net.vertices.size() - 1;
									addEdge(edgePoints, info.startVertexIndex, endVertexIndex, net, pOriginal, origShift, storeAllEdgePoints, smoothingSigma, maxDisplacement);
								}
							}
						}
					}

					// If the end point is on the edge of the image, this is incomplete vertex.
					if (!found && img.isOnEdge(end))
					{
						found = true;
						addEdge(edgePoints, info.startVertexIndex, -1, net, pOriginal, origShift, storeAllEdgePoints, smoothingSigma, maxDisplacement);
					}

					if (!found)
					{
						// Discard single-pixel edges
						if (edgePoints.size() > 1)
							throw ITLException("Found an edge that does not end in an intersection point, is not an incomplete edge, and is not a loop.");
					}
				}
			}
		}

		/**
		Classifies skeleton image for skeleton tracing.
		*/
		//template<typename pixel_t> void classifyForTracing(Image<pixel_t>& img, pixel_t edgeValue = internals::BRANCHING)
		template<typename pixel_t> void classifyForTracing(Image<pixel_t>& img)
		{
			//classifySkeleton(img, true, false, false, edgeValue);
			// Set edges to zero as non-zero edges cause confusion between real and block edges in multithreading and distributed processing.
			setEdges(img, 0);
			classifySkeleton(img, true, false, false);

			// Convert non-CURVE labels to BRANCHING so that all intersection regions have the same value.
			// (We will flood fill them later)
			#pragma omp parallel for if(!omp_in_parallel() && img.pixelCount() > PARALLELIZATION_THRESHOLD)
			for (coord_t n = 0; n < img.pixelCount(); n++)
			{
				pixel_t& p = img(n);
				if (p != 0 && p != internals::CURVE)
					p = internals::BRANCHING;
			}
		}

		/**
		Trace line skeleton into a graph structure using single-threaded processing.
		All pixels in the image are set to zero.
		@param classified Classified skeleton image. The image will be empty at output.
		@param original Original non-skeletonized image used for shape measurements.
		@net The network is inserted to this object. Each skeleton branch will be an edge in the network, and each bifurcation or intersection point will be a vertex.
		*/
		template<typename pixel_t, typename orig_t> void traceLineSkeleton(Image<pixel_t>& classified, const Image<orig_t>* pOriginal, const Vec3d& origShift, bool storeAllEdgePoints, double smoothingSigma, double maxDisplacement, Network& net, size_t& counter, size_t progressMax)
		{
			if (pOriginal)
			{
				classified.mustNotBe(*pOriginal);
				//classified.checkSize(*pOriginal);
			}

			// Find intersection point or fibre end point
			for (coord_t z = 0; z < classified.depth(); z++)
			{
				for (coord_t y = 0; y < classified.height(); y++)
				{
					for (coord_t x = 0; x < classified.width(); x++)
					{
						pixel_t label = classified(x, y, z);
						if(label == internals::BRANCHING)
						{
							// Find fibres that start from this intersection area or end point
							std::vector<internals::BranchStartInfo> branches = internals::addIntersection(classified, Vec3sc((int32_t)x, (int32_t)y, (int32_t)z), net);

							// Follow all branches and process the found intersection areas
							internals::trace(classified, pOriginal, origShift, storeAllEdgePoints, smoothingSigma, maxDisplacement, branches, net);
						}
					}
				}

				showThreadProgress(counter, progressMax);
			}


			// Trace incomplete edges that have no vertices in the image.
			// These may take place in the edges of the image.
			// NOTE: here we will miss loops without intersection points.
			for (coord_t z = 0; z < classified.depth(); z++)
			{
				for (coord_t y = 0; y < classified.height(); y++)
				{
					for (coord_t x = 0; x < classified.width(); x++)
					{
						if (classified.isOnEdge(Vec3c(x, y, z)))
						{
							pixel_t label = classified(x, y, z);
							if (label == internals::CURVE)
							{
								Vec3c p(x, y, z);
								Image<pixel_t> nb(3, 3, 3);
								getNeighbourhood(classified, p, Vec3c(1, 1, 1), nb, BoundaryCondition::Zero);

								coord_t fgc = foregroundCount(nb);
								
								// 1 foreground pixel  => single point, should be traced
								// 2 foreground pixels => end of branch, should be traced
								// 3 foreground pixels => middle of curve, should not be traced
								// 4 or more fg pixels => intersection, should be impossible in this phase of processing
								if(fgc == 1 || fgc == 2)
								{

									internals::BranchStartInfo info;
									info.start = Vec3sc(p);
									info.startVertexIndex = -1;

									internals::trace(classified, pOriginal, origShift, storeAllEdgePoints, smoothingSigma, maxDisplacement, { info }, net);
								}
								else if (fgc >= 4)
								{
									throw ITLException("Not all intersection were traced.");
								}
							}
						}
					}
				}
			}

			// Now trace loops without intersection points.
			// We generate intersection point for each loop in location where it is first encountered,
			// i.e. in the curve point that has minimal z, y, and x coordinates (in that order)
			for (coord_t z = 0; z < classified.depth(); z++)
			{
				for (coord_t y = 0; y < classified.height(); y++)
				{
					for (coord_t x = 0; x < classified.width(); x++)
					{
						pixel_t label = classified(x, y, z);
						if (label == internals::CURVE)
						{
							classified(x, y, z) = internals::BRANCHING;

							// Find fibres that start from this intersection area or end point
							std::vector<internals::BranchStartInfo> branches = internals::addIntersection(classified, Vec3sc((int32_t)x, (int32_t)y, (int32_t)z), net);

							// Follow all branches and process the found intersection areas
							internals::trace(classified, pOriginal, origShift, storeAllEdgePoints, smoothingSigma, maxDisplacement, branches, net);
						}
					}
				}
			}
		}

		inline std::tuple<coord_t, coord_t> calcMinMaxZ(coord_t blockCount, coord_t imageDepth, coord_t idx, bool allowSmallBlocks)
		{
			// Calculate amount of slices single thread should process.
			// The value is capped so that we don't divide the work unnecessarily too much, if the image is small.
			coord_t size = imageDepth / blockCount;

			coord_t minSize = 100;
			if (allowSmallBlocks)
				minSize = 5;

			if (size < minSize)
				size = minSize;

			coord_t minZ = idx * size;

			coord_t maxZ = minZ + size - 1;
			if (maxZ >= imageDepth)
			{
				maxZ = imageDepth - 1;
			}
			else
			{
				if (idx == blockCount - 1)
				{
					// The last thread processes possible "rounding error" slices
					maxZ = imageDepth - 1;
				}
				else
				{
					// If maxZ is near to image end, then adjust size of this block so that this thread processes the remaining slices.
					// This ensures that the last block is not too small.

					if (imageDepth - maxZ < minSize)
						maxZ = imageDepth - 1;

				}
			}

			return std::make_tuple(minZ, maxZ);
		}

		/**
		Traces line skeleton in blocks, and adds the network resulting from each block trace to subNets vector.
		@param origin If processing a block of larger image, set this to the origin of the block in full image coordinates. This value is added to all the vertex coordinates.
		*/
		template<typename pixel_t, typename orig_t> void traceLineSkeletonBlocks(Image<pixel_t>& classified, const Image<orig_t>* pOriginal, bool storeAllEdgePoints, double smoothingSigma, double maxDisplacement, std::vector<Network>& subNets, const Vec3sc& origin = Vec3sc(), coord_t blockCount = 0)
		{
			if (pOriginal)
			{
				classified.mustNotBe(*pOriginal);
				classified.checkSize(*pOriginal);
			}

			size_t counter = 0;
			std::cout << "Tracing skeleton..." << std::endl;

			// Allow small blocks only if the user has specifically requested count of blocks.
			// Otherwise block size is capped to larger value.
			bool allowSmallBlocks = blockCount > 0;
			if (blockCount <= 0)
				blockCount = omp_get_max_threads();

			// Determine real block count (condsidering block size limitations etc.)
			for (coord_t idx = 0; idx < blockCount; idx++)
			{
				coord_t minZ, maxZ;
				std::tie(minZ, maxZ) = calcMinMaxZ(blockCount, classified.depth(), idx, allowSmallBlocks);
				std::cout << "Thread " << idx << " processes range " << Vec2c(minZ, maxZ) << std::endl;
				if (maxZ >= classified.depth() - 1)
				{
					blockCount = idx + 1;
					break;
				}
			}

			//std::cout << "Tracing in " << blockCount << " blocks." << std::endl;

			for (int idx = 0; idx < blockCount; idx++)
				subNets.push_back(Network());

			#pragma omp parallel for if(blockCount > 1)
			for(int idx = 0; idx < blockCount; idx++)
			{
				coord_t minZ, maxZ;
				std::tie(minZ, maxZ) = calcMinMaxZ(blockCount, classified.depth(), idx, allowSmallBlocks);

				// If there's nothing to do for all the threads, the excess threads will just skip processing.
				if (minZ < classified.depth())
				{

					// Get view of part of the image
					Image<pixel_t> block(classified, minZ, maxZ);

					// Trace the region
					Network subNet;
					internals::traceLineSkeleton(block, pOriginal, Vec3d(0, 0, (double)minZ), storeAllEdgePoints, smoothingSigma, maxDisplacement, subNet, counter, classified.depth());

					// Convert vertices to global coordinates (they are in block coordinates)
					for (size_t n = 0; n < subNet.vertices.size(); n++)
					{
						subNet.vertices[n] += Vec3f(origin) + Vec3f(0, 0, (float)minZ);
					}
					for (size_t n = 0; n < subNet.edges.size(); n++)
					{
						//subNet.edges[n].properties.pointOnEdge += origin + Vec3sc(0, 0, (int32_t)minZ);
						subNet.edges[n].properties.edgePoints += (origin + Vec3sc(0, 0, (int32_t)minZ));
					}

					for (size_t n = 0; n < subNet.incompleteVertices.size(); n++)
					{
						for (size_t m = 0; m < subNet.incompleteVertices[n].points.size(); m++)
							subNet.incompleteVertices[n].points[m] += origin + Vec3sc(0, 0, (int32_t)minZ);
					}
					for (size_t n = 0; n < subNet.incompleteEdges.size(); n++)
					{
						for (size_t m = 0; m < subNet.incompleteEdges[n].points.size(); m++)
						{
							subNet.incompleteEdges[n].points[m] += origin + Vec3sc(0, 0, (int32_t)minZ);
							//subNet.incompleteEdges[n].properties.pointOnEdge += origin + Vec3sc(0, 0, (int32_t)minZ);
							subNet.incompleteEdges[n].properties.edgePoints += (origin + Vec3sc(0, 0, (int32_t)minZ));
						}
					}

					subNets[idx] = subNet;
				}
			}

			//// Erase empty networks for clarity
			//for (int idx = 0; idx < blockCount; idx++)
			//{
			//	if (subNets[idx].vertices.size() <= 0 && subNets[idx].incompleteVertices.size() <= 0)
			//	{
			//		subNets.erase(subNets.begin() + n);
			//		n--;
			//	}
			//}
		}


		void insertPoint(const Vec3sc& p, size_t n, Image<std::vector<size_t> >& grid, const Vec3sc& m, const Vec3sc& M);

		inline bool isNeighbour(const Vec3sc& a, const Vec3sc& b)
		{
			return (a - b).abs().max() <= 1;
		}

		coord_t findVertex(const std::vector<IncompleteVertex>& incompleteVertices, const Vec3sc& p, const coord_t invalidIndex);

		/**
		Combines incomplete edges found in the network.
		@param completedVertices List of incomplete vertices that were completed before this call.
		*/
		template<typename orig_t> void combineIncompleteEdges(Network& net, const std::vector<IncompleteVertex>& completedVertices, const Image<orig_t>* pOriginal, bool isFinal, bool storeAllEdgePoints, double smoothingSigma, double maxDisplacement)
		{
			std::vector<IncompleteEdge> incompleteEdges;
			swap(incompleteEdges, net.incompleteEdges);

			// Process all edges where both end points are known.
			{
				std::vector<IncompleteEdge> incompleteEdgesTemp;
				size_t counter = 0;
				for (size_t n = 0; n < incompleteEdges.size(); n++)
				{
					IncompleteEdge& e = incompleteEdges[n];
					if (e.verts[0] >= 0 && e.verts[1] >= 0)
					{
						// This is an edge that starts and/or ends in an vertex that was incomplete.
						// It needs only measuring.

						addEdge<orig_t>(e.points, e.verts[0], e.verts[1], net, pOriginal, Vec3d(), storeAllEdgePoints, smoothingSigma, maxDisplacement, e.properties.area);

						//incompleteEdges.erase(incompleteEdges.begin() + n);
						//n--;
					}
					else
					{
						// End points of this edge are not known so it requires further processing.
						incompleteEdgesTemp.push_back(e);
					}

					showThreadProgress(counter, incompleteEdges.size());
				}
				incompleteEdges = incompleteEdgesTemp;
			}

			// Find neighbouring incomplete edges
			// The forest contains indices into net.incompleteEdges list.
			IndexForest forest;
			forest.initialize(incompleteEdges.size());

			std::cout << "Determine bounding box..." << std::endl;
			Vec3sc m = Vec3sc(std::numeric_limits<int32_t>::max(), std::numeric_limits<int32_t>::max(), std::numeric_limits<int32_t>::max());
			Vec3sc M = Vec3sc(std::numeric_limits<int32_t>::lowest(), std::numeric_limits<int32_t>::lowest(), std::numeric_limits<int32_t>::lowest());
			for (size_t n = 0; n < incompleteEdges.size(); n++)
			{
				IncompleteEdge& e = incompleteEdges[n];

				if (e.verts[0] < 0) // verts[pi] < 0 means that we don't know where the edge is connected in this end.
				{
					AABox<int32_t> b(e.points.front(), e.points.front() + Vec3sc(1, 1, 1));
					b.inflate(1);
					m = min(m, b.minc);
					M = max(M, b.maxc);
				}

				if (e.verts[1] < 0) // verts[pi] < 0 means that we don't know where the edge is connected in this end.
				{
					AABox<int32_t> b(e.points.back(), e.points.back() + Vec3sc(1, 1, 1));
					b.inflate(1);
					m = min(m, b.minc);
					M = max(M, b.maxc);
				}
			}

			Image<std::vector<size_t> > grid(100, 100, 100);

			std::cout << "Divide end points to a grid..." << std::endl;
			size_t counter = 0;
			for (size_t n = 0; n < incompleteEdges.size(); n++)
			{
				IncompleteEdge& e = incompleteEdges[n];

				if (e.verts[0] < 0)
					insertPoint(e.points.front(), n, grid, m, M);

				if (e.verts[1] < 0)
					insertPoint(e.points.back(), n, grid, m, M);

				showThreadProgress(counter, incompleteEdges.size());
			}


			std::cout << "Determine overlaps..." << std::endl;
			counter = 0;
			#pragma omp parallel for if(!omp_in_parallel())
			for (coord_t i = 0; i < grid.pixelCount(); i++)
			{
				auto& list = grid(i);

				for (coord_t ni = 0; ni < (coord_t)list.size(); ni++)
				{
					size_t n = list[ni];
					const IncompleteEdge& v1 = incompleteEdges[n];

					for (size_t mi = ni + 1; mi < list.size(); mi++)
					{
						size_t m = list[mi];
						const IncompleteEdge& v2 = incompleteEdges[m];

						// Only test if the vertices have not been connected yet.
						// There is a race condition between union_sets below and find_sets on this line, but
						// that may only result in the test below evaluating true even though it should be
						// false, and in that case we just do some extra work.
						if (forest.find_set(n) != forest.find_set(m))
						{
							// Test if the incomplete edges are neighbours in some way
							if (isNeighbour(v2.points.front(), v1.points.back()) ||
								isNeighbour(v2.points.front(), v1.points.front()) ||
								isNeighbour(v2.points.back(), v1.points.front()) ||
								isNeighbour(v2.points.back(), v1.points.back()))
							{
								// v1 and v2 are neighbours, so the edges must be combined
								#pragma omp critical(forestInsert)
								{
									forest.union_sets(n, m);
								}
							}
						}
					}
				}

				showThreadProgress(counter, grid.pixelCount());
			}


			// Determine groups of edges that are to be combined
			std::cout << "Process and measure combined edges..." << std::endl;

			std::vector<std::vector<size_t> > combinedEdges(incompleteEdges.size());
			for (coord_t n = 0; n < (coord_t)incompleteEdges.size(); n++)
			{
				size_t base = forest.find_set(n);
				combinedEdges[base].push_back(n);
			}


			//Vec3sc stopPoint(179, 78, 42);

			// Combine edge point groups
			// TODO: This loop is slow-ish compared to other processes that happen in this function. Are findVertex function calls the reason?
			std::vector<IncompleteEdge> newIncompleteEdges;
			counter = 0;
			#pragma omp parallel for if(!omp_in_parallel())
			for (coord_t n = 0; n < (coord_t)combinedEdges.size(); n++)
			{
				std::vector<size_t> indices = combinedEdges[n];
				if (indices.size() > 0)
				{


					//for (size_t i : indices)
					//{
					//	IncompleteEdge& e = incompleteEdges[i];
					//	if (e.points.front() == stopPoint ||
					//		e.points.back() == stopPoint)
					//	{
					//		std::cout << "Stop point encountered" << std::endl;
					//	}
					//}


					// Measure average area
					double w = 0;
					double sum = 0;
					for (size_t i : indices)
					{
						EdgeMeasurements& m = incompleteEdges[i].properties;
						if (!std::isnan(m.area)) // area may be nan if the length of the branch is too short for area to be measured.
						{
							sum += (double)m.area * (double)m.length;
							w += m.length;
						}
					}
					float32_t areaApproximation = (float32_t)(sum / w);

					// Build point list and measure length
					Vec2c finalVerts = incompleteEdges[indices[0]].verts;
					std::vector<Vec3sc> finalPoints = incompleteEdges[indices[0]].points;

					indices.erase(indices.begin());
					while (indices.size() > 0)
					{
						// Find the next edge to be combined with the current one.
						bool found = false;
						for (size_t m = 0; m < indices.size(); m++)
						{
							IncompleteEdge& candidate = incompleteEdges[indices[m]];
							if (candidate.verts[0] < 0 && finalVerts[1] < 0 &&
								isNeighbour(finalPoints.back(), candidate.points.front()))
							{
								// Combine like this: finalPoints -- candidate

								finalPoints.insert(finalPoints.end(), candidate.points.begin(), candidate.points.end());

								//if (candidate.verts[0] >= 0 || finalVerts[1] >= 0)
								//	throw logic_error("Invalid graph (back-front).");

								finalVerts[1] = candidate.verts[1];
								found = true;
							}
							else if (candidate.verts[1] < 0 && finalVerts[1] < 0 &&
								isNeighbour(finalPoints.back(), candidate.points.back()))
							{
								// Combine like this: finalPoints -- candidate reversed
								std::reverse(candidate.points.begin(), candidate.points.end());
								finalPoints.insert(finalPoints.end(), candidate.points.begin(), candidate.points.end());

								//if (candidate.verts[1] >= 0 || finalVerts[1] >= 0)
								//	throw logic_error("Invalid graph (back-back).");

								finalVerts[1] = candidate.verts[0];
								found = true;
							}
							else if (candidate.verts[1] < 0 && finalVerts[0] < 0 &&
								isNeighbour(finalPoints.front(), candidate.points.back()))
							{
								// Combine like this: candidate -- finalPoints

								finalPoints.insert(finalPoints.begin(), candidate.points.begin(), candidate.points.end());

								//if (candidate.verts[1] >= 0 || finalVerts[0] >= 0)
								//	throw logic_error("Invalid graph (front-back).");

								finalVerts[0] = candidate.verts[0];
								found = true;
							}
							else if (candidate.verts[0] < 0 && finalVerts[0] < 0 &&
								isNeighbour(finalPoints.front(), candidate.points.front()))
							{
								// Combine like this: candidate reversed -- finalPoints
								std::reverse(candidate.points.begin(), candidate.points.end());
								finalPoints.insert(finalPoints.begin(), candidate.points.begin(), candidate.points.end());

								//if (candidate.verts[0] >= 0 || finalVerts[0] >= 0)
								//	throw logic_error("Invalid graph (front-front).");

								finalVerts[0] = candidate.verts[1];
								found = true;
							}

							if (found)
							{
								indices.erase(indices.begin() + m);
								break;
							}
						}

						if (!found)
							throw std::logic_error("Continuation not found.");
					}

					// If some end of the edge is still unknown, try connecting it to one of the incomplete vertices.
					// Try also completed vertices if the incomplete ones do not match.
					// NOTE: This does not account for the possibility that there are multiple incomplete vertices that should
					// be connected to this edge. Is that situation even possible?
					//if (finalVerts[0] < 0)
					//	finalVerts[0] = findVertex(net.incompleteVertices, finalPoints.front(), -1);
					//if (finalVerts[0] < 0)
					//	finalVerts[0] = findVertex(completedVertices, finalPoints.front(), -1);

					//if (finalVerts[1] < 0)
					//	finalVerts[1] = findVertex(net.incompleteVertices, finalPoints.back(), -1);
					//if (finalVerts[1] < 0)
					//	finalVerts[1] = findVertex(completedVertices, finalPoints.back(), -1);

					if (finalVerts[0] < 0)
					{
						finalVerts[0] = findVertex(completedVertices, finalPoints.front(), -1);
						if (finalVerts[0] < 0)
						{
							#pragma omp critical(netIncompleteVerticesAccess)
							{
								finalVerts[0] = findVertex(net.incompleteVertices, finalPoints.front(), -1);
							}
						}
					}

					if (finalVerts[1] < 0)
					{
						// Try non-loop first
						finalVerts[1] = findVertex(completedVertices, finalPoints.back(), finalVerts[0]);
						if (finalVerts[1] < 0)
						{
							#pragma omp critical(netIncompleteVerticesAccess)
							{
								finalVerts[1] = findVertex(net.incompleteVertices, finalPoints.back(), finalVerts[0]);
							}
						}

						if (isFinal) // If this is not final combine call, networks to be added might contain node that corresponds to this end without causing a loop.
						{
							// If non-loop possibility is not found, try to find loop. (this should always succeed!)
							if (finalVerts[1] < 0)
							{
								finalVerts[1] = findVertex(completedVertices, finalPoints.back(), -1);
								if (finalVerts[1] < 0)
								{
									#pragma omp critical(netIncompleteVerticesAccess)
									{
										finalVerts[1] = findVertex(net.incompleteVertices, finalPoints.back(), -1);
									}
								}
							}
						}
						else
						{
							// If there is a loop possibility, we must store the possible vertex in the incomplete list again so that we can check
							// the status in the final round.
							if (finalVerts[1] < 0)
							{
								coord_t possibility = findVertex(completedVertices, finalPoints.back(), -1);
								if (possibility >= 0)
								{
									auto obj = std::find_if(completedVertices.begin(), completedVertices.end(), [&](const IncompleteVertex& v) { return v.vertexIndex == possibility; });
									#pragma omp critical(netIncompleteVerticesAccess)
									{
										if (!net.isIncompleteVertex(obj->vertexIndex)) // We may have moved the vertex to the 'alive incomplete' list already.
										{
											net.incompleteVertices.push_back(*obj);
										}
									}
								}
							}
						}
					}





					//if (finalPoints.front() == stopPoint ||
					//	finalPoints.back() == stopPoint)
					//{
					//	std::cout << "stop" << std::endl;
					//}


					#pragma omp critical(addEdge)
					{
						addEdge<orig_t>(finalPoints, finalVerts[0], finalVerts[1], net, pOriginal, Vec3d(), storeAllEdgePoints, smoothingSigma, maxDisplacement, areaApproximation);
					}

#if defined(DEBUG)
					// Sanity check
					std::sort(finalPoints.begin(), finalPoints.end(), vecComparer<int32_t>);
					if (std::unique(finalPoints.begin(), finalPoints.end()) != finalPoints.end())
					{
						std::cout << "final vert 0: " << finalVerts[0] << std::endl;
						std::cout << "final vert 1: " << finalVerts[1] << std::endl;
						for (const auto& p : finalPoints)
							std::cout << p << std::endl;

						throw logic_error("edge contains non-unique points.");
					}
#endif
				}
				showThreadProgress(counter, combinedEdges.size());
			}


			if (isFinal)
			{
				// Here we might have loop edges in the incomplete edges list. Those must be resolved similarly
				// to single-block version.
				std::vector<IncompleteEdge> incompleteEdgesTemp;
				std::cout << "Determine edges that remain incomplete..." << std::endl;
				for (size_t n = 0; n < net.incompleteEdges.size(); n++)
				{
					IncompleteEdge& e = net.incompleteEdges[n];
					if (isNeighbour(e.points.front(), e.points.back()) &&
						e.verts[0] < 0 && e.verts[1] < 0)
					{
						// This is a loop without intersection points.
						// Generate intersection at minimal z, y, x coordinates and
						// generate edge beginning and ending there.
						auto it = std::min_element(e.points.begin(), e.points.end(), vecComparer<int32_t>);
						net.vertices.push_back(Vec3f(*it));
						coord_t vertexIndex = (coord_t)net.vertices.size() - 1;

						// Sort the edge point so that the first point is *it
						std::rotate(e.points.begin(), it, e.points.end());

						// The point corresponding to the vertex must not be in the points list.
						e.points.erase(e.points.begin());

						addEdge<orig_t>(e.points, vertexIndex, vertexIndex, net, pOriginal, Vec3d(), storeAllEdgePoints, smoothingSigma, maxDisplacement, e.properties.area);

						//net.incompleteEdges.erase(net.incompleteEdges.begin() + n);
						//n--;
					}
					else
					{
						incompleteEdgesTemp.push_back(e);
					}
				}
				net.incompleteEdges = incompleteEdgesTemp;
			}

		}


		std::vector<IncompleteVertex> combineIncompleteVertices(Network& net, const Vec3c& imageDimensions, bool isFinal);

		/**
		Combines traced subnetworks into one network.
		@param subNets The traced networks.
		@param net The full network is inserted here.
		@param freeSubnets If set to true, memory allocated to the subnetworks is cleared as soon as the data is added to the full network. If true, all networks in subNets vector will contain no edges and vertices after a call to this function.
		@param isFinal Set to true if this is the last call to combineTracedBlocks and the combined network should be finalized.
		*/
		template<typename orig_t> void combineTracedBlocks(const std::vector<Network>& subNets, Network& net, const Vec3c& imageDimensions, bool freeSubnets, const Image<orig_t>* pOriginal, bool isFinal, bool storeAllEdgePoints, double smoothingSigma, double maxDisplacement)
		{
			// Combine overlapping intersection regions
			std::cout << "Combine calculation blocks..." << std::endl;

			if (subNets.size() <= 0)
				return;

			size_t counter = 0;
			net = subNets[0];
			for (size_t n = 1; n < subNets.size(); n++)
			{
				// Combine net and subNets[n]. Shift vertex indices of subNets[n] before combining.
				Network net2 = subNets[n];
				size_t vertexIndexShift = net.vertices.size();
				for (size_t m = 0; m < net2.edges.size(); m++)
				{
					net2.edges[m].verts += Vec2c(vertexIndexShift, vertexIndexShift);
				}

				for (size_t m = 0; m < net2.incompleteVertices.size(); m++)
				{
					net2.incompleteVertices[m].vertexIndex += vertexIndexShift;
				}

				for (size_t m = 0; m < net2.incompleteEdges.size(); m++)
				{
					if (net2.incompleteEdges[m].verts[0] >= 0)
						net2.incompleteEdges[m].verts[0] += vertexIndexShift;
					if (net2.incompleteEdges[m].verts[1] >= 0)
						net2.incompleteEdges[m].verts[1] += vertexIndexShift;
				}

				net.vertices.insert(net.vertices.end(), net2.vertices.begin(), net2.vertices.end());
				net.edges.insert(net.edges.end(), net2.edges.begin(), net2.edges.end());
				net.incompleteVertices.insert(net.incompleteVertices.end(), net2.incompleteVertices.begin(), net2.incompleteVertices.end());
				net.incompleteEdges.insert(net.incompleteEdges.end(), net2.incompleteEdges.begin(), net2.incompleteEdges.end());

				if (freeSubnets)
				{
					net2.edges.clear();
					net2.edges.shrink_to_fit();
					net2.vertices.clear();
					net2.vertices.shrink_to_fit();
					net2.incompleteVertices.clear();
					net2.incompleteVertices.shrink_to_fit();
					net2.incompleteEdges.clear();
					net2.incompleteEdges.shrink_to_fit();
				}

				//coord_t minZ, maxZ;
				//std::tie(minZ, maxZ) = calcMinMaxZ(subNets.size(), imageDimensions.z, n);
				//Vec3c partialDimensions(imageDimensions.x, imageDimensions.y, maxZ + 1);

				//bool isLast = n == subNets.size() - 1;

				//// Process incomplete vertices
				//std::cout << "Combine incomplete vertices on processing block boundaries..." << std::endl;
				//vector<IncompleteVertex> completedVertices = internals::combineIncompleteVertices(net, partialDimensions, isLast);

				//std::cout << "Combine incomplete edges on processing block boundaries..." << std::endl;
				//internals::combineIncompleteEdges(net, completedVertices, pOriginal, isLast);

				showThreadProgress(counter, subNets.size());
			}

			
			std::cout << "Combine incomplete vertices on processing block boundaries..." << std::endl;
			std::vector<IncompleteVertex> completedVertices = internals::combineIncompleteVertices(net, imageDimensions, isFinal);

			std::cout << "Combine incomplete edges on processing block boundaries..." << std::endl;
			internals::combineIncompleteEdges(net, completedVertices, pOriginal, isFinal, storeAllEdgePoints, smoothingSigma, maxDisplacement);
		}
	
		/**
		Moves incomplete vertices and edges to the complete vertices and edges lists.
		Throws exception if the network contains edges with unknown end points.
		*/
		inline void finalize(Network& net)
		{
			for (const IncompleteEdge& ie : net.incompleteEdges)
			{
				if (ie.verts[0] < 0 || ie.verts[1] < 0)
					throw ITLException("The incomplete edges of the network cannot be completed because one of them has unknown end point: " + toString(ie.verts[0]) + " to " + toString(ie.verts[1]));
			}

			net.incompleteVertices.clear();

			for (const IncompleteEdge& ie : net.incompleteEdges)
			{
				net.edges.push_back(Edge(ie.verts[0], ie.verts[1], ie.properties));
			}

			net.incompleteEdges.clear();
		}
	}

	/**
	Converts line skeleton into a graph.
	Divides the image into blocks and processes one block per thread.
	All pixels in the image are set to zero.
	@param net The graph is inserted into this object.
	@param img Skeletonized image. Will be set to zero at output.
	@param original Original non-skeletonized image used for measurements. This image is not modified.
	@param blockCount Number of threads to use for tracing. Set to zero to determine thread count automatically. Set to one to use single-threaded processing.
	*/
	template<typename pixel_t, typename orig_t> void traceLineSkeleton(Image<pixel_t>& img, const Image<orig_t>* pOriginal, bool storeAllEdgePoints, double smoothingSigma, double maxDisplacement, Network& net, int blockCount = 0)
	{
		if (pOriginal)
			img.mustNotBe(*pOriginal);

		internals::classifyForTracing(img);

		std::vector<Network> subNets;

		if (blockCount != 1)
		{
			internals::traceLineSkeletonBlocks(img, pOriginal, storeAllEdgePoints, smoothingSigma, maxDisplacement, subNets, Vec3sc(0, 0, 0), blockCount);

			internals::combineTracedBlocks(subNets, net, img.dimensions(), true, pOriginal, true, storeAllEdgePoints, smoothingSigma, maxDisplacement);
		}
		else
		{
			size_t counter = 0;
			itl2::internals::traceLineSkeleton(img, pOriginal, Vec3d(0, 0, 0), storeAllEdgePoints, smoothingSigma, maxDisplacement, net, counter, img.depth());
		}

		internals::finalize(net);
	}

	/**
	Converts line skeleton into a graph.
	Divides the image into blocks and processes one block per thread.
	All pixels in the image are set to zero.
	@param net The graph is inserted into this object.
	@param img Skeletonized image. Will be set to zero at output.
	@param blockCount Number of threads to use for tracing. Set to zero to determine thread count automatically. Set to one to use single-threaded processing.
	*/
	template<typename pixel_t> void traceLineSkeleton(Image<pixel_t>& img, bool storeAllEdgePoints, double smoothingSigma, double maxDisplacement, Network& net, int blockCount = 0)
	{
		traceLineSkeleton<pixel_t, pixel_t>(img, nullptr, storeAllEdgePoints, smoothingSigma, maxDisplacement, net, blockCount);
	}
	

	namespace tests
	{
		//void classifySkeleton();
		void traceSkeleton();
		void traceSkeletonRealData();
		void lineLength();
	}
}
