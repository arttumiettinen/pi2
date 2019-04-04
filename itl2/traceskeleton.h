#pragma once

#include "image.h"
#include "floodfill.h"
#include "pointprocess.h"
#include "neighbourhood.h"
#include "math/vec2.h"
#include "math/vec3.h"
#include "utilities.h"
#include "network.h"
#include "box.h"
#include "math/matrix3x3.h"
#include "math/vectoroperations.h"
#include "interpolation.h"
#include "projections.h"

#include <omp.h>

using math::Vec3c;
using math::Vec3f;
using math::Vec3d;
using math::Matrix3x3d;

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
		//constexpr int UNKNOWN = 8;
		

		/**
		Classifies all pixels in the image. Overwrites the image.
		*/
		template<typename pixel_t, bool test(Image<pixel_t>& nb), int value>
		void classify(Image<pixel_t>& img, bool showProgress)
		{
			size_t totalProcessed = 0;
			#pragma omp parallel if(!omp_in_parallel())
			{

				Image<pixel_t> nb(math::Vec3c(3, 3, 3));

				#pragma omp for
				for (coord_t z = 0; z < img.depth(); z++)
				{
					for (coord_t y = 0; y < img.height(); y++)
					{
						for (coord_t x = 0; x < img.width(); x++)
						{
							if (img(x, y, z) == UNCLASSIFIED) // Background pixels do not need to be classified.
							{
								getNeighbourhood(img, math::Vec3c(x, y, z), math::Vec3c(1, 1, 1), nb, BoundaryCondition::Zero);

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
		Set non-background edge voxels to UNKNOWN.
		See also setEdges function.
		*/
		template<typename pixel_t> void handleEdges(Image<pixel_t>& img)
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
								//img(coords) = internals::UNKNOWN;
								img(coords) = internals::BRANCHING;

							coords[skip] = img.dimension(skip) - 1;
							if (img(coords) != 0)
								//img(coords) = internals::UNKNOWN;
								img(coords) = internals::BRANCHING;
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
	template<typename pixel_t> void classifySkeleton(Image<pixel_t>& img, bool curveEnds, bool curveAndIntersectionOnly, bool showProgress)
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

		internals::handleEdges(img);
	}


	namespace internals
	{
		class BranchStartInfo
		{
		public:
			Vec3sc start;
			size_t startVertexIndex;
		};

		template<typename pixel_t> vector<BranchStartInfo> addIntersection(Image<pixel_t>& img, const Vec3sc& pos, Network& network)
		{
			// Find all points in the intersection region
			vector<Vec3sc> filledPoints;
			floodfill<pixel_t>(img, Vec3c(pos), 0, 0, Connectivity::AllNeighbours, 0, &filledPoints);

			if (filledPoints.size() <= 0)
				throw ITLException("Invalid position passed to intersection tracer. There is no intersection at that location.");

			// Calculate average position of the intersection
			bool isOnEdge = false;
			math::Vec3f center(0, 0, 0);
			for (const Vec3sc& p : filledPoints)
			{
				center += math::Vec3f(p);
				if (img.isOnEdge(p))
					isOnEdge = true;
			}
			center /= (float32_t)filledPoints.size();

			// Insert the new vertex to the network
			size_t vertexIndex = network.vertices.size();
			network.vertices.push_back(center);

			// If the vertex is on edge, insert point list to the network, too.
			if(isOnEdge)
				network.incompleteVertices.push_back(IncompleteVertex(vertexIndex, filledPoints));

			// Find all branches that start in point neighbouring the intersection region.
			// Add them to processing list.
			vector<BranchStartInfo> branchStarts;
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
							}
						}
					}
				}
			}

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
		Trace single branch of the skeleton and return the other end point of the branch, and erases it from the image.
		@param count The count of pixels traversed is stored in this variable on output.
		*/
		template<typename pixel_t> tuple<vector<Vec3f>, Vec3sc> traceBranch(Image<pixel_t>& img, Vec3sc start, size_t& count)
		{
			vector<Vec3f> points;
			
			count = 0;
			while (true)
			{
				points.push_back(Vec3f(start));
				img(start) = 0;
				count++;

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
		float32_t lineLength(vector<Vec3f>& points, float32_t* pStraightLength, double sigma = 1);

		/**
		Estimates length of straight line between a and b using same method than lineLength function.
		*/
		//float32_t straightLineLength(const Vec3f& a, const Vec3f& b, double sigma = 1);

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
					Vec3c xc = math::round(x);
					
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
		*/
		template<typename orig_t> EdgeMeasurements measureEdge(const Image<orig_t>* pOriginal, vector<Vec3f> edgePoints)
		{
			EdgeMeasurements result;

			result.pointCount = (float32_t)edgePoints.size();

			if (edgePoints.size() > 0)
			{
				result.distance = (edgePoints[0] - edgePoints[edgePoints.size() - 1]).norm();

				// Calculate length and smoothed point list
				result.length = lineLength(edgePoints, 0);

				result.adjustedStart = edgePoints[0];
				result.adjustedEnd = edgePoints[edgePoints.size() - 1];
			}
			else
			{
				result.distance = numeric_limits<float32_t>::signaling_NaN();
			}

			// Default value for area used if area measurements fail.
			result.area = numeric_limits<float32_t>::signaling_NaN();

			// Calculate area by extracting a slice through each point on the path
			// If there are less than three points, we cannot measure area. (two points would be possible but a special case)
			if (edgePoints.size() >= 3 && pOriginal)
			{
				vector<size_t> areaSamples;
				areaSamples.reserve(edgePoints.size());
				for (size_t n = 1; n < edgePoints.size() - 1; n++)
				{
					Vec3d center = Vec3d(edgePoints[n]);
					Vec3d tangent = Vec3d(edgePoints[n + 1] - edgePoints[n - 1]);
					Image<orig_t> slice(75, 75, 0); // TODO: Hardcoded slice size
					getSlice(*pOriginal, center, tangent, slice);

					//raw::writed(slice, "./skeleton/tracing/slice");

					orig_t M = max(slice);
					threshold(slice, M / 2);

					//raw::writed(slice, "./skeleton/tracing/slice_th");

					vector<Vec3sc> filledPoints;
					floodfill<orig_t>(slice, slice.dimensions() / 2, (orig_t)0, (orig_t)0, Connectivity::AllNeighbours, 0, &filledPoints);

					//raw::writed(slice, "./skeleton/tracing/slice_fill");

					if(filledPoints.size() > 0)
						areaSamples.push_back(filledPoints.size());
				}

				if(areaSamples.size() > 0)
					result.area = (float32_t)mode(areaSamples);
			}

			return result;
		}

		/**
		Trace branches in the given seeds list and continue recursively to all branches that are connected to them.
		Add everything to the network.

		Single pixel long branch connecting to single intersection is missed.
		@param img Classified skeleton image.
		@param orig Original (non-skeletonized) image that is used for branch area and shape measurements.
		*/
		template<typename pixel_t, typename orig_t> void trace(Image<pixel_t>& img, const Image<orig_t>* pOriginal, const vector<BranchStartInfo>& seeds, Network& net)
		{
			// Add branches to the tracing list
			deque<internals::BranchStartInfo> tracingList;
			for (const BranchStartInfo& f : seeds)
				tracingList.push_front(f);

			// Process branches until the tracing list is empty
			while (!tracingList.empty())
			{
				BranchStartInfo info = tracingList.front();
				tracingList.pop_front();

				size_t count = 0;
				tuple<vector<Vec3f>, Vec3sc> result = traceBranch(img, info.start, count);
				vector<Vec3f> edgePoints = get<0>(result);
				Vec3sc end = get<1>(result);
				edgePoints.insert(edgePoints.begin(), net.vertices[info.startVertexIndex]);

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
								vector<BranchStartInfo> newBranches = addIntersection(img, pp, net);
								
								// The current branch is already erased so it is not going to be in the newBranches list.
								for (BranchStartInfo& f : newBranches)
								{
									tracingList.push_front(f);
								}

								// Add edge from start to this intersection. Measure its shape before adding. Process only valid edges.
								size_t endVertexIndex = net.vertices.size() - 1;
								vector<Vec3f> currEdge;
								currEdge.reserve(edgePoints.size() + 1);
								currEdge.insert(currEdge.begin(), edgePoints.begin(), edgePoints.end());
								currEdge.push_back(net.vertices[endVertexIndex]);
								if (*currEdge.begin() != *currEdge.rbegin())
								{
									Edge e(info.startVertexIndex, endVertexIndex, internals::measureEdge(pOriginal, currEdge));
									net.edges.push_back(e);
								}
							}
						}
					}
				}

				// Find known edges where the end point is a start point (to resolve loops)
				for (size_t n = 0; n < tracingList.size(); n++)
				{
					if (tracingList[n].start == end)
					{
						// Found start point corresponding to our end.
						vector<Vec3f> currEdge;
						currEdge.reserve(edgePoints.size() + 1);
						currEdge.insert(currEdge.begin(), edgePoints.begin(), edgePoints.end());
						currEdge.push_back(net.vertices[tracingList[n].startVertexIndex]);
						if (*currEdge.begin() != *currEdge.rbegin())
						{
							Edge e(info.startVertexIndex, tracingList[n].startVertexIndex, internals::measureEdge(pOriginal, currEdge));
							net.edges.push_back(e);
						}

						tracingList.erase(tracingList.begin() + n);
						n--;
					}
				}
			}
		}


		/**
		Trace line skeleton into a graph structure using single-threaded processing.
		All pixels in the image are set to zero.
		@param img Non-classified skeleton image. The image will be empty at output.
		@param original Original non-skeletonized image used for shape measurements.
		@net The network is inserted to this object.
		*/
		template<typename pixel_t, typename orig_t> void traceLineSkeleton(Image<pixel_t>& img, const Image<orig_t>* pOriginal, Network& net, size_t& counter, size_t progressMax)
		{
			if (pOriginal)
			{
				img.mustNotBe(*pOriginal);
				img.checkSize(*pOriginal);
			}

			classifySkeleton(img, true, false, false);

			// Find intersection point or fibre end point
			for (coord_t z = 0; z < img.depth(); z++)
			{
				for (coord_t y = 0; y < img.height(); y++)
				{
					for (coord_t x = 0; x < img.width(); x++)
					{
						pixel_t label = img(x, y, z);
						if (label == internals::ENDPOINT || label == internals::BRANCHING)
						{
							// Find fibres that start from this intersection area or end point
							vector<internals::BranchStartInfo> branches = internals::addIntersection(img, Vec3sc((int32_t)x, (int32_t)y, (int32_t)z), net);

							// Follow all branches and process the found intersection areas
							internals::trace(img, pOriginal, branches, net);
						}
					}
				}

				showThreadProgress(counter, progressMax);
			}
		}

		/**
		Traces line skeleton in blocks, and adds the network resulting from each block trace to subNets vector.
		@param origin If processing a block of larger image, set this to the origin of the block in full image coordinates. This value is added to all the vertex coordinates.
		*/
		template<typename pixel_t, typename orig_t> void traceLineSkeletonBlocks(Image<pixel_t>& img, Image<orig_t>* pOriginal, vector<Network>& subNets, const Vec3sc& origin = Vec3sc())
		{
			if (pOriginal)
			{
				img.mustNotBe(*pOriginal);
				img.checkSize(*pOriginal);
			}

			size_t counter = 0;
			cout << "Tracing skeleton..." << endl;
#pragma omp parallel if(!omp_in_parallel() && img.pixelCount() > PARALLELIZATION_THRESHOLD)
			{
				int idx = omp_get_thread_num();
				int count = omp_get_num_threads();
				if (idx == 0)
				{
#pragma omp critical
					{
						cout << "Tracing in " << count << " blocks." << endl;
					}
				}

				// Calculate amount of slices single thread should process.
				// The value is capped so that we don't divide the work unnecessarily too much, if the image is small.
				coord_t size = img.depth() / count;
				if (size < 20)
					size = 20;

				// If there's nothing to do for all the threads, the excess threads will just skip processing.
				coord_t minZ = idx * size;
				if (minZ < img.depth())
				{

					coord_t maxZ = minZ + size - 1;
					if (maxZ >= img.depth())
						maxZ = img.depth() - 1;

					if (idx == count - 1)
					{
						// The last thread processes possible "rounding error" slices
						maxZ = img.depth() - 1;
					}

					// Get view of part of the image
					Image<pixel_t> block(img, minZ, maxZ);
					Image<orig_t> origBlock;
					if(pOriginal)
						origBlock.init(*pOriginal, minZ, maxZ);

					//raw::writed(block, string("./skeleton/tracing/block") + toString(idx));

					// Trace the region
					Network subNet;
					internals::traceLineSkeleton(block, pOriginal ? &origBlock : 0, subNet, counter, img.depth());

					// Convert vertices to global coordinates (they are in block coordinates)
					for (size_t n = 0; n < subNet.vertices.size(); n++)
						subNet.vertices[n] += math::Vec3f(origin) + math::Vec3f(0, 0, (float)minZ);
					for (size_t n = 0; n < subNet.incompleteVertices.size(); n++)
						for (size_t m = 0; m < subNet.incompleteVertices[n].points.size(); m++)
							subNet.incompleteVertices[n].points[m] += origin + math::Vec3sc(0, 0, (int32_t)minZ);

#pragma omp critical(traceLineSkeleton)
					{
						subNets.push_back(subNet);
					}
				}
			}
		}

		/**
		Combines traced subnetworks into one network.
		@param subNets The traced networks.
		@param net The full network is inserted here.
		@param freeSubnets If set to true, memory allocated to the subnetworks is cleared as soon as the data is added to the full network. If true, all networks in subNets vector will contain no edges and vertices after a call to this function.
		*/
		void combineTracedBlocks(const vector<Network>& subNets, Network& net, bool freeSubnets = true);
	}

	/**
	Converts line skeleton into a graph.
	Divides the image into blocks and processes one block per thread.
	All pixels in the image are set to zero.
	@param net The graph is inserted into this object.
	@param img Skeletonized image. Will be set to zero at output.
	@param original Original non-skeletonized image used for measurements. This image is not modified.
	*/
	template<typename pixel_t, typename orig_t> void traceLineSkeleton(Image<pixel_t>& img, Image<orig_t>* pOriginal, Network& net)
	{
		vector<Network> subNets;
		internals::traceLineSkeletonBlocks(img, pOriginal, subNets, Vec3sc(0, 0, 0));

		internals::combineTracedBlocks(subNets, net, true);
	}

	/**
	Converts line skeleton into a graph.
	Divides the image into blocks and processes one block per thread.
	All pixels in the image are set to zero.
	@param net The graph is inserted into this object.
	@param img Skeletonized image. Will be set to zero at output.
	*/
	template<typename pixel_t> void traceLineSkeleton(Image<pixel_t>& img, Network& net)
	{
		traceLineSkeleton<pixel_t, pixel_t>(img, 0, net);
	}
	

	namespace tests
	{
		//void classifySkeleton();
		void traceSkeleton();
		void traceSkeletonRealData();
		void lineLength();
	}
}
