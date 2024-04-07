
#include "traceskeleton.h"

#include "io/raw.h"
#include "filters.h"
#include "transform.h"
#include "generation.h"

using namespace std;

namespace itl2
{

	namespace internals
	{
		/**
		Smooths dimension dim of points in the given list and anchors them so that they stay near their original locations.
		The distance between original location of a point and new location won't be more than delta.
		*/
		void smoothAndAnchor(vector<Vec3f>& points, double sigma, double delta)
		{

			// Smooth
			vector<Vec3f> smoothed(points.size());

			for (size_t dim = 0; dim < 3; dim++)
			{
				Image<float32_t> c(points.size());
				for (size_t n = 0; n < points.size(); n++)
					c(n) = points[n][dim];

				gaussFilter(c, Vec3d(sigma, 0.0, 0.0), BoundaryCondition::Nearest);

				for (size_t n = 0; n < points.size(); n++)
					smoothed[n][dim] = c(n);
			}

			// Anchor to spheres around pixel centers
			for (size_t n = 0; n < points.size(); n++)
			{
				Vec3f& orig = points[n];
				Vec3f sm = smoothed[n];
				Vec3f dir = sm - orig;
				float32_t l = dir.norm();
				if (l <= delta)
				{
					orig = sm;
				}
				else
				{
					orig = orig + dir / l * (float32_t)delta;
				}
			}


			/*
			for (size_t dim = 0; dim < 3; dim++)
			{
				Image<float32_t> c(points.size());
				for (size_t n = 0; n < points.size(); n++)
					c(n) = points[n][dim];

				gaussFilter(c, Vec3d(sigma, 0.0, 0.0), BoundaryCondition::Nearest);

				for (size_t n = 0; n < points.size(); n++)
				{
					float32_t orig = points[n][dim];
					float32_t smoothed = c(n);
					if (smoothed > orig + delta)
						points[n][dim] = orig + delta;
					else if (smoothed < orig - delta)
						points[n][dim] = orig - delta;
					else
						points[n][dim] = c(n);
				}
			}
			*/
		}

		float32_t lineLength(vector<Vec3f>& points, double sigma, double maxDisplacement)
		{
			if (sigma > 0)
				smoothLine(points, sigma, maxDisplacement);

			float32_t l = 0;
			for (size_t n = 1; n < points.size(); n++)
			{
				l += (points[n] - points[n - 1]).norm();
			}

			return l;
		}



		bool isNeighbour(const vector<Vec3sc>& v1, const Vec3sc& p)
		{
			for (size_t n = 0; n < v1.size(); n++)
			{
				if (isNeighbour(p, v1[n]))
				{
					return true;
				}
			}

			return false;
		}


		/**
		Finds out if there are any neighbouring points in v1 and v2.
		*/
		bool isNeighbour(const vector<Vec3sc>& v1, const vector<Vec3sc>& v2)
		{
			for (size_t n = 0; n < v1.size(); n++)
			{
				if (isNeighbour(v2, v1[n]))
					return true;
				//for (size_t m = 0; m < v2.size(); m++)
				//{
				//	if ((v2[m] - v1[n]).abs().max() <= 1)
				//	{
				//		return true;
				//	}
				//}
			}

			return false;
		}

		bool bordersEdge(const vector<Vec3sc>& points, const Vec3c& dims)
		{
			for (const Vec3sc& p : points)
			{
				if (p.x <= 1 || p.x >= dims.x - 2 ||
					p.y <= 1 || p.y >= dims.y - 2 ||
					p.z <= 1 || p.z >= dims.z - 2)
					return true;
			}
			return false;
		}

		///**
		//Finds out if there are any overlapping points in v1 and v2.
		//*/
		//bool overlaps(const vector<Vec3sc>& v1, const vector<Vec3sc>& v2)
		//{
		//	for (size_t n = 0; n < v1.size(); n++)
		//	{
		//		for (size_t m = 0; m < v2.size(); m++)
		//		{
		//			if(v2[m] == v1[n])
		//			{
		//				return true;
		//			}
		//		}
		//	}

		//	return false;
		//}

		/**
		Finds intersecting incomplete vertex/edge regions and combines their indices in the forest.
		At output the forest will contain all indices of vertices list, and indices of intersecting vertices belong to the same root.

		T must be type with .points member that is a vector containing the points in the region
		*/
		template<typename T>
		void findIntersectingRegions(const vector<T>& vertices, IndexForest& forest)
		{
			// Add all indices to the forest
			cout << "Build forest..." << endl;
			forest.initialize(vertices.size());

			// Calculate bounding boxes for incomplete vertex sets
			cout << "Determine bounding boxes..." << endl;
			vector<AABox<int32_t> > boundingBoxes;
			boundingBoxes.reserve(vertices.size());
			Vec3sc m = Vec3sc(numeric_limits<int32_t>::max(), numeric_limits<int32_t>::max(), numeric_limits<int32_t>::max());
			Vec3sc M = Vec3sc(numeric_limits<int32_t>::lowest(), numeric_limits<int32_t>::lowest(), numeric_limits<int32_t>::lowest());
			for (size_t n = 0; n < vertices.size(); n++)
			{
				AABox<int32_t> b = AABox<int32_t>::boundingBox(vertices[n].points);
				b.inflate(1);
				boundingBoxes.push_back(b);
				m = min(m, b.minc);
				M = max(M, b.maxc);
			}

			Image<vector<size_t> > grid(100, 100, 100);

			cout << "Divide boxes to a grid..." << endl;
			size_t counter = 0;
			for (size_t n = 0; n < boundingBoxes.size(); n++)
			{
				Vec3c mcell((coord_t)floor((double)(boundingBoxes[n].minc.x - m.x) / (double)(M.x - m.x) * (grid.width() - 1)),
					(coord_t)floor((double)(boundingBoxes[n].minc.y - m.y) / (double)(M.y - m.y) * (grid.height() - 1)),
					(coord_t)floor((double)(boundingBoxes[n].minc.z - m.z) / (double)(M.z - m.z) * (grid.depth() - 1)));

				Vec3c Mcell((coord_t)ceil((double)(boundingBoxes[n].maxc.x - m.x) / (double)(M.x - m.x) * (grid.width() - 1)),
					(coord_t)ceil((double)(boundingBoxes[n].maxc.y - m.y) / (double)(M.y - m.y) * (grid.height() - 1)),
					(coord_t)ceil((double)(boundingBoxes[n].maxc.z - m.z) / (double)(M.z - m.z) * (grid.depth() - 1)));

				for (coord_t zc = mcell.z; zc <= Mcell.z; zc++)
				{
					for (coord_t yc = mcell.y; yc <= Mcell.y; yc++)
					{
						for (coord_t xc = mcell.x; xc <= Mcell.x; xc++)
						{
							grid(xc, yc, zc).push_back(n);
						}
					}
				}

				showThreadProgress(counter, boundingBoxes.size());
			}

			cout << "Determine overlaps..." << endl;
			counter = 0;
#pragma omp parallel if(!omp_in_parallel())
			{
				set<tuple<size_t, size_t> > nonOverlapping;
#pragma omp for
				for (coord_t i = 0; i < grid.pixelCount(); i++)
				{
					auto& list = grid(i);

					for (coord_t ni = 0; ni < (coord_t)list.size(); ni++)
					{
						size_t n = list[ni];
						const T& v1 = vertices[n];
						const AABox<int32_t>& b1 = boundingBoxes[n];

						for (size_t mi = ni + 1; mi < list.size(); mi++)
						{
							size_t m = list[mi];
							const T& v2 = vertices[m];
							const AABox<int32_t>& b2 = boundingBoxes[m];

							// Only test if the vertices have not been connected yet.
							// There is a race condition between union_sets below and find_sets on this line, but
							// that may only result in the test below evaluating true even though it should be
							// false, and in that case we just do some extra work.
							if (forest.find_set(n) != forest.find_set(m))
							{
								if (b1.overlapsInclusive(b2)) // Do bounding boxes overlap?
								{
									auto key = make_tuple(std::min(n, m), std::max(n, m));

									// Only test points if they have not been determined to be non-overlapping.
									if (nonOverlapping.find(key) == nonOverlapping.end())
									{
										if (isNeighbour(v1.points, v2.points))
										{
											// v1 and v2 are neighbours, so the intersection regions must be combined
#pragma omp critical(forestInsert)
											{
												forest.union_sets(n, m);
											}
										}
										else
										{
											nonOverlapping.emplace(key);
										}
									}
								}
							}
						}
					}

					showThreadProgress(counter, grid.pixelCount());
				}
			}
		}

		/**
		Combines all incomplete vertices that represent the same vertex.
		After the call net.incompleteVertices is empty.
		The network may contain unnecessary vertices.
		@return list of vertices that were completed. Remaining incomplete vertices cannot be completed without more data (they are at the edges of the image).
		*/
		vector<IncompleteVertex> combineIncompleteVertices(Network& net, const Vec3c& imageDimensions, bool isFinal)
		{
			// Find neighbouring incomplete vertices
			// The forest contains indices into net.incompleteVertices list.
			IndexForest forest;
			findIntersectingRegions(net.incompleteVertices, forest);

			// Move all points to roots
			cout << "Move points to roots..." << endl;
			size_t counter = 0;
			for (coord_t n = 0; n < (coord_t)net.incompleteVertices.size(); n++)
			{
				size_t base = forest.find_set(n);
				if (base != n)
				{
					auto& target = net.incompleteVertices[base].points;
					auto& source = net.incompleteVertices[n].points;
					target.insert(target.end(), source.begin(), source.end());
					source.clear();
					source.shrink_to_fit();
					// Make the vertex to be deleted invalid.
					net.vertices[net.incompleteVertices[n].vertexIndex] = Network::INVALID_VERTEX;
				}

				showThreadProgress(counter, net.incompleteVertices.size());
			}

			// Calculate true position of each incomplete vertex.
			cout << "Calculate new positions for root vertices..." << endl;
			counter = 0;
			#pragma omp parallel for if(!omp_in_parallel()) // Can be parallelized because incompleteVertices[n].vertexIndex is different for each n.
			for (coord_t n = 0; n < (coord_t)net.incompleteVertices.size(); n++)
			{
				auto& v = net.incompleteVertices[n].points;
				if (v.size() > 0)
				{
					// Find unique points
					//sort(v.begin(), v.end(), vecComparer<int32_t>);
					//auto last = unique(v.begin(), v.end());
					////v.erase(last, v.end());
					//testAssert(last == v.end(), "Intersection region contains non-unique points.");

					// Calculate center point of the vertex
					Vec3f center(0, 0, 0);
					for (const auto& p : v)
						center += Vec3f(p);
					center /= (float32_t)v.size();


//if (center == Vec3f(161, 67, 85))
//{
//	cout << "found" << endl;
//}

					net.vertices[net.incompleteVertices[n].vertexIndex] = center;
				}
				showThreadProgress(counter, net.incompleteVertices.size());
			}

			// Replace removed vertices by the new ones.
			//cout << "Replace removed vertices by the new ones..." << endl;
			//counter = 0;
			//#pragma omp parallel for if(!omp_in_parallel()) // Can be parallelized because net.incompleteVertices[n].vertexIndex is different for each n.
			//for (coord_t n = 0; n < (coord_t)net.incompleteVertices.size(); n++)
			//{
			//	if (net.incompleteVertices[n].points.size() <= 0)
			//	{
			//		size_t m = forest.find_set(n);
			//		size_t sourceVertexIndex = net.incompleteVertices[n].vertexIndex;
			//		size_t targetVertexIndex = net.incompleteVertices[m].vertexIndex;

			//		// Replace each occurence of sourceVertexIndex by targetVertexIndex
			//		for (size_t i = 0; i < net.edges.size(); i++)
			//		{
			//			auto& edge = net.edges[i];
			//			if (edge.verts[0] == sourceVertexIndex)
			//				edge.verts[0] = targetVertexIndex;
			//			if (edge.verts[1] == sourceVertexIndex)
			//				edge.verts[1] = targetVertexIndex;
			//		}
			//	}

			//	showThreadProgress(counter, net.incompleteVertices.size());
			//}

			cout << "Create map of vertex indices that have changed..." << endl;
			map<size_t, size_t> vertexIndexMap;
			counter = 0;
			for (coord_t n = 0; n < (coord_t)net.incompleteVertices.size(); n++)
			{
				if (net.incompleteVertices[n].points.size() <= 0)
				{
					size_t m = forest.find_set(n);
					size_t sourceVertexIndex = net.incompleteVertices[n].vertexIndex;
					size_t targetVertexIndex = net.incompleteVertices[m].vertexIndex;
					vertexIndexMap[sourceVertexIndex] = targetVertexIndex;
				}
				showThreadProgress(counter, net.incompleteVertices.size());
			}

			cout << "Change vertex indices in edges list..." << endl;
			counter = 0;
			#pragma omp parallel for if(!omp_in_parallel())
			for (coord_t i = 0; i < (coord_t)net.edges.size(); i++)
			{
				auto& edge = net.edges[i];

				const auto it = vertexIndexMap.find(edge.verts[0]);
				if (it != vertexIndexMap.end())
					edge.verts[0] = it->second;

				const auto it2 = vertexIndexMap.find(edge.verts[1]);
				if (it2 != vertexIndexMap.end())
					edge.verts[1] = it2->second;

				showThreadProgress(counter, net.edges.size());
			}

			// Change indices also in incomplete edge list
			counter = 0;
			#pragma omp parallel for if(!omp_in_parallel())
			for (coord_t i = 0; i < (coord_t)net.incompleteEdges.size(); i++)
			{
				auto& edge = net.incompleteEdges[i];

				const auto it = vertexIndexMap.find(edge.verts[0]);
				if (it != vertexIndexMap.end())
					edge.verts[0] = it->second;

				const auto it2 = vertexIndexMap.find(edge.verts[1]);
				if (it2 != vertexIndexMap.end())
					edge.verts[1] = it2->second;

				showThreadProgress(counter, net.incompleteEdges.size());
			}

			// Remove items corresponding to the completed incomplete vertices,
			// and change vertex indices in remaining incomplete vertices.
			{
				vector<IncompleteVertex> incompleteVerticesTemp;
				for (coord_t n = 0; n < (coord_t)net.incompleteVertices.size(); n++)
				{
					if (net.incompleteVertices[n].points.size() <= 0)
					{
						//net.incompleteVertices.erase(net.incompleteVertices.begin() + n);
						//n--;
					}
					else
					{
						// The vertex remains in the incomplete list.
						incompleteVerticesTemp.push_back(net.incompleteVertices[n]);
					}

					//else
					//{
					//	const auto& it = vertexIndexMap.find(net.incompleteVertices[n].vertexIndex);
					//	if (it != vertexIndexMap.end())
					//	{
					//		cout << "old = " << net.incompleteVertices[n].vertexIndex << ", new = " << it->second << endl;
					//		net.incompleteVertices[n].vertexIndex = it->second;
					//	}
					//}
				}
				net.incompleteVertices = incompleteVerticesTemp;
			}

			// Remove invalid nodes and UPDATE NODE INDICES in ALL arrays.
			net.removeInvalidNodes();

			// Now find all incomplete vertices that are on edge. Those remain as incomplete.
			// Put everything else to another array.
			vector<IncompleteVertex> completedVertices;
			{
				vector<IncompleteVertex> incompleteVerticesTemp;
				for (coord_t n = 0; n < (coord_t)net.incompleteVertices.size(); n++)
				{
					if ((isFinal && !isOnEdge(net.incompleteVertices[n].points, imageDimensions)) ||
						(!isFinal && !bordersEdge(net.incompleteVertices[n].points, imageDimensions)))
					{
						completedVertices.push_back(net.incompleteVertices[n]);
						//net.incompleteVertices.erase(net.incompleteVertices.begin() + n);
						//n--;
					}
					else
					{
						// The vertex remains in the incomplete list.
						incompleteVerticesTemp.push_back(net.incompleteVertices[n]);
					}
				}
				net.incompleteVertices = incompleteVerticesTemp;
			}

			return completedVertices;
		}


		void insertPoint(const Vec3sc& p, size_t n, Image<vector<size_t> >& grid, const Vec3sc& m, const Vec3sc& M)
		{
			Vec3sc minc = p - Vec3sc(1, 1, 1);
			Vec3sc maxc = p + Vec3sc(1, 1, 1);
			Vec3c mcell((coord_t)floor((double)(minc.x - m.x) / (double)(M.x - m.x) * (grid.width() - 1)),
				(coord_t)floor((double)(minc.y - m.y) / (double)(M.y - m.y) * (grid.height() - 1)),
				(coord_t)floor((double)(minc.z - m.z) / (double)(M.z - m.z) * (grid.depth() - 1)));
			Vec3c Mcell((coord_t)ceil((double)(maxc.x - m.x) / (double)(M.x - m.x) * (grid.width() - 1)),
				(coord_t)ceil((double)(maxc.y - m.y) / (double)(M.y - m.y) * (grid.height() - 1)),
				(coord_t)ceil((double)(maxc.z - m.z) / (double)(M.z - m.z) * (grid.depth() - 1)));
			for (coord_t zc = mcell.z; zc <= Mcell.z; zc++)
			{
				for (coord_t yc = mcell.y; yc <= Mcell.y; yc++)
				{
					for (coord_t xc = mcell.x; xc <= Mcell.x; xc++)
					{
						grid(xc, yc, zc).push_back(n);
					}
				}
			}
		}


		/**
		Finds incomplete vertex that neighbours given point.
		Returns vertexIndex of that vertex.
		@param invalidIndex this index is never returned.
		*/
		coord_t findVertex(const vector<IncompleteVertex>& incompleteVertices, const Vec3sc& p, const coord_t invalidIndex)
		{
			// Find all neighbouring indices that are not equal to 'invalidIndex'
			vector<size_t> candidates;
			for (const auto& v : incompleteVertices)
			{
				if (v.vertexIndex != invalidIndex && isNeighbour(v.points, p))
					candidates.push_back(v.vertexIndex);
			}

			// Remove duplicates
			sort(candidates.begin(), candidates.end());
			candidates.erase(unique(candidates.begin(), candidates.end()), candidates.end());

			if (candidates.size() > 0)
				return candidates[0];

			return -1;
		}


	}
	
	
	void smoothLine(vector<Vec3f>& points, double sigma, double delta)
	{
		// Anchor start and end points to their locations,
		// use anchored convolution on other points.
		Vec3f start = points.front();
		Vec3f end = points.back();
		internals::smoothAndAnchor(points, sigma, delta);
		points.front() = start;
		points.back() = end;
	}

	namespace tests
	{
		//void classifySkeleton()
		//{
		//	// NOTE: No asserts!
		
		//	Image<uint8_t> skele;
		//	raw::read(skele, "./skeleton/cavities/cavity_hybrid_skeleton");

		//	classifySkeleton(skele, true, true, true);

		//	raw::writed(skele, "./skeleton/cavities/classified_cavity_hybrid_skeleton");
		//}

		vector<Edge> getEdges(const Network& net, const Vec3f& start, const Vec3f& end)
		{
			vector<Edge> edg;
			for (const Edge& e : net.edges)
			{
				if (e.verts[0] >= 0 && e.verts[1] >= 0)
				{
					Vec3f a = net.vertices[e.verts[0]];
					Vec3f b = net.vertices[e.verts[1]];
					if ((start.equals(a, 0.1f) && end.equals(b, 0.1f)) ||
						(end.equals(a, 0.1f) && start.equals(b, 0.1f)))
					{
						edg.push_back(e);
					}
				}
			}

			return edg;
		}

		bool findEdge(const Network& net, const Vec3f& start, const Vec3f& end, const EdgeMeasurements& props, bool propsMustEqual)
		{
			for (const Edge& e : net.edges)
			{
				if (e.verts[0] >= 0 && e.verts[1] >= 0)
				{
					Vec3f a = net.vertices[e.verts[0]];
					Vec3f b = net.vertices[e.verts[1]];
					if ((start.equals(a, 0.1f) && end.equals(b, 0.1f)) ||
						(end.equals(a, 0.1f) && start.equals(b, 0.1f)))
					{
						if((propsMustEqual && e.properties.equals(props)) || !propsMustEqual)
							return true;
					}
				}
			}

			return false;
		}

		template<typename T> bool findVert(const Vec3<T>& v, const vector<Vec3<T> >& list)
		{
			for (const Vec3<T>& vv : list)
			{
				if (vv.equals(v, (T)0.01))
					return true;
			}
			return false;
		}

		bool equals(const IncompleteVertex& a, const IncompleteVertex& b)
		{
			for (const Vec3sc& v : a.points)
			{
				if (!findVert(v, b.points))
					return false;
			}
			
			for (const Vec3sc& v : b.points)
			{
				if (!findVert(v, a.points))
					return false;
			}

			return true;
		}

		bool findVert(const IncompleteVertex& v, const vector<IncompleteVertex>& list)
		{
			for (const IncompleteVertex& vv : list)
			{
				if (equals(vv, v))
					return true;
			}
			return false;
		}

		bool equals(const vector<Vec3sc>& a, const vector<Vec3sc>& b)
		{
			if (a.size() != b.size())
				return false;

			for (size_t n = 0; n < a.size(); n++)
			{
				if (a[n] != b[n])
					return false;
			}

			return true;
		}

		bool reverseEquals(const vector<Vec3sc>& a, const vector<Vec3sc>& b)
		{
			if (a.size() != b.size())
				return false;

			for (size_t n = 0; n < a.size(); n++)
			{
				if (a[n] != b[b.size() - 1 - n])
					return false;
			}

			return true;
		}

		bool findEdge(const IncompleteEdge& e, const vector<IncompleteEdge>& list)
		{
			for (const auto& vv : list)
			{
				if (equals(vv.points, e.points) || reverseEquals(vv.points, e.points))
					return true;
			}
			return false;
		}


		void assertEdges(const Network& a, const Network& b, const string& aname, const string& bname)
		{
			bool differences = false;
			bool areaDifferencesOnly = true;
			size_t areaDifferenceCount = 0;
			for (const Edge& e : a.edges)
			{
				bool found = false;
				if (testAssert(e.verts[0] >= 0 && e.verts[1] >= 0, "invalid vertex in " + aname))
				{
					Vec3f start = a.vertices[e.verts[0]];
					Vec3f end = a.vertices[e.verts[1]];

					bool fe = findEdge(b, start, end, e.properties, true);
					//if (!testAssert(fe, name + ", a edge from " + toString(start) + " to " + toString(end) + " not found from b. Edge properties: " + toString(e.properties)))
					if (!fe)
					{
						differences = true;

						auto list = getEdges(b, start, end);

						bool areaDifference = false;
						float32_t areaCandidate = 0;
						for (Edge& ee : list)
						{
							Vec3f estart = b.vertices[ee.verts[0]];
							Vec3f eend = b.vertices[ee.verts[1]];
							if ((estart.equals(start, 0.1f) && eend.equals(end, 0.1f) && e.properties.pointCount == ee.properties.pointCount && NumberUtils<float32_t>::equals(e.properties.length, ee.properties.length, 0.1f)) ||
								estart.equals(end, 0.1f) && eend.equals(start, 0.1f) && e.properties.pointCount == ee.properties.pointCount && NumberUtils<float32_t>::equals(e.properties.length, ee.properties.length, 0.1f))
							{
								areaDifference = true;
								areaCandidate = ee.properties.area;
								break;
							}
						}

						if (areaDifference)
						{
							cout << "Area difference, edge " << start << " -> " << end << ". N = " << e.properties.pointCount << ", area in " << aname << " = " << e.properties.area << ", area candidate in " << bname << " = " << areaCandidate << endl;
							areaDifferenceCount++;
						}
						else
						{
							areaDifferencesOnly = false;

							testAssert(false, aname + " edge from " + toString(start) + " to " + toString(end) + " not found from " + bname + ". Edge properties: " + toString(e.properties));

							cout << "Candidates:" << endl;
							if (list.size() > 0)
							{
								for (Edge& e : list)
								{
									Vec3f sv = b.vertices[e.verts[0]];
									Vec3f ev = b.vertices[e.verts[1]];
									cout << sv << " -> " << ev << ", " << e.properties << endl;
								}
							}
							else
							{
								cout << "none" << endl;
							}
						}
					}
				}
			}

			if (differences)
			{
				cout << "DIFFERENCES FOUND" << endl;
				if (areaDifferencesOnly)
					cout << "AREA DIFFERENCES ONLY!" << endl;
				cout << "Area difference count: " << areaDifferenceCount << endl;
			}
		}

		/**
		Checks that two networks have same nodes and edges.
		*/
		void assertNetworks(const Network& a, const Network& b, const string& name)
		{
			// Vertices
			//vector<Vec3f> vertsa = a.vertices;
			//vector<Vec3f> vertsb = b.vertices;
			//sort(vertsa.begin(), vertsa.end(), vecComparer<float>);
			//sort(vertsb.begin(), vertsb.end(), vecComparer<float>);

			testAssert(a.vertices.size() == b.vertices.size(), name + ", vertex count differs.");

			for (const Vec3f& v : a.vertices)
			{
				testAssert(findVert(v, b.vertices), "Vertex " + toString(v) + " not found from b.");
			}

			for (const Vec3f& v : b.vertices)
			{
				testAssert(findVert(v, a.vertices), "Vertex " + toString(v) + " not found from a.");
			}



			// Incomplete edges and vertices
			testAssert(a.incompleteVertices.size() == b.incompleteVertices.size(), name + ", incomplete vertex count differs.");

			for (const auto& v : a.incompleteVertices)
			{
				testAssert(findVert(v, b.incompleteVertices), "Incomplete vertex " + toString(v) + " not found from b.");
			}

			for (const auto& v : b.incompleteVertices)
			{
				testAssert(findVert(v, a.incompleteVertices), "Incomplete vertex " + toString(v) + " not found from a.");
			}


			testAssert(a.incompleteEdges.size() == b.incompleteEdges.size(), name + ", incomplete edge count differs.");


			for (const auto& v : a.incompleteEdges)
			{
				testAssert(findEdge(v, b.incompleteEdges), "Incomplete edge of a from " + toString(v.points.front()) + " to " + toString(v.points.back()) + " not found from b.");
			}

			for (const auto& v : b.incompleteEdges)
			{
				testAssert(findEdge(v, a.incompleteEdges), "Incomplete edge of b from " + toString(v.points.front()) + " to " + toString(v.points.back()) + " not found from a.");
			}


			// Edges
			testAssert(a.edges.size() == b.edges.size(), name + ", different number of edges.");

			assertEdges(a, b, "a", "b");
			assertEdges(b, a, "b", "a");
		}

		void traceSkeleton()
		{
			// NOTE: Not enough asserts!

			throwOnFailedAssertion(true);


			vector<string> inputs = {
				"test1",
				"test1_z",
				"test2",
				"test2_z",
				"test3",
				"test3_z",
				"test4_intersections",		// z-directional skeleton
				"test5",
				"test5_z",
				"test6_single_pixel_edge",		// Single pixel long edge
				"test7",
				"test7_z",
				"test8",
				"test8_z",
				"test9_loop",
				"test9_loop_z",
				"t1-head_bin"	// Real data
			};

			for (const auto& in : inputs)
			{
				cout << "INPUT = " << in << endl;

				//string in = "test4";
				string skeleout = "./skeleton/skele_" + in;


				Image<uint8_t> skele, orig;
				raw::read(skele, "../test_input_data/" + in);
				raw::read(orig, "../test_input_data/" + in);

				//Image<uint8_t> skele0, orig0;
				//raw::read(skele0, in);
				//raw::read(orig0, in);
				//Image<uint8_t> skele(256, 256, 51);
				//Image<uint8_t> orig(skele.dimensions());
				//crop(skele0, skele, Vec3c(0, 0, 0));
				//crop(orig0, orig, Vec3c(0, 0, 0));


				lineSkeleton(skele);
				raw::writed(skele, skeleout);


				

				//// Single-threaded
				Network net0;
				raw::read(skele, skeleout);
				Image<uint8_t> vis0;
				setValue(vis0, skele);
				internals::classifyForTracing(skele);
				raw::writed(skele, skeleout + "_classified");
				size_t counter = 0;
				internals::traceLineSkeleton(skele, &orig, Vec3d(), false, 1.0, 1.0, net0, counter, skele.depth());
				//cout << net0 << endl;
				cout << net0.vertices.size() << " vertices" << endl;
				cout << net0.edges.size() << " edges" << endl;
				
				draw<uint8_t>(vis0, net0, true, 1, 170, true, true, 128);

				raw::writed(vis0, skeleout + "_vis0");


				for (int blockCount = 1; blockCount < 20; blockCount++)
				//int blockCount = 3;
				//int blockCount = 7;
				//int blockCount = 13;
				{
					// Multithreaded, with measurements
					Network net1;
					raw::read(skele, skeleout);
					Image<uint8_t> vis1;
					setValue(vis1, skele);
					traceLineSkeleton(skele, &orig, false, 1.0, 1.0, net1, blockCount);
					//cout << net1 << endl;
					cout << net1.vertices.size() << " vertices" << endl;
					cout << net1.edges.size() << " edges" << endl;
					draw<uint8_t>(vis1, net1, true, 1, 170, true, true, 128);

					raw::writed(vis1, skeleout + "_vis1_" + toString(blockCount) + "_blocks");

					assertNetworks(net0, net1, "Single- and multi-threaded, block count = " + toString(blockCount));

					// This is not a good test as the visualization may differ because of floating point inaccuracy in node positions.
					//testAssert(equals(vis0, vis1), "Difference between visualizations, block count = " + toString(blockCount));


					if (skele.depth() / blockCount <= 20)
						break;
				}
				//for (const auto& v : net1.incompleteEdges[0].points)
				//{
				//	vis(v) = 50;
				//}

				// Multithreaded, no measurements
				//Network net2;
				//raw::read(skele, skeleout);
				//setValue(vis, skele);
				//traceLineSkeleton(skele, net2);
				//cout << net2 << endl;
				//draw<uint8_t>(vis, net2, 2, 170, 128);
				//raw::writed(vis, skeleout + "_vis2");



				//assertNetworks(net0, net2, "Single- and multi-threaded without measurements");
			}
		}

		void traceSkeletonRealData()
		{
			Image<uint8_t> skele, orig;
			raw::read(skele, "../test_input_data/real_skele_200x200x200.raw");
			Network net;
			traceLineSkeleton(skele, false, 1.0, 1.0, net);

			size_t badCount = 0;
			for (const auto& edge : net.edges)
			{
				//float32_t distance = edge.properties.adjustedDistance();
				float32_t distance = (net.vertices[edge.verts[0]] - net.vertices[edge.verts[1]]).norm();
				if (!testAssert(NumberUtils<float32_t>::greaterThanOrEqual(edge.properties.length, distance, 0.0001f), "adjusted distance and length"))
				{
					cout << edge.properties.length << " < " << distance << endl;
					badCount++;
				}
			}

			if (badCount > 0)
				cout << "Found " << badCount << " edges where lengths are not consistent." << endl;
		}

		void lineLength()
		{
			Image<uint8_t> img(10, 20);
			
			vector<Vec3f> points;
			points.push_back(Vec3f(7, 5, 0));
			points.push_back(Vec3f(7, 4, 0));
			points.push_back(Vec3f(6, 3, 0));
			points.push_back(Vec3f(5, 2, 0));
			points.push_back(Vec3f(4, 2, 0));
			points.push_back(Vec3f(3, 2, 0));
			points.push_back(Vec3f(2, 3, 0));
			points.push_back(Vec3f(2, 4, 0));
			points.push_back(Vec3f(1, 5, 0));
			points.push_back(Vec3f(1, 6, 0));
			points.push_back(Vec3f(1, 7, 0));
			points.push_back(Vec3f(1, 8, 0));
			points.push_back(Vec3f(2, 9, 0));
			points.push_back(Vec3f(3, 10, 0));
			points.push_back(Vec3f(4, 11, 0));
			points.push_back(Vec3f(4, 12, 0));
			points.push_back(Vec3f(5, 13, 0));
			points.push_back(Vec3f(6, 14, 0));
			points.push_back(Vec3f(7, 15, 0));
			points.push_back(Vec3f(7, 16, 0));
			points.push_back(Vec3f(7, 17, 0));
			points.push_back(Vec3f(6, 18, 0));
			points.push_back(Vec3f(5, 18, 0));
			points.push_back(Vec3f(4, 18, 0));
			points.push_back(Vec3f(3, 17, 0));
			points.push_back(Vec3f(2, 16, 0));

			for (const auto& p : points)
				img((int)p.x, (int)p.y, (int)p.z) = 255;
			
			raw::writed(img, "./skeleton/line_length/line1");

			float32_t L0straight = (points.back() - points.front()).norm();
			float32_t L0 = internals::lineLength(points,  0);
			float32_t L1straight = (points.back() - points.front()).norm();
			float32_t L1 = internals::lineLength(points, 1.0);
			cout << "Line length (no smoothing) = " << L0 << endl;
			cout << "Straight line length (no smoothing) = " << L0straight << endl;
			cout << "Line length (with smoothing) = " << L1 << endl;
			cout << "Straight line length (no smoothing) = " << L1straight << endl;

			testAssert(L0straight == L1straight, "Straight distance changes.");
			testAssert(L0straight <= L0, "straight line length and smoothed line length 0");
			testAssert(L1straight <= L1, "straight line length and smoothed line length 1");
			

			int S = 10;
			Image<uint8_t> img2(S * img.width(), S * img.height());
			scale(img, img2, false, InterpolationMode::Nearest, BoundaryCondition::Zero);

			for (const auto& p : points)
				img2((int)(S * p.x), (int)(S * p.y), (int)(S * p.z)) = 128;

			raw::writed(img2, "./skeleton/line_length/line1_result_vis");
		}

	}
}
