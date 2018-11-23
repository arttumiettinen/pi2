
#include "traceskeleton.h"

#include "io/raw.h"
#include "hybridskeleton.h"
#include "filters.h"
#include "transform.h"

using namespace math;

namespace itl2
{

	namespace internals
	{
		/**
		Smooths dimension dim of points in the given list and anchors them so that they stay near their original locations.
		The distance between original location of a point and new location won't be more than delta.
		*/
		void smoothAndAnchor(vector<Vec3f>& points, double sigma, float32_t delta)
		{
			
			// Smooth
			vector<Vec3f> smoothed(points.size());

			for (size_t dim = 0; dim < 3; dim++)
			{
				Image<float32_t> c(points.size());
				for (size_t n = 0; n < points.size(); n++)
					c(n) = points[n][dim];

				gaussFilter(c, Vec3d(sigma, 0.0, 0.0), Nearest);

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
					orig = orig + dir / l * delta;
				}
			}
			


			/*
			for (size_t dim = 0; dim < 3; dim++)
			{
				Image<float32_t> c(points.size());
				for (size_t n = 0; n < points.size(); n++)
					c(n) = points[n][dim];

				gaussFilter(c, Vec3d(sigma, 0.0, 0.0), Nearest);

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

		float32_t lineLength(vector<Vec3f>& points, double sigma)
		{
			if (sigma > 0)
			{
				constexpr float32_t delta = 0.5f;
				smoothAndAnchor(points, sigma, delta);
			}

			float32_t l = 0;
			for (size_t n = 1; n < points.size(); n++)
			{
				l += (points[n] - points[n - 1]).norm();
			}

			return l;
		}

		float32_t straightLineLength(const Vec3f& a, const Vec3f& b, double sigma)
		{
			vector<Vec3f> points;
			points.push_back(a);
			points.push_back(b);
			return lineLength(points, sigma);
		}



		/**
		Finds out if there are any neighbouring points in v1 and v2.
		*/
		bool isNeighbour(const vector<Vec3c>& v1, const vector<Vec3c>& v2)
		{
			for (size_t n = 0; n < v1.size(); n++)
			{
				for (size_t m = 0; m < v2.size(); m++)
				{
					if ((v2[m] - v1[n]).abs().max() <= 1)
					{
						return true;
					}
				}
			}

			return false;
		}

		/**
		Finds intersecting incomplete vertex regions and combines their indices in the forest.
		At output the forest will contain all indices of vertices list, and indices of intersecting vertices belong to the same root.
		*/
		void findIntersectingRegions(const vector<IncompleteVertex>& vertices, DisjointSetForest<size_t>& forest)
		{
			// Add all indices to the forest
			cout << "Build forest..." << endl;
			for (size_t n = 0; n < vertices.size(); n++)
				forest.add_element(n);

			// Calculate bounding boxes for incomplete vertex sets
			cout << "Determine bounding boxes..." << endl;
			vector<Box<coord_t> > boundingBoxes;
			boundingBoxes.reserve(vertices.size());
			Vec3c m = Vec3c(numeric_limits<coord_t>::max(), numeric_limits<coord_t>::max(), numeric_limits<coord_t>::max());
			Vec3c M = Vec3c(numeric_limits<coord_t>::lowest(), numeric_limits<coord_t>::lowest(), numeric_limits<coord_t>::lowest());
			for (size_t n = 0; n < vertices.size(); n++)
			{
				Box<coord_t> b = Box<coord_t>::boundingBox(vertices[n].points);
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
#pragma omp parallel for if(!omp_in_parallel())
			for (coord_t i = 0; i < grid.pixelCount(); i++)
			{
				auto& list = grid(i);

				for (coord_t ni = 0; ni < (coord_t)list.size(); ni++)
				{
					size_t n = list[ni];
					const IncompleteVertex& v1 = vertices[n];
					const Box<coord_t>& b1 = boundingBoxes[n];

					for (size_t mi = ni + 1; mi < list.size(); mi++)
					{
						size_t m = list[mi];
						const IncompleteVertex& v2 = vertices[m];
						const Box<coord_t>& b2 = boundingBoxes[m];

						if (b1.overlaps(b2) && isNeighbour(v1.points, v2.points))
						{
							// v1 and v2 are neighbours, so the intersection regions must be combined
#pragma omp critical(forestInsert)
							{
								forest.union_sets(n, m);
							}
						}
					}
				}

				showThreadProgress(counter, grid.pixelCount());
			}
		}

		/**
		Combines all incomplete vertices that represent the same vertex.
		After the call net.incompleteVertices is empty.
		The network may contain unnecessary vertices.
		TODO: This does not adjust edge length values when vertices move. For well-behaved high-resolution skeleton the error should be minimal
		*/
		void combineIncompleteVertices(Network& net)
		{
			// Find neighbouring incomplete vertices
			// The forest contains indices into net.incompleteVertices list.
			DisjointSetForest<size_t> forest;
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
				}

				showThreadProgress(counter, net.incompleteVertices.size());
			}

			// Calculate true position of each incomplete vertex.
			cout << "Calculate new positions for root vertices..." << endl;
			counter = 0;
			#pragma omp parallel for if(!omp_in_parallel()) // Can be parallelized because incompleteVertices[n].vertexIndex is different for each n.
			for (coord_t n = 0; n < (coord_t)net.incompleteVertices.size(); n++)
			{
				if (net.incompleteVertices[n].points.size() > 0)
				{
					math::Vec3f center(0, 0, 0);
					for (const Vec3c& p : net.incompleteVertices[n].points)
						center += math::Vec3f(p);
					center /= (float32_t)net.incompleteVertices[n].points.size();

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

			cout << "Create map of vertex indexes that have changed..." << endl;
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

				const auto& it = vertexIndexMap.find(edge.verts[0]);
				if (it != vertexIndexMap.end())
					edge.verts[0] = it->second;

				const auto& it2 = vertexIndexMap.find(edge.verts[1]);
				if (it2 != vertexIndexMap.end())
					edge.verts[1] = it2->second;

				showThreadProgress(counter, net.edges.size());
			}


			// There are no incomplete vertices anymore
			cout << "Free memory..." << endl;
			net.incompleteVertices.clear();
			net.incompleteVertices.shrink_to_fit();
		}

		//size_t findNextNonRedundantVertex(Network& net, const vector<uint8_t>& degree, size_t first, size_t second)
		//{
		//	size_t previous = first;
		//	size_t curr = second;
		//	while (degree[curr] == 2)
		//	{

		//	}
		//}

		///**
		//Removes vertices that have exactly 2 branches.
		//*/
		//void removeUnncesessaryVertices(Network& net)
		//{
		//	// Build branch count per vertex (== degree of each vertex)
		//	vector<uint8_t> degree(net.vertices.size(), 0);

		//	for (size_t n = 0; n < net.edges.size(); n++)
		//	{
		//		degree[net.edges[n].verts[0]]++;
		//		degree[net.edges[n].verts[1]]++;
		//	}

		//	// Find edge that starts a redundant set of edges
		//	for (size_t n = 0; n < net.edges.size(); n++)
		//	{
		//		Vec2c edge = net.edges[n].verts;
		//		if (degree[edge[0]] == 2 && degree[edge[1]] != 2)
		//		{
		//			// edge[0] vertex is redundant
		//			// connect edge[1] to first nonredunant vertex found
		//			size_t target = findNextNonRedundantVertex(net, degree, edge[1], edge[0]);

		//		}
		//		else if (degree[edge[0]] != 2 && degree[edge[1]] == 2)
		//		{
		//			// edge[1] vertex is redundant

		//		}

		//	}

		//}

		void combineTracedBlocks(const vector<Network>& subNets, Network& net, bool freeSubnets)
		{
			// Combine overlapping intersection regions
			cout << "Combine calculation blocks..." << endl;

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
					net2.edges[m].verts += math::Vec2c(vertexIndexShift, vertexIndexShift);
				for (size_t m = 0; m < net2.incompleteVertices.size(); m++)
					net2.incompleteVertices[m].vertexIndex += vertexIndexShift;

				net.vertices.insert(net.vertices.end(), net2.vertices.begin(), net2.vertices.end());
				net.edges.insert(net.edges.end(), net2.edges.begin(), net2.edges.end());
				net.incompleteVertices.insert(net.incompleteVertices.end(), net2.incompleteVertices.begin(), net2.incompleteVertices.end());

				if (freeSubnets)
				{
					net2.edges.clear();
					net2.edges.shrink_to_fit();
					net2.vertices.clear();
					net2.vertices.shrink_to_fit();
					net2.incompleteVertices.clear();
					net2.incompleteVertices.shrink_to_fit();
				}

				showThreadProgress(counter, subNets.size());
			}

			// Process incomplete vertices
			cout << "Combine incomplete vertices on processing block boundaries..." << endl;
			internals::combineIncompleteVertices(net);

			cout << "Remove straight-through nodes..." << endl;
			net.disconnectStraightThroughNodes(true);

			cout << "Remove isolated nodes..." << endl;
			net.removeIsolatedNodes(true);

			//cout << "Found " << net.vertices.size() << " vertices" << endl;
			//cout << "and   " << net.edges.size() << " edges." << endl;

			// Remove vertices with less than 2 branches
			// TODO: Test first, then continue
			//removeUnnecessaryVertices(net);
		}
	}
	
	
	namespace tests
	{
		//void classifySkeleton()
		//{
		//	// NOTE: No asserts!
		
		//	Image<uint8_t> skele;
		//	raw::readd(skele, "./skeleton/cavities/cavity_hybrid_skeleton");

		//	classifySkeleton(skele, true, true, true);

		//	raw::writed(skele, "./skeleton/cavities/classified_cavity_hybrid_skeleton");
		//}

		void traceSkeleton()
		{
			// NOTE: No asserts!

			string in = "./test4";
			string skeleout = "./skele_test4";
			

			Image<uint8_t> skele, orig;
			raw::readd(skele, in);
			raw::readd(orig, in);

			hybridSkeleton(skele);
			raw::writed(skele, skeleout);

			Network net;
			traceLineSkeleton(skele, orig, net);
			cout << net;
		}

		void lineLength()
		{
			// NOTE: No asserts!

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

			float32_t L0 = internals::lineLength(points, 0);
			float32_t L1 = internals::lineLength(points, 1.0);
			cout << "Length (no smoothing) = " << L0 << endl;
			cout << "Length (with smoothing) = " << L1 << endl;
			cout << "Straight line length = " << (*points.begin() - *points.rbegin()).norm() << endl;
			cout << "Adjusted straight line length = " << internals::straightLineLength(*points.begin(), *points.rbegin()) << endl;

			int S = 10;
			Image<uint8_t> img2(S * img.width(), S * img.height());
			scale(img, img2, NearestNeighbourInterpolator<uint8_t, uint8_t>(Zero));

			for (const auto& p : points)
				img2((int)(S * p.x), (int)(S * p.y), (int)(S * p.z)) = 128;

			raw::writed(img2, "./skeleton/line_length/line1_result_vis");
		}
	}
}