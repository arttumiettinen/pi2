
#include "network.h"
#include "io/raw.h"
#include "traceskeleton.h"

#include <fstream>
#include <iostream>

#include <random>
#include <functional>

using namespace std;


namespace itl2
{
	const Vec3sc EdgeMeasurements::INVALID_POINT = Vec3sc(numeric_limits<int>::max(), numeric_limits<int>::max(), numeric_limits<int>::max());

	const Edge Edge::INVALID = Edge(numeric_limits<size_t>::max(), numeric_limits<size_t>::max(), EdgeMeasurements());

	const Vec3f Network::INVALID_VERTEX = Vec3f(numeric_limits<float32_t>::signaling_NaN(), numeric_limits<float32_t>::signaling_NaN(), numeric_limits<float32_t>::signaling_NaN());

	void Network::degree(vector<size_t>& deg, bool reportProgress) const
	{
		deg.resize(vertices.size());

		size_t counter = 0;
		#pragma omp parallel if(!omp_in_parallel() && edges.size() + incompleteEdges.size() > PARALLELIZATION_THRESHOLD)
		{
			vector<size_t> degPrivate;
			degPrivate.resize(vertices.size());

			#pragma omp for nowait
			for (coord_t n = 0; n < (coord_t)edges.size(); n++)
			{
				degPrivate[edges[n].verts[0]]++;
				degPrivate[edges[n].verts[1]]++;

				//showThreadProgress(counter, edges.size() + incompleteEdges.size(), reportProgress);
			}

			#pragma omp for nowait
			for (coord_t n = 0; n < (coord_t)incompleteEdges.size(); n++)
			{
				if(incompleteEdges[n].verts[0] >= 0)
					degPrivate[incompleteEdges[n].verts[0]]++;

				if (incompleteEdges[n].verts[1] >= 0)
					degPrivate[incompleteEdges[n].verts[1]]++;

				//showThreadProgress(counter, edges.size() + incompleteEdges.size(), reportProgress);
			}

			#pragma omp critical(degree_reduction)
			{
				for (coord_t n = 0; n < (coord_t)deg.size(); n++)
				{
					deg[n] += degPrivate[n];
				}
			}
		}
	}

	void Network::inEdges(size_t n, vector<size_t>& edg) const
	{
		#pragma omp parallel for if(!omp_in_parallel() && edges.size() > PARALLELIZATION_THRESHOLD)
		for (coord_t i = 0; i < (coord_t)edges.size(); i++)
		{
			if (edges[i].verts[1] == n)
			{
				#pragma omp critical(inEdges)
				{
					edg.push_back(i);
				}
			}
		}
	}

	void Network::outEdges(size_t n, vector<size_t>& edg) const
	{
		#pragma omp parallel for if(!omp_in_parallel() && edges.size() > PARALLELIZATION_THRESHOLD)
		for (coord_t i = 0; i < (coord_t)edges.size(); i++)
		{
			if (edges[i].verts[0] == n)
			{
				#pragma omp critical(outEdges)
				{
					edg.push_back(i);
				}
			}
		}
	}

	void Network::inOutEdges(size_t n, vector<size_t>& inEdg, vector<size_t>& outEdg) const
	{
		#pragma omp parallel for if(!omp_in_parallel() && edges.size() > PARALLELIZATION_THRESHOLD)
		for (coord_t i = 0; i < (coord_t)edges.size(); i++)
		{
			if (edges[i].verts[0] == n)
			{
				#pragma omp critical(outEdges)
				{
					outEdg.push_back(i);
				}
			}

			if (edges[i].verts[1] == n)
			{
				#pragma omp critical(inEdges)
				{
					inEdg.push_back(i);
				}
			}
		}
	}

	void Network::neighbours(size_t n, vector<size_t>& edgeIndices) const
	{
		#pragma omp parallel for if(!omp_in_parallel() && edges.size() > PARALLELIZATION_THRESHOLD)
		for (coord_t i = 0; i < (coord_t)edges.size(); i++)
		{
			if (edges[i].verts[0] == n || edges[i].verts[1] == n)
			{
				#pragma omp critical(outEdges)
				{
					edgeIndices.push_back(i);
				}
			}
		}
	}

	//const IncompleteVertex* Network::getIncompleteVertex(size_t vertexIndex) const
	//{
	//	for (const IncompleteVertex& v : incompleteVertices)
	//	{
	//		if (v.vertexIndex == vertexIndex)
	//			return &v;
	//	}

	//	return nullptr;
	//}

	//bool Network::isIncompleteVertex(size_t vertexIndex) const
	//{
	//	return getIncompleteVertex(vertexIndex) != nullptr;
	//}

	bool Network::isIncompleteVertex(size_t vertexIndex) const
	{
		for (const IncompleteVertex& v : incompleteVertices)
		{
			if (v.vertexIndex == vertexIndex)
				return true;
		}

		return false;
	}

	//void Network::disconnectInternal(size_t n)
	//{
	//	vector<size_t> neighbourIndices;
	//	neighbourIndices.reserve(2);

	//	neighbours(n, neighbourIndices);
	//	
	//	for (size_t n = 0; n < neighbourIndices.size(); n++)
	//	{
	//		Edge& edgen = edges[n];

	//		for (size_t m = n + 1; m < neighbourIndices.size(); m++)
	//		{
	//			Edge& edgem = edges[m];
	//			float32_t L = edgen.length + edgem.length;
	//			float32_t A = 0.5f * (edgem.area + edgen.area);

	//		}
	//	}

	//	// Find parents and children
	//	/*
	//	vector<size_t> inEdg, outEdg;
	//	inEdg.reserve(2);
	//	outEdg.reserve(2);
	//	
	//	inOutEdges(n, inEdg, outEdg);

	//	// Add edges from parents to children
	//	for (size_t parentIndex : inEdg)
	//	{
	//	    Edge& parentEdge = edges[parentIndex];
	//		for (size_t childIndex : outEdg)
	//		{
 //   			Edge& childEdge = edges[childIndex];
 //   			
	//			float32_t L = parentEdge.length + childEdge.length;
	//			float32_t A = 0.5f * (parentEdge.area + childEdge.area);
	//			size_t startvi = parentEdge.verts[0];
	//			size_t endvi = childEdge.verts[1];
	//			
	//			// Replace child edge by the new edge
	//			//edges.push_back(Edge(startvi, endvi, L, A));
	//			childEdge.verts = Vec2c(startvi, endvi);
	//			childEdge.length = L;
	//			childEdge.area = A;
	//		}

	//		// Invalidate parent edge
	//		parentEdge = Edge::INVALID;
	//	}
	//	*/

	//	// Invalidate all edges where n is either start or end
	//	//#pragma omp parallel for if(!omp_in_parallel() && edges.size() > PARALLELIZATION_THRESHOLD)
	//	//for (coord_t i = 0; i < (coord_t)edges.size(); i++)
	//	//{
	//	//	if (edges[i].verts[0] == n || edges[i].verts[1] == n)
	//	//	{
	//	//		edges[i] = Edge::INVALID;
	//	//	}
	//	//}
	//}

	//void Network::disconnect(size_t n, bool reportProgress)
	//{
	//	disconnectInternal(n);
	//	clean(reportProgress);
	//}

	void Network::clean(bool reportProgress)
	{
	    vector<Edge> newEdges;
		for (size_t n = 0; n < edges.size(); n++)
		{
			if (edges[n] == Edge::INVALID)
			{
				//edges.erase(edges.begin() + n);
				//n--;
			}
			else
			{
			    newEdges.push_back(edges[n]);
			}

			showProgress(n, edges.size(), reportProgress);
		}
		
		edges = newEdges;
	}

	size_t find(size_t ind, const vector<tuple<size_t, EdgeMeasurements> >& list)
	{
		for (size_t n = 0; n < list.size(); n++)
		{
			if (get<0>(list[n]) == ind)
				return n;
		}

		throw ITLException("Impossible situation: Unable to find correct index from neighbour list.");
	}

	void remove(size_t ind, vector<tuple<size_t, EdgeMeasurements> >& list)
	{
		size_t n = find(ind, list);
		list.erase(list.begin() + n);
	}

	void Network::markStraightThroughNodes(bool reportProgress)
	{
		// For each node, list of [neighbour node index, edge properties]
		vector<vector<tuple<size_t, EdgeMeasurements> > > net(vertices.size());

		// Populate net list from edges
		cout << "Convert to neighbour list format..." << endl;
		for (size_t n = 0; n < edges.size(); n++)
		{
			const Edge& edge = edges[n];
			net[edge.verts[0]].push_back(make_tuple(edge.verts[1], edge.properties));
			net[edge.verts[1]].push_back(make_tuple(edge.verts[0], edge.properties));
			showProgress(n, edges.size(), reportProgress);
		}
		edges.clear();

		// Find nodes with exactly 2 neighbours
		cout << "Erase connections to nodes with 2 neighbours..." << endl;
		for (size_t n = 0; n < net.size(); n++)
		{
			if (net[n].size() == 2)
			{
				// Connect a to b by replacing
				// - n by b in a:s neighbour list and
				// - n by a in b:s neighbour list.
				size_t a = get<0>(net[n][0]);
				size_t b = get<0>(net[n][1]);

				if (a != n && b != n) // Don't adjust connections if any of the neighbours is this node, i.e. there is a loop.
				{
					EdgeMeasurements p1 = get<1>(net[n][0]);
					EdgeMeasurements p2 = get<1>(net[n][1]);

					EdgeMeasurements combined;
					combined.pointCount = p1.pointCount + p2.pointCount;
					combined.length = p1.length + p2.length;
					combined.area = (p1.length * p1.area + p2.length * p2.area) / (p1.length + p2.length);


					// Combine edgePoints lists in correct order
					combined.edgePoints.reserve(p1.edgePoints.size() + p2.edgePoints.size());

					Vec3f startVertex = vertices[a];
					Vec3f endVertex = vertices[b];

					vector<vector<Vec3sc>> pointLists;
					pointLists.push_back(p1.edgePoints);
					pointLists.push_back(p2.edgePoints);
					Vec3f currentEnd = startVertex;
					while (pointLists.size() > 0)
					{
						coord_t minIndex = -1;
						bool isStart = true;
						float32_t minDist = numeric_limits<float32_t>::infinity();

						for (size_t n = 0; n < pointLists.size(); n++)
						{
							const vector<Vec3sc>& list = pointLists[n];
							if (list.size() > 0)
							{
								// Test if one end of the list is closer to the current end point than any previous end.
								float32_t d = (currentEnd - Vec3f(list.front())).norm<float32_t>();
								if (d < minDist)
								{
									minDist = d;
									minIndex = n;
									isStart = true;
								}

								d = (currentEnd - Vec3f(list.back())).norm<float32_t>();
								if (d < minDist)
								{
									minDist = d;
									minIndex = n;
									isStart = false;
								}
							}
						}

						// No non-empty list available.
						if (minIndex < 0)
							break;

						vector<Vec3sc>& list = pointLists[minIndex];
						if (!isStart)
							reverse(list.begin(), list.end());

						combined.edgePoints.insert(combined.edgePoints.end(), list.begin(), list.end());

						pointLists.erase(pointLists.begin() + minIndex);

						currentEnd = Vec3f(combined.edgePoints.back());
					}

					//combined.edgePoints.insert(combined.edgePoints.end(), p1.edgePoints.begin(), p1.edgePoints.end());
					//// Test if order of p2.edgePoints must be reversed
					//if (p2.edgePoints.size() > 0)
					//{
					//	float32_t d1 = (combined.edgePoints.back() - p2.edgePoints.back()).norm<float32_t>();
					//	float32_t d2 = (combined.edgePoints.back() - p2.edgePoints.front()).norm<float32_t>();
					//	float32_t d3 = (combined.edgePoints.front() - p2.edgePoints.back()).norm<float32_t>();
					//	float32_t d4 = (combined.edgePoints.front() - p2.edgePoints.front()).norm<float32_t>();

					//	if(d1 <= d2 && d1 <= d3 && d1 <= d4)
					//		combined.edgePoints.insert(combined.edgePoints.end(), p2.edgePoints.rbegin(), p2.edgePoints.rend());
					//	else if (d2 <= d1 && d2 <= d3 && d2 <= d4)
					//		combined.edgePoints.insert(combined.edgePoints.end(), p2.edgePoints.begin(), p2.edgePoints.end());
					//	else if (d3 <= d1 && d3 <= d2 && d3 <= d4)
					//		combined.edgePoints.insert(combined.edgePoints.begin(), p2.edgePoints.begin(), p2.edgePoints.end());
					//	else
					//		combined.edgePoints.insert(combined.edgePoints.begin(), p2.edgePoints.rbegin(), p2.edgePoints.rend());
					//}


					// Make sure that we don't propagate nan areas
					if (std::isnan(combined.area))
					{
						if (!std::isnan(p1.area))
							combined.area = p1.area;
						else
							combined.area = p2.area;
					}

					size_t nInd;
					nInd = find(n, net[a]);
					net[a][nInd] = make_tuple(b, combined);
					nInd = find(n, net[b]);
					net[b][nInd] = make_tuple(a, combined);

					net[n].clear();

					// The node should be removed (and it does not have any connections), so mark it as invalid.
					vertices[n] = INVALID_VERTEX;
				}
			}
			showProgress(n, net.size(), reportProgress);
		}

		// Convert network back to edge list format
		cout << "Convert to edge list format..." << endl;
		for (size_t n = 0; n < net.size(); n++)
		{
			for (size_t m = 0; m < net[n].size(); m++)
			{
				size_t b = get<0>(net[n][m]);
				EdgeMeasurements p = get<1>(net[n][m]);
				edges.push_back(Edge(n, b, p));

				// Remove the other item corresponding to this edge.
				remove(n, net[b]);
			}
			showProgress(n, net.size(), reportProgress);
		}
	}

	void Network::removeStraightThroughNodes(bool reportProgress)
	{
		// Mark nodes
		markStraightThroughNodes(reportProgress);

		// Remove nodes that were marked as invalid.
		removeInvalidNodes(reportProgress);
	}

	void Network::markIsolatedNodes(bool reportProgress)
	{
		// Remove nodes with no connections
		cout << "Calculate degree..." << endl;
		vector<size_t> deg;
		degree(deg, reportProgress);

		cout << "Remove nodes..." << endl;
		for (size_t n = 0; n < vertices.size(); n++)
		{
			if (deg[n] <= 0)
			{
				vertices[n] = INVALID_VERTEX;
			}

			showProgress(n, vertices.size(), reportProgress);
		}
	}

	void Network::removeInvalidNodes(bool reportProgress)
	{
		// Remove nodes with no connections
		cout << "Calculate degree..." << endl;
		vector<size_t> deg;
		degree(deg, reportProgress);

		cout << "Remove nodes..." << endl;
		vector<Vec3f> newNodes;
		vector<size_t> oldToNewIndex(vertices.size());
		for (size_t n = 0; n < vertices.size(); n++)
		{
			// remove if deg[n] <= 0 && (isnan(vertices[n].x) || isnan(vertices[n].y) || isnan(vertices[n].z))
			//if (deg[n] > 0)
			if(!(deg[n] <= 0 && (std::isnan(vertices[n].x) || std::isnan(vertices[n].y) || std::isnan(vertices[n].z))))
			{
				newNodes.push_back(vertices[n]);
				oldToNewIndex[n] = newNodes.size() - 1;
			}

			showProgress(n, vertices.size(), reportProgress);
		}

		size_t removedCount = vertices.size() - newNodes.size();

		vertices = newNodes;

		cout << "Removed " << removedCount << " nodes." << endl;

		cout << "Update edges..." << endl;
		size_t counter = 0;
		#pragma omp parallel for if(!omp_in_parallel() && edges.size() > PARALLELIZATION_THRESHOLD)
		for (coord_t n = 0; n < (coord_t)edges.size(); n++)
		{
			edges[n].verts[0] = oldToNewIndex[edges[n].verts[0]];
			edges[n].verts[1] = oldToNewIndex[edges[n].verts[1]];
			showThreadProgress(counter, edges.size(), reportProgress);
		}

		if (incompleteEdges.size() > 0)
		{
			cout << "Update incomplete edges..." << endl;
			counter = 0;
			#pragma omp parallel for if(!omp_in_parallel() && incompleteEdges.size() > PARALLELIZATION_THRESHOLD)
			for (coord_t n = 0; n < (coord_t)incompleteEdges.size(); n++)
			{
				Vec2c& verts = incompleteEdges[n].verts;
				if(verts[0] >= 0)
					verts[0] = oldToNewIndex[verts[0]];
				if (verts[1] >= 0)
					verts[1] = oldToNewIndex[verts[1]];

				showThreadProgress(counter, incompleteEdges.size(), reportProgress);
			}
		}

		if (incompleteVertices.size() > 0)
		{
			cout << "Update incomplete vertices..." << endl;
			counter = 0;
			#pragma omp parallel for if(!omp_in_parallel() && incompleteVertices.size() > PARALLELIZATION_THRESHOLD)
			for (coord_t n = 0; n < (coord_t)incompleteVertices.size(); n++)
			{
				incompleteVertices[n].vertexIndex = oldToNewIndex[incompleteVertices[n].vertexIndex];

				showThreadProgress(counter, incompleteVertices.size(), reportProgress);
			}
		}
	}

	void Network::removeEdges(const vector<coord_t>& edgeIndices, bool removeStraightThrough, bool removeIsolated, bool reportProgress)
	{
		for (coord_t n = 0; n < (coord_t)edgeIndices.size(); n++)
		{
			coord_t i = edgeIndices[n];
			if (i < 0 || i >= (coord_t)edges.size())
				throw ITLException("Edge index out of range: " + toString(i));
		}


		vector<size_t> deg;
		degree(deg, reportProgress);

		for (coord_t n = 0; n < (coord_t)edgeIndices.size(); n++)
		{
			coord_t i = edgeIndices[n];

			coord_t v1i = edges[i].verts.x;
			coord_t v2i = edges[i].verts.y;
			edges[i] = Edge::INVALID;

			// Mark empty node for removal
			if (deg[v1i] <= 1)
				vertices[v1i] = INVALID_VERTEX;

			if (deg[v2i] <= 1)
				vertices[v2i] = INVALID_VERTEX;
		}

		// Remove edges marked as INVALID
		cout << "Remove unnecessary edges..." << endl;
		clean(reportProgress);

		// Combine x--x--x sequences into x--x (x=node, --=edge) (mark unnecessary nodes for removal)
		if (removeStraightThrough)
			markStraightThroughNodes(reportProgress);

		// Remove (mark for removal) nodes with degree = 0
		if (removeIsolated)
			markIsolatedNodes(reportProgress);

		// Remove nodes that were marked as invalid.
		removeInvalidNodes(reportProgress);
	}

	void Network::prune(float32_t maxLength, bool removeStraightThrough, bool removeIsolated, bool reportProgress)
	{
		cout << "Calculate degree..." << endl;
		vector<size_t> deg;
		degree(deg, reportProgress);

        cout << "Find edges and nodes to prune..." << endl;
        #pragma omp parallel for if(!omp_in_parallel() && edges.size() >= PARALLELIZATION_THRESHOLD)
		for (coord_t n = 0; n < (coord_t)edges.size(); n++)
		{
			coord_t v1i = edges[n].verts.x;
			coord_t v2i = edges[n].verts.y;

			if (deg[v1i] <= 1 || deg[v2i] <= 1)
			{
				// This edge connects to a leaf node (node without other edges connected to it).
				if (NumberUtils<float32_t>::lessThanOrEqual(edges[n].properties.length, maxLength))
				{
					// The edge is short, too.
					// Remove it.
					edges[n] = Edge::INVALID;

					// Mark empty node for removal
					if (deg[v1i] <= 1)
						vertices[v1i] = INVALID_VERTEX;

					if (deg[v2i] <= 1)
						vertices[v2i] = INVALID_VERTEX;
				}
			}
		}

		// Remove edges marked as INVALID
		cout << "Remove unnecessary edges..." << endl;
		clean(reportProgress);

		// Combine x--x--x sequences into x--x (x=node, --=edge) (mark unnecessary nodes for removal)
		if(removeStraightThrough)
			markStraightThroughNodes(reportProgress);

		// Remove (mark for removal) nodes with degree = 0
		if(removeIsolated)
			markIsolatedNodes(reportProgress);

		// Remove nodes that were marked as invalid.
		removeInvalidNodes(reportProgress);
	}

	void Network::toImage(Image<float32_t>& vertices, Image<uint64_t>& edges, Image<float32_t>* pEdgeMeasurements, Image<int32_t>* pEdgePoints) const
	{
		vertices.ensureSize(3, this->vertices.size());
		edges.ensureSize(2, this->edges.size());

		#pragma omp parallel for if(!omp_in_parallel())
		for (coord_t n = 0; n < (coord_t)this->vertices.size(); n++)
		{
			const Vec3f& p = this->vertices[n];
			vertices(0, n) = p.x;
			vertices(1, n) = p.y;
			vertices(2, n) = p.z;
		}

		#pragma omp parallel for if(!omp_in_parallel())
		for (coord_t n = 0; n < (coord_t)this->edges.size(); n++)
		{
			const Edge& p = this->edges[n];
			edges(0, n) = p.verts.x;
			edges(1, n) = p.verts.y;
		}

		if (pEdgeMeasurements)
		{
			//pEdgeMeasurements->ensureSize(4+3+3, this->edges.size());
			pEdgeMeasurements->ensureSize(3, this->edges.size());
			#pragma omp parallel for if(!omp_in_parallel())
			for (coord_t n = 0; n < (coord_t)this->edges.size(); n++)
			{
				const Edge& p = this->edges[n];
				(*pEdgeMeasurements)(0, n) = p.properties.pointCount;
				(*pEdgeMeasurements)(1, n) = p.properties.length;
				(*pEdgeMeasurements)(2, n) = p.properties.area;
				//(*pEdgeMeasurements)(3, n) = pixelRound<float32_t>(p.properties.pointOnEdge[0]);
				//(*pEdgeMeasurements)(4, n) = pixelRound<float32_t>(p.properties.pointOnEdge[1]);
				//(*pEdgeMeasurements)(5, n) = pixelRound<float32_t>(p.properties.pointOnEdge[2]);
			}
		}

		if (pEdgePoints)
		{
			// Count total number of points
			coord_t totalPointCount = 0;
			for (size_t n = 0; n < this->edges.size(); n++)
				totalPointCount += this->edges[n].properties.edgePoints.size();

			pEdgePoints->ensureSize(this->edges.size() + 3 * totalPointCount);

			// Write point counts
			for (coord_t n = 0; n < (coord_t)this->edges.size(); n++)
				(*pEdgePoints)(n) = (int32_t)this->edges[n].properties.edgePoints.size();

			// Write points
			coord_t pos = (coord_t)this->edges.size();
			for (coord_t n = 0; n < (coord_t)this->edges.size(); n++)
			{
				for (const auto& x : this->edges[n].properties.edgePoints)
				{
					(*pEdgePoints)(pos) = x[0];
					pos++;
					(*pEdgePoints)(pos) = x[1];
					pos++;
					(*pEdgePoints)(pos) = x[2];
					pos++;
				}
			}
		}
	}

	void Network::fromImage(const Image<float32_t>& verts, const Image<uint64_t>& edg, const Image<float32_t>* pEdgeMeasurements, const Image<int32_t>* pEdgePoints)
	{
		vertices.clear();
		edges.clear();
		incompleteVertices.clear();

		vertices.resize(verts.height());
		#pragma omp parallel for if(!omp_in_parallel())
		for (coord_t n = 0; n < verts.height(); n++)
		{
			vertices[n].x = verts(0, n);
			vertices[n].y = verts(1, n);
			vertices[n].z = verts(2, n);
		}

		edges.resize(edg.height());
		#pragma omp parallel for if(!omp_in_parallel())
		for (coord_t n = 0; n < edg.height(); n++)
		{
			size_t v1 = edg(0, n);
			size_t v2 = edg(1, n);

			if (v1 >= vertices.size() || v2 >= vertices.size())
				throw ITLException(string("Loading invalid network from image: vertex index for edge ") + toString(n) + " is greater than or equal to vertex count.");

			edges[n].verts.x = v1;
			edges[n].verts.y = v2;
		}

		if (pEdgeMeasurements)
		{
			pEdgeMeasurements->checkSize(Vec3c(3, edg.height(), 1));
			#pragma omp parallel for if(!omp_in_parallel())
			for (coord_t n = 0; n < edg.height(); n++)
			{
				edges[n].properties.pointCount = (*pEdgeMeasurements)(0, n);
				edges[n].properties.length = (*pEdgeMeasurements)(1, n);
				edges[n].properties.area = (*pEdgeMeasurements)(2, n);
				//edges[n].properties.pointOnEdge = Vec3sc(pixelRound<int>((*pEdgeMeasurements)(3, n)), pixelRound<int>((*pEdgeMeasurements)(4, n)), pixelRound<int>((*pEdgeMeasurements)(5, n)));
			}
		}

		if (pEdgePoints)
		{
			if ((size_t)pEdgePoints->pixelCount() < edges.size())
				throw ITLException("Edge points image is invalid. It contain less pixels than count of edges.");

			// Read point count for each edge
			coord_t totalCount = 0;
			for (size_t n = 0; n < edges.size(); n++)
			{
				auto& points = edges[n].properties.edgePoints;
				int32_t pointCount = (*pEdgePoints)(n);
				points.insert(points.end(), pointCount, Vec3sc());
				totalCount += pointCount;
			}

			size_t trueSize = edges.size() + 3 * totalCount;
			if (pEdgePoints->pixelCount() != trueSize)
				throw ITLException(string("Edge points image does not contain correct number of pixels. Required pixel count is ") + itl2::toString(trueSize) + ", but the image contains " + itl2::toString(pEdgePoints->pixelCount()) + " pixels.");

			// Read points
			coord_t pos = edges.size();
			for (size_t n = 0; n < edges.size(); n++)
			{
				auto& points = edges[n].properties.edgePoints;
				for (size_t m = 0; m < points.size(); m++)
				{
					int32_t x = (*pEdgePoints)(pos);
					pos++;
					int32_t y = (*pEdgePoints)(pos);
					pos++;
					int32_t z = (*pEdgePoints)(pos);
					pos++;
					points[m] = Vec3sc(x, y, z);
				}
			}
		}
	}

	void writeEdge(ofstream& out, const Edge& e)
	{
		const Vec2c& v = e.verts;
		out.write((const char*)&v.components[0], sizeof(Vec2c));
		out.write((const char*)&e.properties.pointCount, sizeof(float32_t));
		out.write((const char*)&e.properties.length, sizeof(float32_t));
		out.write((const char*)&e.properties.area, sizeof(float32_t));
		
		//out.write((const char*)&e.properties.pointOnEdge, sizeof(Vec3sc));
		
		size_t edgePointCount = e.properties.edgePoints.size();
		out.write((const char*)&edgePointCount, sizeof(size_t));
		for(size_t n = 0; n < edgePointCount; n++)
			out.write((const char*)&(e.properties.edgePoints[n]), sizeof(Vec3sc));
	}

	void Network::write(const string& filename, bool append) const
	{
		//FileLock lock(filename);
		
		createFoldersFor(filename);

		ios::openmode mode;
		if (append)
			mode = ios_base::out | ios_base::app | ios_base::binary;
		else
			mode = ios_base::out | ios_base::trunc | ios_base::binary;
		
		ofstream out(filename, mode);

		// If the file does not exist, create it.
		if (!out.is_open())
			out.open(filename, ios_base::out | ios_base::trunc | ios_base::binary);

		if (!out)
			throw ITLException(string("Unable to open file ") + filename);

		// File structure
		// vert count (size_t)
		//	vert x
		//	vert y
		//	vert z
		//	...
		// edge count
		//	vertex index 1
		//	vertex index 2
		//	point count
		//	length
		//  (other properties)
		// ....
		// incomplete vertex count
		//	vertex index
		//	point count
		//		point x
		//		point y
		//		point z
		//		...
		//	...
		// incomplete edge count
		//	vertex index 1
		//	vertex index 2
		//	point count
		//	length
		//  (other properties)
		//  ....
		//  point count
		//	    point x
		//      point y
		//      point z
		//      ...
		//  ...

		size_t count = vertices.size();
		out.write((char*)&count, sizeof(size_t));

		for (size_t n = 0; n < vertices.size(); n++)
		{
			const Vec3f& v = vertices[n];
			out.write((const char*)&v.components[0], sizeof(Vec3f));
		}

		count = edges.size();
		out.write((char*)&count, sizeof(size_t));

		for (size_t n = 0; n < edges.size(); n++)
		{
			writeEdge(out, edges[n]);
		}

		count = incompleteVertices.size();
		out.write((char*)&count, sizeof(size_t));

		for (size_t n = 0; n < incompleteVertices.size(); n++)
		{
			const IncompleteVertex& v = incompleteVertices[n];
			out.write((char*)&v.vertexIndex, sizeof(size_t));
			count = v.points.size();
			out.write((char*)&count, sizeof(size_t));
			for (size_t m = 0; m < v.points.size(); m++)
			{
				const Vec3sc& p = v.points[m];
				out.write((const char*)&p.components[0], sizeof(Vec3sc));
			}
		}

		count = incompleteEdges.size();
		out.write((char*)&count, sizeof(size_t));

		for (size_t n = 0; n < incompleteEdges.size(); n++)
		{
			const IncompleteEdge& v = incompleteEdges[n];
			
			writeEdge(out, v);

			count = v.points.size();
			out.write((char*)&count, sizeof(size_t));
			for (size_t m = 0; m < v.points.size(); m++)
			{
				const Vec3sc& p = v.points[m];
				out.write((const char*)&p.components[0], sizeof(Vec3sc));
			}
		}
	}

	Edge readEdge(ifstream& in)
	{
		Vec2c v;
		EdgeMeasurements p;
		in.read((char*)&v.components[0], sizeof(Vec2c));
		in.read((char*)&p.pointCount, sizeof(float32_t));
		in.read((char*)&p.length, sizeof(float32_t));
		in.read((char*)&p.area, sizeof(float32_t));
		
		//in.read((char*)&p.pointOnEdge, sizeof(Vec3sc));
		size_t edgePointCount = 0;
		in.read((char*)&edgePointCount, sizeof(size_t));
		p.edgePoints.reserve(edgePointCount);
		for (size_t n = 0; n < edgePointCount; n++)
		{
			Vec3sc x;
			in.read((char*)&x, sizeof(Vec3sc));
			p.edgePoints.push_back(x);
		}
		

		//if (v.x < 0 || v.y < 0 || (size_t)v.x >= net.vertices.size() || (size_t)v.y >= net.vertices.size())
		//	throw ITLException(string("Loading invalid network: vertex index for edge ") + toString(n) + " is negative or greater than or equal to vertex count.");

		return Edge(v.x, v.y, p);
	}

	void Network::read(const string& filename, vector<Network>& nets)
	{
		//FileLock lock(filename);

		ifstream in(filename, ios_base::in | ios_base::binary);

		if (!in)
			throw ITLException(string("Unable to open file ") + filename);

		do
		{
			Network net;

			size_t count;

			in.read((char*)&count, sizeof(size_t));
			if (!in)
			{
			    break;
			}
			net.vertices.reserve(count);
			for (size_t n = 0; n < count; n++)
			{
				Vec3f v;
				in.read((char*)&v.components[0], sizeof(Vec3f));
				net.vertices.push_back(v);
			}

			in.read((char*)&count, sizeof(size_t));
			net.edges.reserve(count);
			for (size_t n = 0; n < count; n++)
			{
				net.edges.push_back(readEdge(in));
			}

			in.read((char*)&count, sizeof(size_t));
			net.incompleteVertices.reserve(count);
			for (size_t n = 0; n < count; n++)
			{
				IncompleteVertex v;
				in.read((char*)&v.vertexIndex, sizeof(size_t));

				size_t vcount;
				in.read((char*)&vcount, sizeof(size_t));
				v.points.reserve(vcount);
				for (size_t m = 0; m < vcount; m++)
				{
					Vec3sc p;
					in.read((char*)&p.components[0], sizeof(Vec3sc));
					v.points.push_back(p);
				}

				net.incompleteVertices.push_back(v);
			}

			in.read((char*)&count, sizeof(size_t));
			net.incompleteEdges.reserve(count);
			for (size_t n = 0; n < count; n++)
			{
				Edge e = readEdge(in);

				IncompleteEdge v;
				v.verts = e.verts;
				v.properties = e.properties;
				
				size_t vcount;
				in.read((char*)&vcount, sizeof(size_t));
				v.points.reserve(vcount);
				for (size_t m = 0; m < vcount; m++)
				{
					Vec3sc p;
					in.read((char*)&p.components[0], sizeof(Vec3sc));
					v.points.push_back(p);
				}

				net.incompleteEdges.push_back(v);
			}

			nets.push_back(net);
		} while (!in.eof());
		
	}

	namespace tests
	{
		void disconnections()
		{
			Network net;


			//   |-1  0   1   2   3   4  5
			//----------------------------
			//   |
			// -1|    2
			//   |      \                 /\
			//  0|    1 - 3 - 4 - 5 - 6  7  |
			//   |      /                 \/   
			//  1|    0  
			//   |
			
			net.vertices.push_back(Vec3f(0, 1, 0));
			net.vertices.push_back(Vec3f(0, 0, 0));
			net.vertices.push_back(Vec3f(0, -1, 0));
			net.vertices.push_back(Vec3f(1, 0, 0));
			net.vertices.push_back(Vec3f(2, 0, 0));
			net.vertices.push_back(Vec3f(3, 0, 0));
			net.vertices.push_back(Vec3f(4, 0, 0));

			net.vertices.push_back(Vec3f(5, 0, 0));

			net.edges.push_back(Edge(0, 3, EdgeMeasurements(2, 2, 1)));
			net.edges.push_back(Edge(1, 3, EdgeMeasurements(2, 1, 1)));
			net.edges.push_back(Edge(2, 3, EdgeMeasurements(2, 2, 1)));

			net.edges.push_back(Edge(3, 4, EdgeMeasurements(2, 1, 1)));
			net.edges.push_back(Edge(5, 4, EdgeMeasurements(2, 1, 1)));
			net.edges.push_back(Edge(5, 6, EdgeMeasurements(2, 1, 1)));

			net.edges.push_back(Edge(7, 7, EdgeMeasurements(2, 1, 1)));

			vector<size_t> deg;
			net.degree(deg, true);
			if (testAssert(deg.size() == 8, "Wrong degree array size."))
			{
				testAssert(deg[0] == 1, "deg 0");
				testAssert(deg[1] == 1, "deg 1");
				testAssert(deg[2] == 1, "deg 2");
				testAssert(deg[3] == 4, "deg 3");
				testAssert(deg[4] == 2, "deg 4");
				testAssert(deg[5] == 2, "deg 5");
				testAssert(deg[6] == 1, "deg 6");
				testAssert(deg[7] == 2, "deg 7");
			}

			net.removeStraightThroughNodes(true);
			deg.clear();
			net.degree(deg, true);
			if (testAssert(deg.size() == 6, "Wrong degree array size."))
			{
				testAssert(deg[0] == 1, "deg 0");
				testAssert(deg[1] == 1, "deg 1");
				testAssert(deg[2] == 1, "deg 2");
				testAssert(deg[3] == 4, "deg 3");
				testAssert(deg[4] == 1, "deg 4 (old vertex 6)");
				testAssert(deg[5] == 2, "deg 5 (old vertex 7)");
				//testAssert(deg[4] == 0, "deg 4");
				//testAssert(deg[5] == 0, "deg 5");
				//testAssert(deg[6] == 1, "deg 6");
				//testAssert(deg[7] == 0, "deg 6");
			}

			//net.removeIsolatedNodes(true);
			//testAssert(net.vertices.size() == 5, "size");
			//deg.clear();
			//net.degree(deg, true);
			//testAssert(deg[0] == 1, "deg 0");
			//testAssert(deg[1] == 1, "deg 1");
			//testAssert(deg[2] == 1, "deg 2");
			//testAssert(deg[3] == 4, "deg 3");
			//testAssert(deg[4] == 1, "deg 4");
		}

		void disconnectStraightThroughPerformance()
		{
			// NOTE: No asserts!

			Network net;
			const size_t NODE_COUNT = 10000000;
			net.vertices.reserve(NODE_COUNT);
			for (size_t n = 0; n < NODE_COUNT; n++)
			{
				net.vertices.push_back(Vec3f(0, 0, 0));
			}

			std::mt19937 gen;
			std::uniform_int_distribution<size_t> dist(0, NODE_COUNT - 1);
			auto dice = std::bind(dist, gen);

			const size_t EDGE_COUNT = 10000000;
			net.edges.reserve(EDGE_COUNT);
			for (size_t n = 0; n < EDGE_COUNT; n++)
			{
				size_t start = dice();
				size_t end = dice();
				if(start != end)
					net.edges.push_back(Edge(start, end, EdgeMeasurements(2, 1, 1)));
			}

			net.removeStraightThroughNodes(true);

//			net.removeIsolatedNodes(true);
		}


		void networkio()
		{
			// NOTE: No asserts!

			Network orig;
			orig.vertices.push_back(Vec3f(1.1f, 2.2f, 3.3f));
			orig.vertices.push_back(Vec3f(4.4f, 5.5f, 6.6f));

			orig.edges.push_back(Edge(0, 1, EdgeMeasurements(7, 3.14f, 5)));
			orig.edges.push_back(Edge(1, 0, EdgeMeasurements(7, 2*3.14f, 6)));
			orig.edges.push_back(Edge(1, 1, EdgeMeasurements(7, 3*3.14f, 7)));

			IncompleteVertex iv;
			iv.vertexIndex = 1;
			iv.points.push_back(Vec3sc(1, 2, 3));
			iv.points.push_back(Vec3sc(4, 5, 6));
			orig.incompleteVertices.push_back(iv);

			IncompleteEdge ie;
			ie.verts = Vec2c(1, -1);
			ie.points.push_back(Vec3sc(10, 20, 30));
			ie.points.push_back(Vec3sc(40, 50, 60));
			orig.incompleteEdges.push_back(ie);

			orig.write("./network/original.dat", false);
			orig.write("./network/original.dat", true);

			cout << "Original network" << endl;
			cout << orig << endl;

			vector<Network> read;
			Network::read("./network/original.dat", read);

			for(size_t n = 0; n < read.size(); n++)
			{
				cout << "Read from file (index = " << n << ")" << endl;
				cout << read[n] << endl;

				testAssert(orig.vertices == read[n].vertices, "verts");
				testAssert(orig.edges == read[n].edges, "edges");
				testAssert(orig.incompleteVertices == read[n].incompleteVertices, "incomplete verts");
				testAssert(orig.incompleteEdges == read[n].incompleteEdges, "incomplete edges");
			}
		}

		void pruning()
		{
			Network net;

			// o3
			// |     o4
			// |     |     o5
			// |     |     |
			// o0 -- o1 -- o2

			net.vertices.push_back(Vec3f(0, 0, 0));
			net.vertices.push_back(Vec3f(1, 0, 0));
			net.vertices.push_back(Vec3f(2, 0, 0));

			net.vertices.push_back(Vec3f(0, 3, 0));
			net.vertices.push_back(Vec3f(1, 2, 0));
			net.vertices.push_back(Vec3f(2, 1, 0));

			net.edges.push_back(Edge(0, 1, EdgeMeasurements(0, 1, 0)));
			net.edges.push_back(Edge(1, 2, EdgeMeasurements(0, 1, 0)));

			net.edges.push_back(Edge(0, 3, EdgeMeasurements(0, 3, 0)));
			net.edges.push_back(Edge(1, 4, EdgeMeasurements(0, 2, 0)));
			net.edges.push_back(Edge(2, 5, EdgeMeasurements(0, 1, 0)));

			net.prune(0.5, false, false);

			if (testAssert(net.edges.size() == 5, "edge count"))
			{
				testAssert(net.edges[0].verts == Vec2c(0, 1), "edge 0");
				testAssert(net.edges[1].verts == Vec2c(1, 2), "edge 1");
				testAssert(net.edges[2].verts == Vec2c(0, 3), "edge 2");
				testAssert(net.edges[3].verts == Vec2c(1, 4), "edge 3");
				testAssert(net.edges[4].verts == Vec2c(2, 5), "edge 4");
			}

			net.prune(1.0f, false, false);

			if (testAssert(net.edges.size() == 4, "edge count"))
			{
				testAssert(net.edges[0].verts == Vec2c(0, 1), "edge 0");
				testAssert(net.edges[1].verts == Vec2c(1, 2), "edge 1");
				testAssert(net.edges[2].verts == Vec2c(0, 3), "edge 2");
				testAssert(net.edges[3].verts == Vec2c(1, 4), "edge 3");
				//testAssert(net.edges[4].verts == Vec2c(2, 5), "edge 4");
			}

			net.prune(1.0f, false, false);

			// o2
			// |     o3
			// |     |
			// |     |
			// o0 -- o1

			if (testAssert(net.edges.size() == 3, "edge count"))
			{
				testAssert(net.edges[0].verts == Vec2c(0, 1), "edge 0");
				testAssert(net.edges[1].verts == Vec2c(0, 2), "edge 1");
				testAssert(net.edges[2].verts == Vec2c(1, 3), "edge 2");
			}
		}
	}
}
