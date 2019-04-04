
#include "network.h"
#include "io/raw.h"
#include "traceskeleton.h"

#include <fstream>
#include <iostream>

#include <random>
#include <functional>

using namespace std;
using namespace math;

namespace itl2
{


	const Edge Edge::INVALID = Edge(numeric_limits<size_t>::max(), numeric_limits<size_t>::max(), EdgeMeasurements());

	void Network::degree(vector<size_t>& deg, bool reportProgress) const
	{
		deg.resize(vertices.size());

		size_t counter = 0;
		#pragma omp parallel if(!omp_in_parallel() && vertices.size() > PARALLELIZATION_THRESHOLD)
		{
			vector<size_t> degPrivate;
			degPrivate.resize(vertices.size());

			#pragma omp for nowait
			for (coord_t n = 0; n < (coord_t)edges.size(); n++)
			{
				degPrivate[edges[n].verts[0]]++;
				degPrivate[edges[n].verts[1]]++;

				showThreadProgress(counter, edges.size(), reportProgress);
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
		for (size_t n = 0; n < edges.size(); n++)
		{
			if (edges[n] == Edge::INVALID)
			{
				edges.erase(edges.begin() + n);
				n--;
			}

			showProgress(n, edges.size(), reportProgress);
		}
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

	void Network::disconnectStraightThroughNodes(bool reportProgress)
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

				EdgeMeasurements p1 = get<1>(net[n][0]);
				EdgeMeasurements p2 = get<1>(net[n][1]);

				EdgeMeasurements combined;
				combined.pointCount = p1.pointCount + p2.pointCount;
				combined.length = p1.length + p2.length;
				combined.area = (p1.length * p1.area + p2.length * p2.area) / (p1.length + p2.length);

				combined.distance = (vertices[a] - vertices[b]).norm();
				//combined.adjustedDistance = internals::straightLineLength(vertices[a], vertices[b]);
				if ((p1.adjustedStart - vertices[a]).norm() < (p1.adjustedEnd - vertices[a]).norm())
					combined.adjustedStart = p1.adjustedStart;
				else
					combined.adjustedStart = p1.adjustedEnd;
				if ((p2.adjustedStart - vertices[b]).norm() < (p2.adjustedEnd - vertices[b]).norm())
					combined.adjustedEnd = p2.adjustedStart;
				else
					combined.adjustedEnd = p2.adjustedEnd;


				/*float32_t L1 = get<1>(net[n][0]);
				float32_t L2 = get<1>(net[n][1]);
				float32_t A1 = get<2>(net[n][0]);
				float32_t A2 = get<2>(net[n][1]);
				float32_t L = L1 + L2;
				float32_t A = (L1 * A1 + L2 * A2) / L;*/


				size_t nInd;
				nInd = find(n, net[a]);
				net[a][nInd] = make_tuple(b, combined);
				nInd = find(n, net[b]);
				net[b][nInd] = make_tuple(a, combined);

				net[n].clear();
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
				//float32_t L = get<1>(net[n][m]);
				//float32_t A = get<2>(net[n][m]);
				//edges.push_back(Edge(n, b, L, A));
				edges.push_back(Edge(n, b, p));

				// Remove the other item corresponding to this edge.
				remove(n, net[b]);
			}
			showProgress(n, net.size(), reportProgress);
		}
	}

	void Network::removeIsolatedNodes(bool reportProgress)
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
			if (deg[n] > 0)
			{
				newNodes.push_back(vertices[n]);
				oldToNewIndex[n] = newNodes.size() - 1;
			}

			showProgress(n, vertices.size(), reportProgress);
		}

		vertices = newNodes;

		cout << "Update edges..." << endl;
		size_t counter = 0;
		#pragma omp parallel for if(!omp_in_parallel() && edges.size() > PARALLELIZATION_THRESHOLD)
		for (coord_t n = 0; n < (coord_t)edges.size(); n++)
		{
			edges[n].verts[0] = oldToNewIndex[edges[n].verts[0]];
			edges[n].verts[1] = oldToNewIndex[edges[n].verts[1]];
			showThreadProgress(counter, edges.size(), reportProgress);
		}
	}

	void Network::toImage(Image<float32_t>& vertices, Image<size_t>& edges, Image<float32_t>* pEdgeMeasurements) const
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
			pEdgeMeasurements->ensureSize(4+3+3, this->edges.size());
			#pragma omp parallel for if(!omp_in_parallel())
			for (coord_t n = 0; n < (coord_t)this->edges.size(); n++)
			{
				const Edge& p = this->edges[n];
				(*pEdgeMeasurements)(0, n) = p.properties.pointCount;
				(*pEdgeMeasurements)(1, n) = p.properties.length;
				(*pEdgeMeasurements)(2, n) = p.properties.area;
				(*pEdgeMeasurements)(3, n) = p.properties.distance;
				(*pEdgeMeasurements)(4, n) = p.properties.adjustedStart.x;
				(*pEdgeMeasurements)(5, n) = p.properties.adjustedStart.y;
				(*pEdgeMeasurements)(6, n) = p.properties.adjustedStart.z;
				(*pEdgeMeasurements)(7, n) = p.properties.adjustedEnd.x;
				(*pEdgeMeasurements)(8, n) = p.properties.adjustedEnd.y;
				(*pEdgeMeasurements)(9, n) = p.properties.adjustedEnd.z;
			}
		}
	}

	void Network::fromImage(const Image<float32_t>& verts, const Image<size_t>& edg, const Image<float32_t>* pEdgeMeasurements)
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
			edges[n].verts.x = edg(0, n);
			edges[n].verts.y = edg(1, n);
		}

		if (pEdgeMeasurements)
		{
			#pragma omp parallel for if(!omp_in_parallel())
			for (coord_t n = 0; n < edg.height(); n++)
			{
				edges[n].properties.pointCount = (*pEdgeMeasurements)(0, n);
				edges[n].properties.length = (*pEdgeMeasurements)(1, n);
				edges[n].properties.area = (*pEdgeMeasurements)(2, n);
				edges[n].properties.distance = (*pEdgeMeasurements)(3, n);
				edges[n].properties.adjustedStart = Vec3f((*pEdgeMeasurements)(4, n), (*pEdgeMeasurements)(5, n), (*pEdgeMeasurements)(6, n));
				edges[n].properties.adjustedEnd = Vec3f((*pEdgeMeasurements)(7, n), (*pEdgeMeasurements)(8, n), (*pEdgeMeasurements)(9, n));
			}
		}
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
			const Vec2c& v = edges[n].verts;
			out.write((const char*)&v.components[0], sizeof(Vec2c));
			out.write((const char*)&edges[n].properties.pointCount, sizeof(float32_t));
			out.write((const char*)&edges[n].properties.length, sizeof(float32_t));
			out.write((const char*)&edges[n].properties.area, sizeof(float32_t));
			out.write((const char*)&edges[n].properties.distance, sizeof(float32_t));
			//out.write((const char*)&edges[n].properties.adjustedDistance, sizeof(float32_t));
			out.write((const char*)&edges[n].properties.adjustedStart.x, sizeof(float32_t));
			out.write((const char*)&edges[n].properties.adjustedStart.y, sizeof(float32_t));
			out.write((const char*)&edges[n].properties.adjustedStart.z, sizeof(float32_t));
			out.write((const char*)&edges[n].properties.adjustedEnd.x, sizeof(float32_t));
			out.write((const char*)&edges[n].properties.adjustedEnd.y, sizeof(float32_t));
			out.write((const char*)&edges[n].properties.adjustedEnd.z, sizeof(float32_t));
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
				Vec2c v;
				EdgeMeasurements p;
				in.read((char*)&v.components[0], sizeof(Vec2c));
				in.read((char*)&p.pointCount, sizeof(float32_t));
				in.read((char*)&p.length, sizeof(float32_t));
				in.read((char*)&p.area, sizeof(float32_t));
				in.read((char*)&p.distance, sizeof(float32_t));
				in.read((char*)&p.adjustedStart.x, sizeof(float32_t));
				in.read((char*)&p.adjustedStart.y, sizeof(float32_t));
				in.read((char*)&p.adjustedStart.z, sizeof(float32_t));
				in.read((char*)&p.adjustedEnd.x, sizeof(float32_t));
				in.read((char*)&p.adjustedEnd.y, sizeof(float32_t));
				in.read((char*)&p.adjustedEnd.z, sizeof(float32_t));
				net.edges.push_back(Edge(v.x, v.y, p));
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

			nets.push_back(net);
		} while (!in.eof());
		
	}

	namespace tests
	{
		void disconnections()
		{
			Network net;
			
			net.vertices.push_back(Vec3f(0, 1, 0));
			net.vertices.push_back(Vec3f(0, 0, 0));
			net.vertices.push_back(Vec3f(0, -1, 0));
			net.vertices.push_back(Vec3f(1, 0, 0));
			net.vertices.push_back(Vec3f(2, 0, 0));
			net.vertices.push_back(Vec3f(3, 0, 0));
			net.vertices.push_back(Vec3f(4, 0, 0));

			net.vertices.push_back(Vec3f(5, 0, 0));

			net.edges.push_back(Edge(0, 3, EdgeMeasurements(2, 2, 1, 0, net.vertices[0], net.vertices[3])));
			net.edges.push_back(Edge(1, 3, EdgeMeasurements(2, 1, 1, 0, net.vertices[1], net.vertices[3])));
			net.edges.push_back(Edge(2, 3, EdgeMeasurements(2, 2, 1, 0, net.vertices[2], net.vertices[3])));

			net.edges.push_back(Edge(3, 4, EdgeMeasurements(2, 1, 1, 0, net.vertices[3], net.vertices[4])));
			net.edges.push_back(Edge(5, 4, EdgeMeasurements(2, 1, 1, 0, net.vertices[5], net.vertices[4])));
			net.edges.push_back(Edge(5, 6, EdgeMeasurements(2, 1, 1, 0, net.vertices[5], net.vertices[6])));

			net.edges.push_back(Edge(7, 7, EdgeMeasurements(2, 1, 1, 0, net.vertices[7], net.vertices[7])));

			vector<size_t> deg;
			net.degree(deg, true);
			testAssert(deg[0] == 1, "deg 0");
			testAssert(deg[1] == 1, "deg 1");
			testAssert(deg[2] == 1, "deg 2");
			testAssert(deg[3] == 4, "deg 3");
			testAssert(deg[4] == 2, "deg 4");
			testAssert(deg[5] == 2, "deg 5");
			testAssert(deg[6] == 1, "deg 6");
			testAssert(deg[7] == 2, "deg 7");

			net.disconnectStraightThroughNodes(true);
			deg.clear();
			net.degree(deg, true);
			testAssert(deg[0] == 1, "deg 0");
			testAssert(deg[1] == 1, "deg 1");
			testAssert(deg[2] == 1, "deg 2");
			testAssert(deg[3] == 4, "deg 3");
			testAssert(deg[4] == 0, "deg 4");
			testAssert(deg[5] == 0, "deg 5");
			testAssert(deg[6] == 1, "deg 6");
			testAssert(deg[7] == 0, "deg 6");

			net.removeIsolatedNodes(true);
			testAssert(net.vertices.size() == 5, "size");
			deg.clear();
			net.degree(deg, true);
			testAssert(deg[0] == 1, "deg 0");
			testAssert(deg[1] == 1, "deg 1");
			testAssert(deg[2] == 1, "deg 2");
			testAssert(deg[3] == 4, "deg 3");
			testAssert(deg[4] == 1, "deg 4");
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
					net.edges.push_back(Edge(start, end, EdgeMeasurements(2, 1, 1, 0, net.vertices[start], net.vertices[end])));
			}

			net.disconnectStraightThroughNodes(true);

			net.removeIsolatedNodes(true);
		}


		void networkio()
		{
			// NOTE: No asserts!

			Network orig;
			orig.vertices.push_back(Vec3f(1.1f, 2.2f, 3.3f));
			orig.vertices.push_back(Vec3f(4.4f, 5.5f, 6.6f));

			orig.edges.push_back(Edge(0, 1, EdgeMeasurements(7, 3.14f, 5, 1, orig.vertices[0], orig.vertices[1])));
			orig.edges.push_back(Edge(1, 0, EdgeMeasurements(7, 2*3.14f, 6, 3, orig.vertices[1], orig.vertices[0])));
			orig.edges.push_back(Edge(1, 1, EdgeMeasurements(7, 3*3.14f, 7, 3, orig.vertices[1], orig.vertices[1])));

			IncompleteVertex iv;
			iv.vertexIndex = 1;
			iv.points.push_back(Vec3sc(1, 2, 3));
			iv.points.push_back(Vec3sc(4, 5, 6));
			orig.incompleteVertices.push_back(iv);

			orig.write("./network/original.dat", false);
			orig.write("./network/original.dat", true);

			cout << "Original" << endl;
			cout << orig << endl;

			vector<Network> read;
			Network::read("./network/original.dat", read);

			for(size_t n = 0; n < read.size(); n++)
			{
				cout << "Read " << n << endl;
				cout << read[n] << endl;
			}
		}
	}
}
