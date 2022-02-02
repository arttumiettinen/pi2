#pragma once

#include "math/vec2.h"
#include "math/vec3.h"
#include "math/numberutils.h"
#include "image.h"

#include <vector>


namespace itl2
{
	struct EdgeMeasurements
	{
	public:

		static const Vec3sc INVALID_POINT;

		/**
		Count of pixels in the image that belong to this edge.
		In intersection regions the pixels belonging to the edge might not be neighbours.
		*/
		float32_t pointCount;

		/**
		Length of the edge, calculated using anchored convolution approach.
		*/
		float32_t length;

		/**
		Average cross-sectional area, calculated as in CSA paper.
		*/
		float32_t area;

		/**
		Location of some pixel on the edge.
		This is used later in network drawing algorithms to fill the edge with specific color.
		*/
		//Vec3sc pointOnEdge;

		/**
		Locations of some points on the edge.
		This is needed for fillskeleton-type commands.
		*/
		std::vector<Vec3sc> edgePoints;


		template<typename out_t> out_t get(size_t index) const
		{
			if (index == 0)
				return pixelRound<out_t>(pointCount);
			if (index == 1)
				return pixelRound<out_t>(length);
			if (index == 2)
				return pixelRound<out_t>(area);
			//if (index == 3)
			//	return pixelRound<out_t>(pointOnEdge.x);
			//if (index == 4)
			//	return pixelRound<out_t>(pointOnEdge.y);
			//if (index == 5)
			//	return pixelRound<out_t>(pointOnEdge.z);
			throw ITLException("Invalid property index.");
		}

		EdgeMeasurements() :
			pointCount(0),
			length(0),
			area(0)
			//pointOnEdge(Vec3sc())
			//distance(0),
			//adjustedStart(0, 0, 0),
			//adjustedEnd(0, 0, 0)
			//adjustedDistance(0)
		{

		}

		EdgeMeasurements(float32_t pointCount, float32_t length, float32_t area/*, const Vec3sc& pointOnEdge, float32_t distance, float32_t adjustedDistance const Vec3f& adjustedStart, const Vec3f& adjustedEnd*/) :
			pointCount(pointCount),
			length(length),
			area(area)
			//pointOnEdge(pointOnEdge)
			//distance(distance),
			//adjustedStart(adjustedStart),
			//adjustedEnd(adjustedEnd)
			//adjustedDistance(adjustedDistance)
		{
		}

		friend std::ostream& operator<<(std::ostream& stream, const EdgeMeasurements& e)
		{
			//stream << "N = " << e.pointCount << ", L = " << e.length << ", A = " << e.area << ", d = " << e.distance << ", d'= " << e.adjustedDistance();
			stream << "N = " << e.pointCount << ", L = " << e.length << ", A = " << e.area << ", " << e.edgePoints.size() << " points";
				// << ", p = " << e.pointOnEdge;
			return stream;
		}

		bool operator==(const EdgeMeasurements& r) const
		{
			//return pointCount == r.pointCount && length == r.length && area == r.area && distance == r.distance && /*adjustedDistance == r.adjustedDistance*/ adjustedStart == r.adjustedStart && adjustedEnd == r.adjustedEnd;
			return pointCount == r.pointCount && length == r.length && area == r.area && edgePoints == r.edgePoints;
				// && pointOnEdge == pointOnEdge;
		}

		/**
		Tests if two edges equal. Allows numerical inaccuracy in floating point values.
		*/
		bool equals(const EdgeMeasurements& r) const
		{
			return NumberUtils<float32_t>::equals(pointCount, r.pointCount, 1e-3f) &&
				NumberUtils<float32_t>::equals(length, r.length, 1e-3f) &&
				(NumberUtils<float32_t>::equals(area, r.area, 1e-3f) || (std::isnan(area) && std::isnan(r.area))) && // NOTE: nan == nan is false
				edgePoints == r.edgePoints;
				//pointOnEdge.equals(r.pointOnEdge);
		}
	};


	class Edge
	{
	public:
		/**
		Indices of vertices that this edge connects.
		*/
		Vec2c verts;

		/**
		Measured properties of this edge.
		*/
		EdgeMeasurements properties;


		Edge() :
			verts(INVALID.verts),
			properties(INVALID.properties)
		{

		}

		Edge(size_t vertex1, size_t vertex2, EdgeMeasurements props) :
			verts(vertex1, vertex2),
			properties(props)
		{

		}

		friend std::ostream& operator<<(std::ostream& stream, const Edge& e)
		{
			stream << e.verts << ", " << e.properties;
			return stream;
		}

		static const Edge INVALID;

		bool operator==(const Edge& r) const
		{
			return verts == r.verts && properties == r.properties;
		}
	};

	class IncompleteVertex
	{
	public:
		/**
		The index of the vertex that corresponds to this edge intersection
		*/
		size_t vertexIndex;

		/**
		Points from which the vertex position was calculated.
		*/
		std::vector<Vec3sc> points;

		IncompleteVertex() :
			vertexIndex(0)
		{
		}

		IncompleteVertex(size_t vertexIndex, const std::vector<Vec3sc>& points) :
			vertexIndex(vertexIndex)
		{
			this->points.reserve(points.size());
			this->points.insert(this->points.begin(), points.begin(), points.end());
			this->points.shrink_to_fit();
		}

		friend std::ostream& operator<<(std::ostream& stream, const IncompleteVertex& e)
		{
			stream << e.vertexIndex << " (" << e.points.size() << " points)";
			return stream;
		}

		bool operator==(const IncompleteVertex& r) const
		{
			return vertexIndex == r.vertexIndex && points == r.points;
		}
	};

	class IncompleteEdge : public Edge
	{
	public:
		
		/**
		Points forming this edge.
		*/
		std::vector<Vec3sc> points;

		IncompleteEdge()
		{
		}

		IncompleteEdge(coord_t startVertexIndex, coord_t endVertexIndex, const std::vector<Vec3sc>& points) :
			Edge(startVertexIndex, endVertexIndex, EdgeMeasurements())
		{
			this->points.reserve(points.size());
			this->points.insert(this->points.begin(), points.begin(), points.end());
		}

		IncompleteEdge(coord_t startVertexIndex, coord_t endVertexIndex, const std::vector<Vec3sc>& points, EdgeMeasurements props) :
			Edge(startVertexIndex, endVertexIndex, props)
		{
			this->points.reserve(points.size());
			this->points.insert(this->points.begin(), points.begin(), points.end());
			this->points.shrink_to_fit();
		}

		friend std::ostream& operator<<(std::ostream& stream, const IncompleteEdge& e)
		{
			stream << (const Edge&)e << " (" << e.points.size() << " points)";
			return stream;
		}

		bool operator==(const IncompleteEdge& r) const
		{
			return (Edge::operator == (r)) && (points == r.points);
		}
	};

	class Network
	{
	private:

		/**
		Disconnects node n but does not call clean().
		The network is not in clean state before clean() is called.
		*/
		//void disconnectInternal(size_t n);

		/**
		Sets straight-through nodes to INVALID_VERTEX.
		*/
		void markStraightThroughNodes(bool reportProgress = true);

		/**
		Sets isolated nodes to INVALID_VERTEX.
		*/
		void markIsolatedNodes(bool reportProgress = true);

	public:

		static const Vec3f INVALID_VERTEX;


		/**
		Positions of vertices.
		*/
		std::vector<Vec3f> vertices;

		/**
		Edges between vertices.
		*/
		std::vector<Edge> edges;

		/**
		Incomplete vertices
		*/
		std::vector<IncompleteVertex> incompleteVertices;

		/**
		Incomplete edges
		*/
		std::vector<IncompleteEdge> incompleteEdges;


		/**
		Calculate degree of each node and place the values to given list.
		Degree of a node is number of edges connected to it.
		Loops that begin and end in the same node are counted as two edges (i.e. two connections to the node).
		*/
		void degree(std::vector<size_t>& deg, bool reportProgress) const;

		/**
		Finds parent nodes of node n and adds indices of edges from parent nodes to node n to the given list.
		*/
		void inEdges(size_t n, std::vector<size_t>& edg) const;

		/**
		Finds child nodes of node n and adds indices of edges from n to child nodes to the given list.
		*/
		void outEdges(size_t n, std::vector<size_t>& edg) const;

		/**
		Finds indices of edges where this node is target and indices of edges where this node is source.
		*/
		void inOutEdges(size_t n, std::vector<size_t>& inEdg, std::vector<size_t>& outEdg) const;

		/**
		Finds neighbours of node n, and stores indices of edges leading to neighbours to the given list.
		*/
		void neighbours(size_t n, std::vector<size_t>& edgeIndices) const;

		///**
		//Return incomplete vertex having specified index.
		//*/
		//const IncompleteVertex* getIncompleteVertex(size_t vertexIndex) const;

		//IncompleteVertex* getIncompleteVertex(size_t vertexIndex)
		//{
		//	return const_cast<IncompleteVertex*>(const_cast<const Network*>(this)->getIncompleteVertex(vertexIndex));
		//}

		/**
		Tests if the vertex that has the given index is incomplete.
		*/
		bool isIncompleteVertex(size_t vertexIndex) const;

		/**
		Invalidates all connections where node n is source or target, and adds new connections from parent nodes of n to child nodes of n.
		*/
		//void disconnect(size_t n, bool reportProgress);

		/**
		Remove edges marked as invalid.
		*/
		void clean(bool reportProgress);

		/**
		Removes connections to nodes with exactly 2 edges ("straight-through nodes").
		Retains edges from a node to itself.
		*/
		void removeStraightThroughNodes(bool reportProgress = true);

		/**
		Removes all nodes that have no connections.
		*/
		//void removeIsolatedNodes(bool reportProgress = true);

		/**
		Remove nodes that have no connections and whose coordinates are nan.
		*/
		void removeInvalidNodes(bool reportProgress = true);

		/**
		Prunes the network.
		Removes all edges that end in a node with no other branches connected to it (degree = 1) and whose length is less than maxLength.
		@param removeStraightThrough If set to true, all straight-through nodes (nodes with degree = 2) are removed after pruning. This operation might change the network even if no edges are pruned.
		@param removeIsolated If set to true, all isolated nodes (nodes with degree = 0) will be removed after pruning. This operation might change the network even if no edges are pruned.
		*/
		void prune(float32_t maxLength, bool removeStraightThrough, bool removeIsolated, bool reportProgress = true);

		/**
		Convert network to a pair of images.
		@param net Network to convert.
		@param vertices Image where vertex coordinates will be set. The size of the image is set to Nx3 where N is the number of vertices.
		@param edges Image where vertex indices corresponding to each edge will be set. The size of the image is set to Mx2 where M is the number of edges.
		@param pEdgeMeasurements Pointer to image that stores (pointCount, length, cross-sectional area) of each edge. Set to null if this information is not required.
		@param pEdgePoints Pointer to image that stores locations of one or more points that are on the edge, for each edge. The data is stored such that
		the first edgeCount pixels of the image store count of points for each edge. The remaining pixels store (x, y, z) coordinates of each point and each edge
		sequentially. The final format for two edges that have 1 and 2 points is therefore
		1
		2
		x11
		y11
		z11
		x21
		y21
		z21
		x22
		y22
		z22
		Set this parameter to nullptr if this information is not required.
		*/
		void toImage(Image<float32_t>& vertices, Image<uint64_t>& edges, Image<float32_t>* pEdgeMeasurements, Image<int32_t>* pEdgePoints) const;

		/**
		Replaces this network by one read from given vertex, edge and measurements images.
		*/
		void fromImage(const Image<float32_t>& verts, const Image<uint64_t>& edg, const Image<float32_t>* pEdgeMeasurements, const Image<int32_t>* pEdgePoints);


		/**
		Write network to file.
		*/
		void write(const string& filename, bool append) const;

		/**
		Reads one or more networks from a file and places them to given list.
		*/
		static void read(const string& filename, std::vector<Network>& nets);

		/**
		Converts this object to string.
		*/
		friend std::ostream& operator<<(std::ostream& stream, const Network& net)
		{
			stream << "Vertices:" << std::endl;
			for (size_t n = 0; n < net.vertices.size(); n++)
			{
				stream << net.vertices[n] << std::endl;
			}
			stream << "Edges:" << std::endl;
			for (size_t n = 0; n < net.edges.size(); n++)
			{
				stream << net.edges[n] << std::endl;
			}
			stream << "Incomplete vertices:" << std::endl;
			for (size_t n = 0; n < net.incompleteVertices.size(); n++)
			{
				stream << net.incompleteVertices[n] << std::endl;
			}
			stream << "Incomplete edges:" << std::endl;
			for (size_t n = 0; n < net.incompleteEdges.size(); n++)
			{
				stream << net.incompleteEdges[n];
				if (n < net.incompleteEdges.size() - 1)
					stream << std::endl;
			}
			return stream;
		}

		/**
		Removes all edges whose index is given in the array.
		*/
		void removeEdges(const std::vector<coord_t>& edgeIndices, bool removeStraightThrough, bool removeIsolated, bool reportProgress);
	};

	namespace tests
	{
		void disconnections();
		void disconnectStraightThroughPerformance();
		void networkio();
		void pruning();
	}
}
