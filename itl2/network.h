#pragma once

#include "math/vec2.h"
#include "math/vec3.h"
#include "image.h"

#include <vector>

using std::vector;

namespace itl2
{
	struct EdgeMeasurements
	{
	public:
		/**
		Count of pixels in the image that belong to this edge.
		In intersection regions the pixels belonging to the edge might not be neighbours.
		*/
		float32_t pointCount;

		/**
		Length of the edge.
		*/
		float32_t length;

		/**
		Average cross-sectional area.
		*/
		float32_t area;

		/**
		Euclidean distance between the first and the last pixels on the edge.
		*/
		float32_t distance;

		/**
		Euclidean distance between the first and the last pixels on the edge, calculated comparably to the length of the edge.
		*/
		float32_t adjustedDistance;

		EdgeMeasurements() :
			pointCount(0),
			length(0),
			area(0),
			distance(0),
			adjustedDistance(0)
		{

		}

		EdgeMeasurements(float32_t pointCount, float32_t length, float32_t area, float32_t distance, float32_t adjustedDistance) :
			pointCount(pointCount),
			length(length),
			area(area),
			distance(distance),
			adjustedDistance(adjustedDistance)
		{
		}

		friend ostream& operator<<(ostream& stream, const EdgeMeasurements& e)
		{
			stream << "N = " << e.pointCount << ", L = " << e.length << ", A = " << e.area << ", d = " << e.distance << ", d'= " << e.adjustedDistance;
			return stream;
		}

		bool operator==(const EdgeMeasurements& r) const
		{
			return pointCount == r.pointCount && length == r.length && area == r.area && distance == r.distance && adjustedDistance == r.adjustedDistance;
		}
	};


	class Edge
	{
	public:
		/**
		Indices of vertices that this edge connects.
		*/
		math::Vec2c verts;

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

		friend ostream& operator<<(ostream& stream, const Edge& e)
		{
			stream << e.verts + math::Vec2c(1, 1) << ", " << e.properties;
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
		vector<math::Vec3c> points;

		IncompleteVertex() :
			vertexIndex(0)
		{
		}

		IncompleteVertex(size_t vertexIndex, const vector<math::Vec3c>& points) :
			vertexIndex(vertexIndex)
		{
			this->points.reserve(points.size());
			this->points.insert(this->points.begin(), points.begin(), points.end());
		}

		friend ostream& operator<<(ostream& stream, const IncompleteVertex& e)
		{
			stream << e.vertexIndex << " (" << e.points.size() << " points)";
			return stream;
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

	public:
		/**
		Positions of vertices.
		*/
		vector<math::Vec3f> vertices;

		/**
		Edges between vertices.
		*/
		vector<Edge> edges;

		/**
		Incomplete vertices
		*/
		vector<IncompleteVertex> incompleteVertices;


		/**
		Calculate degree of each node and place the values to given list.
		Degree of a node is number of edges connected to it.
		*/
		void degree(vector<size_t>& deg, bool reportProgress) const;

		/**
		Finds parent nodes of node n and adds indices of edges from parent nodes to node n to the given list.
		*/
		void inEdges(size_t n, vector<size_t>& edg) const;

		/**
		Finds child nodes of node n and adds indices of edges from n to child nodes to the given list.
		*/
		void outEdges(size_t n, vector<size_t>& edg) const;

		/**
		Finds indices of edges where this node is target and indices of edges where this node is source.
		*/
		void inOutEdges(size_t n, vector<size_t>& inEdg, vector<size_t>& outEdg) const;

		/**
		Finds neighbours of node n, and stores indices of edges leading to neighbours to the given list.
		*/
		void neighbours(size_t n, vector<size_t>& edgeIndices) const;

		/**
		Invalidates all connections where node n is source or target, and adds new connections from parent nodes of n to child nodes of n.
		*/
		//void disconnect(size_t n, bool reportProgress);

		/**
		Remove edges marked as invalid.
		*/
		void clean(bool reportProgress);

		/**
		Removes all nodes with exactly 2 edges ("straight-through nodes").
		*/
		void disconnectStraightThroughNodes(bool reportProgress);

		/**
		Removes all nodes that have no connections.
		*/
		void removeIsolatedNodes(bool reportProgress);

		/**
		Convert network to a pair of images.
		@param net Network to convert.
		@param vertices Image where vertex coordinates will be set. The size of the image is set to Nx3 where N is the number of vertices.
		@param edges Image where vertex indices corresponding to each edge will be set. The size of the image is set to Mx2 where M is the number of edges.
		@param pEdgeMeasurements Pointer to image that stores (pointCount, length, cross-sectional area, end distance, adjusted end distance) of each edge. Set to null if this information is not required.
		*/
		void toImage(Image<float32_t>& vertices, Image<size_t>& edges, Image<float32_t>* pEdgeMeasurements) const;

		/**
		Replaces this network by one read from given vertex, edge and measurements images.
		*/
		void fromImage(const Image<float32_t>& verts, const Image<size_t>& edg, const Image<float32_t>* pEdgeMeasurements);


		/**
		Write network to file.
		*/
		void write(const string& filename, bool append) const;

		/**
		Reads one or more networks from a file and places them to given list.
		*/
		static void read(const string& filename, vector<Network>& nets);

		/**
		Converts this object to string.
		*/
		friend ostream& operator<<(ostream& stream, const Network& net)
		{
			stream << "Vertices:" << endl;
			for (size_t n = 0; n < net.vertices.size(); n++)
			{
				stream << net.vertices[n] << endl;
			}
			stream << "Edges:" << endl;
			for (size_t n = 0; n < net.edges.size(); n++)
			{
				stream << net.edges[n] << endl;
			}
			stream << "Incomplete (edge) vertices:" << endl;
			for (size_t n = 0; n < net.incompleteVertices.size(); n++)
			{
				stream << net.incompleteVertices[n];
				if (n < net.incompleteVertices.size() - 1)
					stream << endl;
			}
			return stream;
		}

	};

	namespace tests
	{
		void disconnections();
		void disconnectStraightThroughPerformance();
		void networkio();
	}
}
