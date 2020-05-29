
#include "traceskeletonpoints.h"

#include "image.h"
#include "generation.h"
#include "floodfill.h"
#include "traceskeleton.h"

#include "io/fileutils.h"

#include <map>

using namespace std;

namespace itl2
{

	void getPointsAndLines(const Network& net, std::vector<Vec3f>& points, std::vector<std::vector<size_t>>& lines, double smoothingSigma, double maxDisplacement)
	{
		points.clear();
		map<Vec3f, size_t, decltype(vecComparer<float32_t>)* > indexMap(vecComparer<float32_t>);

		// Add all vertices to the points list
		for (int n = 0; n < net.vertices.size(); n++)
		{
			Vec3f p = net.vertices[n];
			indexMap[p] = points.size();
			points.push_back(p);
		}

		// Add all edge points to the points list and create corresponding index lists.
		lines.clear();
		for (int n = 0; n < net.edges.size(); n++)
		{
			const Edge& edge = net.edges[n];

			vector<size_t> indices;
			indices.reserve(edge.properties.edgePoints.size() + 2);

			for (size_t m = 0; m < edge.properties.edgePoints.size(); m++)
			{
				Vec3f p = Vec3f(edge.properties.edgePoints[m]);
				if (indexMap.count(p) <= 0)
				{
					indexMap[p] = points.size();
					points.push_back(p);
				}

				indices.push_back(indexMap[p]);
			}

			bool reverse = false;

			if (edge.properties.edgePoints.size() > 0)
			{
				float32_t d1 = (Vec3f(edge.properties.edgePoints.front()) - net.vertices[edge.verts[0]]).norm<float32_t>();
				float32_t d2 = (Vec3f(edge.properties.edgePoints.front()) - net.vertices[edge.verts[1]]).norm<float32_t>();
				reverse = d2 < d1;
			}
			
			if (!reverse)
			{
				indices.insert(indices.begin(), indexMap[net.vertices[edge.verts[0]]]);
				indices.push_back(indexMap[net.vertices[edge.verts[1]]]);
			}
			else
			{
				indices.insert(indices.begin(), indexMap[net.vertices[edge.verts[1]]]);
				indices.push_back(indexMap[net.vertices[edge.verts[0]]]);
			}

			lines.push_back(indices);
		}

		// If smoothing is requested, do it here
		if (smoothingSigma > 0)
		{
			// Note that each line has distinct set of points except intersection points in its ends.
			// We leave the end points to the original locations and smooth other points.
			for (const vector<size_t>& indices : lines)
			{
				// First get points on the line
				vector<Vec3f> linePoints;
				linePoints.reserve(indices.size());
				for (size_t i = 0; i < indices.size(); i++)
					linePoints.push_back(points[indices[i]]);

				// Smooth them
				smoothLine(linePoints, smoothingSigma, maxDisplacement);

				// Assign the smoothed points back. We don't change the first and the last point.
				for (size_t i = 1; i < indices.size() - 1; i++)
					points[indices[i]] = linePoints[i];
			}
		}


		// DEBUG
		//for (const vector<size_t>& indices : lines)
		//{
		//	float32_t l = 0;
		//	for (size_t i = 1; i < indices.size(); i++)
		//	{
		//		l += (points[indices[i]] - points[indices[i - 1]]).norm();
		//	}

		//	cout << l << endl;
		//}

	}

	namespace vtk
	{
		namespace internals
		{
			/**
			Helper to write
			FIELD FieldData N
			[Array name] 1 M
			[element 1]
			[element 2]
			...
			[element M]
			[Array name] 1 K
			[element 1]
			...
			*/
			void writeVtkFields(ofstream& out, const std::vector<std::tuple<std::string, std::vector<float32_t>>>& data)
			{
				out << "FIELD FieldData " << data.size() << endl; // TODO: one "out << pointDataName << " 1 " << pPointData->size() << " float" <<endl;" for each array
				for (size_t n = 0; n < data.size(); n++)
				{
					const string& arrayName = get<0>(data[n]);
					const auto& arr = get<1>(data[n]);
					out << arrayName << " 1 " << arr.size() << " float" << endl;

					for (size_t n = 0; n < arr.size(); n++)
					{
						out << arr[n] << endl;
					}
				}
			}
		}
		
		void write(const std::vector<Vec3f>& points, const std::vector<std::vector<size_t>>& lines, const string& filename,
			const std::vector<std::tuple<std::string, std::vector<float32_t>>>* pPointData,
			const std::vector<std::tuple<std::string, std::vector<float32_t>>>* pLineData)
		{
			createFoldersFor(filename);

			ofstream out(filename.c_str(), ios_base::out | ios_base::binary);

			if (!out)
				throw ITLException(std::string("Unable to open ") + filename + string(", ") + getStreamErrorMessage());

			out << "# vtk DataFile Version 3.0\n";
			out << "itl2 vtk output\n";
			out << "ASCII\n";
			out << "DATASET POLYDATA\n";
			out << "POINTS " << points.size() << " float\n";

			for (size_t n = 0; n < points.size(); n++)
			{
				Vec3f p = points[n];
				//out.write((char*)&p.components[0], 3 * sizeof(float32_t));
				out << p.x << " " << p.y << " " << p.z << "\n";
			}

			size_t totalCount = 0;
			for (size_t n = 0; n < lines.size(); n++)
			{
				totalCount += 1 + lines[n].size();
			}

			out << "LINES " << lines.size() << " " << totalCount << "\n";

			for (size_t n = 0; n < lines.size(); n++)
			{
				vector<size_t> indices = lines[n];
				out << indices.size() << " ";
				for (size_t m = 0; m < indices.size(); m++)
				{
					out << indices[m];
					if (m < indices.size() - 1)
						out << " ";
				}
				out << "\n";
			}

			if (pPointData)
			{
				out << "POINT_DATA " << points.size() << endl;
				internals::writeVtkFields(out, *pPointData);
			}

			if (pLineData)
			{
				out << "CELL_DATA " << lines.size() << endl;
				internals::writeVtkFields(out, *pLineData);
			}
		}
		
		void writed(const std::vector<Vec3f>& points, const std::vector<std::vector<size_t>>& lines, const string& filename,
			const std::vector<std::tuple<std::string, std::vector<float32_t>>>* pPointData,
			const std::vector<std::tuple<std::string, std::vector<float32_t>>>* pLineData)
		{
			if (endsWithIgnoreCase(filename, ".vtk"))
				write(points, lines, filename, pPointData, pLineData);
			else
				write(points, lines, filename + ".vtk", pPointData, pLineData);
		}
	}

	namespace tests
	{
		void skeletonToPointsAndLines()
		{
			Image<uint8_t> geometry(100, 100);
			Vec3f p0(1, 1, 0);
			//Vec3f p0(50, 25, 0);
			Vec3f p1(25, 50, 0);
			Vec3f p2(75, 50, 0);
			Vec3f p3(5, 95, 0);
			Vec3f p4(95, 95, 0);

			Vec3f pOut(50, 0, 0);

			Vec3f pEdge1(25, 0, 0);
			Vec3f pEdge1in(25, 5, 0);
			Vec3f pEdge2(75, 0, 0);
			Vec3f pEdge2in(75, 5, 0);


			draw(geometry, Line(p0, p1), (uint8_t)255);
			draw(geometry, Line(p1, p2), (uint8_t)255);
			draw(geometry, Line(p0, p2), (uint8_t)255);
			draw(geometry, Line(p1, p3), (uint8_t)255);
			draw(geometry, Line(p2, p4), (uint8_t)255);

			draw(geometry, Line(pOut, p2), (uint8_t)255);
			draw(geometry, Line(pEdge1, pEdge2), (uint8_t)255);
			draw(geometry, Line(pEdge1, pEdge1in), (uint8_t)255);
			draw(geometry, Line(pEdge2, pEdge2in), (uint8_t)255);

			raw::writed(geometry, "./skeleton_points_and_lines/geometry");
			lineSkeleton(geometry);
			raw::writed(geometry, "./skeleton_points_and_lines/skeleton");


			Network net;
			traceLineSkeleton(geometry, true, 1.5, 1.0, net);

			/*
			vector<coord_t> removeIndices;
			removeIndices.push_back(2);
			net.removeEdges(removeIndices, true, true, true);
			*/

			//raw::read(geometry, "./skeleton_points_and_lines/skeleton");
			vector<Vec3f> points;
			vector<vector<size_t>> lines;
			getPointsAndLines(net, points, lines, 1.0, 1.0);

			// Draw
			for (size_t n = 0; n < lines.size(); n++)
			{
				const vector<size_t>& line = lines[n];
				for (size_t m = 0; m < line.size(); m++)
				{
					geometry(round(points[line[m]])) = pixelRound<uint8_t>(m + 1);
				}
			}

			raw::writed(geometry, "./skeleton_points_and_lines/point_line_vis");


			vector<float32_t> radius;
			radius.reserve(points.size());
			for (size_t n = 0; n < points.size(); n++)
			{
				radius.push_back((n + 1) / 10.0f);
			}

			vector<float32_t> index;
			index.reserve(points.size());
			for (size_t n = 0; n < points.size(); n++)
			{
				index.push_back((float32_t)n);
			}

			vector<float32_t> lineRadius;
			lineRadius.reserve(lines.size());
			for (size_t n = 0; n < lines.size(); n++)
			{
				lineRadius.push_back((n + 1) / 10.0f);
			}

			auto pointDataArrays = vector<tuple<string, vector<float32_t>>>{ make_tuple("radius", radius), make_tuple("index", index) };
			auto lineDataArrays = vector<tuple<string, vector<float32_t>>>{ make_tuple("radius", lineRadius) };
			vtk::writed(points, lines, "../testing/skeleton_points_and_lines/point_line_vis", &pointDataArrays, &lineDataArrays);
		}
	}
}