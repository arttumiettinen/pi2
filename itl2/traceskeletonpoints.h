#pragma once

#include "network.h"
#include <vector>

namespace itl2
{

	/**
	Finds all skeleton pixels and lines connecting them.
	That is, converts traced skeleton to a representation where each pixel on the skeleton is a point and each edge forms a line that connects the points.
	Additional points at possibly non-integer coordinates are inserted for complex intersection regions consisting of many pixels.
	@param net The traced skeleton. The tracing must have been made with storeAllEdgePoints parameter enabled.
	@param points Coordinates of skeleton points will be stored here.
	@param lines Edges will be stored here. Each edge consists of a list of indices into the point list.
	@param smoothingSigma If non-zero, the returned lines will be smoothed using anchored convolution with Gaussian kernel, and this value is the standard deviation of the Gaussian.
	@param maxDisplacement Maximum displacement allowed when smoothing the lines with anchored convolution.
	*/
	void getPointsAndLines(const Network& net, std::vector<Vec3f>& points, std::vector<std::vector<size_t>>& lines, double smoothingSigma = 0, double maxDisplacement = 0.5);
	
	namespace vtk
	{
		/**
		Writes points and lines data to .vtk file.
		Does not append extension to file name.
		*/
		void write(const std::vector<Vec3f>& points, const std::vector<std::vector<size_t>>& lines, const string& filename,
			const string& pointDataName = "", const std::vector<float32_t>* pPointData = nullptr,
			const string& lineDataName = "", const std::vector<float32_t>* pLineData = nullptr);

		/**
		Writes points and lines data to .vtk file.
		Appends extension to file name if it does not contain one.
		*/
		void writed(const std::vector<Vec3f>& points, const std::vector<std::vector<size_t>>& lines, const string& filename,
			const string& pointDataName = "", const std::vector<float32_t>* pPointData = nullptr,
			const string& lineDataName = "", const std::vector<float32_t>* pLineData = nullptr);
	}

	namespace tests
	{
		void skeletonToPointsAndLines();
	}

}
