
#include "registration.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include "io/raw.h"
#include "conversions.h"
#include "transform.h"
#include "projections.h"
#include "pointprocess.h"
#include "filters.h"
#include "inpaint.h"
#include "generation.h"

using namespace std;

namespace itl2
{

	void filterDisplacements(const PointGrid3D<coord_t>& refGrid, Image<Vec3d>& defPoints, Image<float32_t>& accuracy, size_t filterRadius, float32_t threshold)
	{

		// Calculate x, y, and z displacements
		Image<double> u(defPoints.dimensions());
		Image<double> v(defPoints.dimensions());
		Image<double> w(defPoints.dimensions());
		Image<float32_t> gof(defPoints.dimensions());

		for (coord_t z = 0; z < defPoints.depth(); z++)
		{
			for (coord_t y = 0; y < defPoints.height(); y++)
			{
				for (coord_t x = 0; x < defPoints.width(); x++)
				{
					coord_t n = (coord_t)defPoints.getLinearIndex(x, y, z);
					Vec3d U = defPoints(n) - Vec3d(refGrid(x, y, z));
					u(n) = U.x;
					v(n) = U.y;
					w(n) = U.z;
					gof(n) = accuracy(x, y, z);
				}
			}
		}


		// Subtract mean displacement
		double Mu = mean(u);
		double Mv = mean(v);
		double Mw = mean(w);

		subtract(u, Mu);
		subtract(v, Mv);
		subtract(w, Mw);


		// Determine which values are bad by comparing to median
		Image<double> uMed(u.dimensions());
		Image<double> vMed(v.dimensions());
		Image<double> wMed(w.dimensions());

		// Don't account values corresponding to gof = 0 in the calculation of median.
		constexpr double MEDIAN_FLAG = numeric_limits<double>::max();
		for (coord_t n = 0; n < u.pixelCount(); n++)
		{
			if (gof(n) <= 0)
			{
				u(n) = MEDIAN_FLAG;
				v(n) = MEDIAN_FLAG;
				w(n) = MEDIAN_FLAG;
			}
		}

		maskedMedianFilter(u, uMed, filterRadius, MEDIAN_FLAG, NeighbourhoodType::Rectangular, BoundaryCondition::Nearest);
		maskedMedianFilter(v, vMed, filterRadius, MEDIAN_FLAG, NeighbourhoodType::Rectangular, BoundaryCondition::Nearest);
		maskedMedianFilter(w, wMed, filterRadius, MEDIAN_FLAG, NeighbourhoodType::Rectangular, BoundaryCondition::Nearest);


		coord_t badCount = 0;
		constexpr double FLAG = numeric_limits<double>::signaling_NaN();
		for (coord_t n = 0; n < u.pixelCount(); n++)
		{
			bool bad = (fabs(u(n) - uMed(n)) > threshold) ||
						(fabs(v(n) - vMed(n)) > threshold) ||
						(fabs(w(n) - wMed(n)) > threshold) ||
						(gof(n) <= 0);
			if (bad)
			{
				badCount++;
				u(n) = FLAG;
				v(n) = FLAG;
				w(n) = FLAG;
				gof(n) = 0;
			}
		}

		// If there are only bad values, replace everything with zeroes to produce a smooth field without any details.
		if (badCount >= u.pixelCount())
		{
			setValue(u, 0);
			setValue(v, 0);
			setValue(w, 0);
		}


		// Calculate values that will replace the bad ones
		inpaintGarcia(u, FLAG, true);
		inpaintGarcia(v, FLAG, true);
		inpaintGarcia(w, FLAG, true);

		// Add mean displacement back
		add(u, Mu);
		add(v, Mv);
		add(w, Mw);

		// Invert decomposition into components
		for (coord_t z = 0; z < defPoints.depth(); z++)
		{
			for (coord_t y = 0; y < defPoints.height(); y++)
			{
				for (coord_t x = 0; x < defPoints.width(); x++)
				{
					coord_t n = (coord_t)defPoints.getLinearIndex(x, y, z);
					Vec3c X0 = refGrid(x, y, z);
					defPoints(n).x = X0.x + u(n);
					defPoints(n).y = X0.y + v(n);
					defPoints(n).z = X0.z + w(n);
					accuracy(n) = gof(n);
				}
			}
		}
	}
	


	inline void writeBlockMatchResult(const string& filenamePrefix, const PointGrid3D<coord_t>& refPoints, const Image<Vec3d>& defPoints, const Image<float32_t>& gof, double normFact, double normFactStd, double meanDef)
	{
		if (refPoints.pointCount() != defPoints.pixelCount () || gof.pixelCount() != refPoints.pointCount())
			throw ITLException("Point counts in reference, deformed and gof objects must be equal.");

		raw::writed(defPoints, filenamePrefix + "_defpoints");
		raw::writed(gof, filenamePrefix + "_gof");

		ofstream out;
		out.open(filenamePrefix + "_refpoints.txt");
		out << refPoints.xg.first << ", " << refPoints.xg.maximum << ", " << refPoints.xg.step << endl;
		out << refPoints.yg.first << ", " << refPoints.yg.maximum << ", " << refPoints.yg.step << endl;
		out << refPoints.zg.first << ", " << refPoints.zg.maximum << ", " << refPoints.zg.step << endl;
		out << normFact << endl;
		out << normFactStd << endl;
		out << meanDef << endl;

	}

	PointGrid1D<coord_t> readPointGrid1D(std::ifstream& in)
	{
		std::string line;
		std::getline(in, line);

		std::stringstream str;
		str << line;

		std::string val1, val2, val3;
		getline(str, val1, ',');
		getline(str, val2, ',');
		getline(str, val3);

		PointGrid1D<coord_t> g1(fromString<coord_t>(val1), fromString<coord_t>(val2), fromString<coord_t>(val3));
		return g1;
	}

	void readBlockMatchResult(const string& filenamePrefix, PointGrid3D<coord_t>& refPoints, Image<Vec3d>& defPoints, Image<float32_t>& gof, double& normFact, double& normFactStd, double& meanDef)
	{
		ifstream in(filenamePrefix + "_refpoints.txt");

		if (!in.good())
			throw ITLException(string("File not found: ") + filenamePrefix + "_refpoints.txt");

		PointGrid1D<coord_t> g1 = readPointGrid1D(in);
		PointGrid1D<coord_t> g2 = readPointGrid1D(in);
		PointGrid1D<coord_t> g3 = readPointGrid1D(in);
		refPoints = PointGrid3D<coord_t>(g1, g2, g3);

		std::string line;
		getline(in, line);
		normFact = fromString<double>(line);
		getline(in, line);
		normFactStd = fromString<double>(line);
		getline(in, line);
		meanDef = fromString<double>(line);

		raw::read(defPoints, filenamePrefix + "_defpoints_" + toString(g1.pointCount()) + "x" + toString(g2.pointCount()) + "x" + toString(g3.pointCount()) + ".raw");
		raw::read(gof, filenamePrefix + "_gof_" + toString(g1.pointCount()) + "x" + toString(g2.pointCount()) + "x" + toString(g3.pointCount()) + ".raw");
	}


	///*
	//Writes result of block matching to file.
	//*/
	//void writeBlockMatchResult(const std::string& filename, const vector<Vec3c>& refPoints, const vector<Vec3d>& defPoints, const std::vector<double> gof)
	//{
	//	if (refPoints.size() != defPoints.size() || gof.size() != refPoints.size())
	//		throw ITLException("Sizes of reference, deformed and gof point lists must be equal.");

	//	std::ofstream out;
	//	out.open(filename);
	//	for (coord_t n = 0; n < (coord_t)refPoints.size(); n++)
	//	{
	//		Vec3c refPoint = refPoints[n];
	//		Vec3d defPoint = defPoints[n];

	//		out << refPoint.x << ", " << refPoint.y << ", " << refPoint.z << ", " << defPoint.x << ", " << defPoint.y << ", " << defPoint.z << ", " << gof[n] << endl;
	//	}
	//}

	///*
	//Reads file written by writeBlockMatchResult.
	//*/
	//void readBlockMatchResult(const std::string& filename, vector<Vec3d>& refPoints, std::vector<Vec3d>& defPoints, std::vector<double>& gof)
	//{
	//	std::ifstream in(filename);

	//	if (!in.good())
	//		throw ITLException(std::string("Unable to read file ") + filename);

	//	std::vector<double> row;
	//	row.reserve(7);
	//	while (in.good())
	//	{
	//		std::string line;
	//		getline(in, line);

	//		row.clear();
	//		if (!line.empty())
	//		{
	//			std::stringstream str;
	//			str << line;

	//			while (str.good())
	//			{
	//				std::string svalue;
	//				getline(str, svalue, ',');
	//				if (svalue.length() > 0)
	//				{
	//					row.push_back(fromString<double>(svalue));
	//				}
	//			}
	//		}

	//		if (row.size() == 7)
	//		{
	//			Vec3d refPoint(row[0], row[1], row[2]);
	//			Vec3d defPoint(row[3], row[4], row[5]);
	//			double g = row[6];

	//			refPoints.push_back(refPoint);
	//			defPoints.push_back(defPoint);
	//			gof.push_back(g);
	//		}
	//	}
	//}

	namespace tests
	{
		void blockMatch2Match()
		{
			// NOTE: No asserts!

			Image<uint16_t> head16;
			raw::read(head16, "./input_data/t1-head_256x256x129.raw");

			Image<float32_t> reference(head16.dimensions());
			convert(head16, reference);

			Image<uint16_t> deformed16;
			raw::read(deformed16, "./input_data/t1-head_rot_trans_256x256x129.raw");

			Image<float32_t> deformed(deformed16.dimensions());
			convert(deformed16, deformed);


			// Create calculation point grid.
			Vec3d initialShift(0, 0, 0);
			Vec3c step(40, 40, 40);
			Vec3c compRadius(40, 40, 40);

			Vec3c dims = reference.dimensions();
			
			PointGrid3D<coord_t> refPoints(PointGrid1D<coord_t>(0, dims.x, step.x), PointGrid1D<coord_t>(0, dims.y, step.y), PointGrid1D<coord_t>(0, dims.z, step.z));
			Vec3c pointCount = refPoints.pointCounts();
			Image<Vec3d> defPoints(pointCount.x, pointCount.y, pointCount.z);
			Image<float32_t> fitGoodness;

			// Construct initial guess of the deformed points
			for (coord_t zi = 0; zi < pointCount.z; zi++)
			{
				for (coord_t yi = 0; yi < pointCount.y; yi++)
				{
					for (coord_t xi = 0; xi < pointCount.x; xi++)
					{
						defPoints(xi, yi, zi) = Vec3d(refPoints(xi, yi, zi));
					}
				}
			}

			// Perform block matching
			blockMatch(reference, deformed, refPoints, defPoints, fitGoodness, compRadius);

			// Output shifts
			int w = 8;
			cout << "Measurements:" << endl;
			for (coord_t z = 0; z < defPoints.depth(); z++)
			{
				for (coord_t y = 0; y < defPoints.height(); y++)
				{
					for (coord_t x = 0; x < defPoints.width(); x++)
					{
						Vec3d refPoint = Vec3d(refPoints(x, y, z));
						Vec3d defPoint = defPoints(x, y, z);

						Vec3d shift = defPoint - refPoint;

						cout << setw(w) << refPoint.x << ", " << setw(w) << refPoint.y << ", " << setw(w) << refPoint.z << ", " << setw(w) << shift.x << ", " << setw(w) << shift.y << ", " << setw(w) << shift.z << ", " << setw(w) << fitGoodness(x, y, z) << endl;

					}
				}
			}

			writeBlockMatchResult("./registration/blockmatch2_result", refPoints, defPoints, fitGoodness, 0, 1, 0);
		}

		void blockMatch2Pullback()
		{
			// NOTE: No asserts!

			PointGrid3D<coord_t> refPoints;
			Image<Vec3d> defPoints;
			Image<float32_t> fitGoodness;
			double normFact, normFactStd, meanDef;
			readBlockMatchResult("./registration/blockmatch2_result", refPoints, defPoints, fitGoodness, normFact, normFactStd, meanDef);

			Vec3c referenceDimensions(256, 256, 129);

			Image<uint16_t> deformed16;
			raw::read(deformed16, "./input_data/t1-head_rot_trans_256x256x129.raw");

			Image<float32_t> deformed(deformed16.dimensions());
			convert(deformed16, deformed);


			// Calculate reverse deformation

			Image<float32_t> pullback(referenceDimensions);
			
			reverseDeformation(deformed, pullback, refPoints, defPoints);

			raw::writed(pullback, "./registration/blockmatch2_head_pullback");
		}


		/*
		Calculates translation between two 3D images by Maximum Intensity projecting them in two planes
		and by using phase correlation on the projections.
		*/
		template<typename pixel_t> Vec3d mipMatch(const Image<pixel_t>& ref, const Image<pixel_t>& def)
		{
			// XY
			Image<float32_t> refP, defP;
			max(ref, 0, refP);
			max(def, 0, defP);

			raw::writed(refP, "./registration/mipmatch_xy_ref");
			raw::writed(defP, "./registration/mipmatch_xy_def");

			double goodness;
			Vec3c maxXY;
			maxXY.x = refP.width() / 2;
			maxXY.y = refP.height() / 2;
			maxXY.z = 0;
			Vec3d shiftXY = phaseCorrelation(refP, defP, maxXY, goodness);

			raw::writed(refP, "./registration/mipmatch_xy_correlation");

			// YZ
			max(ref, 1, refP);
			max(def, 1, defP);

			raw::writed(refP, "./registration/mipmatch_yz_ref");
			raw::writed(defP, "./registration/mipmatch_yz_def");

			Vec3c maxYZ;
			maxYZ.x = refP.width() / 2;
			maxYZ.y = refP.height() / 2;
			maxYZ.z = 0;
			Vec3d shiftYZ = phaseCorrelation(refP, defP, maxYZ, goodness);

			raw::writed(refP, "./registration/mipmatch_yz_correlation");

			// Total shift
			Vec3d shift(shiftYZ.x, shiftXY.y, -(shiftXY.x - shiftYZ.y) / 2.0);

			return shift;
		}


		void blockMatch1()
		{
			// NOTE: No asserts!

			Image<uint16_t> head16;
			raw::read(head16, "./input_data/t1-head_256x256x129.raw");

			Image<float32_t> head(head16.dimensions());
			convert(head16, head);



			Vec3d shiftGT(-10.75, 8.5, -12.02);

			Image<float32_t> headShifted(head.dimensions());
			itl2::translate(head, headShifted, shiftGT, LinearInterpolator<float32_t, float32_t>(BoundaryCondition::Zero));

			raw::writed(headShifted, "./registration/blockmatch1_head_shifted_gt");


			// Create calculation point grid.
			Vec3d initialShift(0, 0, 0);
			Vec3c step(50, 50, 50);
			Vec3c compRadius(50, 50, 50);

			Vec3c dims = head.dimensions();

			PointGrid3D<coord_t> refPoints(PointGrid1D<coord_t>(0, dims.x, step.x), PointGrid1D<coord_t>(0, dims.y, step.y), PointGrid1D<coord_t>(0, dims.z, step.z));
			Vec3c pointCount = refPoints.pointCounts();
			Image<Vec3d> defPoints(pointCount.x, pointCount.y, pointCount.z);
			Image<float32_t> fitGoodness;

			// Construct initial guess of the deformed points
			for (coord_t zi = 0; zi < pointCount.z; zi++)
			{
				for (coord_t yi = 0; yi < pointCount.y; yi++)
				{
					for (coord_t xi = 0; xi < pointCount.x; xi++)
					{
						defPoints(xi, yi, zi) = Vec3d(refPoints(xi, yi, zi));
					}
				}
			}
			
			// Perform block matching
			blockMatch(head, headShifted, refPoints, defPoints, fitGoodness, compRadius);


			// Output shifts, calculate average shift
			Vec3d avg;
			double weight = 0;

			int w = 8;
			cout << "True:" << endl;
			cout << setw(w) << shiftGT.x << ", " << setw(w) << shiftGT.y << ", " << setw(w) << shiftGT.z << endl;
			cout << "Measurements:" << endl;
			for (coord_t z = 0; z < defPoints.depth(); z++)
			{
				for (coord_t y = 0; y < defPoints.height(); y++)
				{
					for (coord_t x = 0; x < defPoints.width(); x++)
					{
						Vec3d refPoint = Vec3d(refPoints(x, y, z));
						Vec3d defPoint = defPoints(x, y, z);

						Vec3d shift = refPoint - defPoint;

						cout << setw(w) << shift.x << ", " << setw(w) << shift.y << ", " << setw(w) << shift.z << ", " << setw(w) << fitGoodness(x, y, z) << endl;

						avg += shift * fitGoodness(x, y, z);
						weight += fitGoodness(x, y, z);
					}
				}
			}
			avg /= weight;

			cout << "average shift = " << avg << endl;

			Vec3d shiftMIP = mipMatch<float32_t>(head, headShifted);
			cout << "mipMatch shift = " << shiftMIP << endl;


			Image<float32_t> headShifted2;
			translate(head, headShifted2, avg);

			raw::writed(headShifted2, "./registration/blockmatch1_head_shifted_meas_avg_shift");


			// Calculate reverse deformation
			Image<float32_t> headPullback(head.dimensions());
			
			reverseDeformation(headShifted, headPullback, refPoints, defPoints);

			raw::writed(headPullback, "./registration/blockmatch1_head_pullback");
		}


		void mipMatch()
		{
			Image<uint16_t> head16;
			raw::read(head16, "./input_data/t1-head_256x256x129.raw");

			Image<float32_t> head(head16.dimensions());
			convert(head16, head);

			Vec3d shiftGT(-10.75, 8.5, -12.02);

			Image<float32_t> headShifted(head.dimensions());
			itl2::translate(head, headShifted, shiftGT, LinearInterpolator<float32_t, float32_t>(BoundaryCondition::Zero));

			raw::writed(headShifted, "./registration/mipmatch_head_shifted_GT");

			Vec3d shift = -mipMatch<float32_t>(head, headShifted);

			cout << "GT = " << shiftGT << endl;
			cout << "meas = " << shift << endl;

			testAssert((shift - shiftGT).norm() < 1, "MIP shift");
		}

		void pointsToDeformed()
		{
			Image<uint8_t> img1(100, 100, 1);
			setValue(img1, 1);
			draw(img1, Sphere<double>(Vec3d(45.0, 50.0, 0.0), 10.0), (uint8_t)255);
			raw::writed(img1, "./pointstodeformed/img1");

			Image<uint8_t> img2(100, 100, 1);
			setValue(img2, 1);
			draw(img2, Sphere<double>(Vec3d(50.0, 50.0, 0.0), 10.0), (uint8_t)255);
			raw::writed(img2, "./pointstodeformed/img2");


			Vec3c step(50, 50, 50);
			Vec3c compRadius(50, 50, 1);

			Vec3c dims = img1.dimensions();

			PointGrid3D<coord_t> refPoints(PointGrid1D<coord_t>(0, dims.x, step.x), PointGrid1D<coord_t>(0, dims.y, step.y), PointGrid1D<coord_t>(0, dims.z, step.z));
			Vec3c pointCount = refPoints.pointCounts();
			Image<Vec3d> defPoints(pointCount.x, pointCount.y, pointCount.z);
			Image<float32_t> fitGoodness;

			// Construct initial guess of the deformed points
			for (coord_t zi = 0; zi < pointCount.z; zi++)
			{
				for (coord_t yi = 0; yi < pointCount.y; yi++)
				{
					for (coord_t xi = 0; xi < pointCount.x; xi++)
					{
						defPoints(xi, yi, zi) = Vec3d(refPoints(xi, yi, zi));
					}
				}
			}

			blockMatch(img1, img2, refPoints, defPoints, fitGoodness, compRadius);


			for (coord_t z = 0; z < defPoints.depth(); z++)
			{
				for (coord_t y = 0; y < defPoints.height(); y++)
				{
					for (coord_t x = 0; x < defPoints.width(); x++)
					{
						Vec3d refPoint = Vec3d(refPoints(x, y, z));
						Vec3d defPoint = defPoints(x, y, z);

						cout << refPoint << " --> " << defPoint << endl;
					}
				}
			}

			std::vector<Vec3d> points;
			Vec3d p1(45, 50, 0);
			points.push_back(p1);

			itl2::pointsToDeformed(points, refPoints, defPoints);

			cout << "Reference point = " << p1 << endl;
			cout << "Deformed point = " << points[0] << endl;
		}
	}
}
