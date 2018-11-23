#pragma once

#include <iostream>
#include <algorithm>

#include "image.h"
#include "floodfill.h"
#include "utilities.h"

//#include "Eigen/Core"
//#include "Eigen/Eigenvalues"
#include "math/eig3.h"
#include "math/eig2.h"
#include "math/matrix3x3.h"

//#include "convexhull.h"
//#include "marchingcubes.h"
//#include "sphere.h"

#include "math/vec2.h"
#include "math/vec3.h"


namespace itl2
{

	using math::Vec2c;
	using math::Vec3c;
	using math::pixelRound;

	///**
	//Number labels such that the numbering begins from one and is consecutive.
	//*/
	//template<typename T> void renumberLabels(Image<T>& image)
	//{
	//	T maxLabel = 0;
	//	map<T, T> oldtonew;	// Maps color of current image to its color in result.
	//	for(size_t n = 0; n < image.pixelCount(); n++)
	//	{
	//		T oldval = **p;
	//		T newval = oldtonew[oldval];
	//		if(newval == 0)
	//		{
	//			// Take unused color (unused in labels and unused in new labels)
	//			maxLabel++;
	//			newval = maxLabel;
	//			oldtonew[oldval] = newval;
	//		}
	//		**p = newval;
	//	}
	//}

	///**
	//Re-colors particles such that particles in the same disjoint set are colored with the same color.
	//@param image The image containing the particles.
	//@param particleRelations Forest describing the relations of the particles.
	//*/
	//template<typename T> void recolorParticles(Image<T>& image, DisjointSetForest<T>& particleRelations, T backgroundColor = 0)
	//{
	//	CursorPointer<T, LinearCursor<T> > p(image, image.createLinearCursor());
	//	for(; p->hasNext(); p->proceed())
	//	{
	//		T oldColor = **p;
	//		T newColor;

	//		if(oldColor != backgroundColor)
	//			newColor = particleRelations.find_set(oldColor);
	//		else
	//			newColor = 0;

	//		if(oldColor != newColor)
	//			**p = newColor;
	//	}
	//}

	/**
	* Label all particles with distinct colors beginning from one.
	* It is assumed that background pixels are set to zero.
	* Uses simple flood fill algorithm.
	* @param image Image containing the particles.
	* @param particleColor Color of particles. This color will be skipped in labeling. Pass zero to label particles of any color.
	* @param firstLabelValue Value of the first label.
	* @returns The largest label value used in the image.
	*/
	template<typename pixel_t> pixel_t labelParticles(Image<pixel_t>& image, pixel_t particleColor, pixel_t firstLabelValue = 1, Connectivity connectivity = NearestNeighbours, bool showProgress = true)
	{
		size_t counter = 0;

		pixel_t particleNumber = firstLabelValue;
		if(particleNumber == particleColor)
			particleNumber++;

		for (coord_t z = 0; z < image.depth(); z++)
		{
			for (coord_t y = 0; y < image.height(); y++)
			{
				for (coord_t x = 0; x < image.width(); x++)
				{
					pixel_t pixel = image(x, y, z);

					if ((particleColor != 0 && pixel == particleColor) || (particleColor == 0 && pixel != 0))
					{
						floodfill(image, Vec3c(x, y, z), particleNumber, particleNumber, connectivity);

						particleNumber++;
						if (particleNumber == particleColor)
							particleNumber++;
					}

					showThreadProgress(counter, image.pixelCount(), showProgress);
				}
			}
		}

		return particleNumber - 1;
	}

    /**
    Class that specifies one analyzer for particle analysis.
	@param POINT Storage type for point coordinate class. Should be indexable with [0] and have size() member.
    */
    template<class POINT, class PIXELTYPE> class Analyzer
    {
    public:
		/**
			* Virtual destructor
			*/
		virtual ~Analyzer()
		{
		}

        /**
        Get title of this analysis.
        */
        virtual vector<string> getTitles() = 0;

        /**
        Analyze the particles and return the result.
        */
		virtual vector<double> analyze(const PIXELTYPE pixelValue, const vector<POINT>& points) = 0;
    };

	/**
	* Vector for analyzers. Deletes its contents automatically.
	*/
	template<class POINT, typename PIXELTYPE> class AnalyzerVector : public vector<Analyzer<POINT, PIXELTYPE>*>
	{
	private:
		AnalyzerVector(const AnalyzerVector<POINT, PIXELTYPE>& right)
		{
		}
	public:
		AnalyzerVector()
		{
		}

		virtual ~AnalyzerVector()
		{
			for(size_t n = 0; n < vector<Analyzer<POINT, PIXELTYPE>*>::size(); n++)
				delete vector<Analyzer<POINT, PIXELTYPE>*>::operator[](n);
		}

		/**
		Get list of column headers for all the analyzers in this vector.
		*/
		virtual vector<string> getHeaders() const
		{
			vector<string> headers;
			headers.reserve(this->size());
			for(size_t n = 0; n < this->size(); n++)
			{
				vector<string> titles = (*this)[n]->getTitles();
				for(size_t m = 0; m < titles.size(); m++)
					headers.push_back(titles[m]);
			}

			return headers;
		}

		/**
		Perform analysis using all the analyzers in this vector.
		@param particlePoints Points on the particle.
		@param resultLine Results will be added to this array.
		*/
		virtual void analyze(const PIXELTYPE pixel, const vector<POINT>& particlePoints, vector<double>& resultLine)
		{
			size_t s = vector<Analyzer<POINT, PIXELTYPE>*>::size();

			if(particlePoints.size() <= 0)
				throw ITLException("Particle analyzer received no input points.");

			resultLine.reserve(s);
            for(size_t n = 0; n < s; n++)
            {
				vector<double> currentResults = (*this)[n]->analyze(pixel, particlePoints);
				for(size_t m = 0; m < currentResults.size(); m++)
					resultLine.push_back(currentResults[m]);
            }
		}
	};

	namespace internal
	{

		template<typename T> class SpecialColors
		{
		public:
			static T fillColor();
			static T largeColor();
		};

		template<> class SpecialColors<uint8_t>
		{
		public:
			static uint8_t fillColor()
			{
				return 128;
			}

			static uint8_t largeColor()
			{
				return 230;
			}
		};

		template<> class SpecialColors<float32_t>
		{
		public:
			static float32_t fillColor()
			{
				return numeric_limits<float32_t>::max();
			}

			static float32_t largeColor()
			{
				return numeric_limits<float32_t>::lowest();
			}
		};

		template<> class SpecialColors<uint16_t>
		{
		public:
			static uint16_t fillColor()
			{
				return numeric_limits<uint16_t>::max();
			}

			static uint16_t largeColor()
			{
				return numeric_limits<uint16_t>::max() - 1;
			}
		};

		template<> class SpecialColors<uint32_t>
		{
		public:
			static uint32_t fillColor()
			{
				return numeric_limits<uint32_t>::max();
			}

			static uint32_t largeColor()
			{
				return numeric_limits<uint32_t>::max() - 1;
			}
		};

		template<> class SpecialColors<uint64_t>
		{
		public:
			static uint64_t fillColor()
			{
				return numeric_limits<uint64_t>::max();
			}

			static uint64_t largeColor()
			{
				return numeric_limits<uint64_t>::max() - 1;
			}
		};
	}

	/**
    Perform particle analysis, 3D version.
    Has no restrictions on particle count. (uses flood fill algorithm)
    Needs image that can store at least four distinct colours - original particle color, temporary color, big particle color and background color.
    Particles colored with big particle color are skipped. Additionally, particles that are larger than volumeLimit will be colored with big
    particle color. Other particles will be colored with temporary color.
    @param image The image containing the particles.
    @param analyzers List of analyzers to apply to the particles.
    @param results List for result matrix.
    @param volumeLimit Only include particles smaller than this value in the results.
    @param fillColor The analyzed particles will be colored with this color.
    @param largeColor Particles that are skipped because their size is larger than volumeLimit are colored with this color.
    @param backgroundColor Color of background.
    */
    template<typename pixel_t> void analyzeParticles3D(Image<pixel_t>& image, AnalyzerVector<Vec3c, pixel_t>& analyzers, vector<vector<double> >& results, Connectivity connectivity = AllNeighbours, size_t volumeLimit = numeric_limits<size_t>::max(), pixel_t fillColor = internal::SpecialColors<pixel_t>::fillColor(), pixel_t largeColor = internal::SpecialColors<pixel_t>::largeColor(), pixel_t backgroundColor = 0)
    {
		if(image.dimensionality() != 3)
			throw ITLException("This method supports only 3-dimensional images.");

        if(fillColor == largeColor)
            throw ITLException("Fill color and large color must not be equal.");
        if(largeColor == backgroundColor)
            throw ITLException("Background color and large color must not be equal.");

        
        vector<Vec3c> particlePoints;
        particlePoints.reserve(1000);

        size_t counter = 0;
		for (coord_t z = 0; z < image.depth(); z++)
		{
			for (coord_t y = 0; y < image.height(); y++)
			{
				for (coord_t x = 0; x < image.width(); x++)
				{
					pixel_t pixel = image(x, y, z);

					if (pixel != fillColor && pixel != backgroundColor && pixel != largeColor)
					{
						if (floodfill(image, Vec3c(x, y, z), fillColor, largeColor, connectivity, &particlePoints, volumeLimit))
						{
							vector<double> resultLine;
							analyzers.analyze(pixel, particlePoints, resultLine);
							results.push_back(resultLine);
						}
						else
						{
							// The fill was ended by size limit. Mark the filled points as belonging to a large particle.

							for (size_t n = 0; n < particlePoints.size(); n++)
								image(particlePoints[n]) = largeColor;
						}

						particlePoints.clear();
					}


					showThreadProgress(counter, image.pixelCount());
				}
			}
		}

    }

	///**
	//Perform greedy coloring of particles in 3D.
	//Colors each region in image such that its neighbours are colored with different colors.
	//Uses greedy algorithm.
	//Needs image that can store at least four distinct colours - original particle color, temporary color, big particle color and background color.
 //   Particles colored with big particle color are skipped. Additionally, particles that are larger than volumeLimit will be colored with big
 //   particle color. Other particles will be colored with colors such that neighbouring particles are colored with different colors.
 //   @param image The image containing the particles.
 //   @param volumeLimit Only include particles smaller than this value in the results.
 //   @param fillColor The analyzed particles will be colored with this color.
 //   @param largeColor Particles that are skipped because their size is larger than volumeLimit are colored with this color.
 //   @param backgroundColor Color of background.
 //   */
 //   template<typename T> void colorMapGreedy3D(Image<T>& image, Connectivity connectivity = All, size_t volumeLimit = numeric_limits<size_t>::max(), T fillColor = internal::SpecialColors<T>::fillColor(), T largeColor = internal::SpecialColors<T>::largeColor(), T backgroundColor = 0)
 //   {
	//	if(image.dimensionality() != 3)
	//		throw ITLException("This method supports only 3-dimensional images.");

 //       if(fillColor == largeColor)
 //           throw ITLException("Fill color and large color must not be equal.");
 //       if(largeColor == backgroundColor)
 //           throw ITLException("Background color and large color must not be equal.");

 //       CursorPointer<T, LinearCursor<T> > cursor(image, image.createLinearCursor());

	//	set<T> usedColors;

	//	set<T> neighbourColors;

 //       vector<Vec3c> particlePoints;
 //       particlePoints.reserve(1000);

 //       size_t n = 0;
 //       for(; cursor->hasNext(); cursor->proceed())
 //       {
 //           T pixel = **cursor;
 //           if(pixel != fillColor && pixel != backgroundColor && pixel != largeColor && usedColors.find(pixel) == usedColors.end())
 //           {
 //               if(singlethreaded::floodfill3D(image, cursor->getPosition(), fillColor, largeColor, connectivity, &particlePoints, volumeLimit, &neighbourColors))
 //               {
 //                   // Analyze the points only if the volume of the particle is small enough.

	//				// Get color not in use by neighbours
	//				T newColor = 0;
	//				for(; newColor < numeric_limits<T>::max(); newColor++)
	//				{
	//					if(newColor != fillColor && newColor != backgroundColor && newColor != largeColor && neighbourColors.find(newColor) == neighbourColors.end())
	//						break;
	//				}

	//				usedColors.insert(newColor);

	//				// Fill the particle with that color
	//				for(size_t n = 0; n < particlePoints.size(); n++)
 //                       image.setPixel(particlePoints[n], newColor);

 //               }
 //               else
 //               {
 //                   // The fill was ended by size limit. Mark the filled points as belonging to a large particle.

 //                   for(size_t n = 0; n < particlePoints.size(); n++)
 //                       image.setPixel(particlePoints[n], largeColor);
 //               }

 //               particlePoints.clear();
	//			neighbourColors.clear();
 //           }

 //           n++;
	//		utilities::showProgress(n, image.pixelCount());
 //       }

 //   }

 //   /**
 //   Perform particle analysis, 2D version.
 //   Has no restrictions on particle count. (uses flood fill algorithm)
 //   Needs image that can store at least four distinct colours - original particle color, temporary color, big particle color and background color.
 //   Particles colored with big particle color are skipped. Additionally, particles that are larger than volumeLimit will be colored with big
 //   particle color. Other particles will be colored with temporary color.
 //   @param image The image containing the particles.
 //   @param analyzers List of analyzers to apply to the particles.
 //   @param results List for result matrix.
 //   @param volumeLimit Only include particles smaller than this value in the results.
 //   @param fillColor The analyzed particles will be colored with this color.
 //   @param largeColor Particles that are skipped because their size is larger than volumeLimit are colored with this color.
 //   @param backgroundColor Color of background.
 //   */
 //   template<typename T> void analyzeParticles2D(Image<T>& image, AnalyzerVector<Vec2c, T>& analyzers, vector<vector<double> >& results, Connectivity connectivity = All, size_t volumeLimit = numeric_limits<size_t>::max(), T fillColor = 128, T largeColor = 230, T backgroundColor = 0)
 //   {
	//	if(image.dimensionality() != 2)
	//		throw ITLException("This method supports only 2-dimensional images.");

 //       if(fillColor == largeColor)
 //           throw ITLException("Fill color and large color must not be equal.");
 //       if(largeColor == backgroundColor)
 //           throw ITLException("Background color and large color must not be equal.");

 //       CursorPointer<T, LinearCursor<T> > cursor(image, image.createLinearCursor());

 //       vector<Vec2c> particlePoints;
 //       particlePoints.reserve(1000);

 //       size_t n = 0;
 //       for(; cursor->hasNext(); cursor->proceed())
 //       {
 //           T pixel = **cursor;
 //           if(pixel != fillColor && pixel != backgroundColor && pixel != largeColor)
 //           {
 //               //cout << (int)mathutils::round(cursor->getProgress() * 100) << " %, last particle at " << cursor->getPosition() << "\r";

 //               if(singlethreaded::floodfill2D(image, cursor->getPosition2D(), fillColor, largeColor, connectivity, &particlePoints, volumeLimit))
 //               {
 //                   // Analyze the points only if the volume of the particle is small enough.

 //                   vector<double> resultLine;
	//				analyzers.analyze(pixel, particlePoints, resultLine);
 //                   results.push_back(resultLine);
 //               }
 //               else
 //               {
 //                   // The fill was ended by size limit. Mark the filled points as belonging to a large particle.

 //                   for(size_t n = 0; n < particlePoints.size(); n++)
 //                       image.setPixel(particlePoints[n], largeColor);
 //               }

 //               particlePoints.clear();
 //           }

 //           n++;
 //           utilities::showProgress(n, image.pixelCount());
 //       }
 //   }


	/**
	Find in index of element name in headers array.
	*/
	inline size_t getColumnIndex(const string& name, const vector<string> headers)
	{
		vector<string>::const_iterator result = std::find(headers.begin(), headers.end(), name);
		if(result == headers.end())
			throw ITLException("Column not found.");

		return result - headers.begin();
	}

	///**
	//* Fill all particles that have some measurement value above or below some threshold.
	//* @param image The image.
	//* @param results Particle analysis results vector.
	//* @param fillColor Color to use when filling the particles.
	//* @param condi Index of column containing the value to test.
	//* @param condThreshold Threshold for the condition value.
	//* @param higher If set to true, fill if condition value for particle is larger than the threshold value; If set to false, fill if condition value is lower than threshold value.
	//*/
 //   template<typename T> void fillParticles2D(Image<T>& image, const vector<string>& headers, const vector<vector<double> >& results, T fillColor, int condi, double condThreshold, Connectivity conn, bool higher = false)
	//{

	//	size_t xi = getColumnIndex("X [pixel]", headers);
	//	size_t yi = getColumnIndex("Y [pixel]", headers);

	//	// TODO: Make n-dimensional version
	//	// TODO: This can be threaded easily

 //       Vec2c pos;

	//	for(size_t n = 0; n < results.size(); n++)
	//	{
	//		double value = results[n][condi];
	//		bool fill = higher ? value > condThreshold : value <= condThreshold;

	//		if(fill)
	//		{
	//			pos.x = mathutils::round(results[n][xi]);
	//			pos.y = mathutils::round(results[n][yi]);

	//			singlethreaded::floodfill2D<T>(image, pos, fillColor, fillColor, conn);
	//		}

	//		utilities::showProgress(n, results.size());
	//	}
	//}

	///**
	//* Fill all particles with their corresponding measurement value.
	//* @param image The image.
	//* @param results Particle analysis results vector.
	//* @param condi Index of column containing the fill value
	//* @param multiplier Multiplier for fill value.
	//*/
 //   template<typename T> void fillParticles2D(Image<T>& image, const vector<string>& headers, const vector<vector<double> >& results, int condi, double multiplier, Connectivity conn)
	//{
	//	size_t xi = getColumnIndex("X [pixel]", headers);
	//	size_t yi = getColumnIndex("Y [pixel]", headers);

	//	// TODO: This can be threaded easily

 //       Vec2c pos;

	//	for(size_t n = 0; n < results.size(); n++)
	//	{
	//		T value = pixelRound<T>(results[n][condi] * multiplier);

	//		pos.x = mathutils::round(results[n][xi]);
	//		pos.y = mathutils::round(results[n][yi]);

	//		singlethreaded::floodfill2D<T>(image, pos, value, value, conn);

	//		utilities::showProgress(n, results.size());
	//	}
	//}

    /**
	* Fill all particles that have some measurement value above or below some threshold.
	* @param image The image.
	* @param results Particle analysis results vector.
	* @param fillColor Color to use when filling the particles.
	* @param condi Index of column containing the value to test.
	* @param condThreshold Threshold for the condition value.
	* @param higher If set to true, fill if condition value for particle is larger than the threshold value; If set to false, fill if condition value is lower than threshold value.
	*/
    template<typename pixel_t> void fillParticles3D(Image<pixel_t>& image, const vector<string>& headers, const vector<vector<double> >& results, pixel_t fillColor, int condi, double condThreshold, Connectivity conn, bool higher = false)
	{
		size_t xi = getColumnIndex("X [pixel]", headers);
		size_t yi = getColumnIndex("Y [pixel]", headers);
		size_t zi = getColumnIndex("Z [pixel]", headers);

		size_t counter = 0;
		#pragma omp parallel for
		for(coord_t n = 0; n < (coord_t)results.size(); n++)
		{
			double value = results[n][condi];
			bool fill = higher ? value > condThreshold : value <= condThreshold;

			if(fill)
			{
				Vec3c pos;
				pos.x = (coord_t)::round(results[n][xi]);
				pos.y = (coord_t)::round(results[n][yi]);
				pos.z = (coord_t)::round(results[n][zi]);

				floodfill<pixel_t>(image, pos, fillColor, fillColor, conn);
			}

			showThreadProgress(counter, results.size());
		}
	}

	/**
	* Fill all particles with their corresponding measurement value.
	* @param image The image.
	* @param results Particle analysis results vector.
	* @param condi Index of column containing the fill value
	* @param multiplier Multiplier for fill value.
	*/
    template<typename pixel_t> void fillParticles3D(Image<pixel_t>& image, const vector<string>& headers, const vector<vector<double> >& results, int condi, double multiplier, Connectivity conn)
	{
		size_t xi = getColumnIndex("X [pixel]", headers);
		size_t yi = getColumnIndex("Y [pixel]", headers);
		size_t zi = getColumnIndex("Z [pixel]", headers);

        
		size_t counter = 0;
		#pragma omp parallel for
		for(coord_t n = 0; n < (coord_t)results.size(); n++)
		{
			pixel_t value = pixelRound<pixel_t>(results[n][condi] * multiplier);

			Vec3c pos;
			pos.x = ::round(results[n][xi]);
			pos.y = ::round(results[n][yi]);
			pos.z = ::round(results[n][zi]);

			floodfill<pixel_t>(image, pos, value, value, conn);

			showThreadProgress(counter, results.size());
		}
	}

	/**
	* Fill all particles with given color.
	* @param image The image.
	* @param results Particle analysis results vector.
	* @param fillColor Fill color.
	*/
    template<typename pixel_t> void fillParticles3D(Image<pixel_t>& image, const vector<string>& headers, const vector<vector<double> >& results, pixel_t fillColor, Connectivity conn)
	{
		size_t xi = getColumnIndex("X [pixel]", headers);
		size_t yi = getColumnIndex("Y [pixel]", headers);
		size_t zi = getColumnIndex("Z [pixel]", headers);

		size_t counter = 0;
		#pragma omp parallel for
		for(coord_t n = 0; n < (coord_t)results.size(); n++)
		{
			Vec3c pos;
			pos.x = ::round(results[n][xi]);
			pos.y = ::round(results[n][yi]);
			pos.z = ::round(results[n][zi]);

			floodfill<pixel_t>(image, pos, fillColor, fillColor, conn);

			showThreadProgress(counter, results.size());
		}
	}

	///**
	//Gets value of left side of ellipsoid equation at p for ellipsoid located at c that has
	//semi-axis lengths l1, l2 and l3	and whose semi-axis orientations are given by (phiN, thetaN), N=1..3.
	//Points inside the ellipsoid are characterized by return value <= 1.
	//*/
	//inline double getEllipsoidFunctionValue(const Vec3d& p,
	//							const Vec3d& c,
	//							double l1, double l2, double l3,
	//							double phi1, double theta1,
	//							double phi2, double theta2,
	//							double phi3, double theta3)
	//{
	//	// Make the ellipsoid centered to origin
	//	Vec3d pdot = p - c;

	//	// Rotate pdot such that ellipsoid axes are aligned with coordinate axes
	//	Vec3d u1 = tocartesian(1.0, phi1, theta1);
	//	Vec3d u2 = tocartesian(1.0, phi2, theta2);
	//	Vec3d u3 = tocartesian(1.0, phi3, theta3);

	//	Matrix3x3d R(u1.x, u2.x, u3.x,
	//				 u1.y, u2.y, u3.y,
	//				 u1.z, u2.z, u3.z);

	//	R.transpose();
	//	pdot = R * pdot;

	//	//Matrix3x3d Ri;
	//	//R.inverse(Ri);
	//	//pdot = Ri * pdot;

	//	// Use ellipsoid equation
	//	double f = (pdot.x * pdot.x) / (l1 * l1) + (pdot.y * pdot.y) / (l2 * l2) + (pdot.z * pdot.z) / (l3 * l3);

	//	return f;
	//}

	///**
	//Test if point p is in ellipsoid located at c that has semi-axis lengths l1, l2 and l3
	//and whose semi-axis orientations are given by (phiN, thetaN), N=1..3.
	//*/
	//inline bool isInEllipsoid(const Vec3d& p,
	//							const Vec3d& c,
	//							double l1, double l2, double l3,
	//							double phi1, double theta1,
	//							double phi2, double theta2,
	//							double phi3, double theta3)
	//{
	//	double f = getEllipsoidFunctionValue(p, c, l1, l2, l3, phi1, theta1, phi2, theta2, phi3, theta3);
	//	return f <= 1;
	//}

	///**
	//Test if point p is in ellipsoid located at c that has semi-axis lengths l1, l2 and l3
	//and that can be transformed to axis-aligned ellipsoid by rotating with rotation matrix Rinv
	//*/
	//inline bool isInEllipsoid(const Vec3d& p,
	//							const Vec3d& c,
	//							double l1, double l2, double l3,
	//							const Matrix3x3d& Rinv)
	//{
	//	// Make the ellipsoid centered to origin
	//	Vec3d pdot = p - c;

	//	// Rotate pdot such that ellipsoid axes are aligned with coordinate axes
	//	pdot = Rinv * pdot;

	//	// Use ellipsoid equation
	//	double f = (pdot.x * pdot.x) / (l1 * l1) + (pdot.y * pdot.y) / (l2 * l2) + (pdot.z * pdot.z) / (l3 * l3);

	//	return f <= 1;
	//}

	///**
	//Draws single ellipsoid.
	//@param pMask Positions where *pMask == 0 are not filled.
	//@return Count of filled pixels.
	//*/
	//template<typename T, typename Tmask> unsigned long long drawEllipsoid(Image<T>& image, const Vec3d& pos,
	//																		double l1, double l2, double l3,
	//																		double phi1, double theta1,
	//																		double phi2, double theta2,
	//																		double phi3, double theta3,
	//																		T color,
	//																		Image<Tmask>* pMask)
	//{
	//	double maxl = mathutils::max(l1, l2);
	//	maxl = mathutils::max(maxl, l3);

	//	Vec3c minPos = round(pos - Vec3d(maxl, maxl, maxl));
	//	Vec3c maxPos = round(pos + Vec3d(maxl, maxl, maxl));
	//	
	//	Vec3c dimensions((int)image.getDimension(0), (int)image.getDimension(1), (int)image.getDimension(2));
	//	clamp(minPos, Vec3c(0, 0, 0), dimensions);
	//	clamp(maxPos, Vec3c(0, 0, 0), dimensions);

	//	Vec3d u1 = tocartesian(1.0, phi1, theta1);
	//	Vec3d u2 = tocartesian(1.0, phi2, theta2);
	//	Vec3d u3 = tocartesian(1.0, phi3, theta3);

	//	Matrix3x3d R(u1.x, u2.x, u3.x,
	//				 u1.y, u2.y, u3.y,
	//				 u1.z, u2.z, u3.z);

	//	R.transpose();
	//	
	//	unsigned long long filledCount = 0;
	//	#pragma omp parallel for reduction(+:filledCount)
	//	for(coord_t z = minPos.z; z < maxPos.z; z++)
	//	{
	//		for(size_t y = minPos.y; y < maxPos.y; y++)
	//		{
	//			for(size_t x = minPos.x; x < maxPos.x; x++)
	//			{
	//				if(pMask == 0 || pMask->getPixel((coord_t)x, (coord_t)y, (coord_t)z) != 0)
	//				{
	//					if(isInEllipsoid(Vec3d((double)x, (double)y, (double)z), pos, l1, l2, l3, R))
	//					{
	//						image.setPixel((coord_t)x, (coord_t)y, (coord_t)z, color);
	//						filledCount++;
	//					}
	//				}
	//			}
	//		}
	//	}
	//	
	//	//line<T>(image, pos, pos + tocartesian(l1, phi1, theta1), 200);
	//	//line<T>(image, pos, pos + tocartesian(l2, phi2, theta2), 220);
	//	//line<T>(image, pos, pos + tocartesian(l3, phi3, theta3), 240);

	//	return filledCount;
	//}

	///**
	//Draws bounding PCA ellipsoids.
	//@param pMask Pointer to mask image. Only pixels where mask != 0 will be colored.
	//@param pResults Pointer to results array. If not null, count of filled pixels for each bounding ellipsoid is added to the array.
	//*/
	//template<typename T, typename Tmask> void drawEllipsoids(Image<T>& image, const vector<string>& headers, const vector<vector<double> >& results, T ellipsoidColor,
	//										Image<Tmask>* pMask = 0, vector<unsigned long long>* pResults = 0)
	//											/*int cxi = 5, int cyi = 6, int czi = 7,
	//											int l1i = 9, int l2i = 10, int l3i = 11,
	//											int phi1i = 12, int theta1i = 13,
	//											int phi2i = 14, int theta2i = 15,
	//											int phi3i = 16, int theta3i = 17)*/
	//{
	//	size_t cxi = getColumnIndex("CX [pixel]", headers);
	//	size_t cyi = getColumnIndex("CY [pixel]", headers);
	//	size_t czi = getColumnIndex("CZ [pixel]", headers);
	//	size_t l1i = getColumnIndex("l1 [pixel]", headers);
	//	size_t l2i = getColumnIndex("l2 [pixel]", headers);
	//	size_t l3i = getColumnIndex("l3 [pixel]", headers);
	//	size_t phi1i = getColumnIndex("phi1 [rad]", headers);
	//	size_t theta1i = getColumnIndex("theta1 [rad]", headers);
	//	size_t phi2i = getColumnIndex("phi2 [rad]", headers);
	//	size_t theta2i = getColumnIndex("theta2 [rad]", headers);
	//	size_t phi3i = getColumnIndex("phi3 [rad]", headers);
	//	size_t theta3i = getColumnIndex("theta3 [rad]", headers);
	//	size_t d1i = getColumnIndex("d1 [pixel]", headers);
	//	size_t d2i = getColumnIndex("d2 [pixel]", headers);
	//	size_t d3i = getColumnIndex("d3 [pixel]", headers);
	//	size_t boundingScalei = getColumnIndex("bounding scale", headers);

	//	Vec3d pos;
	//	for(size_t n = 0; n < results.size(); n++)
	//	{
	//		pos.x = results[n][cxi];
	//		pos.y = results[n][cyi];
	//		pos.z = results[n][czi];

	//		double l1 = results[n][l1i];
	//		double l2 = results[n][l2i];
	//		double l3 = results[n][l3i];

	//		double phi1 = results[n][phi1i];
	//		double phi2 = results[n][phi2i];
	//		double phi3 = results[n][phi3i];

	//		double theta1 = results[n][theta1i];
	//		double theta2 = results[n][theta2i];
	//		double theta3 = results[n][theta3i];

	//		double d1 = results[n][d1i];
	//		double d2 = results[n][d2i];
	//		double d3 = results[n][d3i];

	//		double boundingScale = results[n][boundingScalei];

	//		// Scale the ellipsoid so that it becomes bounding ellipsoid
	//		l1 = boundingScale * sqrt(l1);
	//		l2 = boundingScale * sqrt(l2);
	//		l3 = boundingScale * sqrt(l3);

	//		unsigned long long count = 0;

	//		// Only draw if the ellipsoid is sane!
	//		if(!isinf(l1) && !isnan(l1) && l1 != 0 &&
	//			!isinf(l2) && !isnan(l2) && l2 != 0 &&
	//			!isinf(l3) && !isnan(l3) && l3 != 0)
	//			count = drawEllipsoid<T, Tmask>(image, pos, l1, l2, l3, phi1, theta1, phi2, theta2, phi3, theta3, ellipsoidColor, pMask);

	//		if(pResults != 0)
	//			pResults->push_back(count);

	//		utilities::showProgress(n, results.size());
	//	}
	//}

	///**
	//Draws single sphere
	//@param pMask Positions where *pMask == 0 are not filled.
	//@return Count of filled pixels.
	//*/
	//template<typename T, typename COORDT, typename Tmask> unsigned long long drawSphere(Image<T>& image, const Sphere<COORDT>& sphere, T color, Image<T>* pMask)
	//{
	//	Vec3c minPos = round(sphere.center - sphere.radius * Vec3d(1, 1, 1));
	//	Vec3c maxPos = round(sphere.center + sphere.radius * Vec3d(1, 1, 1));

	//	Vec3c dimensions((int)image.getDimension(0), (int)image.getDimension(1), (int)image.getDimension(2));
	//	clamp(minPos, Vec3c(0, 0, 0), dimensions);
	//	clamp(maxPos, Vec3c(0, 0, 0), dimensions);

	//	unsigned long long filledCount = 0;
	//	#pragma omp parallel for reduction(+:filledCount)
	//	for(coord_t z = minPos.z; z < maxPos.z; z++)
	//	{
	//		for(size_t y = minPos.y; y < maxPos.y; y++)
	//		{
	//			for(size_t x = minPos.x; x < maxPos.x; x++)
	//			{
	//				if(pMask == 0 || pMask->getPixel((coord_t)x, (coord_t)y, (coord_t)z) != 0)
	//				{
	//					if(sphere.contains(Vec3<COORDT>((COORDT)x, (COORDT)y, (COORDT)z)))
	//					{
	//						image.setPixel((coord_t)x, (coord_t)y, (coord_t)z, color);
	//						filledCount++;
	//					}
	//				}
	//			}
	//		}
	//	}

	//	return filledCount;
	//}

	///**
	//Draws bounding spheres.
	//@param pMask Pointer to mask image. Only pixels where mask != 0 will be colored.
	//@param pResults Pointer to results array. If not null, contains count of filled pixels for each bounding sphere.
	//*/
	//template<typename T, typename Tmask> void drawBoundingSpheres(Image<T>& image, const vector<string>& headers, const vector<vector<double> >& results, T sphereColor,
	//												Image<Tmask>* pMask = 0, vector<unsigned long long>* pResults = 0)
	//{
	//	//int cxi = 25, int cyi = 26, int czi = 27, int radiusi = 28
	//	size_t cxi = getColumnIndex("bcenter X [pixel]", headers);
	//	size_t cyi = getColumnIndex("bcenter Y [pixel]", headers);
	//	size_t czi = getColumnIndex("bcenter Z [pixel]", headers);
	//	size_t radiusi = getColumnIndex("bradius [pixel]", headers);

	//	Vec3d pos;
	//	Sphere<double> sphere;
	//	for(size_t n = 0; n < results.size(); n++)
	//	{
	//		pos.x = results[n][cxi];
	//		pos.y = results[n][cyi];
	//		pos.z = results[n][czi];

	//		double radius = results[n][radiusi];
	//		
	//		sphere.center = pos;
	//		sphere.radius = radius;

	//		unsigned long long count = drawSphere<T, double, Tmask>(image, sphere, sphereColor, pMask);
	//		if(pResults != 0)
	//			pResults->push_back(count);

	//		utilities::showProgress(n, results.size());
	//	}
	//}

    /**
    Contains particle analyzers.
    */
    namespace analyzers
    {
		/**
		Returns the coordinates of the first point on the particle.
		*/
        template<class POINT, typename PIXELTYPE> class PointCoords2D : public Analyzer<POINT, PIXELTYPE>
        {
        public:
            virtual vector<string> getTitles()
            {
				vector<string> labels;
				labels.push_back("X [pixel]");
				labels.push_back("Y [pixel]");
                return labels;
            }

            virtual vector<double> analyze(const PIXELTYPE pixelValue, const vector<POINT >& points)
            {
				vector<double> results;
				results.push_back(points[0][0]);
				results.push_back(points[0][1]);
                return results;
            }
        };

        /**
		Returns the color of the point where the particle was found.
		*/
        template<class POINT, typename PIXELTYPE> class Color : public Analyzer<POINT, PIXELTYPE>
        {
        public:
            virtual vector<string> getTitles()
            {
				vector<string> labels;
				labels.push_back("Color");
                return labels;
            }

            virtual vector<double> analyze(const PIXELTYPE pixelValue, const vector<POINT >& points)
            {
				vector<double> results;
				results.push_back((double)pixelValue);
                return results;
            }
        };

        /**
		Returns user-specified value.
		*/
        template<class POINT, typename PIXELTYPE> class ArbitraryValue : public Analyzer<POINT, PIXELTYPE>
        {
        private:
            const string title;
            double value;
        public:
            ArbitraryValue(const string& title, const double value) : title(title), value(value)
            {
            }

            virtual vector<string> getTitles()
            {
				vector<string> labels;
				labels.push_back(title);
                return labels;
            }

            virtual vector<double> analyze(const PIXELTYPE pixelValue, const vector<POINT >& points)
            {
				vector<double> results;
				results.push_back(value);
                return results;
            }

            void setValue(double newValue)
            {
                value = newValue;
            }
        };

		/**
		Returns the coordinates of the first point on the particle.
		*/
        template<class POINT, typename PIXELTYPE> class PointCoords3D : public Analyzer<POINT, PIXELTYPE>
        {
        public:
            virtual vector<string> getTitles()
            {
				vector<string> labels;
				labels.push_back("X [pixel]");
				labels.push_back("Y [pixel]");
				labels.push_back("Z [pixel]");
                return labels;
            }

            virtual vector<double> analyze(const PIXELTYPE pixelValue, const vector<POINT >& points)
            {
				vector<double> results;
				results.push_back((double)points[0][0]);
				results.push_back((double)points[0][1]);
				results.push_back((double)points[0][2]);
                return results;
            }
        };

		/**
		Calculates the volume of the particle.
		*/
        template<class POINT, typename PIXELTYPE> class Volume : public Analyzer<POINT, PIXELTYPE>
        {
        public:
            virtual vector<string> getTitles()
            {
				vector<string> labels;
				labels.push_back("Volume [pixel]");
                return labels;
            }

            virtual vector<double> analyze(const PIXELTYPE pixelValue, const vector<POINT >& points)
            {
				vector<double> results;
				results.push_back((double)points.size());
				return results;
            }
        };

		///**
		//Calculates the surface area of a 3D object using Marching Cubes.
		//*/
		//template<class POINT, typename PIXELTYPE> class SurfaceArea3D : public Analyzer<POINT, PIXELTYPE>
		//{
		//public:
  //          virtual vector<string> getTitles()
  //          {
		//		vector<string> labels;
		//		labels.push_back("Surface area [pixel^2]");
  //              return labels;
  //          }

  //          virtual vector<double> analyze(const PIXELTYPE pixelValue, const vector<POINT >& points)
  //          {
		//		// Calculate bounding box for the points.
		//		POINT min;
		//		POINT max;
		//		for(size_t n = 0; n < points.size(); n++)
		//		{
		//			POINT p = points[n];
		//			if(p.x < min.x)
		//				min.x = p.x;
		//			if(p.y < min.y)
		//				min.y = p.y;
		//			if(p.z < min.z)
		//				min.z = p.z;

		//			if(p.x > max.x)
		//				max.x = p.x;
		//			if(p.y > max.y)
		//				max.y = p.y;
		//			if(p.z > max.z)
		//				max.z = p.z;
		//		}

		//		coord_t w = (coord_t)round(max.x) - (coord_t)round(min.x) + 2;
		//		coord_t h = (coord_t)round(max.y) - (coord_t)round(min.y) + 2;
		//		coord_t d = (coord_t)round(max.z) - (coord_t)round(min.z) + 2;

		//		// Plot the points into temporary image
		//		Image<uint8_t> block(w, h, d);

		//		for(size_t n = 0; n < points.size(); n++)
		//		{
		//			POINT p = points[n];
		//			p.x++;
		//			p.y++;
		//			p.z++;
		//			block.setPixel(p, 255);
		//		}

		//		// Calculate area of the particle using Marching Cubes
		//		double area = getMarchingCubesArea<uint8_t>(block, 128);

		//		vector<double> results;
		//		results.push_back(area);
		//		return results;
  //          }
		//};

		///**
		//Calculates measures of the convex hull of the particle.
		//*/
  //      template<class POINT, typename PIXELTYPE> class ConvexHull2D : public Analyzer<POINT, PIXELTYPE>
  //      {
		//private:
		//	/**
		//	Calculates area of convex 2D polygon.
		//	*/
		//	double convexPolyArea(const vector<POINT>& v)
		//	{
		//		double area = 0;
		//		for(size_t n = 1; n < v.size()-1; n++)
		//		{
		//			POINT A = v[0];
		//			POINT B = v[n];
		//			POINT C = v[n + 1];
		//			area += abs((A.x-C.x)*(B.y-A.y)-(A.x-B.x)*(C.y-A.y));
		//		}
		//		area *= 0.5;

		//		return area;
		//	}

		//	/**
		//	Calculates edge length of polygon.
		//	*/
		//	double polyPerimeter(const vector<POINT>& v)
		//	{
		//		double p = 0;
		//		for(size_t n = 0; n < v.size() - 1; n++)
		//			p += (v[n + 1] - v[n]).norm();
		//		p += (v[0] - v[v.size() - 1]).norm();
		//		return p;
		//	}

  //      public:
  //          virtual vector<string> getTitles()
  //          {
		//		vector<string> labels;
		//		labels.push_back("Convex area [pixel^2]");
		//		labels.push_back("Convex perimeter [pixel]");
  //              return labels;
  //          }

  //          virtual vector<double> analyze(const PIXELTYPE pixelValue, const vector<POINT >& points)
  //          {
		//		vector<double> results;

		//		// Calculate convex hull from edge points so that area of convex hull is always >= nonconvex area.
		//		vector<POINT > edgePoints;
		//		for(size_t n = 0; n < points.size(); n++)
		//		{
		//			POINT p = points[n];
		//			int ix = (int)p.x;
		//			int iy = (int)p.y;

		//			edgePoints.push_back(POINT(ix, iy));
		//			edgePoints.push_back(POINT(ix+1, iy));
		//			edgePoints.push_back(POINT(ix, iy+1));
		//			edgePoints.push_back(POINT(ix+1, iy+1));
		//		}

		//		vector<POINT> hull;
		//		convexHull2D(edgePoints, hull);
		//		double area = convexPolyArea(hull);
		//		double perimeter = polyPerimeter(hull);

		//		results.push_back(area);
		//		results.push_back(perimeter);
		//		return results;
  //          }
  //      };

		///**
		//Calculates count of holes and their area by using labeling algorithm on the background.
		//*/
  //      template<class POINT, typename PIXELTYPE> class HoleSizeCount2D : public Analyzer<POINT, PIXELTYPE>
  //      {
		//private:
		//	

  //      public:
  //          virtual vector<string> getTitles()
  //          {
		//		vector<string> labels;
		//		labels.push_back("Hole count [1]");
		//		labels.push_back("Total hole area [pixel^2]");
  //              return labels;
  //          }

  //          virtual vector<double> analyze(const PIXELTYPE pixelValue, const vector<POINT >& points)
  //          {
		//		// Calculate bounds
		//		POINT min = points[0];
		//		POINT max = points[0];
		//		for(size_t n = 1; n < points.size(); n++)
		//		{
		//			if(points[n].x < min.x)
		//				min.x = points[n].x;
		//			if(points[n].y < min.y)
		//				min.y = points[n].y;
		//			if(points[n].x > max.x)
		//				max.x = points[n].x;
		//			if(points[n].y > max.y)
		//				max.y = points[n].y;
		//		}

		//		min.x--;
		//		min.y--;
		//		max.x++;
		//		max.y++;

		//		coord_t w = max.x - min.x + 1;
		//		coord_t h = max.y - min.y + 1;

		//		Image<PIXELTYPE> img(w, h);

		//		// Plot the points
		//		for(size_t n = 0; n < points.size(); n++)
		//		{
		//			img.setPixel(Vec2c(points[n].x - min.x, points[n].y - min.y), pixelValue);
		//		}

		//		// Fill outside
		//		floodfill2D(img, Vec2c(0, 0), pixelValue, pixelValue);

		//		// Go through the image and flood fill holes
		//		double holeCount = 0;
		//		double totalHoleArea = 0;
		//		for(coord_t y = 0; y < h; y++)
		//		{
		//			for(coord_t x = 0; x < w; x++)
		//			{
		//				if(img.getPixel(Vec2c(x, y)) != pixelValue)
		//				{
		//					// Hole found, measure area by filling it.
		//					holeCount++;
		//					vector<Vec2c> filledPoints;
		//					floodfill2D(img, Vec2c(x, y), pixelValue, pixelValue, Connectivity::Nearest, &filledPoints);
		//					totalHoleArea += filledPoints.size();
		//				}
		//			}
		//		}

		//		vector<double> results;
		//		results.push_back(holeCount);
		//		results.push_back(totalHoleArea);
		//		return results;
  //          }
  //      };

		/**
		Tests whether the particle touches image edge.
		*/
		template<class POINT, typename PIXELTYPE> class IsOnEdge : public Analyzer<POINT, PIXELTYPE>
		{
		private:
			Vec3c dimensions;
		public:
			/**
			Constructor
			@param dimensions Dimensions of the image.
			*/
			IsOnEdge(const Vec3c& dimensions)
			{
				this->dimensions = dimensions;
			}

			virtual vector<string> getTitles()
			{
				vector<string> titles;
				titles.push_back("Is on edge [0/1]");
				return titles;
			}

			virtual vector<double> analyze(const PIXELTYPE pixelValue, const vector<POINT >& points)
			{
				vector<double> results;
				results.push_back(0.0);

				// Traverse all the points on the particle
				for(size_t n = 0; n < points.size(); n++)
				{
					// Test each dimension of the point for image edge status
					POINT p = points[n];
					for(size_t dimension = 0; dimension < p.size(); dimension++)
					{
						if(p[dimension] <= 0 || p[dimension] >= dimensions[dimension] - 1)
						{
							results[0] = 1.0;
							n = points.size();
							break;
						}
					}
				}

				return results;
			}
		};

		///**
		//Calculates principal components.
		//*/
		//template<class POINT, typename PIXELTYPE> class PCA2D : public Analyzer<POINT, PIXELTYPE>
		//{
		//	virtual vector<string> getTitles()
		//	{
		//		vector<string> titles;
		//		titles.push_back("CX [pixel]");
		//		titles.push_back("CY [pixel]");
		//		titles.push_back("e");
		//		titles.push_back("l1 [pixel]");
		//		titles.push_back("l2 [pixel]");
		//		titles.push_back("theta1 [rad]");
		//		titles.push_back("theta2 [rad]");
		//		titles.push_back("rmax [pixel]");	// Overall maximum radius from center point
		//		titles.push_back("d1 [pixel]");     // Maximum diameter in direction of first principal component.
		//		titles.push_back("d2 [pixel]");     // Maximum diameter in direction of second principal component.
		//		return titles;
		//	}

		//	virtual vector<double> analyze(const PIXELTYPE pixelValue, const vector<POINT >& points)
		//	{
		//		vector<double> results;

		//		// Calculate centroid
		//		double mx = 0.0;
		//		double my = 0.0;

		//		for(size_t n = 0; n < points.size(); n++)
		//		{
		//			mx += points[n][0];
		//			my += points[n][1];
		//		}
		//		mx /= points.size();
		//		my /= points.size();

		//		results.push_back(mx);
		//		results.push_back(my);

		//		// Calculate covariance matrix
		//		double sumxx = 0.0;
		//		double sumyy = 0.0;
		//		double sumxy = 0.0;

		//		for(size_t n = 0; n < points.size(); n++)
		//		{
		//			double currX = points[n][0] - mx;
		//			double currY = points[n][1] - my;

		//			sumxx += currX * currX;
		//			sumyy += currY * currY;

		//			sumxy += currX * currY;
		//		}

		//		double K = (double)points.size();
		//		sumxx /= K;
		//		sumyy /= K;
		//		sumxy /= K;

		//		double S[2][2];

		//		S[0][0] = sumxx;
		//		S[0][1] = sumxy;

		//		S[1][0] = sumxy;
		//		S[1][1] = sumyy;

		//		// Symmetric matrix A => eigenvectors in columns of V, corresponding eigenvalues in d.
		//		double V[2][2];
		//		double d[2];
		//		double lambda1 = 0;
		//		double lambda2 = 0;
		//		if(eigen_decomposition2(S, V, d))
		//		{
		//			lambda1 = d[1];
		//			lambda2 = d[0];
		//		}

		//		if(lambda1 < lambda2)
		//			throw ITLException(string("Eigenvalues are unsorted: ") + toString(lambda1) + string(" and ") + toString(lambda2));

		//		double e = 0.0;
		//		if(abs(lambda1) > 1e-6)
		//			e = sqrt(1.0 - (lambda2 * lambda2) / (lambda1 * lambda1));

		//		results.push_back(e);
		//		results.push_back(lambda1);
		//		results.push_back(lambda2);

		//		// Remember to reverse order of vectors here as was done for eigenvalues above.
		//		results.push_back(atan2(V[1][1], V[0][1]));
		//		results.push_back(atan2(V[1][0], V[0][0]));

		//		// Calculate maximum radius and maximum diameters overall and in directions of eigenvectors.
		//		Vec2d t1dir(V[0][1], V[1][1]);
		//		Vec2d t2dir(V[0][0], V[1][0]);
		//		t1dir.normalize();
		//		t2dir.normalize();
		//		double maxr = 0.0;
		//		double mind1 = 0.0;
		//		double maxd1 = 0.0;
		//		double mind2 = 0.0;
		//		double maxd2 = 0.0;
		//		for(size_t n = 0; n < points.size(); n++)
		//		{
		//			Vec2d currx(points[n][0] - mx, points[n][1] - my);
		//			double r = currx.norm();
		//			if(r > maxr)
		//				maxr = r;

		//			r = t1dir.dot(currx);
		//			if(r > maxd1)
		//				maxd1 = r;
  //                  if(r < mind1)
  //                      mind1 = r;

		//			r = t2dir.dot(currx);
		//			if(r > maxd2)
		//				maxd2 = r;
  //                  if(r < mind2)
  //                      mind2 = r;
		//		}

		//		results.push_back(maxr);
		//		results.push_back(maxd1 - mind1);
		//		results.push_back(maxd2 - mind2);

		//		return results;
		//	}
		//};

		///**
		//Calculates principal components.
		//*/
		//template<class POINT, typename PIXELTYPE> class PCA3D : public Analyzer<POINT, PIXELTYPE>
		//{
		//	virtual vector<string> getTitles()
		//	{
		//		vector<string> titles;
		//		titles.push_back("CX [pixel]");
		//		titles.push_back("CY [pixel]");
		//		titles.push_back("CZ [pixel]");
		//		titles.push_back("e (meridional)");
		//		titles.push_back("l1 [pixel]");
		//		titles.push_back("l2 [pixel]");
		//		titles.push_back("l3 [pixel]");
		//		titles.push_back("phi1 [rad]");
		//		titles.push_back("theta1 [rad]");
		//		titles.push_back("phi2 [rad]");
		//		titles.push_back("theta2 [rad]");
		//		titles.push_back("phi3 [rad]");
		//		titles.push_back("theta3 [rad]");
		//		titles.push_back("rmax [pixel]");	// Overall maximum radius from center point
		//		titles.push_back("d1 [pixel]");     // Maximum diameter in direction of first principal component.
		//		titles.push_back("d2 [pixel]");     // Maximum diameter in direction of second principal component.
		//		titles.push_back("d3 [pixel]");     // Maximum diameter in direction of third principal component.
		//		titles.push_back("bounding scale"); // Semi-axis of bounding ellipsoid are calculated as (bounding scale)*l1, (bounding scale)*l2, etc.
		//		return titles;
		//	}

		//	virtual vector<double> analyze(const PIXELTYPE pixelValue, const vector<POINT >& points)
		//	{
		//		vector<double> results;

		//		// Calculate centroid
		//		double mx = 0.0;
		//		double my = 0.0;
		//		double mz = 0.0;

		//		for(size_t n = 0; n < points.size(); n++)
		//		{
		//			mx += points[n][0];
		//			my += points[n][1];
		//			mz += points[n][2];
		//		}
		//		mx /= points.size();
		//		my /= points.size();
		//		mz /= points.size();

		//		results.push_back(mx);
		//		results.push_back(my);
		//		results.push_back(mz);

		//		// Calculate covariance matrix
		//		double sumxx = 0.0;
		//		double sumyy = 0.0;
		//		double sumzz = 0.0;
		//		double sumxy = 0.0;
		//		double sumxz = 0.0;
		//		double sumyz = 0.0;

		//		for(size_t n = 0; n < points.size(); n++)
		//		{
		//			double currX = points[n][0] - mx;
		//			double currY = points[n][1] - my;
		//			double currZ = points[n][2] - mz;

		//			sumxx += currX * currX;
		//			sumyy += currY * currY;
		//			sumzz += currZ * currZ;

		//			sumxy += currX * currY;
		//			sumxz += currX * currZ;
		//			sumyz += currY * currZ;
		//		}

		//		double K = (double)points.size();
		//		sumxx /= K;
		//		sumyy /= K;
		//		sumzz /= K;
		//		sumxy /= K;
		//		sumxz /= K;
		//		sumyz /= K;

		//		//Matrix3d S;
		//		//SelfAdjointEigenSolver<Matrix3d> solver;

		//		//S(0, 0) = sumxx;
		//		//S(0, 1) = sumxy;
		//		//S(0, 2) = sumxz;

		//		//S(1, 0) = sumxy;
		//		//S(1, 1) = sumyy;
		//		//S(1, 2) = sumyz;

		//		//S(2, 0) = sumxz;
		//		//S(2, 1) = sumyz;
		//		//S(2, 2) = sumzz;

		//		//solver.compute(S);

		//		//double lambda1 = 0;
		//		//double lambda2 = 0;
		//		//double lambda3 = 0;

		//		////Matrix3d V;

		//		//if(solver.info() == Success)
		//		//{
		//		//	lambda1 = solver.eigenvalues()[2];
		//		//	lambda2 = solver.eigenvalues()[1];
		//		//	lambda3 = solver.eigenvalues()[0];

		//		//	//V = solver.eigenvectors();
		//		//}
		//		//else
		//		//{
		//		//	//V.setZero();
		//		//}

		//		double S[3][3];

		//		S[0][0] = sumxx;
		//		S[0][1] = sumxy;
		//		S[0][2] = sumxz;

		//		S[1][0] = sumxy;
		//		S[1][1] = sumyy;
		//		S[1][2] = sumyz;

		//		S[2][0] = sumxz;
		//		S[2][1] = sumyz;
		//		S[2][2] = sumzz;

		//		// Symmetric matrix A => eigenvectors in columns of V, corresponding eigenvalues in d.
		//		// See also analysis of structure tensor - approx the same code is there, too.
		//		double V[3][3];
		//		double d[3];
		//		eigen_decomposition(S, V, d);
		//		double lambda1 = d[2];
		//		double lambda2 = d[1];
		//		double lambda3 = d[0];

		//		if(lambda1 < lambda2 || lambda2 < lambda3 || lambda1 < lambda3)
		//			throw ITLException("Eigenvalues are unsorted.");

		//		double r1, phi1, theta1;
		//		double r2, phi2, theta2;
		//		double r3, phi3, theta3;

		//		double V1mult = V[0][0] < 0 ? -1 : 1;
		//		double V2mult = V[0][1] < 0 ? -1 : 1;
		//		double V3mult = V[0][2] < 0 ? -1 : 1;
		//		math::tospherical(V1mult * V[0][0], V1mult * V[1][0], V1mult * V[2][0], r3, phi3, theta3);
		//		math::tospherical(V2mult * V[0][1], V2mult * V[1][1], V2mult * V[2][1], r2, phi2, theta2);
		//		math::tospherical(V3mult * V[0][2], V3mult * V[1][2], V3mult * V[2][2], r1, phi1, theta1);

		//		double e = 0.0;
		//		if(abs(lambda1) > 1e-6)
		//			e = sqrt(1.0 - (lambda3 * lambda3) / (lambda1 * lambda1));

		//		results.push_back(e);

		//		results.push_back(lambda1);
		//		results.push_back(lambda2);
		//		results.push_back(lambda3);

		//		results.push_back(phi1);
		//		results.push_back(theta1);
		//		
		//		results.push_back(phi2);
		//		results.push_back(theta2);
		//		
		//		results.push_back(phi3);
		//		results.push_back(theta3);


		//		// Calculate maximum radius
		//		double maxr = 0;
		//		Vec3d t3dir(V1mult * V[0][0], V1mult * V[1][0], V1mult * V[2][0]);
		//		Vec3d t2dir(V2mult * V[0][1], V2mult * V[1][1], V2mult * V[2][1]);
		//		Vec3d t1dir(V3mult * V[0][2], V3mult * V[1][2], V3mult * V[2][2]);
		//		double l1 = sqrt(lambda1);
		//		double l2 = sqrt(lambda2);
		//		double l3 = sqrt(lambda3);
		//		Vec3d c(mx, my, mz);
		//		double maxd1 = 0;
		//		double mind1 = 0;
		//		double maxd2 = 0;
		//		double mind2 = 0;
		//		double maxd3 = 0;
		//		double mind3 = 0;
		//		double boundingScale = 0;
		//		for(size_t n = 0; n < points.size(); n++)
		//		{
		//			Vec3d p(points[n][0] - mx, points[n][1] - my, points[n][2] - mz);
		//			double r = p.norm();
		//			maxr = math::max(maxr, r);

		//			r = t1dir.dot(p);
		//			if(r > maxd1)
		//				maxd1 = r;
  //                  if(r < mind1)
  //                      mind1 = r;

		//			r = t2dir.dot(p);
		//			if(r > maxd2)
		//				maxd2 = r;
  //                  if(r < mind2)
  //                      mind2 = r;

		//			r = t3dir.dot(p);
		//			if(r > maxd3)
		//				maxd3 = r;
  //                  if(r < mind3)
  //                      mind3 = r;


		//			double f = getEllipsoidFunctionValue(Vec3d(points[n][0], points[n][1], points[n][2]),
		//						c,
		//						l1, l2, l3,
		//						phi1, theta1,
		//						phi2, theta2,
		//						phi3, theta3);
		//			double scale = sqrt(f);
		//			if(!isnan(scale) && !isinf(scale))
		//				boundingScale = math::max(boundingScale, scale);
		//		}
		//		results.push_back(maxr);
		//		results.push_back(maxd1 - mind1);
		//		results.push_back(maxd2 - mind2);
		//		results.push_back(maxd3 - mind3);
		//		results.push_back(boundingScale);

		//		return results;
		//	}
		//};

		///**
		//Returns measurements of histogram (from another image) at the points of the particle.
		//*/
  //      template<class POINT, typename PIXELTYPE, typename PTYPE2> class Histogram : public Analyzer<POINT, PIXELTYPE>
  //      {
  //      private:
  //          Image<PTYPE2>& anotherImage;

		//	vector<vector<float> > pixelValues;

  //      public:

		//	const vector<vector<float> >& getExtraResults() const
		//	{
		//		return pixelValues;
		//	}

  //          Histogram(Image<PTYPE2>& anotherImage) : anotherImage(anotherImage)
  //          {
  //          }

  //          virtual vector<string> getTitles()
  //          {
		//		vector<string> labels;
		//		labels.push_back("Mean");
		//		labels.push_back("Stdev");
		//		labels.push_back("Min");
		//		labels.push_back("Max");
  //              return labels;
  //          }

  //          virtual vector<double> analyze(const PIXELTYPE pixelValue, const vector<POINT >& points)
  //          {
		//		pixelValues.push_back(vector<float>());
		//		vector<float>& vec = pixelValues[pixelValues.size()-1];
		//		vec.reserve(points.size());

		//		// Calculate statistics
		//		double mean = 0;
		//		double stdev = 0;
		//		double min = numeric_limits<double>::infinity();
		//		double max = -numeric_limits<double>::infinity();

		//		
		//		for(size_t n = 0; n < points.size(); n++)
		//		{
		//			double pix = anotherImage.getPixel(points[n]);
		//			vec.push_back((float)pix);
		//			mean += pix;
		//			if(pix < min)
		//				min = pix;
		//			if(pix > max)
		//				max = pix;
		//		}

		//		mean /= (double)points.size();

		//		// Calculate standard deviation
		//		for(size_t n = 0; n < points.size(); n++)
		//		{
		//			double pix = anotherImage.getPixel(points[n]);
		//			stdev += (pix - mean) * (pix - mean);
		//		}

		//		stdev = sqrt(stdev/(double)points.size());

		//		vector<double> results;
		//		results.push_back(mean);
		//		results.push_back(stdev);
		//		results.push_back(min);
		//		results.push_back(max);
  //              return results;
  //          }
  //      };


		/**
		Calculates axis-aligned bounding box of the particle.
		*/
        template<class POINT, typename PIXELTYPE> class BoundingBox3D : public Analyzer<POINT, PIXELTYPE>
        {
		private:
			

        public:
            virtual vector<string> getTitles()
            {
				vector<string> labels;
				labels.push_back("minx");
				labels.push_back("maxx");
				labels.push_back("miny");
				labels.push_back("maxy");
				labels.push_back("minz");
				labels.push_back("maxz");
                return labels;
            }

            virtual vector<double> analyze(const PIXELTYPE pixelValue, const vector<POINT >& points)
            {
				// Calculate bounds
				POINT min = points[0];
				POINT max = points[0];
				for(size_t n = 1; n < points.size(); n++)
				{
					if(points[n].x < min.x)
						min.x = points[n].x;
					if(points[n].y < min.y)
						min.y = points[n].y;
					if(points[n].z < min.z)
						min.z = points[n].z;
					if(points[n].x > max.x)
						max.x = points[n].x;
					if(points[n].y > max.y)
						max.y = points[n].y;
					if(points[n].z > max.z)
						max.z = points[n].z;
				}

				vector<double> results;
				results.push_back(min.x);
				results.push_back(max.x);
				results.push_back(min.y);
				results.push_back(max.y);
				results.push_back(min.z);
				results.push_back(max.z);
                return results;
			}
		};

		/**
		Calculates axis-aligned bounding box of the particle.
		*/
        template<class POINT, typename PIXELTYPE> class BoundingBox2D : public Analyzer<POINT, PIXELTYPE>
        {
		private:
			

        public:
            virtual vector<string> getTitles()
            {
				vector<string> labels;
				labels.push_back("minx");
				labels.push_back("maxx");
				labels.push_back("miny");
				labels.push_back("maxy");
                return labels;
            }

            virtual vector<double> analyze(const PIXELTYPE pixelValue, const vector<POINT >& points)
            {
				// Calculate bounds
				POINT min = points[0];
				POINT max = points[0];
				for(size_t n = 1; n < points.size(); n++)
				{
					if(points[n].x < min.x)
						min.x = points[n].x;
					if(points[n].y < min.y)
						min.y = points[n].y;
					if(points[n].x > max.x)
						max.x = points[n].x;
					if(points[n].y > max.y)
						max.y = points[n].y;
				}

				vector<double> results;
				results.push_back(min.x);
				results.push_back(max.x);
				results.push_back(min.y);
				results.push_back(max.y);
                return results;
			}
		};


		/**
		Calculates minimum bounding sphere of the particle
		*/
        //template<class POINT, typename PIXELTYPE> class BoundingSphere3D : public Analyzer<POINT, PIXELTYPE>
        //{
		//private:
		//	

        //public:
        //    virtual vector<string> getTitles()
        //    {
		//		vector<string> labels;
		//		labels.push_back("bcenter X [pixel]");
		//		labels.push_back("bcenter Y [pixel]");
		//		labels.push_back("bcenter Z [pixel]");
		//		labels.push_back("bradius [pixel]");
        //       return labels;
        //    }

        //    virtual vector<double> analyze(const PIXELTYPE pixelValue, const vector<POINT >& points)
        //    {
		//		Sphere<double> bs = Sphere<double>::miniball(points);

		//		vector<double> results;
		//		results.push_back(bs.center.x);
		//		results.push_back(bs.center.y);
		//		results.push_back(bs.center.z);
		//		results.push_back(bs.radius);
        //        return results;
		//	}
		//};

		///**
		//Calculates quantities related to skeleton branches.
		//*/
		//template<class POINT, typename PIXELTYPE> class SkeletonBranch3D : public Analyzer<POINT, PIXELTYPE>
		//{
		//private:
		//	

  //      public:
  //          virtual vector<string> getTitles()
  //          {
		//		vector<string> labels;
		//		labels.push_back("pixel count [1]");
		//		labels.push_back("euclidean distance [pixel]");
  //              return labels;
  //          }

  //          virtual vector<double> analyze(const PIXELTYPE pixelValue, const vector<POINT >& points)
  //          {
		//		double distance = numeric_limits<double>::signaling_NaN();

		//		// Find two points that have only one neighbour
		//		vector<size_t> endPointIndices;
		//		for(size_t n = 0; n < points.size(); n++)
		//		{
		//			int neighbourCount = 0;
		//			for(size_t m = 0; m < points.size(); m++)
		//			{
		//				if(n != m)
		//				{
		//					if((points[n] - points[m]).norm() <= sqrt(1*1 + 1*1 + 1*1))
		//					{
		//						// Points m and n are neighbours
		//						neighbourCount++;
		//					}
		//				}
		//			}

		//			if(neighbourCount == 1)
		//				endPointIndices.push_back(n);
		//		}

		//		if(endPointIndices.size() == 2)
		//		{
		//			size_t i1 = endPointIndices[0];
		//			size_t i2 = endPointIndices[1];
		//			distance = (points[i1] - points[i2]).norm();
		//		}


		//		vector<double> results;
		//		results.push_back((double)points.size());
		//		results.push_back(distance);
  //              return results;
		//	}
		//};
    }

	/**
	Creates analyzer vector for 3D particle analysis.
	*/
	template<typename PIXELTYPE> AnalyzerVector<Vec3c, PIXELTYPE>* build3DAnalyzers(vector<size_t> dimensions)
	{
		if(dimensions.size() != 3)
			throw ITLException("This method supports only three dimensions.");

		AnalyzerVector<Vec3c, PIXELTYPE>* analyzers = new AnalyzerVector<Vec3c, PIXELTYPE>();
		analyzers->push_back(new analyzers::PointCoords3D<Vec3c, PIXELTYPE>());
		analyzers->push_back(new analyzers::IsOnEdge<Vec3c, PIXELTYPE>(dimensions));
		analyzers->push_back(new analyzers::Volume<Vec3c, PIXELTYPE>());
		//analyzers->push_back(new analyzers::PCA3D<Vec3c, PIXELTYPE>());
		//analyzers->push_back(new analyzers::SurfaceArea3D<Vec3c, PIXELTYPE>());
		analyzers->push_back(new analyzers::Color<Vec3c, PIXELTYPE>());
		analyzers->push_back(new analyzers::BoundingBox3D<Vec3c, PIXELTYPE>());
		//analyzers->push_back(new analyzers::BoundingSphere3D<Vec3c, PIXELTYPE>());

		return analyzers;
	}

	/**
	Creates analyzer vector for 2D particle analysis.
	*/
	template<typename PIXELTYPE> AnalyzerVector<Vec2c, PIXELTYPE>* build2DAnalyzers(vector<size_t> dimensions)
	{
		if(dimensions.size() != 2)
			throw ITLException("This method supports only two dimensions.");

		AnalyzerVector<Vec2c, PIXELTYPE>* analyzers = new AnalyzerVector<Vec2c, PIXELTYPE>();
		analyzers->push_back(new analyzers::PointCoords2D<Vec2c, PIXELTYPE>());
		analyzers->push_back(new analyzers::IsOnEdge<Vec2c, PIXELTYPE>(dimensions));
		analyzers->push_back(new analyzers::Volume<Vec2c, PIXELTYPE>());
		//analyzers->push_back(new analyzers::PCA2D<Vec2c, PIXELTYPE>());
		//analyzers->push_back(new analyzers::ConvexHull2D<Vec2c, PIXELTYPE>());
		///analyzers->push_back(new analyzers::HoleSizeCount2D<Vec2c, PIXELTYPE>());
		analyzers->push_back(new analyzers::BoundingBox2D<Vec2c, PIXELTYPE>());

		return analyzers;
	}
}

