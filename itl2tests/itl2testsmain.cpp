
#include "math/numberutils.h"
#include "math/matrix.h"
#include "test.h"
#include "raytrace.h"
#include "io/vol.h"
#include "utilities.h"
#include "danielsson.h"
#include "io/io.h"
#include "io/itltiff.h"
#include "sphere.h"
#include "testutils.h"
#include "io/itlpng.h"
#include "neighbourhood.h"
#include "misc.h"
#include "filters.h"
#include "transform.h"
#include "registration.h"
#include "stitching.h"
#include "histogram.h"
#include "floodfill.h"
#include "lineskeleton.h"
#include "hybridskeleton.h"
#include "traceskeleton.h"
#include "structure.h"
#include "particleanalysis.h"
#include "regionremoval.h"

#include <stack>

using namespace itl2;
using namespace math;
using namespace std;


int main()
{
	//test(math::tests::numberUtils, "Number utilities");
	//test(math::tests::matrix3x3, "3x3 matrix");
	//test(math::tests::matrix, "Matrix");

	//test(itl2::tests::image, "Image");
	//test(raw::tests::expandFilename, "Raw filename expansion");
	//test(raw::tests::raw, "Raw reader");
	//test(io::tests::readWrite, "IO read");
	//test(raw::tests::writeBlock, "Block based raw reader & writer");
	//test(raw::tests::writeBlockFast, "Optimized block based raw reader & writer");
	//test(vol::tests::volio, ".vol input/output");
	//test(itl2::png::tests::png, "Png read and write");
	//test(itl2::tiff::tests::readWrite, "Tiff read and write");

	//test(itl2::sequence::tests::match, "Matching");
	//test(itl2::sequence::tests::sequence, "Image sequence");
	//test(itl2::sequence::tests::readWriteBlock, "Image sequence block");
	//test(itl2::sequence::tests::readWriteBlockOptimization, "Image sequence block write optimization");
	//test(itl2::sequence::tests::fileFormats, "Sequence file formats");

	//test(itl2::tests::siddonProject, "Siddon algorithm");
	//test(math::tests::geometry, "Geometry");

	//test(itl2::tests::neighbourhoodTools, "neighbourhood tools");
	//test(itl2::tests::edges, "Set edges");
	//test(itl2::tests::separableOptimization, "filtering with separable optimization");

	//test(itl2::tests::gaussFilters, "Gaussian filters and derivatives");
	//test(itl2::tests::projections, "projections");
	//test(itl2::tests::fourierTransformPair, "Fourier transforms");
	//test(itl2::tests::dctPair, "DCT");
	//test(itl2::tests::bandpass, "Bandpass filtering");
	//test(itl2::tests::projections2, "projections 2");
	//test(itl2::tests::filters, "filtering");

	//test(itl2::tests::bilateral, "bilateral filter");

	//test(itl2::tests::translate, "translation");

	//test(itl2::tests::pointProcess, "point processes");
	//test(itl2::tests::pointProcessComplex, "point processes on complex numbers");
	

	//test(itl2::tests::phaseCorrelation, "phase correlation");
	//test(itl2::tests::modulo, "modulo function");
	
	//test(itl2::tests::phaseCorrelation2, "phase correlation 2 (rotation)");

	//test(itl2::tests::blockMatch1, "block match 1");
	//test(itl2::tests::blockMatch2Match, "block match 2 (match)");
	//test(itl2::tests::blockMatch2Pullback, "block match 2 (pullback)");

	//test(itl2::tests::inpaintNearest, "Inpainting");
	//test(itl2::tests::inpaintGarcia, "Inpainting (Garcia)");
	//test(itl2::tests::inpaintGarcia2, "Inpainting 2 (Garcia)");
	//test(itl2::tests::dmap1, "Distance map");

	//test(itl2::tests::buffers, "Disk mapped buffer");
	//test(itl2::tests::histogram, "Histogram");

	//test(itl2::tests::binning, "Binning");
	//test(itl2::tests::genericTransform, "Generic geometric transform");
	
	//test(itl2::tests::floodfill, "Flood fill");
	//test(itl2::tests::floodfillSanityChecks, "sanity checks of flood fill implementations");
	//test(itl2::tests::floodfillLeaks, "Flood fill leak tests");

	//test(itl2::tests::mipMatch, "MIP Match");

	//test(itl2::tests::hybridSkeleton, "Hybrid skeleton");
	//test(itl2::tests::lineSkeleton, "Line skeleton");

	//test(itl2::tests::traceSkeleton, "trace skeleton");
	//test(itl2::tests::traceSkeletonRealData, "trace skeleton (real data)");
	//test(itl2::tests::networkio, "network I/O");
	//test(itl2::tests::disconnections, "network connect, disconnect, degree, etc.");
	//test(itl2::tests::disconnectStraightThroughPerformance, "network optimization performance");
	//test(itl2::tests::lineLength, "line length calculation");

	//test(itl2::tests::curvature, "curvature");
	//test(itl2::tests::structureTensor, "structure tensor");
	//test(itl2::tests::lineFilter, "line filtering");
	//test(itl2::tests::canny, "Canny edge detection");

	//test(itl2::tests::normalizeZ, "Normalize Z");

	//test(itl2::tests::analyzeParticlesSanity, "Analyze particles sanity checks");
	//test(itl2::tests::analyzeParticlesSanity2, "Analyze particles sanity checks 2");
	//test(itl2::tests::analyzeParticlesVolumeLimit, "Analyze particles volume limit");
	//test(itl2::tests::analyzeParticlesThreading, "Analyze particles threading");
	//test(itl2::tests::analyzeParticlesThreadingBig, "Analyze particles threading, big volumes"); // This is a long test

	//test(itl2::tests::regionRemoval, "Region removal");

	//test(itl2::tests::lineMax, "Line maximum");
	//test(itl2::tests::lineMin, "Line minimum");
	//test(itl2::tests::sphereMaxSpeed, "Sphere max filtering speed");

	//test(itl2::tests::danielssonTableSpeedTest, "Danielsson lookup table calculation speed");
	//test(itl2::tests::fullDanielssonTables, "Danielsson lookup table calculation and caching");
	//test(itl2::tests::danielsson, "Danielsson algorithm for centers of locally maximal spheres");
	


	testReport();

	//cout << "Press return to exit..." << endl;
	//char dummy;
	//cin.getline(&dummy, 1);
	return 0;
}
