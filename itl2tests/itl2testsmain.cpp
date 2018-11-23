
#include "itl2.h"
#include "math/numberutils.h"
#include "test.h"
#include "raytrace.h"
#include "math/conjugategradient.h"

#include "io/itlpng.h"

using namespace itl2;
using namespace math;


int main()
{
	// Enable here the tests you would like to run. Please set working directory to testing/ so that test data is found correctly.
	// (otherwise many tests will fail because of missing data)

	
	//test(math::tests::numberUtils, "Number utilities");
	//test(math::tests::matrix, "Matrix");
	//
	//test(itl2::tests::image, "Image");
	//test(raw::tests::raw, "Raw reader");
	//test(raw::tests::readWriteBlock, "Block based raw reader & writer");
	//
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
	//

	//test(itl2::tests::phaseCorrelation, "phase correlation");
	//
	//test(itl2::tests::phaseCorrelation2, "phase correlation 2 (rotation)");

	//test(itl2::tests::blockMatch1, "block match 1");
	//test(itl2::tests::blockMatch2Match, "block match 2 (match)");
	//test(itl2::tests::blockMatch2Pullback, "block match 2 (pullback)");
	//
	//test(itl2::tests::inpaintNearest, "Inpainting");
	//test(itl2::tests::inpaintGarcia, "Inpainting (Garcia)");
	//test(itl2::tests::inpaintGarcia2, "Inpainting 2 (Garcia)");
	//test(itl2::tests::dmap1, "Distance map");

	//test(itl2::tests::buffers, "Disk mapped buffer");
	//test(itl2::tests::histogram, "Histogram");

	//test(itl2::tests::binning, "Binning");
	//test(itl2::tests::genericTransform, "Generic geometric transform");
	//
	//test(itl2::tests::floodfill, "Flood fill");

	//test(itl2::tests::mipMatch, "MIP Match");

	//test(itl2::tests::hybridSkeleton, "Hybrid skeleton");
	//test(itl2::tests::lineSkeleton, "Line skeleton");
	////test(itl2::tests::cavities, "hybrid skeleton of structure with cavity");
	////test(itl2::tests::classifySkeleton, "classify skeleton");

	//test(itl2::tests::traceSkeleton, "trace skeleton");
	//test(itl2::tests::networkio, "network I/O");
	//test(itl2::tests::disconnections, "network connect, disconnect, degree, etc.");
	//test(itl2::tests::disconnectStraightThroughPerformance, "network optimization performance");
	//test(itl2::tests::lineLength, "line length calculation");

	//test(itl2::tests::curvature, "curvature");
	//test(itl2::tests::structureTensor, "structure tensor");
	//test(itl2::tests::lineFilter, "line filtering");

	//test(itl2::tests::normalizeZ, "Normalize Z");

	//test(itl2::png::tests::png, "Png read and write");

	//test(itl2::sequence::tests::match, "Matching");
	//test(itl2::sequence::tests::sequence, "Image sequence");
	//test(itl2::sequence::tests::readWriteBlock, "Image sequence block");

	//test(itl2::tests::regionRemoval, "Region removal");


	testReport();

	//cout << "Press return to exit..." << endl;
	//char dummy;
	//cin.getline(&dummy, 1);
	return 0;
}
