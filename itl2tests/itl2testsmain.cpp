
#include "math/numberutils.h"
#include "math/matrix.h"
#include "test.h"
#include "raytrace.h"
#include "io/vol.h"
#include "utilities.h"
#include "danielsson.h"
#include "io/io.h"
#include "io/itltiff.h"
#include "io/nrrd.h"
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
#include "surfaceskeleton.h"
#include "surfaceskeleton2.h"
#include "traceskeleton.h"
#include "structure.h"
#include "particleanalysis.h"
#include "regionremoval.h"
#include "fastbilateralfilter.h"
#include "minhash.h"
#include "tomo/fbp.h"
#include "thickmap.h"
#include "fillskeleton.h"
#include "surfacecurvature.h"
#include "autothreshold.h"
#include "traceskeletonpoints.h"
#include "generation.h"
#include "noise.h"
#include "maxima.h"
#include "carpet.h"
#include "montage.h"
#include "math/conjugategradient.h"
#include "tomo/siddonprojections.h"
#include "iteration.h"
#include "progress.h"
#include "imagemetadata.h"


using namespace itl2;
using namespace std;





int main()
{

	//test(itl2::tests::progress, "progress indicator");

	//test(itl2::tests::intermediateTypes, "intermediate type determination");
	//test(itl2::tests::equals, "equals");
	//test(itl2::tests::saturatingArithmetic, "saturating arithmetic");
	//test(itl2::tests::matrix3x3, "3x3 matrix");
	//test(itl2::tests::matrix, "Matrix");
	//test(itl2::tests::solve, "Matrix inverse and solution of group of linear equations");
	//test(itl2::tests::leastSquares, "Least squares solution");
	//test(itl2::tests::conjugateGradient, "Conjugate gradient");
	//test(itl2::tests::cgne, "CGNE");
	//test(itl2::tests::image, "Image");

	
	//test(raw::tests::expandFilename, "Raw filename expansion");
	//test(raw::tests::raw, "Raw reader");
	//test(io::tests::readWrite, "IO read");
	//test(raw::tests::writeBlock, "Block based raw reader & writer");
	//test(raw::tests::writeBlockFast, "Optimized block based raw reader & writer");
	//test(vol::tests::volio, ".vol input/output");
	//test(itl2::png::tests::png, "Png read and write");
	//test(itl2::tiff::tests::readWrite, "Tiff read and write");
	//test(itl2::nrrd::tests::readWrite, "NRRD read and write");

	//test(itl2::sequence::tests::match, "Matching");
	//test(itl2::sequence::tests::sequence, "Image sequence");
	//test(itl2::sequence::tests::fileFormats, "Sequence file formats");
	//test(itl2::sequence::tests::readWriteBlock, "Image sequence block");
	//test(itl2::sequence::tests::readWriteBlockOptimization, "Image sequence block write optimization");
	

	//test(itl2::tests::siddonProject, "Siddon algorithm");
	//test(itl2::tests::geometry, "Geometry");

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

	//test(itl2::tests::broadcast, "Broadcasted point process");
	//test(itl2::tests::bilateral, "bilateral filter");

	//test(itl2::tests::scale, "scaling");
	//test(itl2::tests::translate, "translation");
	//test(itl2::tests::rot90, "90 deg rotations");
	//test(itl2::tests::rotate, "rotations around general axes");
	//test(itl2::tests::reslice, "reslice");
	//test(itl2::tests::crop, "crop and reverse crop");

	//test(itl2::tests::pointProcess, "point processes");
	//test(itl2::tests::pointProcessComplex, "point processes on complex numbers");
	//test(itl2::tests::byteOrder, "byte order swaps");


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
	//test(itl2::tests::histogramIntermediateType, "Intermediate types in histogram");
	//test(itl2::tests::histogram, "Histogram");
	//test(itl2::tests::histogram2d, "Bivariate histogram");

	//test(itl2::tests::binning, "Binning");
	//test(itl2::tests::genericTransform, "Generic geometric transform");

	//test(itl2::tests::floodfill, "Flood fill");
	//test(itl2::tests::floodfillSanityChecks, "sanity checks of flood fill implementations");
	//test(itl2::tests::floodfillLeaks, "Flood fill leak tests");

	//test(itl2::tests::mipMatch, "MIP Match");

	//test(itl2::tests::hashPow, "hash pow function");
	//test(itl2::tests::nbHash, "neighbourhood hash function");
	//test(itl2::tests::minHash, "minHash function");

	//test(itl2::tests::surfaceSkeleton, "Surface skeleton");
	//test(itl2::experimental::tests::surfaceSkeleton2, "Hybrid skeleton 2");
	//test(itl2::tests::lineSkeleton, "Line skeleton");

	//test(itl2::tests::traceSkeleton, "trace skeleton");
	//test(itl2::tests::traceSkeletonRealData, "trace skeleton (real data)");
	//test(itl2::tests::networkio, "network I/O");
	//test(itl2::tests::disconnections, "network connect, disconnect, degree, etc.");
	//test(itl2::tests::disconnectStraightThroughPerformance, "network optimization performance");
	//test(itl2::tests::pruning, "pruning");
	//test(itl2::tests::lineLength, "line length calculation");
	//test(itl2::tests::skeletonToPointsAndLines, "skeleton to point-line form");

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

	//test(itl2::tests::scaleLabels, "label/binary image scaling");
	//test(itl2::tests::fastBilateralSampling, "fast bilateral filtering (sampling approximation)");

	//test(itl2::tests::growPriority, "Meyer's growing algorithm");
	//test(itl2::tests::growAll, "region growing");
	//test(itl2::tests::growComparison, "region growing algorithm comparison");


	//test(itl2::tests::fillSkeleton, "skeleton filling");
	//test(itl2::tests::vectorAngles, "calculation of angle between vectors");

	//test(itl2::tests::surfaceCurvature, "surface curvature");

	test(itl2::tests::recSettings, "Rec settings");
	//test(itl2::tests::paganin, "Paganin method");
	//// NOTE: Data for these tests is not publicly available (yet)
	////test(itl2::tests::fbp, "Filtered backprojection");
	////test(itl2::tests::openCLBackProjection, "OpenCL filtered backprojection");
	////test(itl2::tests::openCLBackProjectionRealBin2, "OpenCL filtered backprojection, real dataset, binning 2");
	////test(itl2::tests::openCLBackProjectionRealBin1, "OpenCL filtered backprojection, real dataset, binning 1"); 

	

	//test(itl2::tests::thickmapsEquality, "equality of different thickness map implementations");	
	//test(itl2::tests::dimred2D, "dimensionality reduction approach to local thickness 2D");
	//test(itl2::tests::dimred3D, "dimensionality reduction approach to local thickness 3D");
	//test(itl2::tests::discretizedCircles, "inclusion of discretized circles");
	//test(itl2::tests::indexItem, "RiStorageItem saving test");
	//test(itl2::tests::readWriteRi, "saving of non-trivially copyable image"); 
	//test(itl2::tests::testDataForFiji, "generate comparison dataset for Fiji implementation");
	//test(itl2::tests::thickmapRounding, "Rounding before thickness map calculation");
	//test(itl2::tests::thickmapScaling, "large thickness map test");

	//test(itl2::tests::autothreshold, "automatic thresholding");
	//test(itl2::tests::localThreshold, "local thresholding");
	//test(itl2::tests::localMaxima, "local maxima search");

	//test(itl2::tests::carpet, "surface finding");
	//test(itl2::tests::ellipsoid, "drawing ellipsoids");

	//test(itl2::tests::montage, "2D montage of 3D stack");

	test(itl2::tests::imagemetadata, "image metadata");

	

	// Experimental tests - these are mostly work in progress and data for them is not available yet

	//test(itl2::tests::createProjection, "Single projection");
	//test(itl2::tests::createProjections, "Projections");
	//test(itl2::tests::createBackprojection, "Back-projection");
	//test(itl2::tests::projectionConsistency, "Consistency");
	//test(itl2::tests::create10Projections, "10 projections");
	//test(itl2::tests::create36Projections, "36 projections");


	//test(tomo::tests::cgneReconstruction1, "CGNE reconstruction");
	//test(tomo::tests::create36ProjectionsSheppLogan, "Shepp-Logan");
	//test(tomo::tests::cgneReconstruction2, "CGNE reconstruction (Shepp-Logan)");
	//test(itl2::tests::simplifyRodSkeleton, "Simplification of rod skeleton");
	//test(itl2::tests::simplifyRodSkeleton2, "Simplification of rod skeleton (real data)");


	//try
	//{
	//	tests::thickmapScaling();
	//}
	//catch (ITLException& e)
	//{
	//	cout << e.message() << endl;
	//}
	//catch (std::bad_alloc& e)
	//{
	//	cout << "Error: Out of memory (" << e.what() << ")" << endl;

	//	return 2;
	//}
	//catch (std::exception& e)
	//{
	//	cout << "Error: " << e.what() << endl;

	//	return 3;
	//}
	//catch (...)
	//{
	//	cout << "Error: Unknown error" << endl;

	//	return 4;
	//}

	

	testReport();

	//cout << "Press return to exit..." << endl;
	//char dummy;
	//cin.getline(&dummy, 1);
	return 0;
}
