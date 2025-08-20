#pragma once
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include "image.h"
#include <vector>
#include <limits>
#include <omp.h>
#include <optional>
#include <timer.h>
#include <deque>
#include <functional>
#include <cmath>

//#include "logger.h"


using std::vector;
using itl2::Image;
using itl2::Vec3c;
using itl2::coord_t;
using std::abs;
using ::fabsf;
using std::deque;
using itl2::ITLException;
using std::reference_wrapper;
using std::ref;


/// <summary>
/// Class used to segment an image with sFCM clustering.
/// </summary>
/// <typeparam name="pixel_t">Type of a pixel in the image.</typeparam>
template <typename pixel_t>
class SFCM {
private:
	Image<pixel_t>& imageData;
	int clusterCount;
	float fuzzyParam;
	float p;
	float q;
	coord_t imageDataSize;
	coord_t nbApothem;
	Vec3c imageDims;
	int iterMax;
	float stopParam;
	vector<float> centroids;
	deque<vector<float>> centroidsHistory;
	vector<std::unique_ptr<Image<float>>> membershipMatrix;
	vector<std::unique_ptr<Image<float>>> spatialMembershipMatrix;


public:


	/// <summary>
	/// Constructor.
	/// </summary>
	/// <param name="imageData">Image to segment.</param>
	/// <param name="clusterCount">How many classes to segment the image to.</param>
	/// <param name="fuzzyParam">Fuzziness parameter of FCM.</param>
	/// <param name="nbApothem">Apothem/radius of how many pixels to include in spatial neighbourhood calculations.</param>
	/// <param name="iterMax">Maxmimum number of iterations to run the clustering algorithm.</param>
	/// <param name="stopParam">Convergence tolerance threshold value.</param>
	/// <param name="p">Weight parameter of pixel value.</param>
	/// <param name="q">Weight parameter of spatial part.</param>
	SFCM(Image<pixel_t>& imageData, int clusterCount, float fuzzyParam, coord_t nbApothem, int iterMax, float stopParam, float p, float q)
		: imageData(imageData),
		centroids(clusterCount)
	{
		this->clusterCount = clusterCount;
		this->fuzzyParam = fuzzyParam;
		this->iterMax = iterMax;
		this->stopParam = stopParam;
		this->p = p;
		this->q = q;
		this->imageDataSize = imageData.pixelCount();
		this->nbApothem = nbApothem;
		this->imageDims = imageData.dimensions();

		membershipMatrix.reserve(clusterCount);
		for (int i = 0; i < clusterCount; i++) {
			membershipMatrix.push_back(std::make_unique<Image<float>>(imageData.dimensions()));
		}

		initCentroids();

		return;
	}


	/// <summary>
	/// Main loop for clustering.
	/// </summary>
	void cluster() {
		int iteration = 1;
		float change = std::numeric_limits<float>::max();
		const int historySize = 3;
		Timer timerFull;
		Timer timerIteration;
		std::cout << std::fixed << std::setprecision(6);
		timerFull.start();
		if (fuzzyParam != 1) { // if fuzzy Cmeans
			if (q != 0) { // if spatial (sFCM)
				spatialMembershipMatrix.reserve(clusterCount);
				for (int i = 0; i < clusterCount; i++) {
					spatialMembershipMatrix.push_back(std::make_unique<Image<float>>(imageData.dimensions()));
				}
				while (iteration < iterMax && change > stopParam) {
					timerIteration.start();
					updateMembershipMatrixCMeans();
					updateSpatialMembershipMatrixFCM();
					change = updateCentroidsSFCM();
					iteration++;
					double duration = timerIteration.getTime();
					double durationFull = timerFull.getTime();
					std::cout << "iteration: " << iteration << ", change: " << change << ", time taken: " << (duration) << "ms/" << (durationFull) << "ms" << std::endl;

					if (isOscillating()) {
						std::cout << "Oscillation detected! Terminating early since won't converge further." << std::endl;
						break;
					}
					centroidsHistory.push_back(centroids);
					if (centroidsHistory.size() > historySize) {
						centroidsHistory.pop_front();
					}
				}
			}
			else { // if non-spatial (FCM)
				while (iteration < iterMax && change > stopParam) {
					timerIteration.start();
					updateMembershipMatrixCMeans();
					change = updateCentroidsFCM();
					iteration++;
					double duration = timerIteration.getTime();
					double durationFull = timerFull.getTime();
					std::cout << "iteration: " << iteration << ", change: " << change << ", time taken: " << (duration) << "ms/" << (durationFull) << "ms" << std::endl;

					if (isOscillating()) {
						std::cout << "Oscillation detected! Terminating early since won't converge further." << std::endl;
						break;
					}
					centroidsHistory.push_back(centroids);
					if (centroidsHistory.size() > historySize) {
						centroidsHistory.pop_front();
					}
				}
			}
		}
		else { // if hard Cmeans
			if (q != 0) { // if spatial (sHCM)
				spatialMembershipMatrix.reserve(clusterCount);
				for (int i = 0; i < clusterCount; ++i) {
					spatialMembershipMatrix.push_back(std::make_unique<Image<float>>(imageData.dimensions()));
				}

				while (iteration < iterMax && change > stopParam) {
					timerIteration.start();
					updateMembershipMatrixCMeans();
					updateSpatialMembershipMatrixKMeans();
					change = updateCentroidsSKMeans();
					iteration++;
					double duration = timerIteration.getTime();
					double durationFull = timerFull.getTime();
					std::cout << "iteration: " << iteration << ", change: " << change << ", time taken: " << (duration) << "ms/" << (durationFull) << "ms" << std::endl;

					if (isOscillating()) {
						std::cout << "Oscillation detected! Terminating early since won't converge further." << std::endl;
						break;
					}
					centroidsHistory.push_back(centroids);
					if (centroidsHistory.size() > historySize) {
						centroidsHistory.pop_front();
					}
				}
			}
			else { // if non-spatial (Kmeans)
				while (iteration < iterMax && change > stopParam) {
					timerIteration.start();
					updateMembershipMatrixKMeans();
					change = updateCentroidsKMeans();
					iteration++;
					double duration = timerIteration.getTime();
					double durationFull = timerFull.getTime();
					std::cout << "iteration: " << iteration << ", change: " << change << ", time taken: " << (duration) << "ms/" << (durationFull) << "ms" << std::endl;

					if (isOscillating()) {
						std::cout << "Oscillation detected! Terminating early since won't converge further." << std::endl;
						break;
					}
					centroidsHistory.push_back(centroids);
					if (centroidsHistory.size() > historySize) {
						centroidsHistory.pop_front();
					}
				}
			}
		}
		timerFull.stop();
		timerIteration.stop();
		std::cout << std::setprecision(6);
		std::cout.unsetf(std::ios::fixed);

		/*
		if (logging) {
			log_value("iterationCount", iteration);
			log_value("sfcmTotalTimeTaken", timerFull.getTime());
			log_value("change", change);
		}
		*/

		return;
	}


	/// <summary>
	/// Applies clustering result as hard segmentation to the image.
	/// </summary>
	void applySegmentation() {
		if (fuzzyParam != 1) { // if fuzzy Cmeans
			if (q != 0) { // if spatial (sFCM)

#pragma omp parallel for if(imageDataSize > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
				for (coord_t i = 0; i < imageDataSize; i++) {
					float largestMembersip = 0.0f;
					int largestMembershipIndex;
					for (int j = 0; j < clusterCount; j++) {
						if ((*spatialMembershipMatrix[j])(i) > largestMembersip) {
							largestMembersip = (*spatialMembershipMatrix[j])(i);
							largestMembershipIndex = j;
						}
					}
					imageData(i) = pixelRound<pixel_t>(largestMembershipIndex);
				}
			}
			else { // if non-spatial (FCM)

#pragma omp parallel for if(imageDataSize > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
				for (coord_t i = 0; i < imageDataSize; i++) {
					float largestMembersip = 0.0f;
					int largestMembershipIndex;
					for (int j = 0; j < clusterCount; j++) {
						if ((*membershipMatrix[j])(i) > largestMembersip) {
							largestMembersip = (*membershipMatrix[j])(i);
							largestMembershipIndex = j;
						}
					}
					imageData(i) = pixelRound<pixel_t>(largestMembershipIndex);
				}
			}
		}
		else { // if hard Cmeans
			if (q != 0) { // if spatial (sHCM)

#pragma omp parallel for if(imageDataSize > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
				for (coord_t i = 0; i < imageDataSize; i++) {
					for (int j = 0; j < clusterCount; j++) {
						if ((*spatialMembershipMatrix[j])(i) == 1.0f) {
							imageData(i) = pixelRound<pixel_t>(j);
							break;
						}
					}
				}
			}
			else { // if non-spatial (Kmeans)

#pragma omp parallel for if(imageDataSize > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
				for (coord_t i = 0; i < imageDataSize; i++) {
					for (int j = 0; j < clusterCount; j++) {
						if ((*membershipMatrix[j])(i) == 1.0f) {
							imageData(i) = pixelRound<pixel_t>(j);
							break;
						}
					}
				}
			}
			return;
		}


	}


	/// <summary>
	/// Used to transfer ownership of membershipMatrix outside of the class scope if needed.
	/// </summary>
	/// <returns>membershipMatrix vector</returns>
	vector<std::unique_ptr<Image<float>>> extractMembershipMatrix() {
		return std::move(membershipMatrix);
	}


	/// <summary>
	/// Used to transfer ownership of spatialMembershipMatrix outside of the class scope if needed.
	/// </summary>
	/// <returns>spatialmembershipMatrix vector</returns>
	vector<std::unique_ptr<Image<float>>> extractSpatialMembershipMatrix() {
		return std::move(spatialMembershipMatrix);
	}


private:


	/// <summary>
	/// Checks if the clustering is stuck oscillating between same centroids.
	/// </summary>
	/// <returns>true if stuck, false if not</returns>
	bool isOscillating() {
		if (centroidsHistory.size() == 0) {
			return false;
		}
		for (size_t i = 0; i < centroidsHistory.size(); ++i) {
			int countCentroidSame = 0;
			for (size_t j = 0; j < centroids.size(); ++j) {
				if (fabsf(centroidsHistory[i][j] - centroids[j]) <= stopParam) {
					countCentroidSame++;
				}
			}
			if (countCentroidSame == clusterCount) {
				return true;
			}
		}
		return false;
	}


	/// <summary>
	/// Initializes cluster centroids evenly between maximum and minimum values in the image.
	/// </summary>
	void initCentroids() {
		float min = static_cast<float>(itl2::min(imageData));
		float max = static_cast<float>(itl2::max(imageData));
		float spacing = (max - min) / (clusterCount - 1);
		for (int j = 0; j < clusterCount; j++) {
			float newCentroid = min + j * spacing;
			centroids[j] = newCentroid;
		}
		return;
	}


	/// <summary>
	/// Calculates new centroid values from last iteration spatial membership values.
	/// </summary>
	/// <returns>Largest change in cluster center value</returns>
	float updateCentroidsSFCM() {

		float largestChange = 0.0f;
		vector<float> sumA(clusterCount, 0.0f);
		vector<float> sumB(clusterCount, 0.0f);

#pragma omp parallel if(imageDataSize > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		{

			vector<float> threadA(clusterCount, 0.0f);
			vector<float> threadB(clusterCount, 0.0f);

#pragma omp for nowait
			for (coord_t i = 0; i < imageDataSize; i++) {
				for (int j = 0; j < clusterCount; j++) {
					float fuzzyWeightTerm = (fuzzyParam == 2.0f)
						? (*spatialMembershipMatrix[j])(i) * (*spatialMembershipMatrix[j])(i)
						: powf((*spatialMembershipMatrix[j])(i), fuzzyParam);
					threadA[j] += fuzzyWeightTerm * imageData(i);
					threadB[j] += fuzzyWeightTerm;
				}
			}

#pragma omp critical
			{
				for (int j = 0; j < clusterCount; j++) {
					sumA[j] += threadA[j];
					sumB[j] += threadB[j];
				}
			}
		}

		for (int j = 0; j < clusterCount; j++) {
			float newCentroid = sumA[j] / sumB[j];
			float newChange = fabsf(centroids[j] - newCentroid);

			if (newChange > largestChange) {
				largestChange = newChange;
			}

			std::cout << "old centroid: " << centroids[j] << "; new centroid: " << newCentroid << std::endl;
			centroids[j] = newCentroid;
		}

		return largestChange;
	}


	/// <summary>
	/// Calculates new centroid values from last iteration membership values.
	/// </summary>
	/// <returns>Largest change in cluster center value</returns>
	float updateCentroidsFCM() {

		float largestChange = 0.0f;
		vector<float> sumA(clusterCount, 0.0f);
		vector<float> sumB(clusterCount, 0.0f);

#pragma omp parallel if(imageDataSize > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		{

			vector<float> threadA(clusterCount, 0.0f);
			vector<float> threadB(clusterCount, 0.0f);

#pragma omp for nowait
			for (coord_t i = 0; i < imageDataSize; i++) {
				for (int j = 0; j < clusterCount; j++) {
					float fuzzyWeightTerm = (fuzzyParam == 2.0f)
						? (*membershipMatrix[j])(i) * (*membershipMatrix[j])(i)
						: powf((*membershipMatrix[j])(i), fuzzyParam);
					threadA[j] += fuzzyWeightTerm * imageData(i);
					threadB[j] += fuzzyWeightTerm;
				}
			}

#pragma omp critical
			{
				for (int j = 0; j < clusterCount; j++) {
					sumA[j] += threadA[j];
					sumB[j] += threadB[j];
				}
			}
		}

		for (int j = 0; j < clusterCount; j++) {
			float newCentroid = sumA[j] / sumB[j];
			float newChange = fabsf(centroids[j] - newCentroid);

			if (newChange > largestChange) {
				largestChange = newChange;
			}

			std::cout << "old centroid: " << centroids[j] << "; new centroid: " << newCentroid << std::endl;
			centroids[j] = newCentroid;
		}

		return largestChange;
	}


	/// <summary>
	/// Calculates new centroid values from last iteration spatial membership values.
	/// </summary>
	/// <returns>Largest change in cluster center value</returns>
	float updateCentroidsSKMeans() {

		float largestChange = 0.0f;
		vector<float> sumA(clusterCount, 0.0f);
		vector<float> sumB(clusterCount, 0.0f);

#pragma omp parallel if(imageDataSize > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		{

			vector<float> threadA(clusterCount, 0.0f);
			vector<float> threadB(clusterCount, 0.0f);

#pragma omp for nowait
			for (coord_t i = 0; i < imageDataSize; i++) {
				for (int j = 0; j < clusterCount; j++) {
					float weightTerm = (*spatialMembershipMatrix[j])(i);
					threadA[j] += weightTerm * imageData(i);
					threadB[j] += weightTerm;
				}
			}

#pragma omp critical
			{
				for (int j = 0; j < clusterCount; j++) {
					sumA[j] += threadA[j];
					sumB[j] += threadB[j];
				}
			}
		}

		for (int j = 0; j < clusterCount; j++) {
			float newCentroid = sumA[j] / sumB[j];
			float newChange = fabsf(centroids[j] - newCentroid);

			if (newChange > largestChange) {
				largestChange = newChange;
			}

			std::cout << "old centroid: " << centroids[j] << "; new centroid: " << newCentroid << std::endl;
			centroids[j] = newCentroid;
		}

		return largestChange;
	}


	/// <summary>
	/// Calculates new centroid values from last iteration membership values.
	/// </summary>
	/// <returns>Largest change in cluster center value</returns>
	float updateCentroidsKMeans() {

		float largestChange = 0.0f;
		vector<float> sumA(clusterCount, 0.0f);
		vector<float> sumB(clusterCount, 0.0f);

#pragma omp parallel if(imageDataSize > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		{

			vector<float> threadA(clusterCount, 0.0f);
			vector<float> threadB(clusterCount, 0.0f);

#pragma omp for nowait
			for (coord_t i = 0; i < imageDataSize; i++) {
				for (int j = 0; j < clusterCount; j++) {
					float weightTerm = (*membershipMatrix[j])(i);
					threadA[j] += weightTerm * imageData(i);
					threadB[j] += weightTerm;
				}
			}

#pragma omp critical
			{
				for (int j = 0; j < clusterCount; j++) {
					sumA[j] += threadA[j];
					sumB[j] += threadB[j];
				}
			}
		}

		for (int j = 0; j < clusterCount; j++) {
			float newCentroid = sumA[j] / sumB[j];
			float newChange = fabsf(centroids[j] - newCentroid);

			if (newChange > largestChange) {
				largestChange = newChange;
			}

			std::cout << "old centroid: " << centroids[j] << "; new centroid: " << newCentroid << std::endl;
			centroids[j] = newCentroid;
		}

		return largestChange;
	}


	/// <summary>
	/// Calculates new membership values 0 or 1.
	/// </summary>
	void updateMembershipMatrixKMeans() {

#pragma omp parallel for if(imageDataSize > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		for (coord_t i = 0; i < imageDataSize; i++) {
			float minDifference = std::numeric_limits<float>::max();;
			int minDifferenceIndex = 0;
			for (int j = 0; j < clusterCount; j++) {
				float differenceIJ = fabsf(centroids[j] - imageData(i));
				if (differenceIJ < minDifference) {
					minDifference = differenceIJ;
					minDifferenceIndex = j;
				}
			}
			for (int j = 0; j < clusterCount; j++) {
				if (j == minDifferenceIndex) {
					(*membershipMatrix[j])(i) = 1.0f;
				}
				else {
					(*membershipMatrix[j])(i) = 0.0f;
				}
			}
		}
		return;
	}


	/// <summary>
	/// Calculates new membership values between 0 and 1.
	/// </summary>
	void updateMembershipMatrixCMeans() {

		const float epsilon = 1e-10f; // Prevents division by zero when a pixel is exactly the same as a centroid.

#pragma omp parallel for if(imageDataSize > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		for (coord_t i = 0; i < imageDataSize; i++) {
			for (int j = 0; j < clusterCount; j++) {
				float sumKC = 0.0f;
				float differenceIJ = fabsf(centroids[j] - imageData(i)) + epsilon;
				for (int k = 0; k < clusterCount; k++) {
					float differenceIK = fabsf(centroids[k] - imageData(i)) + epsilon;
					float differenceFraction = differenceIJ / differenceIK;

					sumKC += (fuzzyParam == 2.0f) ? differenceFraction * differenceFraction
						: powf(differenceFraction, 2.0f / (fuzzyParam - 1.0f));

				}
				(*membershipMatrix[j])(i) = (1.0f / sumKC);
			}
		}
		return;
	}


	/// <summary>
	/// Calculates new spatial membership values between 0 and 1.
	/// </summary>
	void updateSpatialMembershipMatrixFCM() {

#pragma omp parallel for if(imageDataSize > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		for (coord_t i = 0; i < imageDataSize; i++) {
			for (int j = 0; j < clusterCount; j++) {
				float sumKC = 0.0;
				for (int k = 0; k < clusterCount; k++) {

					float memberShip1 = (fuzzyParam == 1.0f) ? (*membershipMatrix[k])(i)
						: (fuzzyParam == 2.0f) ? (*membershipMatrix[k])(i) * (*membershipMatrix[k])(i)
						: powf((*membershipMatrix[k])(i), p);

					float nbSum1 = calcNeighbourhoodLegacy(i, k);
					float neighbourhood1 = (fuzzyParam == 1.0f) ? nbSum1
						: (fuzzyParam == 2.0f) ? nbSum1 * nbSum1
						: powf(nbSum1, q);

					sumKC += memberShip1 * neighbourhood1;
				}
				float memberShip2 = (fuzzyParam == 1.0f) ? (*membershipMatrix[j])(i)
					: (fuzzyParam == 2.0f) ? (*membershipMatrix[j])(i) * (*membershipMatrix[j])(i)
					: powf((*membershipMatrix[j])(i), p);

				float nbSum2 = calcNeighbourhoodLegacy(i, j);
				float neighbourhood2 = (fuzzyParam == 1.0f) ? nbSum2
					: (fuzzyParam == 2.0f) ? nbSum2 * nbSum2
					: powf(nbSum2, q);

				(*spatialMembershipMatrix[j])(i) = (memberShip2 * neighbourhood2) / sumKC;
			}
		}
		return;
	}


	/// <summary>
	/// Calculates new spatial membership values 0 or 1.
	/// </summary>
	void updateSpatialMembershipMatrixKMeans() {

#pragma omp parallel for if(imageDataSize > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		for (coord_t i = 0; i < imageDataSize; i++) {
			float largestMembersip = 0.0f;
			int largestMembershipIndex = 0;
			for (int j = 0; j < clusterCount; j++) {
				float sumKC = 0.0;
				for (int k = 0; k < clusterCount; k++) {

					float memberShip1 = (fuzzyParam == 1.0f) ? (*membershipMatrix[k])(i)
						: (fuzzyParam == 2.0f) ? (*membershipMatrix[k])(i) * (*membershipMatrix[k])(i)
						: powf((*membershipMatrix[k])(i), p);

					float nbSum1 = calcNeighbourhoodLegacy(i, k);
					float neighbourhood1 = (fuzzyParam == 1.0f) ? nbSum1
						: (fuzzyParam == 2.0f) ? nbSum1 * nbSum1
						: powf(nbSum1, q);

					sumKC += memberShip1 * neighbourhood1;
				}
				float memberShip2 = (fuzzyParam == 1.0f) ? (*membershipMatrix[j])(i)
					: (fuzzyParam == 2.0f) ? (*membershipMatrix[j])(i) * (*membershipMatrix[j])(i)
					: powf((*membershipMatrix[j])(i), p);

				float nbSum2 = calcNeighbourhoodLegacy(i, j);
				float neighbourhood2 = (fuzzyParam == 1.0f) ? nbSum2
					: (fuzzyParam == 2.0f) ? nbSum2 * nbSum2
					: powf(nbSum2, q);

				float membershipValue = (memberShip2 * neighbourhood2) / sumKC;

				if (membershipValue > largestMembersip) {
					largestMembersip = membershipValue;
					largestMembershipIndex = j;
				}
			}
			for (int j = 0; j < clusterCount; j++) {
				if (j == largestMembershipIndex) {
					(*spatialMembershipMatrix[j])(i) = 1.0f;
				}
				else {
					(*spatialMembershipMatrix[j])(i) = 0.0f;
				}
			}
		}
		return;
	}


	/// <summary>
	/// Calculates sum of neighbouring pixel values.
	/// </summary>
	/// <param name="i">Center pixel of neighbourhood</param>
	/// <param name="j">Cluster to calculate sum for.</param>
	/// <returns></returns>
	float calcNeighbourhoodLegacy(coord_t i, int j) {
		Vec3c pointCoords = imageData.getCoords(i);

		coord_t startingIndX = (pointCoords.x < nbApothem) ? 0 : (pointCoords.x - nbApothem);
		coord_t tempEdge = imageDims.x - 1;
		coord_t endingIndX = (tempEdge < (pointCoords.x + nbApothem)) ? (tempEdge) : (pointCoords.x + nbApothem);

		coord_t startingIndY = (pointCoords.y < nbApothem) ? 0 : (pointCoords.y - nbApothem);
		tempEdge = imageDims.y - 1;
		coord_t endingIndY = (tempEdge < (pointCoords.y + nbApothem)) ? (tempEdge) : (pointCoords.y + nbApothem);

		coord_t startingIndZ = (pointCoords.z < nbApothem) ? 0 : (pointCoords.z - nbApothem);
		tempEdge = imageDims.z - 1;
		coord_t endingIndZ = (tempEdge < (pointCoords.z + nbApothem)) ? (tempEdge) : (pointCoords.z + nbApothem);

		float neighbourhoodSum = 0.0;
		for (coord_t iz = startingIndZ; iz <= endingIndZ; iz++) {
			for (coord_t iy = startingIndY; iy <= endingIndY; iy++) {
				for (coord_t ix = startingIndX; ix <= endingIndX; ix++) {
					neighbourhoodSum += (*membershipMatrix[j])(ix, iy, iz);
				}
			}
		}

		return neighbourhoodSum;
	}


};


/// <summary>
/// Class used to segment image with sFCM clustering. This version takes multiple separate images with each having a different feature. For example R,G,B.
/// </summary>
/// <typeparam name="pixel_t">Type of a pixel in the image.</typeparam>
template <typename pixel_t>
class multiImagesFCM {
private:
	vector<reference_wrapper<Image<pixel_t>>> imageDataRefs;
	int clusterCount;
	int featureCount;
	float fuzzyParam;
	float p;
	float q;
	coord_t imageDataSize;
	coord_t nbApothem;
	Vec3c imageDims;
	int iterMax;
	float stopParam;
	vector<vector<float>> centroids; //feature<cluster<float>>
	deque<vector<vector<float>>> centroidsHistory;
	vector<std::unique_ptr<Image<float>>> membershipMatrix;
	vector<std::unique_ptr<Image<float>>> spatialMembershipMatrix;


public:


	/// <summary>
	/// Constructor.
	/// </summary>
	/// <param name="imageDataRefs">References to the images to segment.</param>
	/// <param name="clusterCount">How many classes to segment the image to.</param>
	/// <param name="fuzzyParam">Fuzziness parameter of FCM.</param>
	/// <param name="nbApothem">Apothem/radius of how many pixels to include in spatial neighbourhood calculations.</param>
	/// <param name="iterMax">Maxmimum number of iterations to run the clustering algorithm.</param>
	/// <param name="stopParam">Convergence tolerance threshold value.</param>
	/// <param name="p">Weight parameter of pixel value.</param>
	/// <param name="q">Weight parameter of spatial part.</param>
	multiImagesFCM(vector<reference_wrapper<Image<pixel_t>>> imageDataRefs, int clusterCount, float fuzzyParam, coord_t nbApothem, int iterMax, float stopParam, float p, float q)
		: imageDataRefs(std::move(imageDataRefs)),
		centroids(this->imageDataRefs.size(), vector<float>(clusterCount))
	{
		this->clusterCount = clusterCount;
		this->featureCount = static_cast<int>(this->imageDataRefs.size());
		this->fuzzyParam = fuzzyParam;
		this->p = p;
		this->q = q;
		this->imageDataSize = this->imageDataRefs[0].get().pixelCount();
		this->nbApothem = nbApothem;
		this->imageDims = this->imageDataRefs[0].get().dimensions();
		this->iterMax = iterMax;
		this->stopParam = stopParam;

		membershipMatrix.reserve(clusterCount);
		for (int i = 0; i < clusterCount; i++) {
			membershipMatrix.push_back(std::make_unique<Image<float>>(imageDims));
		}

		initCentroids();

		return;
	}


	/// <summary>
	/// Main loop for clustering.
	/// </summary>
	void cluster() {
		int iteration = 1;
		float change = std::numeric_limits<float>::max();
		const int historySize = 3;
		Timer timerFull;
		Timer timerIteration;
		std::cout << std::fixed << std::setprecision(6);
		timerFull.start();
		if (fuzzyParam != 1) { // if fuzzy Cmeans
			if (q != 0) { // if spatial (sFCM)
				spatialMembershipMatrix.reserve(clusterCount);
				for (int i = 0; i < clusterCount; i++) {
					spatialMembershipMatrix.push_back(std::make_unique<Image<float>>(imageDims));
				}
				while (iteration < iterMax && change > stopParam) {
					timerIteration.start();
					updateMembershipMatrixCMeans();
					updateSpatialMembershipMatrixFCM();
					change = updateCentroidsSFCM();
					iteration++;
					double duration = timerIteration.getTime();
					double durationFull = timerFull.getTime();
					std::cout << "iteration: " << iteration << ", change: " << change << ", time taken: " << (duration) << "ms/" << (durationFull) << "ms" << std::endl;

					if (isOscillating()) {
						std::cout << "Oscillation detected! Terminating early since won't converge further." << std::endl;
						break;
					}
					centroidsHistory.push_back(centroids);
					if (centroidsHistory.size() > historySize) {
						centroidsHistory.pop_front();
					}
				}
			}
			else { // if non-spatial (FCM)
				while (iteration < iterMax && change > stopParam) {
					timerIteration.start();
					updateMembershipMatrixCMeans();
					change = updateCentroidsFCM();
					iteration++;
					double duration = timerIteration.getTime();
					double durationFull = timerFull.getTime();
					std::cout << "iteration: " << iteration << ", change: " << change << ", time taken: " << (duration) << "ms/" << (durationFull) << "ms" << std::endl;

					if (isOscillating()) {
						std::cout << "Oscillation detected! Terminating early since won't converge further." << std::endl;
						break;
					}
					centroidsHistory.push_back(centroids);
					if (centroidsHistory.size() > historySize) {
						centroidsHistory.pop_front();
					}
				}
			}
		}
		else { // if hard Cmeans
			if (q != 0) { // if spatial (sHCM)
				// spatialMembershipMatrix = vector<Image<float>>(clusterCount, Image<float>(imageData.dimensions()));
				spatialMembershipMatrix.reserve(clusterCount);
				for (int i = 0; i < clusterCount; ++i) {
					spatialMembershipMatrix.push_back(std::make_unique<Image<float>>(imageDataRefs[0].get().dimensions()));
				}

				while (iteration < iterMax && change > stopParam) {
					timerIteration.start();
					updateMembershipMatrixCMeans();
					updateSpatialMembershipMatrixKMeans();
					change = updateCentroidsSKMeans();
					iteration++;
					double duration = timerIteration.getTime();
					double durationFull = timerFull.getTime();
					std::cout << "iteration: " << iteration << ", change: " << change << ", time taken: " << (duration) << "ms/" << (durationFull) << "ms" << std::endl;

					if (isOscillating()) {
						std::cout << "Oscillation detected! Terminating early since won't converge further." << std::endl;
						break;
					}
					centroidsHistory.push_back(centroids);
					if (centroidsHistory.size() > historySize) {
						centroidsHistory.pop_front();
					}
				}
			}
			else { // if non-spatial (Kmeans)
				while (iteration < iterMax && change > stopParam) {
					timerIteration.start();
					updateMembershipMatrixKMeans();
					change = updateCentroidsKMeans();
					iteration++;
					double duration = timerIteration.getTime();
					double durationFull = timerFull.getTime();
					std::cout << "iteration: " << iteration << ", change: " << change << ", time taken: " << (duration) << "ms/" << (durationFull) << "ms" << std::endl;

					if (isOscillating()) {
						std::cout << "Oscillation detected! Terminating early since won't converge further." << std::endl;
						break;
					}
					centroidsHistory.push_back(centroids);
					if (centroidsHistory.size() > historySize) {
						centroidsHistory.pop_front();
					}
				}
			}
		}
		timerFull.stop();
		timerIteration.stop();
		std::cout << std::setprecision(6);
		std::cout.unsetf(std::ios::fixed);
		return;
	}


	/// <summary>
	/// Applies clustering result as hard segmentation to the image.
	/// </summary>
	void applySegmentation() {
		if (fuzzyParam != 1) { // if fuzzy Cmeans
			if (q != 0) { // if spatial (sFCM)

#pragma omp parallel for if(imageDataSize > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
				for (coord_t i = 0; i < imageDataSize; i++) {
					float largestMembersip = 0.0f;
					int largestMembershipIndex;
					for (int j = 0; j < clusterCount; j++) {
						if ((*spatialMembershipMatrix[j])(i) > largestMembersip) {
							largestMembersip = (*spatialMembershipMatrix[j])(i);
							largestMembershipIndex = j;
						}
					}
					imageDataRefs[0].get()(i) = pixelRound<pixel_t>(largestMembershipIndex);
				}
			}
			else { // if non-spatial (FCM)

#pragma omp parallel for if(imageDataSize > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
				for (coord_t i = 0; i < imageDataSize; i++) {
					float largestMembersip = 0.0f;
					int largestMembershipIndex = -1;
					for (int j = 0; j < clusterCount; j++) {
						if ((*membershipMatrix[j])(i) > largestMembersip) {
							largestMembersip = (*membershipMatrix[j])(i);
							largestMembershipIndex = j;
						}
					}
					imageDataRefs[0].get()(i) = pixelRound<pixel_t>(largestMembershipIndex);
				}
			}
		}
		else { // if hard Cmeans
			if (q != 0) { // if spatial (sHCM)

#pragma omp parallel for if(imageDataSize > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
				for (coord_t i = 0; i < imageDataSize; i++) {
					for (int j = 0; j < clusterCount; j++) {
						if ((*spatialMembershipMatrix[j])(i) == 1.0f) {
							imageDataRefs[0].get()(i) = pixelRound<pixel_t>(j);
							break;
						}
					}
				}
			}
			else { // if non-spatial (Kmeans)

#pragma omp parallel for if(imageDataSize > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
				for (coord_t i = 0; i < imageDataSize; i++) {
					for (int j = 0; j < clusterCount; j++) {
						if ((*membershipMatrix[j])(i) == 1.0f) {
							imageDataRefs[0].get()(i) = pixelRound<pixel_t>(j);
							break;
						}
					}
				}
			}
			return;
		}


	}


	/// <summary>
	/// Used to transfer ownership of membershipMatrix outside of the class scope if needed.
	/// </summary>
	/// <returns>membershipMatrix vector</returns>
	vector<std::unique_ptr<Image<float>>> extractMembershipMatrix() {
		return std::move(membershipMatrix);
	}


	/// <summary>
	/// Used to transfer ownership of spatialMembershipMatrix outside of the class scope if needed.
	/// </summary>
	/// <returns>spatialmembershipMatrix vector</returns>
	vector<std::unique_ptr<Image<float>>> extractSpatialMembershipMatrix() {
		return std::move(spatialMembershipMatrix);
	}


private:


	/// <summary>
	/// Checks if the clustering is stuck oscillating between same centroids.
	/// </summary>
	/// <returns>true if stuck, false if not</returns>
	bool isOscillating() {
		if (centroidsHistory.size() == 0) {
			return false;
		}
		for (size_t i = 0; i < centroidsHistory.size(); i++) {
			int countCentroidSame = 0;
			for (size_t j = 0; j < clusterCount; j++) {
				int countFeatureSame = 0;
				for (size_t f = 0; f < featureCount; f++) {
					if (fabsf(centroidsHistory[i][f][j] - centroids[f][j]) <= stopParam) {
						countFeatureSame++;
					}
				}
				if (countFeatureSame == featureCount) {
					countCentroidSame++;
				}
			}
			if (countCentroidSame == clusterCount) {
				return true;
			}
		}
		return false;
	}


	/// <summary>
	/// Initializes cluster centroids evenly between maximum and minimum values in the image.
	/// </summary>
	void initCentroids() {
		std::cout << "featureCount " << featureCount << std::endl;
		std::cout << "clusterCount " << clusterCount << std::endl;
		for (int f = 0; f < featureCount; ++f) {
			float min = static_cast<float>(itl2::min(imageDataRefs[f].get()));
			float max = static_cast<float>(itl2::max(imageDataRefs[f].get()));
			float spacing = (max - min) / (clusterCount - 1);
			for (int j = 0; j < clusterCount; j++) {
				float newCentroid = min + j * spacing;
				centroids[f][j] = newCentroid;
			}
		}
		return;
	}


	/// <summary>
	/// Calculates new centroid values from last iteration spatial membership values.
	/// </summary>
	/// <returns>Largest change in cluster center value</returns>
	float updateCentroidsSFCM() {

		float largestChange = 0.0f;
		vector<vector<float>> sumA(featureCount, vector<float>(clusterCount, 0.0f));
		vector<float> sumB(clusterCount, 0.0f);

#pragma omp parallel if(imageDataSize > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		{

			vector<vector<float>> threadA(featureCount, vector<float>(clusterCount, 0.0f));
			vector<float> threadB(clusterCount, 0.0f);

#pragma omp for nowait
			for (coord_t i = 0; i < imageDataSize; i++) {

				vector<float> fuzzyWeightTerms(clusterCount);

				for (int j = 0; j < clusterCount; j++) {

					fuzzyWeightTerms[j] = (fuzzyParam == 2.0f)
						? (*spatialMembershipMatrix[j])(i) * (*spatialMembershipMatrix[j])(i)
						: powf((*spatialMembershipMatrix[j])(i), fuzzyParam);

					threadB[j] += fuzzyWeightTerms[j];
				}

				for (int f = 0; f < featureCount; f++) {

					for (int j = 0; j < clusterCount; j++) {

						threadA[f][j] += fuzzyWeightTerms[j] * imageDataRefs[f].get()(i);
					}
				}
			}

#pragma omp critical
			{
				for (int f = 0; f < featureCount; f++) {

					for (int j = 0; j < clusterCount; j++) {

						sumA[f][j] += threadA[f][j];
					}
				}
				for (int j = 0; j < clusterCount; j++) {

					sumB += threadB[j];
				}
			}
		}

		for (int f = 0; f < featureCount; f++) {

			for (int j = 0; j < clusterCount; j++) {

				float newCentroid = sumA[f][j] / sumB[j];
				float newChange = fabsf(centroids[f][j] - newCentroid);

				if (newChange > largestChange) {
					largestChange = newChange;
				}

				std::cout << "old centroid: " << centroids[f][j] << "; new centroid: " << newCentroid << std::endl;
				centroids[f][j] = newCentroid;
			}
		}

		return largestChange;
	}


	/// <summary>
	/// Calculates new centroid values from last iteration membership values.
	/// </summary>
	/// <returns>Largest change in cluster center value</returns>
	float updateCentroidsFCM() {

		float largestChange = 0.0f;
		vector<vector<float>> sumA(featureCount, vector<float>(clusterCount, 0.0f));
		vector<float> sumB(clusterCount, 0.0f);

#pragma omp parallel if(imageDataSize > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		{

			vector<vector<float>> threadA(featureCount, vector<float>(clusterCount, 0.0f));
			vector<float> threadB(clusterCount, 0.0f);

#pragma omp for nowait
			for (coord_t i = 0; i < imageDataSize; i++) {

				vector<float> fuzzyWeightTerms(clusterCount);

				for (int j = 0; j < clusterCount; j++) {

					fuzzyWeightTerms[j] = (fuzzyParam == 2.0f)
						? (*membershipMatrix[j])(i) * (*membershipMatrix[j])(i)
						: powf((*membershipMatrix[j])(i), fuzzyParam);

					threadB[j] += fuzzyWeightTerms[j];
				}

				for (int f = 0; f < featureCount; f++) {

					for (int j = 0; j < clusterCount; j++) {

						threadA[f][j] += fuzzyWeightTerms[j] * imageDataRefs[f].get()(i);
					}
				}
			}

#pragma omp critical
			{
				for (int f = 0; f < featureCount; f++) {

					for (int j = 0; j < clusterCount; j++) {

						sumA[f][j] += threadA[f][j];
					}
				}
				for (int j = 0; j < clusterCount; j++) {

					sumB += threadB[j];
				}
			}
		}

		for (int f = 0; f < featureCount; f++) {

			for (int j = 0; j < clusterCount; j++) {

				float newCentroid = sumA[f][j] / sumB[j];
				float newChange = fabsf(centroids[f][j] - newCentroid);

				if (newChange > largestChange) {
					largestChange = newChange;
				}

				std::cout << "old centroid: " << centroids[f][j] << "; new centroid: " << newCentroid << std::endl;
				centroids[f][j] = newCentroid;
			}
		}

		return largestChange;
	}


	/// <summary>
	/// Calculates new centroid values from last iteration spatial membership values.
	/// </summary>
	/// <returns>Largest change in cluster center value</returns>
	float updateCentroidsSKMeans() {

		float largestChange = 0.0f;
		vector<vector<float>> sumA(featureCount, vector<float>(clusterCount, 0.0f));
		vector<float> sumB(clusterCount, 0.0f);

#pragma omp parallel if(imageDataSize > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		{

			vector<vector<float>> threadA(featureCount, vector<float>(clusterCount, 0.0f));
			vector<float> threadB(clusterCount, 0.0f);

#pragma omp for nowait
			for (coord_t i = 0; i < imageDataSize; i++) {

				vector<float> fuzzyWeightTerms(clusterCount);

				for (int j = 0; j < clusterCount; j++) {

					fuzzyWeightTerms[j] = (*spatialMembershipMatrix[j])(i);
					threadB[j] += fuzzyWeightTerms[j];
				}

				for (int f = 0; f < featureCount; f++) {

					for (int j = 0; j < clusterCount; j++) {

						threadA[f][j] += fuzzyWeightTerms[j] * imageDataRefs[f].get()(i);
					}
				}
			}

#pragma omp critical
			{
				for (int f = 0; f < featureCount; f++) {

					for (int j = 0; j < clusterCount; j++) {

						sumA[f][j] += threadA[f][j];
					}
				}
				for (int j = 0; j < clusterCount; j++) {

					sumB += threadB[j];
				}
			}
		}

		for (int f = 0; f < featureCount; f++) {

			for (int j = 0; j < clusterCount; j++) {

				float newCentroid = sumA[f][j] / sumB[j];
				float newChange = fabsf(centroids[f][j] - newCentroid);

				if (newChange > largestChange) {
					largestChange = newChange;
				}

				std::cout << "old centroid: " << centroids[f][j] << "; new centroid: " << newCentroid << std::endl;
				centroids[f][j] = newCentroid;
			}
		}

		return largestChange;
	}


	/// <summary>
/// Calculates new centroid values from last iteration spatial membership values.
/// </summary>
/// <returns>Largest change in cluster center value</returns>
	float updateCentroidsKMeans() {

		float largestChange = 0.0f;
		vector<vector<float>> sumA(featureCount, vector<float>(clusterCount, 0.0f));
		vector<float> sumB(clusterCount, 0.0f);

#pragma omp parallel if(imageDataSize > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		{

			vector<vector<float>> threadA(featureCount, vector<float>(clusterCount, 0.0f));
			vector<float> threadB(clusterCount, 0.0f);

#pragma omp for nowait
			for (coord_t i = 0; i < imageDataSize; i++) {

				vector<float> fuzzyWeightTerms(clusterCount);

				for (int j = 0; j < clusterCount; j++) {

					fuzzyWeightTerms[j] = (*membershipMatrix[j])(i);
					threadB[j] += fuzzyWeightTerms[j];
				}

				for (int f = 0; f < featureCount; f++) {

					for (int j = 0; j < clusterCount; j++) {

						threadA[f][j] += fuzzyWeightTerms[j] * imageDataRefs[f].get()(i);
					}
				}
			}

#pragma omp critical
			{
				for (int f = 0; f < featureCount; f++) {

					for (int j = 0; j < clusterCount; j++) {

						sumA[f][j] += threadA[f][j];
					}
				}
				for (int j = 0; j < clusterCount; j++) {

					sumB += threadB[j];
				}
			}
		}

		for (int f = 0; f < featureCount; f++) {

			for (int j = 0; j < clusterCount; j++) {

				float newCentroid = sumA[f][j] / sumB[j];
				float newChange = fabsf(centroids[f][j] - newCentroid);

				if (newChange > largestChange) {
					largestChange = newChange;
				}

				std::cout << "old centroid: " << centroids[f][j] << "; new centroid: " << newCentroid << std::endl;
				centroids[f][j] = newCentroid;
			}
		}

		return largestChange;
	}


	/// <summary>
	/// Calculates new membership values 0 or 1.
	/// </summary>
	void updateMembershipMatrixKMeans() {

#pragma omp parallel for if(imageDataSize > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		for (coord_t i = 0; i < imageDataSize; i++) {
			float minDifference = std::numeric_limits<float>::max();
			int minDifferenceIndex = 0;
			for (int j = 0; j < clusterCount; j++) {

				float squaredDistance = 0.0f;
				for (int f = 0; f < featureCount; f++) {
					float difference = centroids[f][j] - imageDataRefs[f].get()(i);
					squaredDistance += difference * difference;
				}

				if (squaredDistance < minDifference) {
					minDifference = squaredDistance;
					minDifferenceIndex = j;
				}

			}
			for (int j = 0; j < clusterCount; j++) {
				if (j == minDifferenceIndex) {
					(*membershipMatrix[j])(i) = 1.0f;
				}
				else {
					(*membershipMatrix[j])(i) = 0.0f;
				}
			}
		}
		return;
	}


	/// <summary>
	/// Calculates new membership values between 0 and 1.
	/// </summary>
	void updateMembershipMatrixCMeans() {
		const float epsilon = 1e-10f; // Prevents division by zero when a pixel is exactly the same as a centroid.
#pragma omp parallel for if(imageDataSize > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		for (coord_t i = 0; i < imageDataSize; i++) {
			for (int j = 0; j < clusterCount; j++) {
				float sumKC = 0.0f;
				float distanceIJSquared = 0.0f;
				for (int f = 0; f < featureCount; f++) {
					float difference = centroids[f][j] - imageDataRefs[f].get()(i);
					distanceIJSquared += difference * difference;
				}
				float distanceIJ = std::sqrt(distanceIJSquared) + epsilon;
				for (int k = 0; k < clusterCount; k++) {

					float distanceIKSquared = 0.0f;
					for (int f = 0; f < featureCount; f++) {
						float difference = centroids[f][k] - imageDataRefs[f].get()(i);
						distanceIKSquared += difference * difference;
					}
					float distanceIK = std::sqrt(distanceIKSquared) + epsilon;
					float distanceFraction = distanceIJ / distanceIK;

					sumKC += (fuzzyParam == 2.0f) ? distanceFraction * distanceFraction
						: powf(distanceFraction, 2.0f / (fuzzyParam - 1.0f));

				}
				(*membershipMatrix[j])(i) = (1.0f / sumKC);
			}
		}
		return;
	}


	/// <summary>
	/// Calculates new spatial membership values between 0 and 1.
	/// </summary>
	void updateSpatialMembershipMatrixFCM() {

#pragma omp parallel for if(imageDataSize > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		for (coord_t i = 0; i < imageDataSize; i++) {
			for (int j = 0; j < clusterCount; j++) {
				float sumKC = 0.0;
				for (int k = 0; k < clusterCount; k++) {
					float memberShip1 = (p == 1.0f) ? (*membershipMatrix[k])(i)
						: (p == 2.0f) ? (*membershipMatrix[k])(i) * (*membershipMatrix[k])(i)
						: powf((*membershipMatrix[k])(i), p);

					float nbSum1 = calcNeighbourhoodLegacy(i, k);
					float neighbourhood1 = (q == 1.0f) ? nbSum1
						: (q == 2.0f) ? nbSum1 * nbSum1
						: powf(nbSum1, q);

					sumKC += memberShip1 * neighbourhood1;
				}

				float memberShip2 = (p == 1.0f) ? (*membershipMatrix[j])(i)
					: (p == 2.0f) ? (*membershipMatrix[j])(i) * (*membershipMatrix[j])(i)
					: powf((*membershipMatrix[j])(i), p);

				float nbSum2 = calcNeighbourhoodLegacy(i, j);
				float neighbourhood2 = (q == 1.0f) ? nbSum2
					: (q == 2.0f) ? nbSum2 * nbSum2
					: powf(nbSum2, q);

				(*spatialMembershipMatrix[j])(i) = (memberShip2 * neighbourhood2) / sumKC;
			}
		}
		return;
	}


	/// <summary>
	/// Calculates new spatial membership values 0 or 1.
	/// </summary>
	void updateSpatialMembershipMatrixKMeans() {

#pragma omp parallel for if(imageDataSize > PARALLELIZATION_THRESHOLD && !omp_in_parallel())
		for (coord_t i = 0; i < imageDataSize; i++) {
			float largestMembersip = 0.0f;
			int largestMembershipIndex = 0;
			for (int j = 0; j < clusterCount; j++) {
				float sumKC = 0.0;
				for (int k = 0; k < clusterCount; k++) {

					float memberShip1 = (p == 1.0f) ? (*membershipMatrix[k])(i)
						: (p == 2.0f) ? (*membershipMatrix[k])(i) * (*membershipMatrix[k])(i)
						: powf((*membershipMatrix[k])(i), p);

					float nbSum1 = calcNeighbourhoodLegacy(i, k);
					float neighbourhood1 = (q == 1.0f) ? nbSum1
						: (q == 2.0f) ? nbSum1 * nbSum1
						: powf(nbSum1, q);

					sumKC += memberShip1 * neighbourhood1;
				}
				float memberShip2 = (p == 1.0f) ? (*membershipMatrix[j])(i)
					: (p == 2.0f) ? (*membershipMatrix[j])(i) * (*membershipMatrix[j])(i)
					: powf((*membershipMatrix[j])(i), p);

				float nbSum2 = calcNeighbourhoodLegacy(i, j);
				float neighbourhood2 = (q == 1.0f) ? nbSum2
					: (q == 2.0f) ? nbSum2 * nbSum2
					: powf(nbSum2, q);

				float membershipValue = (memberShip2 * neighbourhood2) / sumKC;

				if (membershipValue > largestMembersip) {
					largestMembersip = membershipValue;
					largestMembershipIndex = j;
				}
			}
			for (int j = 0; j < clusterCount; j++) {
				if (j == largestMembershipIndex) {
					(*spatialMembershipMatrix[j])(i) = 1.0f;
				}
				else {
					(*spatialMembershipMatrix[j])(i) = 0.0f;
				}
			}
		}
		return;
	}


	/// <summary>
	/// Calculates sum of neighbouring pixel values.
	/// </summary>
	/// <param name="i">Center pixel of neighbourhood</param>
	/// <param name="j">Cluster to calculate sum for.</param>
	/// <returns></returns>
	float calcNeighbourhoodLegacy(coord_t i, int j) {
		Vec3c pointCoords = imageDataRefs[0].get().getCoords(i);

		coord_t startingIndX = (pointCoords.x < nbApothem) ? 0 : (pointCoords.x - nbApothem);
		coord_t tempEdge = imageDims.x - 1;
		coord_t endingIndX = (tempEdge < (pointCoords.x + nbApothem)) ? (tempEdge) : (pointCoords.x + nbApothem);

		coord_t startingIndY = (pointCoords.y < nbApothem) ? 0 : (pointCoords.y - nbApothem);
		tempEdge = imageDims.y - 1;
		coord_t endingIndY = (tempEdge < (pointCoords.y + nbApothem)) ? (tempEdge) : (pointCoords.y + nbApothem);

		coord_t startingIndZ = (pointCoords.z < nbApothem) ? 0 : (pointCoords.z - nbApothem);
		tempEdge = imageDims.z - 1;
		coord_t endingIndZ = (tempEdge < (pointCoords.z + nbApothem)) ? (tempEdge) : (pointCoords.z + nbApothem);

		float neighbourhoodSum = 0.0;
		for (coord_t iz = startingIndZ; iz <= endingIndZ; iz++) {
			for (coord_t iy = startingIndY; iy <= endingIndY; iy++) {
				for (coord_t ix = startingIndX; ix <= endingIndX; ix++) {
					neighbourhoodSum += (*membershipMatrix[j])(ix, iy, iz);
				}
			}
		}

		return neighbourhoodSum;
	}


};


/// <summary>
/// Calculates a crude approximate for how much RAM is required to run the algorithm.
/// </summary>
/// <typeparam name="pixel_t">Type of a pixel in the image.</typeparam>
/// <param name="pixelCount">Amount of pixels in the image.</param>
/// <param name="clusterCount">How many classes to segment the image to.</param>
/// <param name="pixel">A pixel of the image, only used to get pixel type.</param>
/// <param name="featureCount">How many features = how many images.</param>
template <typename pixel_t>
static void calcRequiredMemory(double pixelCount, int clusterCount, pixel_t pixel, int featureCount) {
	double sizeImage = (pixelCount * sizeof(pixel_t) * featureCount) / (1024.0 * 1024.0 * 1024.0);
	double sizeMembershipMatrix = (clusterCount * pixelCount * sizeof(float)) / (1024.0 * 1024.0 * 1024.0);
	std::cout << std::fixed << std::setprecision(3);
	std::cout << "Estimated RAM required: " << std::endl;
	std::cout << "Image size = " << sizeImage << "GB" << std::endl;
	std::cout << "Membership Matrix size = " << sizeMembershipMatrix << "GB" << std::endl;
	std::cout << "Spatial Membership Matrix size = " << sizeMembershipMatrix << "GB" << std::endl;
	std::cout << std::setprecision(6);
	std::cout.unsetf(std::ios::fixed);
	return;
}


/// <summary>
/// Checks if given parameters are legal/defined.
/// </summary>
/// <param name="clusterCount">How many classes to segment the image to.</param>
/// <param name="fuzzyParam">Fuzziness parameter of FCM.</param>
/// <param name="neighbourhoodSize">Apothem/radius of how many pixels to include in spatial neighbourhood calculations.</param>
/// <param name="iterMax">Maxmimum number of iterations to run the clustering algorithm.</param>
/// <param name="stopParam">Convergence tolerance threshold value.</param>
/// <param name="p">Weight parameter of pixel value.</param>
/// <param name="q">Weight parameter of spatial part.</param>
/// <param name="forceParams">Used to force acceptance of non legal parameters for debug purposes.</param>
static void checkParams(int clusterCount, float fuzzyParam, coord_t nbApothem, int iterMax, float stopParam, float p, float q, bool forceParams) {

	if (clusterCount < 2) {
		throw ITLException("clusterCount must be at least 2. " + std::to_string(clusterCount) + " was given.");
	}

	if (fuzzyParam < 1 && !forceParams) {
		throw ITLException("fuzzyParam must be at least 1. " + std::to_string(fuzzyParam) + " was given.");
	}

	if (nbApothem < 0 && !q == 0.0f) {
		throw ITLException("neighbourhoodSize must be a positive integer. (or 0, but use q=0.0f instead) " + std::to_string(nbApothem) + " was given.");
	}

	if (iterMax < 1) {
		throw ITLException("iterMax must be a positive integer. " + std::to_string(iterMax) + " was given.");
	}

	if (stopParam < 0.0f) {
		throw ITLException("stopParam must be a positive float. " + std::to_string(stopParam) + " was given.");
	}

	if (p < 0.0f && !forceParams) {
		throw ITLException("p must be a positive float or 0.0f. " + std::to_string(p) + " was given.");
	}

	if (q < 0.0f && !forceParams) {
		throw ITLException("q must be a positive float or 0.0f. " + std::to_string(q) + " was given.");
	}

	return;
}


/// <summary>
/// Main method to call to run sFCM/FCM/sHCM/K-Means segmentation on an image with one dimensional pixel feature.
/// </summary>
/// <typeparam name="pixel_t">Type of a pixel in the image.</typeparam>
/// <param name="imageData">Image to segment.</param>
/// <param name="clusterCount">How many classes to segment the image to.</param>
/// <param name="fuzzyParam">Fuzziness parameter of FCM.</param>
/// <param name="nbApothem">Apothem/radius of how many pixels to include in spatial neighbourhood calculations.</param>
/// <param name="iterMax">Maxmimum number of iterations to run the clustering algorithm.</param>
/// <param name="stopParam">Convergence tolerance threshold value.</param>
/// <param name="p">Weight parameter of pixel value.</param>
/// <param name="q">Weight parameter of spatial part.</param>
/// <param name="externalMembershipMatrix">membershipMatrix is saved to this if given.</param>
/// <param name="externalSpatialMembershipMatrix">spatialMembershipMatrix is saved to this if given.</param>
/// <param name="forceParams">Used to force acceptance of non legal parameters for debug purposes.</param>
template <typename pixel_t>
void doSFCM(Image<pixel_t>& imageData, int clusterCount, float fuzzyParam = 2.0f, coord_t nbApothem = 2, int iterMax = 100, float stopParam = 1e-3f, float p = 1.0f, float q = 1.0f, vector<std::unique_ptr<Image<float>>>* externalMembershipMatrix = nullptr, vector<std::unique_ptr<Image<float>>>* externalSpatialMembershipMatrix = nullptr, bool forceParams = false) {

	checkParams(clusterCount, fuzzyParam, nbApothem, iterMax, stopParam, p, q, forceParams);
	double pixelCount = static_cast<double>(imageData.pixelCount());
	calcRequiredMemory(pixelCount, clusterCount, imageData(0), 1);

	/*
	if (logging) {
		log_value("clusterCount", clusterCount);
		log_value("fuzzyParam", fuzzyParam);
		log_value("nbApothem", nbApothem);
		log_value("stopParam", stopParam);
		log_value("p", p);
		log_value("q", q);
	}
	*/

	SFCM sfcm(imageData, clusterCount, fuzzyParam, nbApothem, iterMax, stopParam, p, q);
	sfcm.cluster();
	sfcm.applySegmentation();

	if (externalMembershipMatrix) {
		*externalMembershipMatrix = std::move(sfcm.extractMembershipMatrix());
	}
	if (externalSpatialMembershipMatrix) {
		*externalSpatialMembershipMatrix = std::move(sfcm.extractSpatialMembershipMatrix());
	}

	return;
}


/// <summary>
/// Main method to call to run sFCM/FCM/sHCM/K-Means segmentation on an image with n features given as separate images inside a pointer vector.
/// Resulting segmentation is saved to the first image.
/// </summary>
/// <typeparam name="pixel_t">Type of a pixel in the image.</typeparam>
/// <param name="imageDataRefs">References to Images to segment.</param>
/// <param name="clusterCount">How many classes to segment the image to.</param>
/// <param name="fuzzyParam">Fuzziness parameter of FCM.</param>
/// <param name="neighbourhoodSize">Apothem/radius of how many pixels to include in spatial neighbourhood calculations.</param>
/// <param name="iterMax">Maxmimum number of iterations to run the clustering algorithm.</param>
/// <param name="stopParam">Convergence tolerance threshold value.</param>
/// <param name="p">Weight parameter of pixel value.</param>
/// <param name="q">Weight parameter of spatial part.</param>
/// <param name="externalMembershipMatrix">membershipMatrix is saved to this if given.</param>
/// <param name="externalSpatialMembershipMatrix">spatialMembershipMatrix is saved to this if given.</param>
/// <param name="forceParams">Used to force acceptance of non legal parameters for debug purposes.</param>
template <typename pixel_t>
void doSFCMVec(vector<reference_wrapper<Image<pixel_t>>> imageDataRefs, int clusterCount, float fuzzyParam = 2.0f, coord_t nbApothem = 2, int iterMax = 100, float stopParam = 1e-3f, float p = 1.0f, float q = 1.0f, vector<std::unique_ptr<Image<float>>>* externalMembershipMatrix = nullptr, vector<std::unique_ptr<Image<float>>>* externalSpatialMembershipMatrix = nullptr, bool forceParams = false) {
	checkParams(clusterCount, fuzzyParam, nbApothem, iterMax, stopParam, p, q, forceParams);
	double pixelCount = static_cast<double>(imageDataRefs[0].get().pixelCount());
	calcRequiredMemory(pixelCount, clusterCount, imageDataRefs[0].get()(0), static_cast<int>(imageDataRefs.size()));
	multiImagesFCM sfcm(imageDataRefs, clusterCount, fuzzyParam, nbApothem, iterMax, stopParam, p, q);
	sfcm.cluster();
	sfcm.applySegmentation();

	if (externalMembershipMatrix) {
		*externalMembershipMatrix = std::move(sfcm.extractMembershipMatrix());
	}
	if (externalSpatialMembershipMatrix) {
		*externalSpatialMembershipMatrix = std::move(sfcm.extractSpatialMembershipMatrix());
	}

	return;
}


/// <summary>
/// Main method to call to run sFCM/FCM/sHCM/K-Means segmentation on an image with 2 features given as separate images.
/// Resulting segmentation is saved to the first image.
/// </summary>
/// <typeparam name="pixel_t">Type of a pixel in the image.</typeparam>
/// <param name="imageData1">Image with feature 1. Result is applied to this image.</param>
/// <param name="imageData2">Image with feature 2.</param>
/// <param name="clusterCount">How many classes to segment the image to.</param>
/// <param name="fuzzyParam">Fuzziness parameter of FCM.</param>
/// <param name="neighbourhoodSize">Apothem/radius of how many pixels to include in spatial neighbourhood calculations.</param>
/// <param name="iterMax">Maxmimum number of iterations to run the clustering algorithm.</param>
/// <param name="stopParam">Convergence tolerance threshold value.</param>
/// <param name="p">Weight parameter of pixel value.</param>
/// <param name="q">Weight parameter of spatial part.</param>
/// <param name="externalMembershipMatrix">membershipMatrix is saved to this if given.</param>
/// <param name="externalSpatialMembershipMatrix">spatialMembershipMatrix is saved to this if given.</param>
/// <param name="forceParams">Used to force acceptance of non legal parameters for debug purposes.</param>
template <typename pixel_t>
void doSFCM(Image<pixel_t>& imageData1, Image<pixel_t>& imageData2, int clusterCount, float fuzzyParam = 2.0f, coord_t nbApothem = 2, int iterMax = 100, float stopParam = 1e-3f, float p = 1.0f, float q = 1.0f, vector<std::unique_ptr<Image<float>>>* externalMembershipMatrix = nullptr, vector<std::unique_ptr<Image<float>>>* externalSpatialMembershipMatrix = nullptr, bool forceParams = false) {

	vector<reference_wrapper<Image<pixel_t>>> imageDataRefs;
	imageDataRefs.reserve(2);
	imageDataRefs.push_back(ref(imageData1));
	imageDataRefs.push_back(ref(imageData2));

	doSFCMVec(imageDataRefs, clusterCount, fuzzyParam, nbApothem, iterMax, stopParam, p, q, externalMembershipMatrix, externalSpatialMembershipMatrix, forceParams);

	return;
}


/// <summary>
/// Main method to call to run sFCM/FCM/sHCM/K-Means segmentation on an image with 3 features given as separate images.
/// Resulting segmentation is saved to the first image.
/// </summary>
/// <typeparam name="pixel_t">Type of a pixel in the image.</typeparam>
/// <param name="imageData1">Image with feature 1. Result is applied to this image.</param>
/// <param name="imageData2">Image with feature 2.</param>
/// <param name="imageData3">Image with feature 3.</param>
/// <param name="clusterCount">How many classes to segment the image to.</param>
/// <param name="fuzzyParam">Fuzziness parameter of FCM.</param>
/// <param name="neighbourhoodSize">Apothem/radius of how many pixels to include in spatial neighbourhood calculations.</param>
/// <param name="iterMax">Maxmimum number of iterations to run the clustering algorithm.</param>
/// <param name="stopParam">Convergence tolerance threshold value.</param>
/// <param name="p">Weight parameter of pixel value.</param>
/// <param name="q">Weight parameter of spatial part.</param>
/// <param name="externalMembershipMatrix">membershipMatrix is saved to this if given.</param>
/// <param name="externalSpatialMembershipMatrix">spatialMembershipMatrix is saved to this if given.</param>
/// <param name="forceParams">Used to force acceptance of non legal parameters for debug purposes.</param>
template <typename pixel_t>
void doSFCM(Image<pixel_t>& imageData1, Image<pixel_t>& imageData2, Image<pixel_t>& imageData3, int clusterCount, float fuzzyParam = 2.0f, coord_t nbApothem = 2, int iterMax = 100, float stopParam = 1e-3f, float p = 1.0f, float q = 1.0f, vector<std::unique_ptr<Image<float>>>* externalMembershipMatrix = nullptr, vector<std::unique_ptr<Image<float>>>* externalSpatialMembershipMatrix = nullptr, bool forceParams = false) {

	vector<reference_wrapper<Image<pixel_t>>> imageDataRefs;
	imageDataRefs.reserve(3);
	imageDataRefs.push_back(ref(imageData1));
	imageDataRefs.push_back(ref(imageData2));
	imageDataRefs.push_back(ref(imageData3));

	doSFCMVec(imageDataRefs, clusterCount, fuzzyParam, nbApothem, iterMax, stopParam, p, q, externalMembershipMatrix, externalSpatialMembershipMatrix, forceParams);

	return;
}


/// <summary>
/// Main method to call to run sFCM/FCM/sHCM/K-Means segmentation on an image with 4 features given as separate images.
/// Resulting segmentation is saved to the first image.
/// </summary>
/// <typeparam name="pixel_t">Type of a pixel in the image.</typeparam>
/// <param name="imageData1">Image with feature 1. Result is applied to this image.</param>
/// <param name="imageData2">Image with feature 2.</param>
/// <param name="imageData3">Image with feature 3.</param>
/// <param name="imageData4">Image with feature 4.</param>
/// <param name="clusterCount">How many classes to segment the image to.</param>
/// <param name="fuzzyParam">Fuzziness parameter of FCM.</param>
/// <param name="neighbourhoodSize">Apothem/radius of how many pixels to include in spatial neighbourhood calculations.</param>
/// <param name="iterMax">Maxmimum number of iterations to run the clustering algorithm.</param>
/// <param name="stopParam">Convergence tolerance threshold value.</param>
/// <param name="p">Weight parameter of pixel value.</param>
/// <param name="q">Weight parameter of spatial part.</param>
/// <param name="externalMembershipMatrix">membershipMatrix is saved to this if given.</param>
/// <param name="externalSpatialMembershipMatrix">spatialMembershipMatrix is saved to this if given.</param>
/// <param name="forceParams">Used to force acceptance of non legal parameters for debug purposes.</param>
template <typename pixel_t>
void doSFCM(Image<pixel_t>& imageData1, Image<pixel_t>& imageData2, Image<pixel_t>& imageData3, Image<pixel_t>& imageData4, int clusterCount, float fuzzyParam = 2.0f, coord_t nbApothem = 2, int iterMax = 100, float stopParam = 1e-3f, float p = 1.0f, float q = 1.0f, vector<std::unique_ptr<Image<float>>>* externalMembershipMatrix = nullptr, vector<std::unique_ptr<Image<float>>>* externalSpatialMembershipMatrix = nullptr, bool forceParams = false) {

	vector<reference_wrapper<Image<pixel_t>>> imageDataRefs;
	imageDataRefs.reserve(4);
	imageDataRefs.push_back(ref(imageData1));
	imageDataRefs.push_back(ref(imageData2));
	imageDataRefs.push_back(ref(imageData3));
	imageDataRefs.push_back(ref(imageData4));

	doSFCMVec(imageDataRefs, clusterCount, fuzzyParam, nbApothem, iterMax, stopParam, p, q, externalMembershipMatrix, externalSpatialMembershipMatrix, forceParams);

	return;
}

