#pragma once
#include <vector>
#include "image.h"
#include <algorithm>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <limits>
#include <cmath>

//#include "logger.h"


using std::vector;
using std::unordered_map;
using std::unordered_set;
using itl2::Image;
using itl2::coord_t;


/// <summary>
/// Holds the results of the segmentation metrics.
/// </summary>
struct segmentationMetrics {
    double segValue;
    double accuracy;
    double sensitivity;
    double specificity;
    double dice;
};


/// <summary>
/// Matches labels of a segmented image to the groundtruth reference image based on unmatched the confusion matrix.
/// Matching is done so that the best match is used.
/// </summary>
/// <typeparam name="pixel_t">Type of the pixel value in groundtruth reference image.</typeparam>
/// <typeparam name="pixel_t2">Type of the pixel value in segmented image.</typeparam>
/// <param name="confMatrix">Confusion matrix</param>
/// <param name="gtLabels">Values of the class labels in the groundtruth reference image.</param>
/// <param name="segLabels">Values of the class labels in the segmented image.</param>
/// <returns>The best match for labels. (Which column to move to index which index in the confusion matrix.)</returns>
template <typename pixel_t, typename pixel_t2>
vector<int> matchLabelsBruteForce(const vector<vector<coord_t>>& confMatrix, const vector<pixel_t>& gtLabels, const vector<pixel_t2>& segLabels) {

    int classCount = gtLabels.size();
    vector<int> bestPermutation(classCount);
    vector<int> currentPermutation(classCount);

    for (int i = 0; i < classCount; i++)
        currentPermutation[i] = i;

    coord_t bestScore = -1;

    do {
        coord_t score = 0;
        for (int i = 0; i < classCount; i++) {
            score += confMatrix[i][currentPermutation[i]];
        }

        if (score > bestScore) {
            bestScore = score;
            bestPermutation = currentPermutation;
        }
    } while (std::next_permutation(currentPermutation.begin(), currentPermutation.end())); // Goes through all the different permutations.

    return bestPermutation;
}


/// <summary>
/// Calculates the confusion matrix given a groundtruth reference image and a segmented image.
/// </summary>
/// <typeparam name="pixel_t">Type of the pixel value in groundtruth reference image.</typeparam>
/// <typeparam name="pixel_t2">Type of the pixel value in segmented image.</typeparam>
/// <param name="gtImage">Groundtruth reference image.</param>
/// <param name="segImage">Segmented Image</param>
/// <param name="confMatrix">Confusion matrix to save the calculated comfusion matrix into.</param>
/// <param name="uniqueGT">Vector to hold all different unique class labels in groundtruth reference image.</param>
/// <param name="uniqueSeg">Vector to hold all different unique class labels in segmented image.</param>
template <typename pixel_t, typename pixel_t2>
void computeConfusionMatrix(const Image<pixel_t>& gtImage, const Image<pixel_t2>& segImage, vector<vector<coord_t>>& confMatrix, vector<pixel_t>& uniqueGT, vector<pixel_t2>& uniqueSeg) {

    unordered_set<pixel_t> setGT;
    unordered_set<pixel_t2> setSeg;
    for (size_t i = 0; i < gtImage.pixelCount(); i++) {
        setGT.insert(gtImage(i));
        setSeg.insert(segImage(i));
    }
    uniqueGT.assign(setGT.begin(), setGT.end());
    uniqueSeg.assign(setSeg.begin(), setSeg.end());
    sort(uniqueGT.begin(), uniqueGT.end());
    sort(uniqueSeg.begin(), uniqueSeg.end());

    if (uniqueGT.size() != uniqueSeg.size()) {
        throw ITLException("Groundtruth image has " + std::to_string(uniqueGT.size()) + " unique classes/values, while segmented image has " + std::to_string(uniqueSeg.size()));
    }

    int classCount = uniqueGT.size();

    confMatrix.assign(classCount, vector<coord_t>(classCount, 0));

    // label = value -> index (maps value to index)
    unordered_map<pixel_t, int> gtLabelToIndex;
    for (int i = 0; i < classCount; i++)
        gtLabelToIndex[uniqueGT[i]] = i;

    unordered_map<pixel_t2, int> segLabelToIndex;
    for (int j = 0; j < classCount; j++)
        segLabelToIndex[uniqueSeg[j]] = j;

    for (size_t i = 0; i < gtImage.pixelCount(); i++) {
        confMatrix[gtLabelToIndex[gtImage(i)]][segLabelToIndex[segImage(i)]] += 1;
    }
}


/// <summary>
/// Calculates the accuracy segmentation metric from the given cnfusion matrix.
/// </summary>
/// <typeparam name="pixel_t">Type of the pixel value in groundtruth reference image.</typeparam>
/// <param name="confMatrix">Confusion matrix</param>
/// <param name="uniqueGT">Vector holding all different unique class labels.</param>
/// <param name="pixelCount">Amount of voxels in the image.</param>
/// <returns>Accuracy results as an unordered map. Class label gives the accuracy result of the given label.</returns>
template <typename pixel_t>
unordered_map<pixel_t, double> accuracyPerLabel(const vector<vector<coord_t>>& confMatrix, const vector<pixel_t>& gtLabels, const coord_t pixelCount) {

    unordered_map<pixel_t, double> accuracyPerLabel;
    int classCount = gtLabels.size();

    for (int i = 0; i < classCount; ++i) {

        coord_t TP = confMatrix[i][i];

        coord_t FN = 0;
        for (int j = 0; j < classCount; j++) {
            if (j != i) {
                FN += confMatrix[i][j];
            }
        }

        coord_t FP = 0;
        for (int k = 0; k < classCount; k++) {
            if (k != i) {
                FP += confMatrix[k][i];
            }
        }

        coord_t TN = pixelCount - TP - FP - FN;

        double acc = static_cast<double>(TP + TN) / static_cast<double>(pixelCount);

        /*
        double acc = (TP + TN == 0) ? 1.0 // The only way to get 100% incorrect is to know what was 100% correct. This shouldn't happen though, since labels should have been already matched...
            : static_cast<double>(TP + TN) / static_cast<double>(pixelCount);
        */

        accuracyPerLabel[gtLabels[i]] = acc;

    }

    return accuracyPerLabel;
}


/// <summary>
/// Calculates the sensitivity segmentation metric from the given cnfusion matrix.
/// </summary>
/// <typeparam name="pixel_t">Type of the pixel value in groundtruth reference image.</typeparam>
/// <param name="confMatrix">Confusion matrix</param>
/// <param name="uniqueGT">Vector holding all different unique class labels.</param>
/// <param name="pixelCount">Amount of voxels in the image.</param>
/// <returns>Sensitivity results as an unordered map. Class label gives the sensitivity result of the given label.</returns>
template <typename pixel_t>
unordered_map<pixel_t, double> sensitivityPerLabel(const vector<vector<coord_t>>& confMatrix, const vector<pixel_t>& gtLabels) {

    unordered_map<pixel_t, double> sensitivityPerLabel;
    int classCount = gtLabels.size();

    for (int i = 0; i < classCount; i++) {

        coord_t TP = confMatrix[i][i];

        coord_t FN = 0;
        for (int j = 0; j < classCount; j++) {
            if (j != i) {
                FN += confMatrix[i][j];
            }
        }

        double sens = static_cast<double>(TP) / static_cast<double>(TP + FN);

        /*
        double sens = (TP + FN == 0) ? 1.0 // The only way to get 100% incorrect is to know what was 100% correct. This shouldn't happen though, since labels should have been already matched...
            : static_cast<double>(TP) / static_cast<double>(TP + FN);
        */

        sensitivityPerLabel[gtLabels[i]] = sens;
    }

    return sensitivityPerLabel;
}


/// <summary>
/// Calculates the specificity segmentation metric from the given cnfusion matrix.
/// </summary>
/// <typeparam name="pixel_t">Type of the pixel value in groundtruth reference image.</typeparam>
/// <param name="confMatrix">Confusion matrix</param>
/// <param name="uniqueGT">Vector holding all different unique class labels.</param>
/// <param name="pixelCount">Amount of voxels in the image.</param>
/// <returns>Specificity results as an unordered map. Class label gives the specificity result of the given label.</returns>
template <typename pixel_t>
unordered_map<pixel_t, double> specificityPerLabel(const vector<vector<coord_t>>& confMatrix, const vector<pixel_t>& gtLabels, const coord_t pixelCount) {

    unordered_map<pixel_t, double> specificityPerLabel;
    int classCount = gtLabels.size();

    for (int i = 0; i < classCount; i++) {

        coord_t TP = confMatrix[i][i];

        coord_t FN = 0;
        for (int j = 0; j < classCount; j++) {
            if (j != i) {
                FN += confMatrix[i][j];
            }
        }

        coord_t FP = 0;
        for (int k = 0; k < classCount; k++) {
            if (k != i) {
                FP += confMatrix[k][i];
            }
        }

        coord_t TN = pixelCount - TP - FP - FN;

        double spec = static_cast<double>(TN) / static_cast<double>(TN + FP);

        /*
        double spec = (TN + FP == 0) ? 1.0 // The only way to get 100% incorrect is to know what was 100% correct. This shouldn't happen though, since labels should have been already matched...
            : static_cast<double>(TN) / static_cast<double>(TN + FP);
        */

        specificityPerLabel[gtLabels[i]] = spec;
    }

    return specificityPerLabel;
}


/// <summary>
/// Calculates the dice coefficent segmentation metric from the given cnfusion matrix.
/// </summary>
/// <typeparam name="pixel_t">Type of the pixel value in groundtruth reference image.</typeparam>
/// <param name="confMatrix">Confusion matrix</param>
/// <param name="uniqueGT">Vector holding all different unique class labels.</param>
/// <param name="pixelCount">Amount of voxels in the image.</param>
/// <returns>Dice coefficent results as an unordered map. Class label gives the dice coefficent result of the given label.</returns>
template <typename pixel_t>
unordered_map<pixel_t, double> diceCoefficientPerLabel(const vector<vector<coord_t>>& confMatrix, const vector<pixel_t>& gtLabels, const coord_t pixelCount) {

    unordered_map<pixel_t, double> dicePerLabel;
    int classCount = gtLabels.size();

    for (int i = 0; i < classCount; i++) {

        coord_t TP = confMatrix[i][i];

        coord_t FN = 0;
        for (int j = 0; j < classCount; j++) {
            if (j != i) {
                FN += confMatrix[i][j];
            }
        }

        coord_t FP = 0;
        for (int k = 0; k < classCount; k++) {
            if (k != i) {
                FP += confMatrix[k][i];
            }
        }

        coord_t denominator = 2 * TP + FN + FP;

        double dice = (2.0 * TP) / denominator;

        /*
        double dice = (denominator == 0) ? 1.0 // The only way to get 100% incorrect is to know what was 100% correct. This shouldn't happen though, since labels should have been already matched...
            : (2.0 * TP) / denominator;
        */

        dicePerLabel[gtLabels[i]] = dice;
    }

    return dicePerLabel;
}


/// <summary>
/// Calculates confusion matrix segmentation metrics (Accuracy, Sensitivity, Specificity, Dice Coefficent).
/// </summary>
/// <typeparam name="pixel_t">Type of the pixel value in groundtruth reference image.</typeparam>
/// <typeparam name="pixel_t2">Type of the pixel value in segmented image.</typeparam>
/// <param name="gtImage">Groundtruth reference image.</param>
/// <param name="segImage">Segmented image.</param>
/// <returns>The results of the segmentation metrics as an unordered map: results[labelValue].metric </returns>
template <typename pixel_t, typename pixel_t2>
unordered_map<pixel_t, segmentationMetrics> calculateSegmentationMetrics(const Image<pixel_t>& gtImage, const Image<pixel_t2>& segImage) {

    vector<vector<coord_t>> confMatrix;
    vector<pixel_t> uniqueGT;
    vector<pixel_t2> uniqueSeg;
    computeConfusionMatrix(gtImage, segImage, confMatrix, uniqueGT, uniqueSeg);
    int classCount = confMatrix.size();

    /*
    cout << "Confusion Matrix (rows: GT, cols: Seg):" << endl;
    for (int i = 0; i < classCount; i++) {
        for (int j = 0; j < classCount; j++) {
            cout << confMatrix[i][j] << "\t";
        }
        cout << endl;
    }
    */

    vector<int> bestMatch = matchLabelsBruteForce(confMatrix, uniqueGT, uniqueSeg);

    vector<vector<coord_t>> matchedConfMatrix(classCount, vector<coord_t>(classCount, 0));

    for (int i = 0; i < classCount; i++) {
        int segCol = bestMatch[i]; // Which column to move to index i.
        for (int row = 0; row < classCount; row++) {
            matchedConfMatrix[row][i] = confMatrix[row][segCol];
        }
    }

    /*
    cout << "Confusion Matrix after matching labels (rows: GT, cols: Seg):" << endl;
    for (int i = 0; i < classCount; i++) {
        for (int j = 0; j < classCount; j++) {
            cout << matchedConfMatrix[i][j] << "\t";
        }
        cout << endl;
    }
    */

    coord_t pixelCount = gtImage.pixelCount();
    unordered_map<pixel_t, double> acc = accuracyPerLabel(matchedConfMatrix, uniqueGT, pixelCount);
    unordered_map<pixel_t, double> sens = sensitivityPerLabel(matchedConfMatrix, uniqueGT);
    unordered_map<pixel_t, double> spec = specificityPerLabel(matchedConfMatrix, uniqueGT, pixelCount);
    unordered_map<pixel_t, double> dice = diceCoefficientPerLabel(matchedConfMatrix, uniqueGT, pixelCount);

    unordered_map<pixel_t, segmentationMetrics> results;

    for (int i = 0; i < classCount; i++) {

        std::cout << "  GT Label:     " << uniqueGT[i] << "\n";
        std::cout << "  Seg Label:    " << uniqueSeg[bestMatch[i]] << "\n";
        results[uniqueGT[i]].segValue = uniqueSeg[bestMatch[i]];
        std::cout << "  Accuracy:    " << acc[uniqueGT[i]] << "\n";
        results[uniqueGT[i]].accuracy = acc[uniqueGT[i]];
        std::cout << "  Sensitivity: " << sens[uniqueGT[i]] << "\n";
        results[uniqueGT[i]].sensitivity = sens[uniqueGT[i]];
        std::cout << "  Specificity: " << spec[uniqueGT[i]] << "\n";
        results[uniqueGT[i]].specificity = spec[uniqueGT[i]];
        std::cout << "  Dice:        " << dice[uniqueGT[i]] << "\n";
        results[uniqueGT[i]].dice = dice[uniqueGT[i]];

        /*
        if (logging) {
            log_label_result(uniqueGT[i], results[uniqueGT[i]].segValue, results[uniqueGT[i]].accuracy, results[uniqueGT[i]].sensitivity, results[uniqueGT[i]].specificity, results[uniqueGT[i]].dice);
        }
        */
    }

    /*
    if (logging) {
        log_run_separator();
    }
    */

    return results;
}

