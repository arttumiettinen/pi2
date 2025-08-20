#pragma once
#include <fstream>
#include <string>
#include <iomanip>

using std::string;

inline bool logging = false;


/// <summary>
/// Used to set logging on or off.
/// </summary>
/// <param name="value">true/false for logging.</param>
inline void setLogging(bool value) {
    logging = value;
}


/// <summary>
/// Logs a value to the logfile.
/// </summary>
/// <param name="key">Name of the value to log.</param>
/// <param name="value">Value of the value to log.</param>
/// <param name="filename">Filename to log to.</param>
inline void log_value(const string& key, double value, const string& filename = "logFile.txt") {
    std::ofstream file(filename, std::ios::app);
    if (file) {
        file << key << "=" << std::fixed << std::setprecision(10) << value << "\n";
    }
}


/// <summary>
/// Adds a "---" line to the logfile.
/// </summary>
/// <param name="filename">Filename to log to.</param>
inline void log_run_separator(const string& filename = "logFile.txt") {
    std::ofstream file(filename, std::ios::app);
    if (file) {
        file << "---\n";
    }
}


/// <summary>
/// Logs segmentation metrics in specific format.
/// </summary>
/// <param name="gtLabel">Class label in groundtruth reference image.</param>
/// <param name="segLabel">Class label in segmented image.</param>
/// <param name="accuracy">Accuracy metric.</param>
/// <param name="sensitivity">Sensitivity metric.</param>
/// <param name="specificity">Specificity metric.</param>
/// <param name="dice">Dice coefficent metric.</param>
/// <param name="filename">Filename to log to.</param>
inline void log_label_result(int gtLabel, int segLabel, double accuracy, double sensitivity, double specificity, double dice, const string& filename = "logFile.txt") {
    std::ofstream file(filename, std::ios::app);
    if (file) {
        file << "LabelResult {\n";
        file << "  GtLabel=" << gtLabel << "\n";
        file << "  segLabel=" << segLabel << "\n";
        file << "  Accuracy=" << accuracy << "\n";
        file << "  Sensitivity=" << sensitivity << "\n";
        file << "  Specificity=" << specificity << "\n";
        file << "  Dice=" << dice << "\n";
        file << "}\n";
    }
}
