#pragma once
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>
#include <cmath>

///
/// Converts the argument to string in an OS-independent manner
/// The point here is that NaN and inf may be printed differently, this
/// takes care of that.
///
void printToStream(std::ostream &out, double in) {
    if (!std::isnan(in)) {
        out << in;
    } else {
        out << "nan";
    }
}

///
/// Writes the contents of the vector 'data' to the textfile filename
/// The output format should be load-able in MATLAB and Numpy using the
/// load-command.
///
template <typename T>
void writeToFile(const std::string &filename, const T &data) {

    std::ofstream file(filename.c_str());
    // Set highest possible precision, this way we are sure we are
    file << std::setprecision(std::numeric_limits<long double>::digits10 + 1);

    // Loop over vector and write output to file
    for (int i = 0; i < data.size(); ++i) {
        printToStream(file, data[i]);
        file << " ";
    }
    file << std::endl;

    // File closes automatically at end of scope!
}

template <typename T>
void writeMatrixToFile(const std::string &filename, const T &data) {
    std::ofstream file(filename.c_str());
    // Set highest possible precision, this way we are sure we are
    file << std::setprecision(std::numeric_limits<long double>::digits10 + 1);

    // Loop over vector and write output to file
    for (int i = 0; i < data.rows(); ++i) {
        for (int j = 0; j < data.cols(); ++j) {
            printToStream(file, data(i, j));
            file << " ";
        }
        file << std::endl;
    }
    file << std::endl;

    // File closes automatically at end of scope!
}
