#include "tools.hpp"

#include <cmath>
#include <numeric>

using namespace std;

double mean(const unique_ptr<vector<double>> &arr, const int size) {
    return reduce(arr->begin(), arr->end()) / size;
}

double stddev(const unique_ptr<vector<double>> &arr, const int size) {
    const double avg = mean(arr, size);

    auto varianceOp = [avg, valSize = arr->size()](double accumulator, double val) {
        return accumulator += (pow(val - avg, 2) / static_cast<double>(valSize));
    };

    return sqrt(accumulate(arr->begin(), arr->end(), 0.0, varianceOp));
}