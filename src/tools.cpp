#include "tools.hpp"

#include <cmath>
#include <numeric>

#include "const.hpp"
#include "part.hpp"

using namespace std;

#include <vector>
#include <algorithm>

template <typename T>
void init_f64(const int n, T &__restrict__ a, const double val) {
    if constexpr (std::is_same_v<T, std::vector<data_t>>) {
        for (int i = 0; i < n; ++i) {
            a[i].x = val;
            a[i].y = val;
            a[i].z = val;
        }
    } else if constexpr (std::is_same_v<T, std::vector<double>>) {
        std::fill(a.begin(), a.begin() + n, val);
    }
}

template void init_f64<std::vector<data_t>>(int, std::vector<data_t> &__restrict__, const double);
template void init_f64<std::vector<double>>(int, std::vector<double> &__restrict__, const double);


double mean(const int n, const unique_ptr<vector<double>> &__restrict__ arr) {
    return reduce(arr->begin(), arr->end()) / n;
}

double stddev(const int n, const unique_ptr<vector<double>> &__restrict__ arr) {
    const double avg = mean(n, arr);

    auto varianceOp = [avg, valn = arr->size()](double accumulator, double val) {
        return accumulator += (pow(val - avg, 2) / static_cast<double>(valn));
    };

    return sqrt(accumulate(arr->begin(), arr->end(), 0.0, varianceOp));
}

double fn_sign(double a, double s) {
    return (s < 0.0) ? a : -a;
}
