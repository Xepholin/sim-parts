#include "tools.hpp"

#include <cmath>
#include <numeric>

#include "const.hpp"
#include "part.hpp"

using namespace std;

template <typename T>
void empty_vec(const int n, unique_ptr<T> &__restrict__ vec) {
    if constexpr (is_same_v<T, array<data_t, TOTAL_PARTS>>) {
        for (int i = 0; i < n; ++i) {
            (*vec)[i].x = 0.0;
            (*vec)[i].y = 0.0;
            (*vec)[i].z = 0.0;
        }
    } else if constexpr (is_same_v<T, dataArray_t>) {
        fill(vec->x.begin(), vec->x.begin() + n, 0.0);
        fill(vec->y.begin(), vec->y.begin() + n, 0.0);
        fill(vec->z.begin(), vec->z.begin() + n, 0.0);
    }
}

template void empty_vec<array<data_t, TOTAL_PARTS>>(int, unique_ptr<array<data_t, TOTAL_PARTS>> &__restrict__);
template void empty_vec<dataArray_t>(int, unique_ptr<dataArray_t> &__restrict__);

double mean(const unique_ptr<vector<double>> &__restrict__ arr, const int size) {
    return reduce(arr->begin(), arr->end()) / size;
}

double stddev(const unique_ptr<vector<double>> &__restrict__ arr, const int size) {
    const double avg = mean(arr, size);

    auto varianceOp = [avg, valSize = arr->size()](double accumulator, double val) {
        return accumulator += (pow(val - avg, 2) / static_cast<double>(valSize));
    };

    return sqrt(accumulate(arr->begin(), arr->end(), 0.0, varianceOp));
}
