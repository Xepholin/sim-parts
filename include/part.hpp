#pragma once

#include <array>
#include <memory>
#include <string>

#include "const.hpp"

using namespace std;

#pragma pack()
struct alignas(32) data_t {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

#pragma pack()
struct alignas(32) dataArray_t {
    array<double, TOTAL_PARTS> x;
    array<double, TOTAL_PARTS> y;
    array<double, TOTAL_PARTS> z;
};

const array<array<double, 3>, 27> vec_translation = {
    array<double, 3>{0.0, -l, -l},
    array<double, 3>{0.0, -l, 0.0},
    array<double, 3>{0.0, -l, l},
    array<double, 3>{0.0, 0.0, -l},
    array<double, 3>{0.0, 0.0, 0.0},
    array<double, 3>{0.0, 0.0, l},
    array<double, 3>{0.0, l, -l},
    array<double, 3>{0.0, l, 0.0},
    array<double, 3>{0.0, l, l},
    array<double, 3>{-l, -l, -l},
    array<double, 3>{-l, -l, 0.0},
    array<double, 3>{-l, -l, l},
    array<double, 3>{-l, 0.0, -l},
    array<double, 3>{-l, 0.0, 0.0},
    array<double, 3>{-l, 0.0, l},
    array<double, 3>{-l, l, -l},
    array<double, 3>{-l, l, 0.0},
    array<double, 3>{-l, l, l},
    array<double, 3>{l, -l, -l},
    array<double, 3>{l, -l, 0.0},
    array<double, 3>{l, -l, l},
    array<double, 3>{l, 0.0, -l},
    array<double, 3>{l, 0.0, 0.0},
    array<double, 3>{l, 0.0, l},
    array<double, 3>{l, l, -l},
    array<double, 3>{l, l, 0.0},
    array<double, 3>{l, l, l}};

template <typename T>
void fill_vec(std::unique_ptr<T> &__restrict__ parts, const int n, const std::string &path);

template <typename T>
double energy_forces(const unique_ptr<T> &__restrict__ parts, unique_ptr<T> &__restrict__ forces, const int n_sym, const int n);
