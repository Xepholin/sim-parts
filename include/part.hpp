#pragma once

#include <array>
#include <memory>
#include <string>

#include "part.hpp"

#define TOTAL_PARTS 1000
#define LOCAL_PARTS 999

using namespace std;

const double l = 42.0;

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

const array<array<double, 3>, 27> vec_trans = {
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
    array<double, 3>{l, l, l}
};

void fill_vec_aos(unique_ptr<array<data_t, TOTAL_PARTS>> &__restrict__ parts, const int n, const string path);
double energy_forces_aos(const unique_ptr<array<data_t, TOTAL_PARTS>> &__restrict__ parts, unique_ptr<array<data_t, TOTAL_PARTS>> &__restrict__ forces, const int n);
double energy_forces_periode_aos(const unique_ptr<array<data_t, TOTAL_PARTS>> &__restrict__ parts, unique_ptr<array<data_t, TOTAL_PARTS>> &__restrict__ forces, const int n_sym, const double l, const int n);

void fill_vec_soa(unique_ptr<dataArray_t> &__restrict__ parts, const int n, const string path);
double energy_forces_soa(const unique_ptr<dataArray_t> &__restrict__ parts, unique_ptr<dataArray_t> &__restrict__ forces, const int n);
double energy_forces_periode_soa(const unique_ptr<dataArray_t> &__restrict__ parts, unique_ptr<dataArray_t> &__restrict__forces, const int n_sym, const double l, const int n);