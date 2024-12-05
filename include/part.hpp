#pragma once

#include <string>
#include <array>

#include "part.hpp"

#define TOTAL_PARTS 1000
#define LOCAL_PARTS 999

const size_t fsize = 3* (TOTAL_PARTS * TOTAL_PARTS);

using namespace std;

struct Particule {
    double x = 0;
    double y = 0;
    double z = 0;

    Particule() = default;
    ~Particule() = default;
};

void *fill_vec(array<Particule, TOTAL_PARTS> *parts, int n, std::string const path);
double compute_lj(const array<Particule, TOTAL_PARTS> * parts, int n);
void compute_forces(array<double, fsize> *forces, const array<Particule, TOTAL_PARTS> * parts, int n);