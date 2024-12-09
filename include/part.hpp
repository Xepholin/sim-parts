#pragma once

#include <array>
#include <memory>
#include <string>

#include "part.hpp"

#define TOTAL_PARTS 1000
#define LOCAL_PARTS 999

const size_t fsize = 3 * (TOTAL_PARTS * TOTAL_PARTS);

using namespace std;

#pragma pack()
struct Particule {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

void fill_vec(unique_ptr<array<Particule, TOTAL_PARTS>> &parts, const int n, const std::string path);
double energy_forces(const unique_ptr<array<Particule, TOTAL_PARTS>> &parts, unique_ptr<Particule> &forces, const int n);
double energy_forces_periode(const unique_ptr<array<Particule, TOTAL_PARTS>> &parts, unique_ptr<Particule> &forces, const int n_sym, const double lx, const double ly, const double lz, const int n);