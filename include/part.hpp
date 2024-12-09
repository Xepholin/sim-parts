#pragma once

#include <array>
#include <memory>
#include <string>

#include "part.hpp"

#define TOTAL_PARTS 1000
#define LOCAL_PARTS 999

using namespace std;

#pragma pack()
struct Particule {
    array<double, TOTAL_PARTS> x;
    array<double, TOTAL_PARTS> y;
    array<double, TOTAL_PARTS> z;
};

void fill_vec(unique_ptr<Particule> &parts, const int n, const std::string path);
double energy_forces(const unique_ptr<Particule> &parts, unique_ptr<array<double, 3>> &forces, const int n);
double energy_forces_periode(const unique_ptr<Particule> &parts, unique_ptr<array<double, 3>> &forces, const int n_sym, const double lx, const double ly, const double lz, const int n);