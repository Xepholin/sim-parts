#pragma once

#include <string>
#include <vector>

#include "part.hpp"

using namespace std;

struct Particule {
    double x = 0;
    double y = 0;
    double z = 0;

    Particule() = default;
    ~Particule() = default;
};

void *fill_vec(vector<Particule> *parts, int n, std::string const path);
double compute_lj(const vector<Particule> * parts, int n);
void compute_forces(vector<double> *forces, const vector<Particule> * parts, int n);