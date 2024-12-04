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
void compute_force(vector<double> *forces, const vector<Particule> *parts, int n);