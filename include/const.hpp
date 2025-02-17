#pragma once

#include <array>
#include <cmath>
#include <string>

using namespace std;

static string DEFAULT_OUTPUT = "../output.pdb";

const double r_star_sq = 9.0;
const double epsilon_star = 0.2;

const double r_cut = 10.0;
const double r_cut_sq = 100.0;
const double l = 42.0;

const double conv_force = 0.0001 * 4.186;
const double const_r = 0.00199;

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
