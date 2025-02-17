#pragma once

#include <array>
#include <cmath>
#include <string>

#define TOTAL_PARTS 1000
#define LOCAL_PARTS 999

using namespace std;

static string DEFAULT_PARTS_FILE = "../particule.xyz";
static string OUTPUT_PARTS_FILE = "../output.pdb";

const double n_sym = 27.0;

const double r_star_sq = 9.0;
const double epsilon_star = 0.2;

const double r_cut = 10.0;
const double r_cut_sq = 100.0;
const double l = 42.0;

const double conv_force = 0.0001 * 4.186;
const double n_dl = 3.0 * TOTAL_PARTS - 3.0;
const double const_r = 0.00199;

const double n_neighbor = (double)TOTAL_PARTS / (l * l * l);
const int n_max_neighbor = (int)(n_neighbor * 4.0 * M_PI * (r_cut_sq * r_cut)) * 2.0;

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
