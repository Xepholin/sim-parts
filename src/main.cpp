#include <algorithm>
#include <array>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <memory>
#include <ratio>
#include <vector>

#include "part.hpp"
#include "tools.hpp"

#define NB_META 33
#define NB_WARM 100
#define NB_ITER 33

static string DEFAULT_PARTS_FILE = "../particule.xyz";

using namespace std;

int main() {
    alignas(alignof(Particule)) auto parts = make_unique<array<Particule, TOTAL_PARTS>>();
    alignas(alignof(Particule)) auto forces = make_unique<Particule>();
    alignas(alignof(Particule)) auto forces_per = make_unique<Particule>();

    double lj = 0.0;
    double lj_sym = 0.0;

    fill_vec(parts, TOTAL_PARTS, DEFAULT_PARTS_FILE);

    cout << "           AoS           " << endl;

    auto t1 = chrono::high_resolution_clock::now();
    // for (int iter = 0; iter < NB_ITER; ++iter) {
    lj = energy_forces(parts, forces, TOTAL_PARTS);
    // }
    auto t2 = chrono::high_resolution_clock::now();

    double fsum = forces->x + forces->y + forces->z;
    double elapsed = chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / NB_ITER;

    cout << fixed << "Lennard Jones: " << lj << endl;
    cout << scientific << "Somme des forces: " << fsum << endl;
    cout << "Temps: " << elapsed * 1e-9 << " s\n" << endl;

    t1 = chrono::high_resolution_clock::now();
    // for (int iter = 0; iter < NB_ITER; ++iter) {
    lj_sym = energy_forces_periode(parts, forces_per, 27.0, 42.0, 42.0, 42.0, TOTAL_PARTS);
    // }
    t2 = chrono::high_resolution_clock::now();

    double fsum_per = forces_per->x + forces_per->y + forces_per->z;
    elapsed = chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / NB_ITER;

    cout << fixed << "Lennard Jones sym 27: " << lj_sym << endl;
    cout << scientific << "Somme des forces sym 27: " << fsum_per << endl;
    cout << "Temps: " << elapsed * 1e-9 << " s\n" << endl;

    return 0;
}

// int main() {
//     alignas(16) auto parts = make_unique<array<Particule, TOTAL_PARTS>>();
//     alignas(16) auto forces = make_unique<array<double, fsize>>();

//     auto elapsed = make_unique<vector<double>>();
//     auto ljs = make_unique<vector<double>>();
//     auto fsums = make_unique<vector<double>>();

//     double lj = 0.0;

//     for (int meta = 0; meta < NB_META; ++meta) {
//         if (!meta) {
//             for (int warm = 0; warm < NB_WARM; ++warm) {
//                 fill_vec(parts, TOTAL_PARTS, DEFAULT_PARTS_FILE);
//                 lj = compute_lj(parts, TOTAL_PARTS);
//                 compute_forces(forces, parts, TOTAL_PARTS);
//             }
//         } else {
//             fill_vec(parts, TOTAL_PARTS, DEFAULT_PARTS_FILE);

//             auto t1 = chrono::high_resolution_clock::now();

//             for (int iter = 0; iter < NB_ITER; ++iter) {
//                 lj = compute_lj(parts, TOTAL_PARTS);
//             }

//             auto t2 = chrono::high_resolution_clock::now();
//             compute_forces(forces, parts, TOTAL_PARTS);

//             elapsed->push_back(chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / NB_ITER);

//             ljs->push_back(lj);

//             double fsum = 0.0;

//             for (size_t i = 0; i < fsize; ++i) {
//                 fsum += (*forces)[i];
//             }

//             fsums->push_back(fsum);
//         }
//     }

//     cout << fixed << "Lennard Jones moy.: " << mean(ljs, NB_META - 1) << endl;
//     cout << scientific << "Somme moy. des forces: " << mean(fsums, NB_META - 1) << endl;
//     cout << "Temps moy.: " << mean(elapsed, NB_META - 1) * 1e-9 << " s" << endl;

//     return 0;
// }
