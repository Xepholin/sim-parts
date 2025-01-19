#include <algorithm>
#include <array>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <memory>
#include <ratio>
#include <vector>

#include "const.hpp"
#include "part.hpp"
#include "tools.hpp"

#define NB_META 33
#define NB_WARM 100
#define NB_ITER 1

static string DEFAULT_PARTS_FILE = "../particule.xyz";

using namespace std;

int main() {
    auto parts_aos = make_unique<array<data_t, TOTAL_PARTS>>();
    auto forces_aos = make_unique<array<data_t, TOTAL_PARTS>>();

    auto parts_soa = make_unique<dataArray_t>();
    auto forces_soa = make_unique<dataArray_t>();

    double lj_aos = 0.0;
    double lj_soa = 0.0;

    empty_vec(TOTAL_PARTS, parts_aos);
    empty_vec(TOTAL_PARTS, forces_aos);
    empty_vec(TOTAL_PARTS, parts_soa);
    empty_vec(TOTAL_PARTS, forces_soa);

    fill_vec(parts_aos, TOTAL_PARTS, DEFAULT_PARTS_FILE);

    cout << "           AoS           " << endl;

    auto t1 = chrono::high_resolution_clock::now();
    for (int iter = 0; iter < NB_ITER; ++iter) {
        lj_aos = energy_forces(parts_aos, forces_aos, 27, TOTAL_PARTS);
    }
    auto t2 = chrono::high_resolution_clock::now();

    double elapsed = chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / NB_ITER;

    cout << fixed << "Lennard Jones sym 27: " << lj_aos << endl;
    cout << "Temps: " << elapsed * 1e-9 << " s\n"
         << endl;

    fill_vec(parts_soa, TOTAL_PARTS, DEFAULT_PARTS_FILE);

    cout << "           SoA           " << endl;

    t1 = chrono::high_resolution_clock::now();
    for (int iter = 0; iter < NB_ITER; ++iter) {
        lj_soa = energy_forces(parts_soa, forces_soa, 27, TOTAL_PARTS);
    }
    t2 = chrono::high_resolution_clock::now();

    elapsed = chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / NB_ITER;

    cout << fixed << "Lennard Jones sym 27: " << lj_soa << endl;
    cout << "Temps: " << elapsed * 1e-9 << " s\n"
         << endl;

    return 0;
}
