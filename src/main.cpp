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

using namespace std;

int main() {
    Simulator sim(TOTAL_PARTS, n_sym, 0.0, 3000.0, 1.0, 300.0, 0.01);
    // double lj = 0.0;

    sim.fillVec(DEFAULT_PARTS_FILE);

    // auto t1 = chrono::high_resolution_clock::now();
    // for (int iter = 0; iter < NB_ITER; ++iter) {
    //     lj = sim.computeEnergyForces();
    // }
    // auto t2 = chrono::high_resolution_clock::now();

    // double elapsed = chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / NB_ITER;

    cout << "Somme des forces sym " << sim.get_n_sym() << ": " << scientific << sim.sumForces() << endl;
    cout << fixed << "Lennard Jones sym " << sim.get_n_sym() << ": " << sim.computeEnergyForces() << endl;
    // cout << "Temps: " << elapsed * 1e-9 << " s\n"
    //      << endl;

    sim.fillMoment();

    sim.start(true, 20, 5);
    return 0;
}
