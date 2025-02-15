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

using namespace std;

int main() {
    Simulator sim(TOTAL_PARTS, n_sym, 0.0, 4000.0, 1, 300.0, 0.01, true);
    streamsize ss = std::cout.precision();

    cout << "\n=================================\n"
         << "      Début de la simulation\n"
         << "=================================\n"
         << endl;

    sim.fillVec(DEFAULT_PARTS_FILE);
    sim.verletList();

    cout << "Somme des forces sym " << sim.get_n_sym() << ": " << scientific << sim.sumForces() << endl;
    cout << fixed << "Lennard Jones sym " << sim.get_n_sym() << ": " << sim.computeEnergyForces() << endl;

    cout << setprecision(2)
         << "Temps initial: " << sim.get_t_min()
         << "\nTemps final: " << sim.get_t_max()
         << "\ndt: " << sim.get_dt() << "\nListe de Verlet: "
         << sim.get_verlet() << '\n'
         << endl;

    cout << setprecision(ss);

    sim.fillMoment();

    cout << "État initial" << endl;
    cout << "Temps: " << sim.get_t_min()
         << " | Énergie cinétique = " << sim.computeKineticEnergy()
         << " | Température = " << sim.computeKineticTemperature() << '\n'
         << endl;

    auto t1 = chrono::high_resolution_clock::now();

    sim.start(true, 15, 10);

    auto t2 = chrono::high_resolution_clock::now();

    double elapsed = chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();

    cout << "\nÉtat final" << endl;
    cout << "Temps: " << sim.get_t_max()
         << " | Énergie cinétique = " << sim.computeKineticEnergy()
         << " | Température = " << sim.computeKineticTemperature() << endl;

    cout << "\nDurée de la simulation: " << elapsed * 1e-9 << "s\n"
         << endl;

    return 0;
}
