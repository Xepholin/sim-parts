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

int main(int argc, char *argv[]) {
    if (argc != 10 && argc != 11) {
        cerr << "Usage: " << argv[0] << " n_particles n_sym t_min t_max dt T0 gamma use_verlet input_file [output_file]\n";
        return 1;
    }

    const int n_particles = atoi(argv[1]);
    const int n_sym = atoi(argv[2]);
    const double t_min = atof(argv[3]);
    const double t_max = atof(argv[4]);
    const double dt = atof(argv[5]);
    const double T0 = atof(argv[6]);
    const double gamma = atof(argv[7]);
    const bool use_verlet = atoi(argv[8]) != 0;
    const string input = argv[9];
    string output;

    if (argc == 11) {
        output = argv[10];
    } else {
        output = DEFAULT_OUTPUT;
    }

    Simulator sim(n_particles, n_sym, t_min, t_max, dt, T0, gamma, use_verlet);
    streamsize ss = std::cout.precision();

    cout << "\n=================================\n"
         << "      Début de la simulation\n"
         << "=================================\n"
         << endl;

    cout << "Lecture dans " << input << " ...\n" << endl;

    sim.fillVec(input);
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

    cout << "Écriture dans " << output << " ...\n" << endl;

    auto t1 = chrono::high_resolution_clock::now();

    sim.start(true, 15, 10, output);

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
