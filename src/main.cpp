#include <iostream>
#include <vector>

#include "part.hpp"

using namespace std;

#define TOTAL_PARTS 1000
#define LOCAL_PARTS 999

static string DEFAULT_PARTS_FILE = "../particule.xyz";

int main() {
    size_t fsize = ((TOTAL_PARTS * (TOTAL_PARTS - 1)) / 2);

    vector<Particule> *parts = new vector<Particule>;
    vector<double> *energies = new vector<double>;
    vector<double> *forces = new vector<double>;

    double lj = 0.0;

    parts->reserve(TOTAL_PARTS);
    energies->reserve(fsize);
    forces->reserve(3*TOTAL_PARTS * TOTAL_PARTS);

    fill_vec(parts, TOTAL_PARTS, DEFAULT_PARTS_FILE);
    lj =  compute_lj(parts, TOTAL_PARTS);
    compute_forces(forces, parts, TOTAL_PARTS);

    cout << "Lennard-Jones: " << lj << endl;

    double fsum = 0.0;

    for (size_t i = 0; i < 3*TOTAL_PARTS * TOTAL_PARTS; ++i) {
        fsum += (*forces)[i];
    }

    cout << "Somme des forces: " << fsum << endl;
    

    return 0;
}