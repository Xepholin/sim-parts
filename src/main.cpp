#include <iostream>
#include <array>

#include "part.hpp"

using namespace std;

static string DEFAULT_PARTS_FILE = "../particule.xyz";

int main() {
    array<Particule, TOTAL_PARTS> *parts = new array<Particule, TOTAL_PARTS>;
    array<double, fsize> *forces = new array<double, fsize>;

    double lj = 0.0;

    fill_vec(parts, TOTAL_PARTS, DEFAULT_PARTS_FILE);
    lj =  compute_lj(parts, TOTAL_PARTS);
    compute_forces(forces, parts, TOTAL_PARTS);

    cout << "Lennard-Jones: " << lj << endl;

    double fsum = 0.0;

    for (size_t i = 0; i < fsize; ++i) {
        fsum += (*forces)[i];
    }

    cout << "Somme des forces: " << fsum << endl;
    

    return 0;
}