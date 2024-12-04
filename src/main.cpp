#include <iostream>
#include <vector>

#include "part.hpp"

using namespace std;

#define TOTAL_PARTS 1000
#define LOCAL_PARTS 999

static string DEFAULT_PARTS_FILE = "../particule.xyz";

int main() {
    size_t fsize = 3 * ((TOTAL_PARTS * (TOTAL_PARTS - 1)) / 2);

    vector<Particule> *parts = new vector<Particule>;
    vector<double> *forces = new vector<double>;

    parts->reserve(TOTAL_PARTS);
    forces->reserve(fsize);

    fill_vec(parts, TOTAL_PARTS, DEFAULT_PARTS_FILE);
    compute_force(forces, parts, TOTAL_PARTS);

    double sum = 0.0;

    for (size_t i = 0; i < fsize; ++i) {
        sum += (*forces)[i];
    }

    cout << "Somme des forces: " << sum << endl;

    return 0;
}