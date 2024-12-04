#include <iostream>
#include <vector>

#include "part.hpp"
#include "utils.hpp"

using namespace std;

#define NB_PART 1000

static string DEFAULT_PARTS_FILE = "../particule.xyz";

int main() {
    vector<Particule> *parts = new vector<Particule>(NB_PART);

    fill_vec(parts, NB_PART, DEFAULT_PARTS_FILE);

    cout << (*parts)[999].x << " " << (*parts)[999].y << " " << (*parts)[999].z << endl;

    return 0;
}