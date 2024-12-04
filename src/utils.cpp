#include "utils.hpp"

#include <fstream>
#include <iostream>
#include <sstream>

void fill_vec(vector<Particule> *parts, int n, std::string const path) {
    string line;
    ifstream fs(path);

    if (fs.is_open()) {
        for (int i = 0; i < n; ++i)   {
            getline(fs, line);
            istringstream ss(line);
            ss >> (*parts)[i].x >> (*parts)[i].y >> (*parts)[i].z;
        }
    } else {
        cout << "Unable to open file";
    }
}