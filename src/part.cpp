#include "part.hpp"

#include <fstream>
#include <iostream>
#include <set>
#include <sstream>

const double r_star = 3.0;
const double epsilon_star = 0.2;

void *fill_vec(vector<Particule> *parts, int n, std::string const path) {
    string line;
    ifstream fs(path);

    if (fs.is_open()) {
        for (int i = 0; i < n; ++i) {
            getline(fs, line);
            istringstream ss(line);
            ss >> (*parts)[i].x >> (*parts)[i].y >> (*parts)[i].z;
        }
    } else {
        cout << "Unable to open file";
    }

    return parts;
}

void compute_force(vector<double> *forces, const vector<Particule> *parts, int n)   {
    double u_ij = 0.0;

    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;

    double ix = 0.0;
    double iy = 0.0;
    double iz = 0.0;

    double jx = 0.0;
    double jy = 0.0;
    double jz = 0.0;

    double r_ij = 0.0;

    double r = 0.0;

    double r2 = 0.0;
    double r4 = 0.0;
    double r8 = 0.0;
    double r14 = 0.0;
    
    for (int i = 0; i < n; ++i) {
        ix = (*parts)[i].x;
        iy = (*parts)[i].y;
        iz = (*parts)[i].z;

        r_ij = 0;
        r_ij = 0;
        r_ij = 0;

        for (int j = i + 1; j < n; ++j) {
            jx = (*parts)[j].x;
            jy = (*parts)[j].y;
            jz = (*parts)[j].z;

            r_ij += (ix - jx) * (ix - jx);
            r_ij += (iy - jy) * (iy - jy);
            r_ij += (iz - jz) * (iz - jz);

            r = r_star / r_ij;

            r2 = r * r;
            r4 = r2 * r2;
            r8 = r4 * r4;
            r14 = r8 * r4 * r2;

            u_ij = -48.0 * epsilon_star * (r14 - r8);

            fx = u_ij * (ix - jx);
            fy = u_ij * (iy - jy);
            fz = u_ij * (iz - jz);

            forces->push_back(fx);
            forces->push_back(fy);
            forces->push_back(fz);
        }
    }
}