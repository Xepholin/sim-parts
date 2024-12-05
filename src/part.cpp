#include "part.hpp"

#include <math.h>

#include <fstream>
#include <iostream>
#include <set>
#include <sstream>

const double r_star_sq = 9.0;
const double epsilon_star = 0.2;

void *fill_vec(array<Particule, TOTAL_PARTS> *parts, int n, const string path) {
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

double compute_lj(const array<Particule, TOTAL_PARTS> *parts, const int n) {
    double lj = 0.0;

    double ix = 0.0;
    double iy = 0.0;
    double iz = 0.0;

    double jx = 0.0;
    double jy = 0.0;
    double jz = 0.0;

    double r_ij_sq = 0.0;

    double rx = 0.0;
    double ry = 0.0;
    double rz = 0.0;

    double r2 = 0.0;
    double r4 = 0.0;
    double r6 = 0.0;
    double r12 = 0.0;

    for (int i = 0; i < n; ++i) {
        ix = (*parts)[i].x;
        iy = (*parts)[i].y;
        iz = (*parts)[i].z;

        for (int j = i + 1; j < n; ++j) {
            jx = (*parts)[j].x;
            jy = (*parts)[j].y;
            jz = (*parts)[j].z;

            rx = (ix - jx);
            ry = (iy - jy);
            rz = (iz - jz);

            r_ij_sq = rx * rx + ry * ry + rz * rz;

            r2 = r_star_sq / r_ij_sq;

            r4 = r2 * r2;
            r6 = r4 * r2;
            r12 = r6 * r6;

            lj += (r12 - 2.0 * r6);
        }
    }

    return 4.0 * epsilon_star * lj;
}

void compute_forces(array<double, fsize> *forces, const array<Particule, TOTAL_PARTS> *parts, const int n) {
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

    double r_ij_sq = 0.0;

    double rx = 0.0;
    double ry = 0.0;
    double rz = 0.0;

    double r2 = 0.0;

    double r4 = 0.0;
    double r8 = 0.0;
    double r14 = 0.0;

    for (int i = 0; i < n; ++i) {
        ix = (*parts)[i].x;
        iy = (*parts)[i].y;
        iz = (*parts)[i].z;

        for (int j = 0; j < n; ++j) {
            if (i == j) {
                (*forces)[i * n + j] = 0.0;
                (*forces)[i * n + j + 1] = 0.0;
                (*forces)[i * n + j + 2] = 0.0;

                continue;
            }

            jx = (*parts)[j].x;
            jy = (*parts)[j].y;
            jz = (*parts)[j].z;

            rx = (ix - jx);
            ry = (iy - jy);
            rz = (iz - jz);

            r_ij_sq = rx * rx + ry * ry + rz * rz;

            r2 = r_star_sq / r_ij_sq;

            r4 = r2 * r2;
            r8 = r4 * r4;
            r14 = r8 * r4 * r2;

            u_ij = -48.0 * epsilon_star * (r14 - r8);

            fx = u_ij * (ix - jx);
            fy = u_ij * (iy - jy);
            fz = u_ij * (iz - jz);

            (*forces)[i * n + j] = fx;
            (*forces)[i * n + j + 1] = fy;
            (*forces)[i * n + j + 2] = fz;
        }
    }
}