#include "part.hpp"

#include <math.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>

const double r_star_sq = 9.0;
const double epsilon_star = 0.2;
const double r_cut = 100.0;

void fill_vec(unique_ptr<array<Particule, TOTAL_PARTS>> &parts, const int n, const std::string path) {
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
}

double energy_forces(const unique_ptr<array<Particule, TOTAL_PARTS>> &parts, unique_ptr<Particule> &forces, const int n) {
    double lj = 0.0;
    double u_ij = 0.0;
    unique_ptr<Particule> neg_forces = make_unique<Particule>();

    for (int i = 0; i < n; ++i) {
        const double ix = (*parts)[i].x;
        const double iy = (*parts)[i].y;
        const double iz = (*parts)[i].z;

        for (int j = i + 1; j < n; ++j) {
            const double jx = (*parts)[j].x;
            const double jy = (*parts)[j].y;
            const double jz = (*parts)[j].z;

            const double rx = (ix - jx);
            const double ry = (iy - jy);
            const double rz = (iz - jz);

            const double r_ij_sq = rx * rx + ry * ry + rz * rz;

            const double r2 = r_star_sq / r_ij_sq;

            const double r4 = r2 * r2;
            const double r6 = r4 * r2;
            const double r8 = r4 * r4;
            const double r12 = r6 * r6;
            const double r14 = r12 * r2;

            u_ij = -48.0 * epsilon_star * (r14 - r8);

            const double fx = u_ij * rx;
            const double fy = u_ij * ry;
            const double fz = u_ij * rz;

            forces->x += fx;
            forces->y += fy;
            forces->z += fz;

            forces->x += -fx;
            forces->y += -fy;
            forces->z += -fz;

            lj += r12 - 2.0 * r6;
        }
    }

    // forces->x += neg_forces->x;
    // forces->y += neg_forces->y;
    // forces->z += neg_forces->z;

    return 4.0 * epsilon_star * lj;
}

double energy_forces_periode(const unique_ptr<array<Particule, TOTAL_PARTS>> &parts, unique_ptr<Particule> &forces, const int n_sym, const double lx, const double ly, const double lz, const int n) {
    double lj = 0.0;
    double u_ij = 0.0;

    array<double, 27 * 3> trans;
    array<double, 3> temp{0.0, 42.0, -42.0};

    size_t index = 0;

    for (double x : temp) {
        for (double y : temp) {
            for (double z : temp) {
                trans[index++] = x;
                trans[index++] = y;
                trans[index++] = z;
            }
        }
    }

    for (int i_sym = 0; i_sym < n_sym * 3; i_sym += 3) {
        const double ix_sym = trans[i_sym];
        const double iy_sym = trans[i_sym + 1];
        const double iz_sym = trans[i_sym + 2];

        for (int i = 0; i < n; ++i) {
            double ix = (*parts)[i].x;
            double iy = (*parts)[i].y;
            double iz = (*parts)[i].z;

            ix = ix - int(ix / lx) * lx;
            iy = iy - int(iy / ly) * ly;
            iz = iz - int(iz / lz) * lz;

            for (int j = 0; j < n; ++j) {
                if (i == j) {
                    continue;
                }

                double jx = (*parts)[j].x;
                double jy = (*parts)[j].y;
                double jz = (*parts)[j].z;

                jx = jx - int(jx / lx) * lx;
                jy = jy - int(jy / ly) * ly;
                jz = jz - int(jz / lz) * lz;

                const double jx_loc = jx + ix_sym;
                const double jy_loc = jy + iy_sym;
                const double jz_loc = jz + iz_sym;

                const double rx = (ix - jx_loc);
                const double ry = (iy - jy_loc);
                const double rz = (iz - jz_loc);

                const double r_ij_sq = rx * rx + ry * ry + rz * rz;

                if (r_ij_sq < r_cut) {
                    const double r2 = r_star_sq / r_ij_sq;

                    const double r4 = r2 * r2;
                    const double r6 = r4 * r2;
                    const double r8 = r4 * r4;
                    const double r12 = r6 * r6;
                    const double r14 = r12 * r2;

                    u_ij = -48.0 * epsilon_star * (r14 - r8);

                    const double fx = u_ij * rx;
                    const double fy = u_ij * ry;
                    const double fz = u_ij * rz;

                    forces->x += fx;
                    forces->y += fy;
                    forces->z += fz;

                    forces->x += -fx;
                    forces->y += -fy;
                    forces->z += -fz;

                    lj += (r12 - 2.0 * r6);
                } else {
                    continue;
                }
            }
        }
    }

    return ((4.0 * epsilon_star) / 2.0) * lj;
}