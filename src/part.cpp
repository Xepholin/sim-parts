#include "part.hpp"

#include <math.h>

#include <fstream>
#include <iostream>
#include <set>
#include <sstream>

const double r_star_sq = 9.0;
const double epsilon_star = 0.2;
const double r_cut = 100.0;

void fill_vec(unique_ptr<Particule> &parts, const int n, const std::string path) {
    string line;
    ifstream fs(path);

    if (fs.is_open()) {
        for (int i = 0; i < n; ++i) {
            getline(fs, line);
            istringstream ss(line);
            ss >> parts->x[i] >> parts->y[i] >> parts->z[i];
        }
    } else {
        cout << "Unable to open file";
    }
}

double energy_forces(const unique_ptr<Particule> &parts, unique_ptr<array<double, 3>> &forces, const int n) {
    double lj = 0.0;
    double u_ij = 0.0;

    for (int i = 0; i < n; ++i) {
        const double ix = parts->x[i];
        const double iy = parts->y[i];
        const double iz = parts->z[i];

        for (int j = i + 1; j < n; ++j) {
            const double  jx = parts->x[j];
            const double  jy = parts->y[j];
            const double  jz = parts->z[j];

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

            (*forces)[0] += fx;
            (*forces)[1] += fy;
            (*forces)[2] += fz;

            (*forces)[0] += -fx;
            (*forces)[1] += -fy;
            (*forces)[2] += -fz;

            lj += (r12 - 2.0 * r6);
        }
    }

    return 4.0 * epsilon_star * lj;
}

double energy_forces_periode(const unique_ptr<Particule> &parts, unique_ptr<array<double, 3>> &forces, const int n_sym, const double lx, const double ly, const double lz, const int n) {
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
            double ix = parts->x[i];
            double iy = parts->y[i];
            double iz = parts->z[i];

            ix = ix - int(ix / lx) * lx;
            iy = iy - int(iy / ly) * ly;
            iz = iz - int(iz / lz) * lz;

            for (int j = 0; j < n; ++j) {
                if (i == j) {
                    continue;
                }

                double jx = parts->x[j];
                double jy = parts->y[j];
                double jz = parts->z[j];

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

                    (*forces)[0] += fx;
                    (*forces)[1] += fy;
                    (*forces)[2] += fz;

                    (*forces)[0] += -fx;
                    (*forces)[1] += -fy;
                    (*forces)[2] += -fz;

                    lj += (r12 - 2.0 * r6);
                } else {
                    continue;
                }
            }
        }
    }

    return ((4.0 * epsilon_star) / 2.0) * lj;
}