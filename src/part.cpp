#include "part.hpp"

#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>

const double r_star_sq = 9.0;
const double epsilon_star = 0.2;
const double r_cut = 100.0;

void fill_vec_aos(unique_ptr<array<data_t, TOTAL_PARTS>> &__restrict__ parts, const int n, const std::string path) {
    string line;
    ifstream fs(path);

    if (fs.is_open()) {
        for (int i = 0; i < n; ++i) {
            getline(fs, line);
            istringstream ss(line);
            ss >> (*parts)[i].x >> (*parts)[i].y >> (*parts)[i].z;
        }
    } else {
        cerr << "File could not be opened!\n";
        cerr << "Error code: " << strerror(errno);
    }

    fs.close();
}

double energy_forces_aos(const unique_ptr<array<data_t, TOTAL_PARTS>> &__restrict__ parts, unique_ptr<array<data_t, TOTAL_PARTS>> &__restrict__ forces, const int n) {
    double lj = 0.0;
    double u_ij = 0.0;

    double sum = 0.0;

    for (int i = 0; i < n; ++i) {
        const double ix = (*parts)[i].x;
        const double iy = (*parts)[i].y;
        const double iz = (*parts)[i].z;

        for (int j = i + 1; j < n; ++j) {
            const double jx = (*parts)[j].x;
            const double jy = (*parts)[j].y;
            const double jz = (*parts)[j].z;

            const double rx = ix - jx;
            const double ry = iy - jy;
            const double rz = iz - jz;

            const double r_ij_sq = rx * rx + ry * ry + rz * rz;

            const double r2 = r_star_sq * (1.0 / r_ij_sq);

            const double r4 = r2 * r2;

            const double r6 = r4 * r2;
            const double r12 = r6 * r4 * r2;

            const double r8 = r6 * r2;
            const double r14 = r8 * r4 * r2;

            u_ij = (r14 - r8) * -48.0 * epsilon_star;

            const double fx = u_ij * rx;
            const double fy = u_ij * ry;
            const double fz = u_ij * rz;

            (*forces)[i].x += fx;
            (*forces)[i].y += fy;
            (*forces)[i].z += fz;

            (*forces)[j].x += -fx;
            (*forces)[j].y += -fy;
            (*forces)[j].z += -fz;

            lj += r12 - 2.0 * r6;
        }

        sum += (*forces)[i].x + (*forces)[i].y + (*forces)[i].z;
    }

    cout << scientific << "Somme des forces: " << sum << endl;

    return 4.0 * epsilon_star * lj;
}

double energy_forces_periode_aos(const unique_ptr<array<data_t, TOTAL_PARTS>> &__restrict__ parts, unique_ptr<array<data_t, TOTAL_PARTS>> &__restrict__ forces, const int n_sym, const double l, const int n) {
    double lj = 0.0;
    double u_ij = 0.0;

    double sum = 0.0;

    for (int i_sym = 0; i_sym < n_sym; ++i_sym) {
        const double ix_sym = vec_trans[i_sym][0];
        const double iy_sym = vec_trans[i_sym][1];
        const double iz_sym = vec_trans[i_sym][2];

        for (int i = 0; i < n; ++i) {
            double ix = (*parts)[i].x;
            double iy = (*parts)[i].y;
            double iz = (*parts)[i].z;

            ix = ix - int(ix * (1.0 / l)) * l;
            iy = iy - int(iy * (1.0 / l)) * l;
            iz = iz - int(iz * (1.0 / l)) * l;

            for (int j = i + 1; j < n; ++j) {
                if (i == j) {
                    continue;
                }

                double jx = (*parts)[j].x;
                double jy = (*parts)[j].y;
                double jz = (*parts)[j].z;

                jx = jx - int(jx * (1.0 / l)) * l;
                jy = jy - int(jy * (1.0 / l)) * l;
                jz = jz - int(jz * (1.0 / l)) * l;

                const double jx_loc = jx + ix_sym;
                const double jy_loc = jy + iy_sym;
                const double jz_loc = jz + iz_sym;

                const double rx = ix - jx_loc;
                const double ry = iy - jy_loc;
                const double rz = iz - jz_loc;

                const double r_ij_sq = rx * rx + ry * ry + rz * rz;

                if (r_ij_sq > r_cut) {
                    continue;
                }

                const double r2 = r_star_sq * (1.0 / r_ij_sq);

                const double r4 = r2 * r2;

                const double r6 = r4 * r2;
                const double r12 = r6 * r4 * r2;

                const double r8 = r6 * r2;
                const double r14 = r4 * r4 * r2 * r2 * r2;

                u_ij = -48.0 * epsilon_star * (r14 - r8);

                const double fx = u_ij * rx;
                const double fy = u_ij * ry;
                const double fz = u_ij * rz;

                (*forces)[i].x += fx;
                (*forces)[i].y += fy;
                (*forces)[i].z += fz;

                (*forces)[j].x += -fx;
                (*forces)[j].y += -fy;
                (*forces)[j].z += -fz;

                lj += r12 - 2.0 * r6;
            }

            sum += (*forces)[i].x + (*forces)[i].y + (*forces)[i].z;
        }
    }

    cout << scientific << "Somme des forces sym 27: " << sum << endl;

    return 2.0 * ((4.0 * epsilon_star) * (1.0 / 2.0)) * lj;
}

void fill_vec_soa(unique_ptr<dataArray_t> &__restrict__ parts, const int n, const string path) {
    string line;
    ifstream fs(path);

    if (fs.is_open()) {
        for (int i = 0; i < n; ++i) {
            getline(fs, line);
            istringstream ss(line);
            ss >> parts->x[i] >> parts->y[i] >> parts->z[i];
        }
    } else {
        cerr << "File could not be opened!\n";
        cerr << "Error code: " << strerror(errno);
    }

    fs.close();
}

double energy_forces_soa(const unique_ptr<dataArray_t> &__restrict__ parts, unique_ptr<dataArray_t> &__restrict__ forces, const int n) {
    double lj = 0.0;
    double u_ij = 0.0;

    double sum = 0.0;

    for (int i = 0; i < n; ++i) {
        const double ix = parts->x[i];
        const double iy = parts->y[i];
        const double iz = parts->z[i];

        for (int j = i + 1; j < n; ++j) {
            const double jx = parts->x[j];
            const double jy = parts->y[j];
            const double jz = parts->z[j];

            const double rx = ix - jx;
            const double ry = iy - jy;
            const double rz = iz - jz;

            const double r_ij_sq = rx * rx + ry * ry + rz * rz;

            const double r2 = r_star_sq * (1.0 / r_ij_sq);

            const double r4 = r2 * r2;

            const double r6 = r4 * r2;
            const double r12 = r6 * r4 * r2;

            const double r8 = r6 * r2;
            const double r14 = r8 * r4 * r2;

            u_ij = (r14 - r8) * -48.0 * epsilon_star;

            const double fx = u_ij * rx;
            const double fy = u_ij * ry;
            const double fz = u_ij * rz;

            forces->x[i] += fx;
            forces->y[i] += fy;
            forces->z[i] += fz;

            forces->x[j] += -fx;
            forces->y[j] += -fy;
            forces->z[j] += -fz;

            lj += r12 - 2.0 * r6;
        }

        sum += forces->x[i] + forces->y[i] + forces->z[i];
    }

    cout << scientific << "Somme des forces: " << sum << endl;

    return 4.0 * epsilon_star * lj;
}

double energy_forces_periode_soa(const unique_ptr<dataArray_t> &__restrict__ parts, unique_ptr<dataArray_t> &__restrict__ forces, const int n_sym, const double l, const int n) {
    double lj = 0.0;
    double u_ij = 0.0;

    double sum = 0.0;

    for (int i_sym = 0; i_sym < n_sym; ++i_sym) {
        const double ix_sym = vec_trans[i_sym][0];
        const double iy_sym = vec_trans[i_sym][1];
        const double iz_sym = vec_trans[i_sym][2];

        for (int i = 0; i < n; ++i) {
            double ix = parts->x[i];
            double iy = parts->y[i];
            double iz = parts->z[i];

            ix = ix - int(ix * (1.0 / l)) * l;
            iy = iy - int(iy * (1.0 / l)) * l;
            iz = iz - int(iz * (1.0 / l)) * l;

            for (int j = i + 1; j < n; ++j) {
                if (i == j) {
                    continue;
                }

                double jx = parts->x[j];
                double jy = parts->y[j];
                double jz = parts->z[j];

                jx = jx - int(jx * (1.0 / l)) * l;
                jy = jy - int(jy * (1.0 / l)) * l;
                jz = jz - int(jz * (1.0 / l)) * l;

                const double jx_loc = jx + ix_sym;
                const double jy_loc = jy + iy_sym;
                const double jz_loc = jz + iz_sym;

                const double rx = ix - jx_loc;
                const double ry = iy - jy_loc;
                const double rz = iz - jz_loc;

                const double r_ij_sq = rx * rx + ry * ry + rz * rz;

                if (r_ij_sq > r_cut) {
                    continue;
                }

                const double r2 = r_star_sq * (1.0 / r_ij_sq);

                const double r4 = r2 * r2;

                const double r6 = r4 * r2;
                const double r12 = r6 * r4 * r2;

                const double r8 = r6 * r2;
                const double r14 = r4 * r4 * r2 * r2 * r2;

                u_ij = -48.0 * epsilon_star * (r14 - r8);

                const double fx = u_ij * rx;
                const double fy = u_ij * ry;
                const double fz = u_ij * rz;

                forces->x[i] += fx;
                forces->y[i] += fy;
                forces->z[i] += fz;

                forces->x[j] += -fx;
                forces->y[j] += -fy;
                forces->z[j] += -fz;

                lj += r12 - 2.0 * r6;
            }

            sum += forces->x[i] + forces->y[i] + forces->z[i];
        }
    }

    cout << scientific << "Somme des forces sym 27: " << sum << endl;

    return 2.0 * ((4.0 * epsilon_star) * (1.0 / 2.0)) * lj;
}