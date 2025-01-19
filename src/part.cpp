#include "part.hpp"

#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>

template <typename T>
void fill_vec(unique_ptr<T> &__restrict__ parts, const int n, const string &path) {
    string line;
    ifstream fs(path);

    if (fs.is_open()) {
        for (int i = 0; i < n; ++i) {
            getline(fs, line);
            istringstream ss(line);
            if constexpr (is_same_v<T, array<data_t, TOTAL_PARTS>>) {
                ss >> (*parts)[i].x >> (*parts)[i].y >> (*parts)[i].z;
            } else if constexpr (is_same_v<T, dataArray_t>) {
                ss >> parts->x[i] >> parts->y[i] >> parts->z[i];
            }
        }
    } else {
        cerr << "File could not be opened!\n";
        cerr << "Error code: " << strerror(errno) << '\n';
    }

    fs.close();
}

template void fill_vec<array<data_t, TOTAL_PARTS>>(unique_ptr<array<data_t, TOTAL_PARTS>> &__restrict__, const int n, const string &path);
template void fill_vec<dataArray_t>(unique_ptr<dataArray_t> &__restrict__, const int n, const string &path);

template <typename T>
double energy_forces(const unique_ptr<T> &__restrict__ parts, unique_ptr<T> &__restrict__ forces, const int n_sym, const int n) {
    double lj = 0.0;
    double u_ij = 0.0;
    double sum = 0.0;

    for (int i_sym = 0; i_sym < n_sym; ++i_sym) {
        const double ix_sym = vec_translation[i_sym][0];
        const double iy_sym = vec_translation[i_sym][1];
        const double iz_sym = vec_translation[i_sym][2];

        for (int i = 0; i < n; ++i) {
            double ix, iy, iz;

            if constexpr (std::is_same_v<T, std::array<data_t, TOTAL_PARTS>>) {
                ix = (*parts)[i].x;
                iy = (*parts)[i].y;
                iz = (*parts)[i].z;
            } else if constexpr (std::is_same_v<T, dataArray_t>) {
                ix = parts->x[i];
                iy = parts->y[i];
                iz = parts->z[i];
            }

            if (l != 0.0) {
                ix = ix - int(ix * (1.0 / l)) * l;
                iy = iy - int(iy * (1.0 / l)) * l;
                iz = iz - int(iz * (1.0 / l)) * l;
            }

            for (int j = i + 1; j < n; ++j) {
                if (i == j) {
                    continue;
                }

                double jx, jy, jz;

                if constexpr (std::is_same_v<T, std::array<data_t, TOTAL_PARTS>>) {
                    jx = (*parts)[j].x;
                    jy = (*parts)[j].y;
                    jz = (*parts)[j].z;
                } else if constexpr (std::is_same_v<T, dataArray_t>) {
                    jx = parts->x[j];
                    jy = parts->y[j];
                    jz = parts->z[j];
                }

                if (l != 0.0) {
                    jx = jx - int(jx * (1.0 / l)) * l;
                    jy = jy - int(jy * (1.0 / l)) * l;
                    jz = jz - int(jz * (1.0 / l)) * l;
                }

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

                if constexpr (std::is_same_v<T, std::array<data_t, TOTAL_PARTS>>) {
                    (*forces)[i].x += fx;
                    (*forces)[i].y += fy;
                    (*forces)[i].z += fz;
                    (*forces)[j].x += -fx;
                    (*forces)[j].y += -fy;
                    (*forces)[j].z += -fz;
                } else if constexpr (std::is_same_v<T, dataArray_t>) {
                    forces->x[i] += fx;
                    forces->y[i] += fy;
                    forces->z[i] += fz;
                    forces->x[j] += -fx;
                    forces->y[j] += -fy;
                    forces->z[j] += -fz;
                }

                lj += r12 - 2.0 * r6;
            }

            if constexpr (std::is_same_v<T, std::array<data_t, TOTAL_PARTS>>) {
                sum += (*forces)[i].x + (*forces)[i].y + (*forces)[i].z;
            } else if constexpr (std::is_same_v<T, dataArray_t>) {
                sum += forces->x[i] + forces->y[i] + forces->z[i];
            }
        }
    }

    cout << scientific << "Somme des forces sym " << n_sym << ": " << sum << endl;

    return 2.0 * ((4.0 * epsilon_star) * (1.0 / 2.0)) * lj;
}

template double energy_forces<array<data_t, TOTAL_PARTS>>(const unique_ptr<array<data_t, TOTAL_PARTS>> &__restrict__, unique_ptr<array<data_t, TOTAL_PARTS>> &__restrict__, const int n_sym, const int n);
template double energy_forces<dataArray_t>(const unique_ptr<dataArray_t> &__restrict__, unique_ptr<dataArray_t> &__restrict__, const int n_sym, const int n);
