#include "part.hpp"

#include <unistd.h>

#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <set>
#include <sstream>

#include "const.hpp"
#include "tools.hpp"

Simulator::Simulator(const int n, const int n_sym, const double t_min, const double t_max, const double dt, const double T0, const double gamma, const bool verlet) {
    this->n = n;
    this->n_sym = n_sym;
    this->t_min = t_min;
    this->t_max = t_max;
    this->dt = dt;
    this->T0 = T0;
    this->gamma = gamma;

    this->verlet = verlet;

    parts = make_unique<array<data_t, TOTAL_PARTS>>();
    forces = make_unique<array<data_t, TOTAL_PARTS>>();
    moments = make_unique<array<data_t, TOTAL_PARTS>>();
    velocities = make_unique<array<data_t, TOTAL_PARTS>>();
    masses = make_unique<array<double, TOTAL_PARTS>>();
    neighbor = make_unique<array<array<int, n_max_neighbor * 2>, TOTAL_PARTS>>();

    fill(parts->begin(), parts->end(), data_t{0.0, 0.0, 0.0});
    fill(forces->begin(), forces->end(), data_t{0.0, 0.0, 0.0});
    fill(moments->begin(), moments->end(), data_t{0.0, 0.0, 0.0});
    fill(masses->begin(), masses->end(), 18.0);

    for (auto& row : *neighbor) {
        std::fill(row.begin(), row.end(), 0);
    }
}

Simulator::~Simulator() = default;

double Simulator::get_n() {
    return n;
}

double Simulator::get_n_sym() {
    return n_sym;
}

double Simulator::get_t_min() {
    return t_min;
}

double Simulator::get_t_max() {
    return t_max;
}

double Simulator::get_dt() {
    return dt;
}

string Simulator::get_verlet() {
    if (verlet) {
        return "oui";
    } else {
        return "non";
    }
}

void Simulator::fillVec(const string& path) {
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
        cerr << "Error code: " << strerror(errno) << '\n';
    }

    fs.close();
}

double Simulator::computeEnergyForces() {
    std::fill(forces->begin(), forces->end(), data_t{0.0, 0.0, 0.0});
    double lj = 0.0;
    double f_ij = 0.0;

    for (int i_sym = 0; i_sym < n_sym; ++i_sym) {
        const double ix_sym = vec_translation[i_sym][0];
        const double iy_sym = vec_translation[i_sym][1];
        const double iz_sym = vec_translation[i_sym][2];

        for (int i = 0; i < n; ++i) {
            double ix = (*parts)[i].x - ix_sym;
            double iy = (*parts)[i].y - iy_sym;
            double iz = (*parts)[i].z - iz_sym;

            for (int j = i + 1; j < n; ++j) {
                double jx = (*parts)[j].x;
                double jy = (*parts)[j].y;
                double jz = (*parts)[j].z;

                double dx = ix - jx;
                double dy = iy - jy;
                double dz = iz - jz;

                const double r_ij_sq = dx * dx + dy * dy + dz * dz;

                if (r_ij_sq > r_cut_sq || r_ij_sq < 1e-6) {
                    continue;
                }

                const double r2 = r_star_sq * (1.0 / r_ij_sq);
                const double r4 = r2 * r2;
                const double r6 = r4 * r2;
                const double r12 = r6 * r6;
                const double r8 = r6 * r2;
                const double r14 = r8 * r6;

                f_ij = -48.0 * epsilon_star * (r14 - 0.5 * r8);

                const double fx = f_ij * dx;
                const double fy = f_ij * dy;
                const double fz = f_ij * dz;

                (*forces)[i].x += fx;
                (*forces)[i].y += fy;
                (*forces)[i].z += fz;

                (*forces)[j].x -= fx;
                (*forces)[j].y -= fy;
                (*forces)[j].z -= fz;

                lj += r12 - r6;
            }
        }
    }

    return 4.0 * epsilon_star * lj;
}

void Simulator::verletList() {
    for (int i = 0; i < n; ++i) {
        (*neighbor)[i][0] = 0;
    }

    for (int i = 0; i < n; ++i) {
        int ni = 0;

        for (int j = i + 1; j < n; ++j) {
            for (int i_sym = 0; i_sym < n_sym; ++i_sym) {
                double xj_loc = (*parts)[j].x + vec_translation[i_sym][0];
                double yj_loc = (*parts)[j].y + vec_translation[i_sym][1];
                double zj_loc = (*parts)[j].z + vec_translation[i_sym][2];

                double dx = (*parts)[i].x - xj_loc;
                double dy = (*parts)[i].y - yj_loc;
                double dz = (*parts)[i].z - zj_loc;
                double r_ij_sq = dx * dx + dy * dy + dz * dz;

                if (r_ij_sq < r_cut_sq) {
                    ni += 1;
                    int mi = 2 * ni - 1;
                    (*neighbor)[i][mi] = j;
                    (*neighbor)[i][mi + 1] = i_sym;
                }
            }
        }
        (*neighbor)[i][0] = ni;
    }
}

double Simulator::computeEnergyForcesVerlet() {
    std::fill(forces->begin(), forces->end(), data_t{0.0, 0.0, 0.0});
    double lj = 0.0;
    double f_ij = 0.0;

    for (int i = 0; i < n; ++i) {
        for (int ni = 1; ni <= (*neighbor)[i][0]; ++ni) {
            int j = (*neighbor)[i][2 * ni - 1];
            int i_sym = (*neighbor)[i][2 * ni];

            double xj_loc = (*parts)[j].x + vec_translation[i_sym][0];
            double yj_loc = (*parts)[j].y + vec_translation[i_sym][1];
            double zj_loc = (*parts)[j].z + vec_translation[i_sym][2];

            double dx = (*parts)[i].x - xj_loc;
            double dy = (*parts)[i].y - yj_loc;
            double dz = (*parts)[i].z - zj_loc;

            double r_ij_sq = dx * dx + dy * dy + dz * dz;

            if (r_ij_sq > r_cut_sq || r_ij_sq < 1e-6) {
                continue;
            }

            const double r2 = r_star_sq * (1.0 / r_ij_sq);
            const double r4 = r2 * r2;
            const double r6 = r4 * r2;
            const double r12 = r6 * r6;
            const double r8 = r6 * r2;
            const double r14 = r8 * r6;

            f_ij = -48.0 * epsilon_star * (r14 - 0.5 * r8);

            const double fx = f_ij * dx;
            const double fy = f_ij * dy;
            const double fz = f_ij * dz;

            (*forces)[i].x += fx;
            (*forces)[i].y += fy;
            (*forces)[i].z += fz;

            (*forces)[j].x -= fx;
            (*forces)[j].y -= fy;
            (*forces)[j].z -= fz;

            lj += r12 - r6;
        }
    }

    return 4.0 * epsilon_star * lj;
}

double Simulator::sumForces() {
    double sum_forces = 0.0;

    for (int i = 0; i < TOTAL_PARTS; ++i) {
        sum_forces += (*forces)[i].x + (*forces)[i].y + (*forces)[i].z;
    }

    return sum_forces;
}

void Simulator::velocityVerlet() {
    double dt_half = 0.5 * dt;
    double l_half = 0.5 * l;

    for (int i = 0; i < n; ++i) {
        (*moments)[i].x -= dt_half * conv_force * (*forces)[i].x;
        (*moments)[i].y -= dt_half * conv_force * (*forces)[i].y;
        (*moments)[i].z -= dt_half * conv_force * (*forces)[i].z;
    }

    for (int i = 0; i < n; ++i) {
        double inv_mass = 1.0 / (*masses)[i];
        (*parts)[i].x += dt * (*moments)[i].x * inv_mass;
        (*parts)[i].y += dt * (*moments)[i].y * inv_mass;
        (*parts)[i].z += dt * (*moments)[i].z * inv_mass;
    }

    if (l != 0.0) {
        for (int i = 0; i < n; ++i) {
            (*parts)[i].x -= l * int((*parts)[i].x / l_half);
            (*parts)[i].y -= l * int((*parts)[i].y / l_half);
            (*parts)[i].z -= l * int((*parts)[i].z / l_half);
        }
    }

    if (verlet) {
        computeEnergyForcesVerlet();
    } else {
        computeEnergyForces();
    }

    for (int i = 0; i < n; ++i) {
        (*moments)[i].x -= dt_half * conv_force * (*forces)[i].x;
        (*moments)[i].y -= dt_half * conv_force * (*forces)[i].y;
        (*moments)[i].z -= dt_half * conv_force * (*forces)[i].z;
    }
}

double Simulator::computeKineticEnergy() {
    double sum = 0.0;

    for (int i = 0; i < n; ++i) {
        double px = (*moments)[i].x;
        double py = (*moments)[i].y;
        double pz = (*moments)[i].z;
        double inv_mass = 1.0 / (*masses)[i];

        sum += (px * px + py * py + pz * pz) * inv_mass;
    }

    return sum / (2.0 * conv_force);
}

double Simulator::computeKineticTemperature() {
    const double energy = computeKineticEnergy();

    return energy / (n_dl * const_r);
}

void Simulator::correctionRatio() {
    const double energy_kinetic_init = computeKineticEnergy();
    const double ratio = sqrt((n_dl * const_r * T0) / energy_kinetic_init);

    for (int i = 0; i < n; ++i) {
        auto& p = (*moments)[i];
        p.x *= ratio;
        p.y *= ratio;
        p.z *= ratio;
    }
}

void Simulator::correctionCenter() {
    double Px = 0.0;
    double Py = 0.0;
    double Pz = 0.0;

    for (int i = 0; i < n; ++i) {
        Px += (*moments)[i].x;
        Py += (*moments)[i].y;
        Pz += (*moments)[i].z;
    }

    Px /= n;
    Py /= n;
    Pz /= n;

    for (int i = 0; i < n; ++i) {
        (*moments)[i].x -= Px;
        (*moments)[i].y -= Py;
        (*moments)[i].z -= Pz;
    }
}

void Simulator::fillMoment() {
    // srand(random_device()());
    srand(42);

    for (int i = 0; i < n; ++i) {
        double c = 0.0;
        double s = 0.0;

        c = (double)rand() / (double)RAND_MAX;
        s = (double)rand() / (double)RAND_MAX;
        double px = fn_sign(c, 0.5 - s);

        c = (double)rand() / (double)RAND_MAX;
        s = (double)rand() / (double)RAND_MAX;
        double py = fn_sign(c, 0.5 - s);

        c = (double)rand() / (double)RAND_MAX;
        s = (double)rand() / (double)RAND_MAX;
        double pz = fn_sign(c, 0.5 - s);

        (*moments)[i].x = px;
        (*moments)[i].y = py;
        (*moments)[i].z = pz;
    }

    correctionRatio();
    correctionCenter();
    correctionRatio();
}

void Simulator::correctionMoment() {
    const double T = computeKineticTemperature();
    const double lambda = gamma * ((T0 / T) - 1.0);

    for (int i = 0; i < n; ++i) {
        (*moments)[i].x += lambda * (*moments)[i].x;
        (*moments)[i].y += lambda * (*moments)[i].y;
        (*moments)[i].z += lambda * (*moments)[i].z;
    }
}

void Simulator::start(const bool out, const int correction_step, const int save_step) {
    ofstream outfile(OUTPUT_PARTS_FILE, ios::trunc);

    if (!outfile.is_open()) {
        cerr << "Failed to open file " << OUTPUT_PARTS_FILE << " for writing!" << endl;
        return;
    }

    double ang = 0.0;
    double temp = 0.0;
    double t = get_t_min();
    int count = 1;

    for (; t <= t_max + 1e-10; t += dt) {
        velocityVerlet();

        ang = computeKineticEnergy();
        temp = computeKineticTemperature();

        if (count % correction_step == 0) {
            correctionMoment();
        }

        if (verlet && count % (int)l == 0) {
            verletList();
        }

        if (count % 100 == 0) {
            cout << "Itération: " << count << " | Énergie: " << ang << " | Température: " << temp << '\n';
        }

        if (out && count % save_step == 0) {
            outfile << n << endl;
            outfile << "Time: " << t << endl;

            for (int i = 0; i < n; ++i) {
                outfile << "C"
                        << " " << fixed << setprecision(6)
                        << parts->at(i).x << " "
                        << parts->at(i).y << " "
                        << parts->at(i).z << endl;
            }
        }

        count++;
    }

    outfile.close();
}
