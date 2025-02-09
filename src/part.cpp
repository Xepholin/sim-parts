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

Simulator::Simulator(const int n, const int n_sym, const double t_min, const double t_max, const double dt, const double T0, const double gamma) {
    this->n = n;
    this->n_sym = n_sym;
    this->t_min = t_min;
    this->t_max = t_max;
    this->dt = dt;
    this->T0 = T0;
    this->gamma = gamma;

    parts = make_unique<array<data_t, TOTAL_PARTS>>();
    forces = make_unique<array<data_t, TOTAL_PARTS>>();
    moment = make_unique<array<data_t, TOTAL_PARTS>>();
    masses = make_unique<array<double, TOTAL_PARTS>>();

    fill(parts->begin(), parts->end(), data_t{0.0, 0.0, 0.0});
    fill(forces->begin(), forces->end(), data_t{0.0, 0.0, 0.0});
    fill(moment->begin(), moment->end(), data_t{0.0, 0.0, 0.0});
    fill(masses->begin(), masses->end(), 18.0);
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

void Simulator::fillVec(const string &path) {
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
    double u_ij = 0.0;

    for (int i_sym = 0; i_sym < n_sym; ++i_sym) {
        const double ix_sym = vec_translation[i_sym][0];
        const double iy_sym = vec_translation[i_sym][1];
        const double iz_sym = vec_translation[i_sym][2];

        for (int i = 0; i < n; ++i) {
            double ix, iy, iz;

            ix = (*parts)[i].x;
            iy = (*parts)[i].y;
            iz = (*parts)[i].z;

            if (l != 0.0) {
                ix = ix - int(ix * (1.0 / l)) * l;
                iy = iy - int(iy * (1.0 / l)) * l;
                iz = iz - int(iz * (1.0 / l)) * l;
            }

            for (int j = 0; j < n; ++j) {
                if (i == j)
                    continue;

                double jx, jy, jz;

                jx = (*parts)[j].x;
                jy = (*parts)[j].y;
                jz = (*parts)[j].z;

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

                if ((r_ij_sq > r_cut) || (r_ij_sq < 1e-1)) {
                    continue;
                }

                const double r2 = r_star_sq * (1.0 / r_ij_sq);
                const double r4 = r2 * r2;
                const double r6 = r4 * r2;
                const double r12 = r6 * r4 * r2;
                const double r8 = r6 * r2;
                const double r14 = r4 * r4 * r2 * r2 * r2;

                u_ij = 48.0 * epsilon_star * (r14 - r8);

                const double fx = u_ij * rx;
                const double fy = u_ij * ry;
                const double fz = u_ij * rz;
                
                (*forces)[j].x += -fx;
                (*forces)[j].y += -fy;
                (*forces)[j].z += -fz;

                (*forces)[j].x += fx;
                (*forces)[j].y += fy;
                (*forces)[j].z += fz;

                lj += r12 - 2 * r6;
            }
        }
    }

    return ((4.0 * epsilon_star) * (1.0 / 2.0)) * lj;
}

double Simulator::sumForces() {
    double sum_forces = 0.0;
    for (int i = 0; i < TOTAL_PARTS; ++i) {
        sum_forces += (*forces)[i].x + (*forces)[i].y + (*forces)[i].z;
    }

    return sum_forces;
}

void Simulator::velocityVerlet() {
    const double dt_half = 0.5 * dt;

    for (int i = 0; i < n; ++i) {
        double vx = (*moment)[i].x * (1 / (*masses)[i]);
        double vy = (*moment)[i].y * (1 / (*masses)[i]);
        double vz = (*moment)[i].z * (1 / (*masses)[i]);

        vx += (*forces)[i].x * dt_half * (1 / (*masses)[i]);
        vy += (*forces)[i].y * dt_half * (1 / (*masses)[i]);
        vz += (*forces)[i].z * dt_half * (1 / (*masses)[i]);

        (*parts)[i].x += vx * dt;
        (*parts)[i].y += vy * dt;
        (*parts)[i].z += vz * dt;
    }

    computeEnergyForces();

    for (int i = 0; i < n; ++i) {
        double vx = (*moment)[i].x * (1 / (*masses)[i]);
        double vy = (*moment)[i].y * (1 / (*masses)[i]);
        double vz = (*moment)[i].z * (1 / (*masses)[i]);

        vx += (*forces)[i].x * dt_half * (1 / (*masses)[i]);
        vy += (*forces)[i].y * dt_half * (1 / (*masses)[i]);
        vz += (*forces)[i].z * dt_half * (1 / (*masses)[i]);

        (*moment)[i].x = vx * (*masses)[i];
        (*moment)[i].y = vy * (*masses)[i];
        (*moment)[i].z = vz * (*masses)[i];
    }

    // double dt_half = 0.5 * dt;

    // // Velocity half-step
    // for (int i = 0; i < n; i++)
    // {
    //     (*moment)[i].x -= dt_half * (*forces)[i].x * (1 / (*masses)[i]);
    //     (*moment)[i].y -= dt_half * (*forces)[i].y * (1 / (*masses)[i]);
    //     (*moment)[i].z -= dt_half * (*forces)[i].z * (1 / (*masses)[i]);
    // }

    // // Position update
    // for (int i = 0; i < n; i++)
    // {
    //     (*parts)[i].x += dt * (*moment)[i].x;
    //     (*parts)[i].y += dt * (*moment)[i].y;
    //     (*parts)[i].z += dt * (*moment)[i].z;
    // }

    // // Compute new forces
    // computeEnergyForces();

    // // Velocity full-step
    // for (int i = 0; i < n; i++)
    // {
    //     (*moment)[i].x -= dt_half * (*forces)[i].x * (1 / (*masses)[i]);
    //     (*moment)[i].y -= dt_half * (*forces)[i].y * (1 / (*masses)[i]);
    //     (*moment)[i].z -= dt_half * (*forces)[i].z * (1 / (*masses)[i]);
    // }
}

double Simulator::computeKineticEnergy() {
    double sum = 0.0;

    for (int i = 0; i < n; ++i) {
        const double moment_x = (*moment)[i].x * (*moment)[i].x;
        const double moment_y = (*moment)[i].y * (*moment)[i].y;
        const double moment_z = (*moment)[i].z * (*moment)[i].z;
        sum += (moment_x + moment_y + moment_z) / (*masses)[i];
    }

    return (1.0 / (2.0 * conv_force)) * sum;
}

double Simulator::computeKineticTemperature() {
    const double energy = computeKineticEnergy();

    return (1.0 / (n_dl * const_r)) * energy;
}

void Simulator::correctionRatio() {
    const double energy_kinetic_init = computeKineticEnergy();
    const double ratio = sqrt((n_dl * const_r * T0) / energy_kinetic_init);

    for (int i = 0; i < n; ++i) {
        (*moment)[i].x *= ratio;
        (*moment)[i].y *= ratio;
        (*moment)[i].z *= ratio;
    }
}

void Simulator::correctionCenter() {
    double Px = 0.0;
    double Py = 0.0;
    double Pz = 0.0;

    for (int i = 0; i < n; ++i) {
        Px += (*moment)[i].x;
        Py += (*moment)[i].y;
        Pz += (*moment)[i].z;
    }

    for (int i = 0; i < n; ++i) {
        (*moment)[i].x -= Px / (double)TOTAL_PARTS;
        (*moment)[i].y -= Py / (double)TOTAL_PARTS;
        (*moment)[i].z -= Pz / (double)TOTAL_PARTS;
    }
}

void Simulator::fillMoment() {
    // srand(random_device()());
    srand(42);

    double px = 0.0;
    double py = 0.0;
    double pz = 0.0;

    for (int i = 0; i < n; ++i) {
        double c = 0.0;
        double s = 0.0;

        c = (double)rand() / (double)RAND_MAX;
        s = (double)rand() / (double)RAND_MAX;
        px = fn_sign(1.0, 0.5 - s) * c;

        c = (double)rand() / (double)RAND_MAX;
        s = (double)rand() / (double)RAND_MAX;
        py = fn_sign(1.0, 0.5 - s) * c;

        c = (double)rand() / (double)RAND_MAX;
        s = (double)rand() / (double)RAND_MAX;
        pz = fn_sign(1.0, 0.5 - s) * c;

        (*moment)[i].x = px;
        (*moment)[i].y = py;
        (*moment)[i].z = pz;
    }

    correctionRatio();
    correctionCenter();
    correctionRatio();
}

void Simulator::correctionMoment() {
    const double T = computeKineticTemperature();
    const double correction = gamma * ((T0 / T) - 1.0);

    for (int i = 0; i < n; ++i) {
        (*moment)[i].x += correction * (*moment)[i].x;
        (*moment)[i].y += correction * (*moment)[i].y;
        (*moment)[i].z += correction * (*moment)[i].z;
    }
}

void Simulator::start(const bool out, const int correction_step, const int save_step) {
    ofstream outfile(OUTPUT_PARTS_FILE, ios::trunc);

    if (!outfile.is_open()) {
        cerr << "Failed to open file " << OUTPUT_PARTS_FILE << " for writing!" << endl;
        return;
    }

    outfile.close();

    outfile.open(OUTPUT_PARTS_FILE, ios::app);

    if (!outfile.is_open()) {
        cerr << "Failed to open file " << OUTPUT_PARTS_FILE << " for appending!" << endl;
        return;
    }

    double ang = 0.0;
    double temp = 0.0;
    double t = get_t_min();
    int count = 0;
    const double dt_mod = 1.0 / dt;

    cout << "Début de la simulation" << endl;
    cout << "Temps initial: " << get_t_min() << " Temps final: " << get_t_max() << " dt: " << get_dt() << endl;

    ang = computeKineticEnergy();
    temp = computeKineticTemperature();

    cout << "Etat initial" << endl;
    cout << "time: " << t << " : Moment cinétique = " << ang << ", Température = " << temp << endl;

    for (; t < t_max; t += dt) {
        velocityVerlet();

        if (count % correction_step == 0) {
            correctionMoment();
        }

        if (out && (int)(t * dt_mod) % save_step == 0) {
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

    ang = computeKineticEnergy();
    temp = computeKineticTemperature();

    cout << "Etat Final" << endl;
    cout << "time: " << t << " : Moment cinétique = " << ang << ", Température = " << temp << endl;

    outfile.close();
}