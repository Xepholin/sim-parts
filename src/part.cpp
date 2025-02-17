#include "part.hpp"

#include <omp.h>
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

    this->n_dl = 3.0 * n - 3.0;
    this->n_neighbor = (double)n / (l * l * l);
    this->n_max_neighbor = (int)(n_neighbor * 4.0 * M_PI * (r_cut_sq * r_cut)) * 2.0;

    parts.resize(n, data_t{0.0, 0.0, 0.0});
    forces.resize(n, data_t{0.0, 0.0, 0.0});
    moments.resize(n, data_t{0.0, 0.0, 0.0});
    velocities.resize(n, data_t{0.0, 0.0, 0.0});
    masses.resize(n, 18.0);

    neighbor.resize(n, std::vector<int>(n_max_neighbor * 2, 0));
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
            ss >> parts[i].x >> parts[i].y >> parts[i].z;
        }
    } else {
        cerr << "File could not be opened!\n";
        cerr << "Error code: " << strerror(errno) << '\n';
    }

    fs.close();
}

double Simulator::computeEnergyForces() {
    std::fill(forces.begin(), forces.end(), data_t{0.0, 0.0, 0.0});
    double lj = 0.0;

#pragma omp parallel for reduction(+ : lj)
    for (int i_sym = 0; i_sym < n_sym; ++i_sym) {
        const double ix_sym = vec_translation[i_sym][0];
        const double iy_sym = vec_translation[i_sym][1];
        const double iz_sym = vec_translation[i_sym][2];

        for (int i = 0; i < n; ++i) {
            double ix = parts[i].x - ix_sym;
            double iy = parts[i].y - iy_sym;
            double iz = parts[i].z - iz_sym;

            for (int j = i + 1; j < n; ++j) {
                double dx = ix - parts[j].x;
                double dy = iy - parts[j].y;
                double dz = iz - parts[j].z;
                double r_ij_sq = dx * dx + dy * dy + dz * dz;

                if (r_ij_sq > r_cut_sq || r_ij_sq < 1e-6) continue;

                const double r2 = r_star_sq * (1.0 / r_ij_sq);
                const double r6 = r2 * r2 * r2;
                const double r12 = r6 * r6;
                const double r8 = r6 * r2;
                const double r14 = r8 * r6;

                double f_ij = -48.0 * epsilon_star * (r14 - 0.5 * r8);

                const double fx = f_ij * dx;
                const double fy = f_ij * dy;
                const double fz = f_ij * dz;

#pragma omp atomic
                forces[i].x += fx;
#pragma omp atomic
                forces[i].y += fy;
#pragma omp atomic
                forces[i].z += fz;

#pragma omp atomic
                forces[j].x -= fx;
#pragma omp atomic
                forces[j].y -= fy;
#pragma omp atomic
                forces[j].z -= fz;

                lj += r12 - r6;
            }
        }
    }

    return 4.0 * epsilon_star * lj;
}

void Simulator::verletList() {
#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        neighbor[i][0] = 0;
    }

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        int ni = 0;
        for (int j = i + 1; j < n; ++j) {
            for (int i_sym = 0; i_sym < n_sym; ++i_sym) {
                double xj_loc = parts[j].x + vec_translation[i_sym][0];
                double yj_loc = parts[j].y + vec_translation[i_sym][1];
                double zj_loc = parts[j].z + vec_translation[i_sym][2];

                double dx = parts[i].x - xj_loc;
                double dy = parts[i].y - yj_loc;
                double dz = parts[i].z - zj_loc;
                double r_ij_sq = dx * dx + dy * dy + dz * dz;

                if (r_ij_sq < r_cut_sq) {
#pragma omp atomic
                    ni += 1;

                    int mi = 2 * ni - 1;
                    neighbor[i][mi] = j;
                    neighbor[i][mi + 1] = i_sym;
                }
            }
        }
        neighbor[i][0] = ni;
    }
}

double Simulator::computeEnergyForcesVerlet() {
    std::fill(forces.begin(), forces.end(), data_t{0.0, 0.0, 0.0});
    double lj = 0.0;

#pragma omp parallel for reduction(+ : lj)
    for (int i = 0; i < n; ++i) {
        double ix = parts[i].x;
        double iy = parts[i].y;
        double iz = parts[i].z;

        for (int ni = 1; ni <= neighbor[i][0]; ++ni) {
            int j = neighbor[i][2 * ni - 1];
            int i_sym = neighbor[i][2 * ni];

            double xj_loc = parts[j].x + vec_translation[i_sym][0];
            double yj_loc = parts[j].y + vec_translation[i_sym][1];
            double zj_loc = parts[j].z + vec_translation[i_sym][2];

            double dx = ix - xj_loc;
            double dy = iy - yj_loc;
            double dz = iz - zj_loc;
            double r_ij_sq = dx * dx + dy * dy + dz * dz;

            if (r_ij_sq > r_cut_sq || r_ij_sq < 1e-6) {
                continue;
            }

            const double r2 = r_star_sq * (1.0 / r_ij_sq);
            const double r6 = r2 * r2 * r2;
            const double r12 = r6 * r6;
            const double r8 = r6 * r2;
            const double r14 = r8 * r6;

            double f_ij = -48.0 * epsilon_star * (r14 - 0.5 * r8);

            const double fx = f_ij * dx;
            const double fy = f_ij * dy;
            const double fz = f_ij * dz;

#pragma omp atomic
            forces[i].x += fx;
#pragma omp atomic
            forces[i].y += fy;
#pragma omp atomic
            forces[i].z += fz;

#pragma omp atomic
            forces[j].x -= fx;
#pragma omp atomic
            forces[j].y -= fy;
#pragma omp atomic
            forces[j].z -= fz;

            lj += r12 - r6;
        }
    }

    return 4.0 * epsilon_star * lj;
}

double Simulator::sumForces() {
    double sum_forces = 0.0;

#pragma omp parallel for reduction(+ : sum_forces)
    for (int i = 0; i < n; ++i) {
        sum_forces += forces[i].x + forces[i].y + forces[i].z;
    }

    return sum_forces;
}

void Simulator::velocityVerlet() {
    double dt_half = 0.5 * dt;
    double l_half = 0.5 * l;

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        moments[i].x -= dt_half * conv_force * forces[i].x;
        moments[i].y -= dt_half * conv_force * forces[i].y;
        moments[i].z -= dt_half * conv_force * forces[i].z;
    }

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        double inv_mass = 1.0 / masses[i];
        parts[i].x += dt * moments[i].x * inv_mass;
        parts[i].y += dt * moments[i].y * inv_mass;
        parts[i].z += dt * moments[i].z * inv_mass;
    }

    if (l != 0.0) {
#pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            parts[i].x -= l * int(parts[i].x / l_half);
            parts[i].y -= l * int(parts[i].y / l_half);
            parts[i].z -= l * int(parts[i].z / l_half);
        }
    }

    if (verlet) {
        computeEnergyForcesVerlet();
    } else {
        computeEnergyForces();
    }

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        moments[i].x -= dt_half * conv_force * forces[i].x;
        moments[i].y -= dt_half * conv_force * forces[i].y;
        moments[i].z -= dt_half * conv_force * forces[i].z;
    }
}

double Simulator::computeKineticEnergy() {
    double sum = 0.0;

#pragma omp parallel for reduction(+ : sum)
    for (int i = 0; i < n; ++i) {
        double px = moments[i].x;
        double py = moments[i].y;
        double pz = moments[i].z;
        double inv_mass = 1.0 / masses[i];

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

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        auto& p = moments[i];
        p.x *= ratio;
        p.y *= ratio;
        p.z *= ratio;
    }
}

void Simulator::correctionCenter() {
    double Px = 0.0;
    double Py = 0.0;
    double Pz = 0.0;

#pragma omp parallel for reduction(+ : Px, Py, Pz)
    for (int i = 0; i < n; ++i) {
        Px += moments[i].x;
        Py += moments[i].y;
        Pz += moments[i].z;
    }

    Px /= n;
    Py /= n;
    Pz /= n;

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        moments[i].x -= Px;
        moments[i].y -= Py;
        moments[i].z -= Pz;
    }
}

void Simulator::fillMoment() {
    // srand(random_device()());
    srand(42);

#pragma omp parallel for
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

        moments[i].x = px;
        moments[i].y = py;
        moments[i].z = pz;
    }

    correctionRatio();
    correctionCenter();
    correctionRatio();
}

void Simulator::correctionMoment() {
    const double T = computeKineticTemperature();
    const double lambda = gamma * ((T0 / T) - 1.0);

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        moments[i].x += lambda * moments[i].x;
        moments[i].y += lambda * moments[i].y;
        moments[i].z += lambda * moments[i].z;
    }
}

void Simulator::start(const bool out, const int correction_step, const int save_step, const string output) {
    ofstream outfile(output, ios::trunc);

    if (!outfile.is_open()) {
        cerr << "Failed to open file " << output << " for writing!" << endl;
        return;
    }

    outfile << "MODEL 1\n";

    for (int i = 0; i < n; ++i) {
        outfile << "ATOM  "
                << setw(5) << i + 1
                << "  C   MOL     1    "
                << fixed << setprecision(3)
                << setw(8) << parts[i].x
                << setw(8) << parts[i].y
                << setw(8) << parts[i].z
                << "  1.00  0.00           C  \n";
    }

    outfile << "ENDMDL\n";

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

        if (verlet && count % 20 == 0) {
            verletList();
        }

        if (count % (int)(t_max / 5) == 0) {
            cout << "Itération: " << count << " | Énergie: " << ang << " | Température: " << temp << '\n';
        }

        if (out && count % save_step == 0) {
            outfile << "MODEL " << count / save_step << "\n";

            for (int i = 0; i < n; ++i) {
                outfile << "ATOM  "
                        << setw(5) << i + 1
                        << "  C   MOL     1    "
                        << fixed << setprecision(3)
                        << setw(8) << parts[i].x
                        << setw(8) << parts[i].y
                        << setw(8) << parts[i].z
                        << "  1.00  0.00           C  \n";
            }

            outfile << "ENDMDL\n";
        }

        count++;
    }

    outfile.close();
}
