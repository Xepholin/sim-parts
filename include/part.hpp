#pragma once

#include <array>
#include <memory>
#include <string>

#include "const.hpp"

using namespace std;

#pragma pack()
struct alignas(64) data_t {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

class Simulator {
public:
    Simulator(const int n, const int n_sym, const double t_min, const double t_max, const double dt, const double T0, const double gamma, const bool verlet);
    ~Simulator();
    
    double get_n();
    double get_n_sym();
    double get_t_min();
    double get_t_max();
    double get_dt();
    string get_verlet();

    void fillVec(const string &path);
    void verletList();
    double computeEnergyForces();
    double computeEnergyForcesVerlet();
    double sumForces();
    void velocityVerlet();
    double computeKineticEnergy();
    double computeKineticTemperature();
    void correctionRatio();
    void correctionCenter();
    void fillMoment();
    void correctionMoment();
    void start(const bool out, const int correction_step, const int save_step);

private:
    int n = 0;
    int n_sym = 0;
    double t_min = 0.0;
    double t_max = 0.0;
    double dt = 0.0;
    double T0 = 0.0;
    double gamma = 0.0;

    bool verlet = false;

    unique_ptr<array<data_t, TOTAL_PARTS>> parts;
    unique_ptr<array<data_t, TOTAL_PARTS>> forces;
    unique_ptr<array<data_t, TOTAL_PARTS>> moments;
    unique_ptr<array<data_t, TOTAL_PARTS>> velocities;
    unique_ptr<array<double, TOTAL_PARTS>> masses;
    unique_ptr<array<array<int, n_max_neighbor * 2>, TOTAL_PARTS>> neighbor;
};
