#pragma once

#include <array>
#include <memory>
#include <string>

#include "const.hpp"

using namespace std;

#pragma pack()
struct alignas(32) data_t {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

class Simulator {
public:
    Simulator(const int n, const int n_sym, const double t_min, const double t_max, const double dt, const double T0, const double gamma);
    ~Simulator();
    
    double get_n();
    double get_n_sym();
    double get_t_min();
    double get_t_max();
    double get_dt();

    void fillVec(const string &path);
    double computeEnergyForces();
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

    unique_ptr<array<data_t, TOTAL_PARTS>> parts;
    unique_ptr<array<data_t, TOTAL_PARTS>> forces;
    unique_ptr<array<data_t, TOTAL_PARTS>> moment;
    unique_ptr<array<double, TOTAL_PARTS>> masses;
};
