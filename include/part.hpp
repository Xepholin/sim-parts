#pragma once

struct Particule {
    double x = 0;
    double y = 0;
    double z = 0;

    Particule(double new_x = 0, double new_y = 0, double new_z = 0);
    ~Particule() = default;
};