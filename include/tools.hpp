#pragma once

#include <memory>
#include <vector>

using namespace std;

template <typename T>
void init_f64(const int n, unique_ptr<T> &__restrict__ a, const double val);

double mean(const int n, const unique_ptr<vector<double>> &arr);
double stddev(const int n, const unique_ptr<vector<double>> &arr);
double fn_sign(double a, double s);
