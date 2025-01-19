#pragma once

#include <memory>
#include <vector>

using namespace std;

template <typename T>
void empty_vec(const int n, unique_ptr<T> &__restrict__ vec);

double mean(const unique_ptr<vector<double>> &arr, int size);
double stddev(const unique_ptr<vector<double>> &arr, int size);
