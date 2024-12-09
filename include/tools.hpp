#pragma once

#include <memory>
#include <vector>

using namespace std;

double mean(const unique_ptr<vector<double>> &arr, int size);
double stddev(const unique_ptr<vector<double>> &arr, int size);