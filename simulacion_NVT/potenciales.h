#pragma once
#include <iostream>
#include <random>
#include <fstream>

void lennard_jones(const int N, const float rc, const float rc2, const float L, const float V, double* E_p, double* dphi, double* d2phi, const double rx[], const double ry[], const double rz[], double ax[], double ay[], double az[]);
