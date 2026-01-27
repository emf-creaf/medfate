#include <vector>

#ifndef NUMERICAL_SOLVING_C_H
#define NUMERICAL_SOLVING_C_H

void tridiagonalSolving_c(const std::vector<double> &a, const std::vector<double> &b, const std::vector<double> &c, const std::vector<double> &d,
                          std::vector<double> &e, std::vector<double> &f, std::vector<double> &sol);
#endif
