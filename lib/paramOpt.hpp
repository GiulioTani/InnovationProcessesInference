#pragma once
#include "lib/bookprob.hpp"
#include <vector>
#include <unordered_map>

namespace popt
{
    std::pair<double, double> param_opt(const bookprob::book &corpus, double alpha0 = -100, double theta0 = -100);
    std::pair<double, double> param_opt(const std::vector<int> &mults, double alpha0 = -100, double theta0 = -100);
    void gradient_ascent(int N, int T, const std::vector<int> &nconk, double *theta, double *alpha);
    double der_alpha(int T, const std::unordered_map<int, int> &Counter, double theta, double alpha);
    double der_theta(int N, int T, double theta, double alpha);
} // namespace popt