// CP2D -- Constrained Probability Poisson-Dirichlet
// Copyright (C) 2023  Giulio Tani Raffaelli
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#pragma once
#include "lib/bookprob.hpp"
#include "alglib/interpolation.h"
#include <vector>
#include <unordered_map>

namespace popt
{
    std::pair<double, double> param_opt_fit(const bookprob::book &corpus, double alpha0 = -100, double theta0 = -100);
    std::pair<double, double> param_opt_fit(std::vector<double> &Dt, double alpha0 = -100, double theta0 = -100);
    void function_cx_1_func(const alglib::real_1d_array &c, const alglib::real_1d_array &x, double &func, void *ptr);
    void function_cx_1_grad(const alglib::real_1d_array &c, const alglib::real_1d_array &x, double &func, alglib::real_1d_array &grad, void *ptr);
    std::pair<double, double> param_opt_ML(const bookprob::book &corpus, double alpha0 = -100, double theta0 = -100);
    std::pair<double, double> param_opt_ML(const std::vector<int> &mults, double alpha0 = -100, double theta0 = -100);
    void gradient_ascent(int N, int T, const std::vector<int> &nconk, double *theta, double *alpha);
    double der_alpha(int T, const std::unordered_map<int, int> &Counter, double theta, double alpha);
    double der_theta(int N, int T, double theta, double alpha);
    std::pair<double, double> param_opt_fit2D(const bookprob::book &corpus, double alpha0 = -100, double theta0 = -100);
    std::pair<double, double> param_opt_fit2D(std::vector<double> &Dt, double alpha0 = -100, double theta0 = -100);
    void function_cx_2_func(const alglib::real_1d_array &c, const alglib::real_1d_array &x, double &func, void *ptr);
    void function_cx_2_grad(const alglib::real_1d_array &c, const alglib::real_1d_array &x, double &func, alglib::real_1d_array &grad, void *ptr);
} // namespace popt