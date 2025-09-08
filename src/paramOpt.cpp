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

/**
 * @file paramOpt.cpp
 * @author Giulio Tani Raffaelli (tani@cs.cas.cz)
 * @brief Optimizes the parameters alpha and theta following a maximum
 * likelyhood principle assuming that every line in the input files (or stdin)
 * is an element of a PD Process.
 * @version 1.0
 * @date 2023-06-27
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "lib/paramOpt.hpp"   //defines namespace popt
#include "lib/bookprob.hpp"   //defines namespace book
#include "lib/supportlib.hpp" //defines namespace supp
#include "math.h"
#include <csignal>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <unordered_map>
#include <vector>
#include <unistd.h>

#define ST_THETA 1000
#define ST_ALPHA 0.3
#define eps_alpha_inf 1e-7 // 1e-6
#define eps_theta_inf 1e-7 // 1e-6
#define _incr_alpha 1e-4
#define _incr_theta 1
#define TOLLERANZA 50

namespace popt
{
  std::pair<double, double> param_opt(const bookprob::book &corpus, double alpha0,
                                      double theta0)
  {
    return param_opt(corpus.get_entry_t(), alpha0, theta0);
  }

  std::pair<double, double> param_opt(const std::vector<int> &t,
                                      double alpha0, double theta0)
  {
    std::vector<double> Dt(t.size(), 0), _t;
    for (size_t i = 0; i < t.size(); i++)
      Dt[i] = i + 1;
    _t.insert(_t.end(), t.begin(), t.end());
    /** */
    alglib::real_1d_array y;
    alglib::real_2d_array x;
    y.attach_to_ptr(Dt.size(), Dt.data());
    x.attach_to_ptr(_t.size(), 1, _t.data());
    /** */
    alglib::real_1d_array c;
    double par[2] = {(alpha0 == -100 ? ST_ALPHA : alpha0), (theta0 == -100 ? t.size() : theta0)};
    c.setcontent(2, par);
    /** */
    double epsx = 0.000001;
    size_t maxits = 0;
    alglib::lsfitstate state;
    alglib::lsfitreport rep;

    alglib::lsfitcreatefg(x, y, c, state);
    alglib::lsfitsetcond(state, epsx, maxits);
    alglib::lsfitfit(state, function_cx_1_func, function_cx_1_grad);
    alglib::lsfitresults(state, c, rep);
    return {c[0], c[1]};
  }

  void function_cx_1_func(const alglib::real_1d_array &c, const alglib::real_1d_array &x, double &func, void *ptr)
  {
    // this callback calculates f(c,x)=exp(-c0*sqr(x0))
    // where x is a position on X-axis and c is adjustable parameter
    func = c[1] / c[0] * (pow(1 + x[0] / c[1], c[0]) - 1);
  }
  void function_cx_1_grad(const alglib::real_1d_array &c, const alglib::real_1d_array &x, double &func, alglib::real_1d_array &grad, void *ptr)
  {
    // this callback calculates f(c,x)=exp(-c0*sqr(x0)) and gradient G={df/dc[i]}
    // where x is a position on X-axis and c is adjustable parameter.
    // IMPORTANT: gradient is calculated with respect to C, not to X
    func = c[1] / c[0] * (pow(1 + x[0] / c[1], c[0]) - 1);
    grad[0] = (c[1] * pow(1 + x[0] / c[1], c[0]) * log10(1 + x[0] / c[1]) - func) / c[0];
    grad[1] = func / c[0] - x[0] / c[1] * pow(1 + x[0] / c[1], c[0] - 1);
  }
}