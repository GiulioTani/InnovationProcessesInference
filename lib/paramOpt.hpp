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