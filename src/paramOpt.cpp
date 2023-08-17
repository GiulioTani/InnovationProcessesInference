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
    double alpha = (alpha0 == -100 ? ST_ALPHA : alpha0), theta;
    auto occurr = corpus.get_dictionary();
    int N = 0, T = (int)occurr->size();
    std::vector<int> mults;
    mults.reserve(T);
    for (auto &w : *occurr)
    {
      N += w.second;
      mults.push_back(w.second);
    }
    theta = (theta0 == -100 ? T : theta0);
    gradient_ascent(N, T, mults, &theta, &alpha);
    return {alpha, theta};
  }

  std::pair<double, double> param_opt(const std::vector<int> &mults,
                                      double alpha0, double theta0)
  {
    double alpha = (alpha0 == -100 ? ST_ALPHA : alpha0), theta;
    int N = 0, T = (int)mults.size();
    for (auto &w : mults)
      N += w;
    theta = (theta0 == -100 ? T : theta0);
    gradient_ascent(N, T, mults, &theta, &alpha);
    return {alpha, theta};
  }

  /**
 * @brief Stores in alpha and theta the
 *   \f[
 *       arg \max_{\alpha,theta}  \frac{(\theta |
 * \alpha)_K}{(\theta)_N}\prod_{k=1}^K (1-\alpha)_{n_k-1} \f] Where \f$ N \f$ is
 * the total number of elements and \f$ K \f$ is the total number of differents
 * elements.
 *
 * @param[in] N Total number of elements.
 * @param[in] T Total number of differents elements.
 * @param[in] nconk Number of occurrences of each different word.
 * @param[in,out] theta Theta parameter of the PD Process.
 * @param[in,out] alpha Alpha parameter of the PD Process.
 */
  void gradient_ascent(int N, int T, const std::vector<int> &nconk, double *theta,
                       double *alpha)
  {
    double deriv_alpha, old_deriv_alpha;
    double deriv_theta, old_deriv_theta;
    double alpha_new = *alpha, theta_new = *theta; // ST_THETA;//
    double Vtheta = 0, Valpha = 0, smorzamento = 0.9;
    int count = 0, tolleroAlpha = TOLLERANZA, tolleroTheta = TOLLERANZA;
    double incr_alpha = _incr_alpha, incr_theta = _incr_theta;
    std::unordered_map<int, int> Counter;

    for (auto nk : nconk)
      if (!Counter.emplace(nk, 1).second)
        Counter[nk]++;

    deriv_alpha = der_alpha(T, Counter, *theta, *alpha);
    deriv_theta = der_theta(N, T, *theta, *alpha);
    *alpha = 0;
    old_deriv_alpha = deriv_alpha;
    *theta = 0;
    old_deriv_theta = deriv_theta;
    while (abs((alpha_new - *alpha) / (*alpha)) > eps_alpha_inf ||
           abs((theta_new - *theta) / (*theta)) > eps_theta_inf)
    {
      if (alpha_new > 0 && alpha_new < 1)
      {
        *alpha = alpha_new;
      }
      else
      {
        if (!tolleroAlpha)
        { // se i parametri escono dall'intervallo di valori
          // possibili mi fermo
          if (isatty(0))
            std::cerr << "Convergence of alpha failed." << std::endl;
          if (alpha_new < 0)
            *alpha = 0.001;
          else
            *alpha = 0.999;
          break;
        }
        else
        {
          if (alpha_new > 1)
          {
            *alpha = 0.9;
            theta_new *= 2;
          }
          else
          {
            *alpha = 0.1;
            theta_new /= 2;
          }
          incr_alpha /= 2;
          tolleroAlpha--;
        }
      }
      if (theta_new > -*alpha)
      {
        *theta = theta_new;
      }
      else
      {
        if (!tolleroTheta)
        { // se i parametri escono dall'intervallo di valori
          // possibili mi fermo
          if (isatty(0))
            std::cerr << "Convergence of theta failed." << smorzamento
                      << std::endl;
          *theta = 0.01 - *alpha; // ST_THETA; //
          break;
        }
        else
        {
          //*alpha *= 0.9;
          *theta = 1;
          smorzamento /= 2;
          incr_theta /= 2;
          tolleroTheta--;
        }
      }
      /* calcolo il passo su alpha */
      old_deriv_alpha = deriv_alpha;
      deriv_alpha = der_alpha(T, Counter, *theta, *alpha);
      Valpha = old_deriv_alpha * deriv_alpha < 0
                   ? -0.5 * Valpha
                   : smorzamento * Valpha + deriv_alpha;
      alpha_new =
          *alpha +
          incr_alpha * Valpha; // evito il flipper facendo solo mezzo passo

      /* calcolo il passo su theta */
      old_deriv_theta = deriv_theta;
      deriv_theta = der_theta(N, T, *theta, *alpha);
      Vtheta = old_deriv_theta * deriv_theta < 0
                   ? -0.5 * Vtheta
                   : smorzamento * Vtheta + deriv_theta;
      theta_new =
          *theta +
          incr_theta * Vtheta; // evito il flipper facendo solo mezzo passo
      count++;
      if (!(count % 100000))
      {
        incr_alpha /= double(count / 100000 + 1) / (count / 100000);
        incr_theta /= double(count / 100000 + 1) / (count / 100000);
      }
      // if (tolleroTheta<5)std::cerr<<count <<" A"<<*alpha<<" "<<Valpha<<"
      // "<<deriv_alpha<<" T"<<*theta<<" "<<Vtheta<<" "<<deriv_theta<<std::endl; if
      // (count>400020)throw std::runtime_error("Not optimized");
    }
    // if (isatty(0))  std::cerr << "(" << count << "-" << N << "-" << T << ")\t";
    // std::cerr << count << std::endl;
    if (count == 1)
    {
      std::cerr << "Not optimized" << std::endl;
      throw std::runtime_error("Not optimized");
    }
  }

  double digamma_imp_large(double x)
  {
    static const double P[] = {
        0.083333333333333333333333333333333333333333333333333,
        -0.0083333333333333333333333333333333333333333333333333,
        0.003968253968253968253968253968253968253968253968254};
    x -= 1;
    double result = log(x);
    result += 1 / (2 * x);
    double z = 1 / (x * x);
    result -= ((P[2] * z + P[1]) * z + P[0]) * z;
    return result;
  }

  double digamma(double x)
  {
    double result = 0;

    while (x < 10)
    {
      result -= 1 / x;
      x += 1;
    }
    result += digamma_imp_large(x);
    return result;
  }

  /**
 * @brief Computes the derivative of log(P) with respect to alpha. Using the
 * formula: \f[ \partial_{\alpha}P \propto \sum_{i=0}^{T-1} \frac{i}{\theta +
 * \alpha i} - \sum_{k=1}^K \sum_{\nu\in [1,n_k-1]} \frac{ 1}{\nu-\alpha} \f]
 *
 *    For optimization purposes (avoiding nested cycles) the formula isn't used
 * as it is but the terms with equal \f$ n_k \f$ are grouped together and all
 * the sums are computed in a row. This trick reduces the execution time by
 * roughly two thirds.
 *
 * @param T Total number of differents elements.
 * @param Counter A map pairing the number of occurrences of a word and the
 * number of words sharing the same number of occurrences.
 * @param maxnk Nuber of occurrences of the most common word(s).
 * @param theta Present value of the theta parameter of the PD Process.
 * @param alpha Present value of the alpha parameter of the PD Process.
 * @return double The value of the derivative.
 */
  double der_alpha(int T, const std::unordered_map<int, int> &Counter,
                   double theta, double alpha)
  {
    double sumleft = 0, sumright = 0;
    sumleft = (theta * (digamma(theta / alpha) - digamma(theta / alpha + T)) +
               T * alpha) /
              (alpha * alpha);
    for (auto &pa : Counter)
      sumright += (digamma(pa.first - alpha) - digamma(1 - alpha)) * pa.second;
    return sumleft - sumright;
  }

  /**
 * @brief Computes the derivative of log(P) with respect to theta. Using the
 * formula: \f[ \partial_{\theta}P\propto \sum_{i=0}^{T-1} \frac{1}{\theta +
 * \alpha i} - \sum_{i=0}^{N-1} \frac{1}{\theta + i} \f]
 *
 * @param N Total number of elements.
 * @param T Total number of differents elements.
 * @param theta Present value of the theta parameter of the PD Process.
 * @param alpha Present value of the alpha parameter of the PD Process.
 * @return double The value of the derivative.
 */
  double der_theta(int N, int T, double theta, double alpha)
  {
    double sumleft = 0, sumright = 0;
    sumleft = (digamma(theta / alpha + T) - digamma(theta / alpha)) / alpha;
    sumright = digamma(theta + N) - digamma(theta);
    return sumleft - sumright;
  }
} // namespace popt
