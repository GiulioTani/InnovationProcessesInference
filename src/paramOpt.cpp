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
  std::pair<double, double> param_opt_fit(const bookprob::book &corpus, double alpha0,
                                          double theta0)
  {
    std::vector<double> Dt(corpus.get_N(), 0);
    size_t i = 0;
    const std::vector<int> &times = corpus.get_entry_t();
    for (size_t D = 0; D < times.size(); D++)
      for (; i < times[D]; i++)
        Dt[i] = D;
    for (; i < corpus.get_N(); i++)
      Dt[i] = times.size();
    return param_opt_fit(Dt, alpha0, theta0);
  }

  std::pair<double, double> param_opt_fit(std::vector<double> &Dt,
                                          double alpha0, double theta0)
  {
    std::vector<double> t(Dt.size(), 0);
    for (size_t i = 0; i < t.size(); i++)
      t[i] = i;
    /** */
    alglib::real_1d_array y;
    alglib::real_2d_array x;
    y.attach_to_ptr(Dt.size(), Dt.data());
    x.attach_to_ptr(t.size(), 1, t.data());
    /** */
    alglib::real_1d_array c;
    double par[2] = {(alpha0 == -100 ? ST_ALPHA : alpha0), (theta0 == -100 ? t.size() : theta0)};
    c.setcontent(2, par);
    /** */
    double epsx = 0.000001;
    size_t maxits = 0;
    alglib::lsfitstate state;
    alglib::lsfitreport rep;
    alglib::real_1d_array bndl = "[0, 1]";
    alglib::real_1d_array bndu = "[0.99999999, +inf]";

    alglib::lsfitcreatefg(x, y, c, state);
    alglib::lsfitsetcond(state, epsx, maxits);
    alglib::lsfitsetbc(state, bndl, bndu);
    alglib::lsfitfit(state, function_cx_1_func, function_cx_1_grad);
    alglib::lsfitresults(state, c, rep);
    return {c[0], c[1]};
  }

  std::pair<double, double> old_param_opt_fit(std::vector<double> &Dt,
                                              double alpha0, double theta0)
  {
    std::vector<double> t(Dt.size(), 0);
    for (size_t i = 0; i < t.size(); i++)
      t[i] = i + 1;
    /** */
    alglib::real_1d_array y;
    alglib::real_2d_array x;
    y.attach_to_ptr(Dt.size(), Dt.data());
    x.attach_to_ptr(t.size(), 1, t.data());
    /** */
    alglib::real_1d_array c;
    double par[2] = {(alpha0 == -100 ? ST_ALPHA : alpha0), (theta0 == -100 ? t.size() : theta0)};
    c.setcontent(2, par);
    /** */
    double epsx = 0.000001;
    size_t maxits = 0;
    alglib::lsfitstate state;
    alglib::lsfitreport rep;
    alglib::real_1d_array bndl = "[0, 1]";
    alglib::real_1d_array bndu = "[0.99999999, +inf]";
    alglib::real_1d_array scale;
    double scales[2] = {1, theta0};
    scale.setcontent(2, scales);
    alglib::lsfitcreatefg(x, y, c, state);
    alglib::lsfitsetcond(state, epsx, maxits);
    alglib::lsfitsetbc(state, bndl, bndu);
    alglib::lsfitsetscale(state, scale);
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

  std::pair<double, double> param_opt_ML(const bookprob::book &corpus, double alpha0,
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

  std::pair<double, double> param_opt_ML(const std::vector<int> &mults,
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
            *alpha = 0.01;
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

  std::pair<double, double> param_opt_fit2D(const bookprob::book &corpus, double alpha0, double theta0)
  {
    std::vector<double> Dt(corpus.get_N(), 0);
    size_t i = 0;
    const std::vector<int> &times = corpus.get_entry_t();
    for (size_t D = 0; D < times.size(); D++)
      for (; i < times[D]; i++)
        Dt[i] = D;
    for (; i < corpus.get_N(); i++)
      Dt[i] = times.size();
    return param_opt_fit2D(Dt, alpha0, theta0);
  }
  std::pair<double, double> param_opt_fit2D(std::vector<double> &Dt, double alpha0, double theta0)
  {
    std::vector<double> DN(2 * Dt.size(), 0), delta(Dt.size(), 0);
    for (size_t i = 0; i < Dt.size(); i++)
    {
      DN[2 * i] = Dt[i];
      DN[2 * i + 1] = i;
    }
    for (size_t i = 1; i < Dt.size(); i++)
      delta[i] = Dt[i] - Dt[i - 1];
    /** */
    alglib::real_1d_array y;
    alglib::real_2d_array x;
    y.attach_to_ptr(delta.size(), delta.data());
    x.attach_to_ptr(DN.size(), 2, DN.data());
    /** */
    alglib::real_1d_array c;
    double par[2] = {(alpha0 == -100 ? ST_ALPHA : alpha0), (theta0 == -100 ? Dt.back() : theta0)};
    c.setcontent(2, par);
    /** */
    double epsx = 0.000001;
    size_t maxits = 0;
    alglib::lsfitstate state;
    alglib::lsfitreport rep;
    alglib::real_1d_array bndl = "[0, 1]";
    alglib::real_1d_array bndu = "[0.99999999, +inf]";
    alglib::real_1d_array scale;
    double scales[2] = {1, theta0};
    scale.setcontent(2, scales);

    try
    {
      alglib::lsfitcreatefg(x, y, c, Dt.size(), 2, 2, state);
      alglib::lsfitsetcond(state, epsx, maxits);
      alglib::lsfitsetbc(state, bndl, bndu);
      alglib::lsfitsetscale(state, scale);
      alglib::lsfitfit(state, function_cx_2_func, function_cx_2_grad);
    }
    catch (alglib::ap_error alglib_exception)
    {
      printf("ALGLIB exception with message '%s'\n", alglib_exception.msg.c_str());
      throw alglib_exception;
    }
    alglib::lsfitresults(state, c, rep);
    return {c[0], c[1]};
  }
  void function_cx_2_func(const alglib::real_1d_array &c, const alglib::real_1d_array &x, double &func, void *ptr)
  {
    func = (c[1] + x[0] * c[0]) / (c[1] + x[1]);
  }
  void function_cx_2_grad(const alglib::real_1d_array &c, const alglib::real_1d_array &x, double &func, alglib::real_1d_array &grad, void *ptr)
  {
    func = (c[1] + x[0] * c[0]) / (c[1] + x[1]);
    grad[0] = x[0] / (c[1] + x[1]);
    grad[1] = (1 - func) / (c[1] + x[1]);
  }
}