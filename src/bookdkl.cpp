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

#include "lib/bookdkl.hpp"
#include <math.h>
#include <iostream>
#include <array>
#define CROSSENTROPY
namespace bookprob
{
    dkl_book::dkl_book(const dkl_book &oth) : dkl_book::book(oth)
    {
        cache_valid = oth.cache_valid;
        if (cache_valid)
        {
            frequency = oth.frequency;
            corr = oth.corr;
            factP0 = oth.factP0;
        }
    }
    dkl_book::dkl_book(dkl_book &&oth) : dkl_book::book(oth)
    {
        cache_valid = oth.cache_valid;
        if (cache_valid)
        {
            frequency = oth.frequency;
            corr = oth.corr;
            factP0 = oth.factP0;
        }
        oth.cache_valid = false;
    }
    void dkl_book::append(const book &other)
    {
        book::join(other);
        cache_valid = false;
    };
    void dkl_book::build_cache() const
    {
        frequency.clear();
        frequency.reserve(order.size());
        for (auto w : order)
        {
            frequency.push_back({w.first, double(w.second) / N});
        }
        corr = alpha / N;
        factP0 = (theta + alpha * order.size()) / N;
        cache_valid = true;
    };

    std::pair<double, double> dkl_book::log_prob(const book &other) const
    {
        const dkl_book &other_ = other;
        double P = 0;
        double rest = 0;
        if (other_.order.size() == 0)
            throw not_init("Other order " + other_.prefix);
        if (order.size() == 0)
            throw not_init("order");
        if (positions.size() == 0)
            throw not_init("positions");
        if (P0 == NULL)
            throw not_init("P0");
        if (alpha == -100)
            throw not_init("alpha");
        if (theta == -100)
            throw not_init("theta");
        if (other_.prefix == "" || prefix == "")
            throw not_init("prefix");
        if (!cache_valid)
            build_cache();
        if (!other_.cache_valid)
            other_.build_cache();
        init_tmpP0(other);
        for (auto &p : other_.frequency)
        {
            if (positions.find(p.first) != positions.end())
#ifdef CROSSENTROPY
                P += p.second * log(frequency[positions.at(p.first)].second - corr); // if the word was already in my vocabulary
#else
                P += p.second * (log(p.second) - log(frequency[positions.at(p.first)].second - corr)); // if the word was already in my vocabulary
#endif // CROSSENTROPY
            else
            {
                try
                {
#ifdef CROSSENTROPY
                    P += p.second * log(factP0 * P0->prob.at(p.first) / nowWeight);
#else
                    P += p.second * (log(p.second) - log(factP0 * P0->prob.at(p.first) / nowWeight));
#endif // CROSSENTROPY
                    rest += p.second;
                }
                catch (std::out_of_range &e)
                {
                    std::cerr << prefix << " | " << other_.prefix << " | "
                              << " \"" << p.first << "\"" << std::endl;
                    throw word_miss(std::to_string(p.first));
                }
            }
        }
#ifdef CROSSENTROPY
        P += log(N / (theta + N));
        return {P / log(10), rest}; // to have the logarithm in base 10
    };
#else
        P -= log(N / (theta + N));                            // CROSSENTROPY
        return {P / log(10), -rest}; // to have the logarithm in base 10
    };
#endif  
    size_t dkl_book::retsize() const
    {
        return sizeof(double) + sizeof(double);
    };

    void dkl_book::log_prob_to_chararr(void *dest) const
    {
        constexpr std::array<unsigned char, 8> bb{0, 0, 0, 0, 0, 0, 240, 255};
        constexpr double dd = 0;
        memcpy(dest, &bb, sizeof(bb));
        memcpy((char *)dest + sizeof(bb), &dd, sizeof(dd));
    };

    void dkl_book::log_prob_to_chararr(const book &other, void *dest) const
    {
        auto tmp = log_prob(other);
        memcpy(dest, &tmp.first, sizeof(tmp.first));
        memcpy((char *)dest + sizeof(tmp.first), &tmp.second, sizeof(tmp.second));
    };
} // namespace bookprob
