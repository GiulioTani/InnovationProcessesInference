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

namespace bookprob
{
    class dkl_book : public book
    {
        mutable bool cache_valid;
        mutable std::vector<std::pair<hash_type, double>> frequency /** The order of the words in the book with their frequency.*/;
        mutable double corr, factP0;

        void build_cache() const;

    public:
        dkl_book(const dkl_book &oth);
        dkl_book(const book &oth) : book(oth) { cache_valid = false; };
        dkl_book(dkl_book &&oth);
        dkl_book(const std::string &path, const glob_prob *inP0 = NULL, double Alpha = -100, double Theta = -100, std::string pref = "") : book(path, inP0, Alpha, Theta, pref) { cache_valid = false; };
        dkl_book(const std::vector<std::string> &list, const glob_prob *inP0 = NULL, double Alpha = -100, double Theta = -100, std::string pref = "") : book(list, inP0, Alpha, Theta, pref) { cache_valid = false; };
        dkl_book(const std::vector<hash_type> &list, const glob_prob *inP0 = NULL, double Alpha = -100, double Theta = -100, std::string pref = "") : book(list, inP0, Alpha, Theta, pref) { cache_valid = false; };
        dkl_book(const std::vector<std::pair<hash_type, int>> &ord, const glob_prob *inP0 = NULL, double Alpha = -100, double Theta = -100, std::string pref = "") : book(ord, inP0, Alpha, Theta, pref) { cache_valid = false; };
        dkl_book join(const dkl_book &other) const { return book::join(other); };
        void append(const book &other);
        std::pair<double, double> log_prob(const book &other) const;
        size_t retsize() const;
        virtual void log_prob_to_chararr(const book &other, void* dest) const;
        void log_prob_to_chararr(void *dest) const;
    };
}