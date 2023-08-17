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
#include "lib/supportlib.hpp"
#include "lib/datatypes.hpp"
#include <map>
#include <mutex>
#include <vector>
#include <set>

namespace asp
{
    extern dt::book_id_t authCount /**< Used to count the number of virtual authors defined.*/;
    extern int authorDimension, ngramSize;
    extern std::mutex autnum_lock;
    extern std::map<dt::auth_id_t, std::map<dt::book_id_t, std::unique_ptr<bp::book>>> precomputed;

    /**
     * @brief Splits an author in shorter virtual authors
     * 
     * @param shelf_long Contains the books
     * @param autNum     the author to split
     * @return std::map<int,bookprob::book> 
     */
    std::map<dt::book_id_t, std::unique_ptr<bp::book>> authorSplitter(supp::sequence &shelf_long, dt::auth_id_t autNum, dt::book_id_t bookExcl = 0, long F = 0);
    std::map<dt::book_id_t, std::unique_ptr<bp::book>> authorSplitter(supp::sequence &shelf_long, dt::auth_id_t autNum, std::set<dt::book_id_t> bookExcl, long F = 0);
    std::vector<std::string> split(std::string seq, int N);
    size_t get_fullLength(supp::sequence &shelf_long, dt::auth_id_t autNum, std::set<dt::book_id_t> bookExcl);
    dt::book_id_t parts_num(size_t num_books, size_t num_excl, size_t fullLength);
}