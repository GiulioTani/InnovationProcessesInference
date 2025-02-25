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
#include <stdint-gcc.h>
#include <vector>
#include <unordered_map>
#include <memory>
namespace dt
{
    typedef uint32_t diff_tok_t;
    typedef uint16_t auth_id_t;
    typedef uint32_t book_id_t;
    typedef uint32_t frag_id_t;
    typedef uint32_t slice_id_t;
    typedef size_t hash_type;

    typedef std::unordered_map<auth_id_t, std::unordered_map<book_id_t, std::vector<hash_type>>> sequence;

    /**
     * @struct TASK supportlib.hpp
     * @brief Contains information on the book to compute.
     *
     */
    struct TASK
    {
        auth_id_t aut1 /** Author to compare with. */, aut2 /** Author of the book to compare. */;
        book_id_t book /** Number of the book to compare */;
        bool operator==(struct TASK oth) { return aut1 == oth.aut1 && aut2 == oth.aut2 && book == oth.book; };
    };

    /**
     * @typedef typedef struct TASK task
     * @brief Defines as a type the struct TASK.
     *
     */
    typedef struct TASK task;

    struct FRAG_LABEL
    {
        task task_id;
        book_id_t aut_frag;
        frag_id_t frag;
        bool operator==(struct FRAG_LABEL oth) { return task_id == oth.task_id && aut_frag == oth.aut_frag && frag == oth.frag; };
    };

    typedef struct FRAG_LABEL frag_label;

    typedef std::vector<std::pair<hash_type, std::pair<double, uint64_t>>> contrib_type;
}