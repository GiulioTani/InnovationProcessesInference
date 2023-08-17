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
#include "lib/datatypes.hpp"

/**
 * @brief Contains all the information needed to load the results relative to one text.
 * 
 *  char *fname: base file name for the probabilities.
 *  unsigned short int bucket: number of the bucket or 9999 if using a single file.
 *  long N: number of probabilities associated to the book.
 *  long offs: offset from the beginning of the file.
 *  long ncids: number of author parts.
 *  dt::auth_id_t *cids: ids associated to the authors.
 *  char *goods: flags activating the single author parts.
 *  long fn: number of fragments for the book.
 *  dt::diff_tok_t *lengths: actual length of every fragment of the book.
 *  unsigned int F: set length for the fragments.
 */
struct dataLoad
{
    char *fname; //base file name for the probabilities.
    unsigned short int bucket; //number of the bucket or 9999 if using a single file.
    long N; //number of probabilities associated to the book.
    long offs; //offset from the beginning of the file.
    long ncids; //number of author parts.
    dt::auth_id_t *cids; //ids associated to the authors.
    char *goods; //flags activating the single author parts.
    long fn; //number of fragments for the book.
    dt::diff_tok_t *lengths; //actual length of every fragment of the book.
    unsigned int F; //set length for the fragments.
};
extern "C"
{
    int load(struct dataLoad data, double **results);
    int attribute(double *results, long fn, double delta, char **processed);
    void clean(void *data);
}
