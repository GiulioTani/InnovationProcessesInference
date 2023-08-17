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

#include "lib/authorSplitter.hpp"
#include <random>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#define _FLEXIBLE_AUTHOR_SIZE_
namespace asp
{
    dt::book_id_t authCount = 0 /**< Used to count the number of virtual authors defined.*/;
    int authorDimension = -2000000000, ngramSize = -2000000000;
    std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
    std::mutex autnum_lock, orders_lock, exclusion_lock, precomputed_lock;
    std::unordered_map<dt::auth_id_t, std::vector<dt::book_id_t>> orders;
    std::unordered_map<dt::auth_id_t, std::map<dt::book_id_t, std::set<dt::book_id_t>>> excluded;
    std::map<dt::auth_id_t, std::map<dt::book_id_t, std::unique_ptr<bp::book>>> precomputed;

    std::vector<std::string> split(std::string seq)
    {
        std::vector<std::string> ret;
        ret.reserve(seq.size() / 5);
        int j = 0;
        for (j = 0; seq[j] == ' '; j++)
            ;
        for (auto i = j; i<(int)seq.size(); i = j> 0 ? j + 1 : seq.size())
        {
            j = seq.find(" ", i);
            auto tmp = seq.substr(i, j - i);
            if (tmp.size() > 0)
                ret.push_back(tmp);
        }
        if (ret.size() == 0)
            ret.push_back("");
        return ret;
    }

    std::vector<std::string> split(std::string seq, int N)
    {
        if (!N)
            return split(seq);
        std::vector<std::string> ret;
        ret.reserve(seq.size() / 10);
        int j = 0;
        for (j = 0; seq[j] == ' ' && seq[j + 1] == ' '; j++)
            ;
        for (auto i = j; i < (int)seq.size() - N + 1; i++)
        {
            j = seq.find(' ', i + 1);
            if (j < 0 || j > i + N - 2)
            {
                ret.push_back(seq.substr(i, N));
                if (ret.back().back() == ' ')
                    ret.back().back() = '_';
                if (ret.back().front() == ' ')
                    ret.back().front() = '_';
            }
        }
        if (ret.size() == 0)
            ret.push_back("");
        return ret;
    }

    unsigned get_delta()
    {
        return authorDimension > 10000 ? 1000 : authorDimension / 10 /**Tolerance on slice size.*/;
    }

    size_t get_fullLength(supp::sequence &shelf_long, dt::auth_id_t autNum, std::set<dt::book_id_t> bookExcl)
    {
        size_t fullLength = 0 /**Total equivalent size of the author.*/;
        for (auto &book : shelf_long.at(autNum))
            if (bookExcl.find(book.first) == bookExcl.end())
                fullLength += book.second.size();
        return fullLength;
    }
    dt::book_id_t parts_num(size_t num_books, size_t num_excl, size_t fullLength)
    {
        unsigned delta = get_delta();
        if (authorDimension < 0)
            throw std::runtime_error("Author Dimension not initialised.");
        if (fullLength <= num_books - num_excl)
        {
            return 0;
        }
        if (fullLength <= unsigned(authorDimension) + delta || authorDimension == 0)
        {
            return 1;
        }
        else
        {
            // determining the number of parts to split into
            dt::book_id_t parts = fullLength / authorDimension;
#ifdef _FLEXIBLE_AUTHOR_SIZE_
            // choses between more smaller parts or fewer bigger to avoid wasting data
            if (fullLength > unsigned(authorDimension) * std::sqrt(parts * (parts + 1)))
                parts++;
#endif
            return parts;
        }
    }
    std::map<dt::book_id_t, std::unique_ptr<bp::book>> authorSplitter(supp::sequence &shelf_long, dt::auth_id_t autNum, dt::book_id_t bookExcl, long F)
    {
        return authorSplitter(shelf_long, autNum, std::set<dt::book_id_t>({bookExcl}), F);
    }
    std::map<dt::book_id_t, std::unique_ptr<bp::book>> authorSplitter(supp::sequence &shelf_long, dt::auth_id_t autNum, std::set<dt::book_id_t> bookExcl, long F)
    {
        size_t fullLength = 0 /**Total equivalent size of the author.*/;
        double slicesize /**Actual size of the slice, may be different from authorDimension*/, virtualSize = 0 /**Current equivalent size of the slice part.*/;
        dt::book_id_t Anum, parts /**Number of parts.*/;
        unsigned delta = get_delta();
        dt::book_id_t nextbook = 0 /**Index of the next book in order.*/;
        int nexttoken = 0 /**Position of the next character of the present book.*/;
        std::vector<bookprob::hash_type> newAut;
        std::vector<std::unique_ptr<bp::book>> proposedAuthors;
        std::map<dt::book_id_t, std::unique_ptr<bp::book>> returnedAuthors;

        if (!bookExcl.size())
        {
            std::lock_guard<std::mutex> prec_guard(precomputed_lock);
            if (precomputed.find(autNum) != precomputed.end())
            {
                for (auto &pc : precomputed[autNum])
                {
                    if (F == 1)
                        returnedAuthors.emplace(pc.first, std::make_unique<bp::dkl_book>(*pc.second));
                    else
                        returnedAuthors.emplace(pc.first, std::make_unique<bp::book>(*pc.second));
                }
                return returnedAuthors;
            }
        }
        fullLength = get_fullLength(shelf_long, autNum, bookExcl);
        parts = parts_num(shelf_long.at(autNum).size(), bookExcl.size(), fullLength);

        if (parts == 0)
        {
            returnedAuthors.clear();
            return returnedAuthors;
        }
        if (parts == 1)
        {
            newAut.clear();
            for (auto &book : shelf_long.at(autNum))
            {
                if (bookExcl.find(book.first) != bookExcl.end())
                    continue;
                newAut.insert(newAut.end(), book.second.begin(), book.second.end());
            }
            if (F == 1)
                proposedAuthors.push_back(std::make_unique<bp::dkl_book>(bp::dkl_book(newAut)));
            else
                proposedAuthors.push_back(std::make_unique<bp::book>(bp::book(newAut)));
        }
        else
        {
            std::vector<dt::book_id_t> order;
            orders_lock.lock();
            if (orders.find(autNum) != orders.end())
                order = orders[autNum];
            else
            { // shuffling the order of the books
                orders_lock.unlock();
                for (auto &book : shelf_long.at(autNum))
                    order.push_back(book.first);
                std::shuffle(order.begin(), order.end(), generator);
                orders_lock.lock();
                orders[autNum] = order;
            }
            orders_lock.unlock();
#ifdef _FLEXIBLE_AUTHOR_SIZE_
            slicesize = fullLength / parts;
#else
            slicesize = authorDimension;
#endif
            // split into parts
            for (dt::book_id_t i = 0; i < parts; i++)
            {
                for (newAut.clear(), newAut.reserve(slicesize + delta); virtualSize < slicesize - delta && nextbook < order.size(); nextbook++)
                {
                    while (nextbook < order.size() && bookExcl.find(order[nextbook]) != bookExcl.end())
                        nextbook++; // skipping excluded books
                    if (nextbook >= order.size())
                        break;
                    if (virtualSize + (shelf_long.at(autNum).at(order[nextbook]).size() - nexttoken) <= slicesize + delta)
                    { // next book fits in
                        newAut.insert(newAut.end(), shelf_long.at(autNum).at(order[nextbook]).begin() + nexttoken, shelf_long.at(autNum).at(order[nextbook]).end());
                        virtualSize += (shelf_long.at(autNum).at(order[nextbook]).size() - nexttoken);
                        nexttoken = 0;
                    }
                    else
                    {                                                // slice extra long book
                        int missingPart = (slicesize - virtualSize); /**True size of the missing part*/
                        ;
                        newAut.insert(newAut.end(), shelf_long.at(autNum).at(order[nextbook]).begin() + nexttoken, shelf_long.at(autNum).at(order[nextbook]).begin() + nexttoken + missingPart);
                        nexttoken += missingPart;
                        break;
                    }
                }
                if (!newAut.size())
                    continue;
                if (F == 1)
                    proposedAuthors.push_back(std::make_unique<bp::dkl_book>(bp::dkl_book(newAut)));
                else
                    proposedAuthors.push_back(std::make_unique<bp::book>(bp::book(newAut)));
                virtualSize = 0;
            }
        }
        Anum = 0;
        exclusion_lock.lock();
        if (excluded.find(autNum) != excluded.end())
            for (auto ex : excluded.at(autNum))
                if (ex.second == bookExcl)
                    Anum = ex.first;
        exclusion_lock.unlock();
        if (!Anum)
        {
            autnum_lock.lock();
            Anum = authCount + 1;
            authCount += proposedAuthors.size();
            autnum_lock.unlock();
            exclusion_lock.lock();
            excluded[autNum][Anum] = bookExcl;
            exclusion_lock.unlock();
        }
        for (auto &&au : proposedAuthors)
        {
            au->set_prefix("A" + std::to_string(Anum) + "_");
            returnedAuthors.emplace(Anum++, std::move(au));
        }
        if (!bookExcl.size())
        {
            std::lock_guard<std::mutex> prec_guard(precomputed_lock);
            for (auto &pc : returnedAuthors)
            {
                if (F == 1)
                    precomputed[autNum].emplace(pc.first, std::make_unique<bp::dkl_book>(*pc.second));
                else
                    precomputed[autNum].emplace(pc.first, std::make_unique<bp::book>(*pc.second));
            }
        }
        return returnedAuthors;
    }
}