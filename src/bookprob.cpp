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

#include "lib/bookprob.hpp"
#include <stdio.h>
#include <iostream>
#include <unordered_set>
#include <fstream>
#include <filesystem>
#include <math.h>
#include <unistd.h>

#define _P0_NORMALIZATION_ 1 /* -1 = uniform, 0 = fixed, 1 = normalize on author, 2 = normalize on author and fragment*/

namespace bookprob
{
    std::unordered_map<hash_type, std::string> missing_words;
    std::mutex missing_words_lock;
    book::book(const book &oth)
    {
        alpha = oth.alpha;
        theta = oth.theta;
        P0 = oth.P0;
        N = oth.N;
        prefix = oth.prefix;
        order = oth.order;
        myP0 = oth.myP0;
        positions = oth.positions;
        fixWeight = -100;
    };
    book::book(book &&oth)
    {
        alpha = oth.alpha;
        oth.alpha = -100;
        theta = oth.theta;
        oth.theta = -100;
        P0 = oth.P0;
        oth.P0 = nullptr;
        N = oth.N;
        oth.N = 0;
        prefix = std::move(oth.prefix);
        order = std::move(oth.order);
        myP0 = std::move(oth.myP0);
        positions = std::move(oth.positions);
        fixWeight = -100;
    };
    book::book(const std::string &path, const glob_prob *inP0, double Alpha, double Theta, std::string pref) : P0(inP0), prefix(pref)
    {
        initpar(Alpha, Theta);
        if (path == "")
            throw std::invalid_argument("Missing dictionary file name.");
        else if (std::filesystem::path(path).extension() == ".seq")
            read_from_seq(path);
        else if (std::filesystem::path(path).extension() == ".wnt")
            read_from_wnt(path);
        else
            throw std::invalid_argument("Dictionary file type unknown.");
        fixWeight = -100;
    };

    book::book(const std::vector<std::string> &list, const glob_prob *inP0, double Alpha, double Theta, std::string pref) : P0(inP0), prefix(pref)
    {
        initpar(Alpha, Theta);
        read_from_list(list);
        fixWeight = -100;
    };

    book::book(const std::vector<hash_type> &list, const glob_prob *inP0, double Alpha, double Theta, std::string pref) : P0(inP0), prefix(pref)
    {
        initpar(Alpha, Theta);
        read_from_list(list);
        fixWeight = -100;
    };

    book::book(const std::vector<std::pair<hash_type, int>> &ord, const glob_prob *inP0, double Alpha, double Theta, std::string pref) : P0(inP0), prefix(pref)
    {
        N = 0;
        initpar(Alpha, Theta);
        order = ord;
        for (auto i = 0; i < (int)order.size(); i++)
        {
            if (!positions.emplace(order[i].first, i).second)
                throw word_doubled(std::to_string(order[i].first), "order");
            N += order[i].second;
        }
        fixWeight = -100;
    };

    book book::join(const book &other) const
    {
        std::string newpref = prefix + other.prefix;
        std::vector<std::pair<hash_type, int>> neworder;
        neworder = order;
        for (auto &p : other.order)
        {
            if (positions.find(p.first) != positions.end())
                neworder[positions.at(p.first)].second += p.second;
            else
                neworder.push_back(p);
        }
        return book(neworder, P0, alpha, theta, newpref);
    };

    void book::append(const book &other)
    {
        for (auto &p : other.order)
        {
            if (!positions.emplace(p.first, (int)order.size()).second)
                order[positions.at(p.first)].second += p.second;
            else
                order.push_back(p);
            N += p.second;
        }
        fixWeight = -100;
    };

    void book::copy_init(const book &other)
    {
        alpha = other.alpha;
        theta = other.theta;
        P0 = other.P0;
    };

    void book::clear()
    {
        P0 = NULL;
        myP0.prob.clear();
        myP0.count = 0;
        tmpP0.clear();
        order.clear();
        positions.clear();
        alpha = theta = -100;
        prefix.clear();
        fixWeight = -100;
        N = 0;
    };

    void book::add_P0(const glob_prob &inP0)
    {
        myP0 = inP0;
        P0 = &myP0;
    };

    void book::add_par(double Alpha, double Theta)
    {
        if (Alpha >= 0 && Alpha < 1 && Theta > -Alpha)
        {
            alpha = Alpha;
            theta = Theta;
        }
        else
        {
            std::ostringstream err_str;
            err_str << "Wrong values for parameters: alpha=" << Alpha << ", theta=" << Theta << ".";
            throw std::invalid_argument(err_str.str());
        }
    };

    std::pair<double, dt::diff_tok_t> book::log_prob() const
    {
        double P = 0;
        if (positions.size() == 0)
            throw not_init("positions");
        if (order.size() == 0)
            throw not_init("order");
        if (P0 == NULL)
            throw not_init("P0");
        if (alpha == -100)
            throw not_init("alpha");
        if (theta == -100)
            throw not_init("theta");
        if (prefix == "")
            throw not_init("prefix");
        init_tmpP0();
        P += log(alpha) * (int)order.size() + lgamma(theta / alpha + (int)order.size()) - lgamma(theta / alpha);
        P -= lgamma(theta + N) - lgamma(theta);
        for (auto i = 0; i < (int)order.size(); i++)
            P += lgamma(order[i].second - alpha) - lgamma(1 - alpha);
        P += logP_words();
        return {P / log(10), dt::diff_tok_t(tmpP0.size())}; // to have the logarithm in base 10
    };

    std::pair<double, dt::diff_tok_t> book::log_prob(const book &other) const
    {
        double P = 0;
        if (other.order.size() == 0)
            throw not_init("Other order " + other.prefix);
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
        if (other.prefix == "" || prefix == "")
            throw not_init("prefix");
        init_tmpP0(other);
        for (auto &p : other.order)
        {
            if (positions.find(p.first) != positions.end())
                P += lgamma(p.second + order[positions.at(p.first)].second - alpha) - lgamma(order[positions.at(p.first)].second - alpha); // if the word was already in my vocabulary
            else
            {
                if (p.second > 1)
                    P += lgamma(p.second - alpha) - lgamma(1 - alpha); // if the word is new
                try
                {
                    tmpP0.push_back(P0->prob.at(p.first)); // getting counts of all the words in my dictionary
                }
                catch (std::out_of_range &e)
                {
                    std::cerr << prefix << " | " << other.prefix << " | "
                              << " \"" << p.first << "\"" << std::endl;
                    throw word_miss(std::to_string(p.first));
                }
            }
        }
        P += log(alpha) * (int)tmpP0.size() + lgamma(theta / alpha + (int)order.size() + (int)tmpP0.size()) - lgamma(theta / alpha + (int)order.size());
        P -= lgamma(theta + N + other.N) - lgamma(theta + N);
        P += logP_words();
        return {P / log(10), dt::diff_tok_t(tmpP0.size())}; // to have the logarithm in base 10
    };

    void book::init_tmpP0(const book &other) const
    {
        tmpP0 = std::vector<int>();
        tmpP0.reserve(std::size_t(other.order.size()/2));
#if _P0_NORMALIZATION_ > 0
        if (fixWeight < 0)
        {
            fixWeight = P0->count;
            for (auto &p : order)
            {
                try
                {
                    fixWeight -= P0->prob.at(p.first);
                }
                catch (std::out_of_range &e)
                {
                    std::cerr << prefix << " | " << other.prefix << " | " << p.first << std::endl;
                    std::lock_guard<std::mutex> lk(missing_words_lock);
                    missing_words.emplace(p.first, "");
                }
            }
        }
        nowWeight = fixWeight;
#else
        nowWeight = P0->count;
#endif
    };

    void book::init_tmpP0() const
    {
        tmpP0 = std::vector<int>(order.size(), 0);
#if _P0_NORMALIZATION_ > 0
        if (fixWeight < 0)
        {
            fixWeight = P0->count;
            for (auto &p : order)
                fixWeight -= P0->prob.at(p.first);
        }
        nowWeight = fixWeight;
#else
        nowWeight = P0->count;
#endif
        for (auto i = 0; i < (int)order.size(); i++)
        {
            try
            {
                tmpP0[i] = P0->prob.at(order[i].first); // getting counts of all the words in my dictionary
            }
            catch (std::out_of_range &e)
            {
                std::cerr << i << " " << order[i].first << std::endl;
                throw word_miss(std::to_string(order[i].first));
            }
        }
    };

    double book::logP_words() const
    {
#if _P0_NORMALIZATION_ < 0
        return tmpP0.size()*log(1./nowWeight);
        throw std::runtime_error("Not quitted.");
#endif
        double p = 0;
        for (auto molt : tmpP0)
        {
            p += log((double)molt / nowWeight);
#if _P0_NORMALIZATION_ > 1
            nowWeight -= molt;
#endif
        }
        return p;
    };

    size_t book::retsize() const
    {
        return sizeof(double) + sizeof(dt::diff_tok_t);
    };

    void book::log_prob_to_chararr(const book &other, void *dest) const
    {
        auto tmp = log_prob(other);
        memcpy(dest, &tmp.first, sizeof(tmp.first));
        memcpy((char*)dest + sizeof(tmp.first), &tmp.second, sizeof(tmp.second));
    };


    void book::log_prob_to_chararr(void* dest) const
    {
        constexpr std::array<unsigned char, 8> bb{0, 0, 0, 0, 0, 0, 240, 255};
        constexpr dt::diff_tok_t dd=0;
        memcpy(dest, &bb, sizeof(bb));
        memcpy((char*)dest + sizeof(bb), &dd, sizeof(dd));
    }

    void book::initpar(double Alpha, double Theta)
    {
        if (Alpha != -100 && Theta != -100)
            add_par(Alpha, Theta);
        else
        {
            alpha = -100;
            theta = -100;
        }
    };

    book &book::operator=(book &&oth)
    {
        if ((&oth == this))
            return *this;

        alpha = oth.alpha;
        oth.alpha = -100;
        theta = oth.theta;
        oth.theta = -100;
        N = oth.N;
        oth.N = 0;
        P0 = oth.P0;
        oth.P0 = nullptr;
        fixWeight = oth.fixWeight;
        oth.fixWeight = -100;
        prefix = std::move(oth.prefix);
        order = std::move(oth.order);
        myP0 = std::move(oth.myP0);
        positions = std::move(oth.positions);
        return *this;
    };
}

/* INPUT FROM FILE */
namespace bookprob
{
    void book::read_from_wnt(std::string file)
    {
        std::ifstream fpinput;
        std::string line, word;
        std::unordered_set<std::string> seen;
        int nk;
        double tk;
        std::stringstream ss;
        N = 0;
        fpinput.open(file);
        if (fpinput.fail())
        {
            fprintf(stderr, "Dictionary file  %s doesn't exist ... exiting\n", file.data());
            throw std::runtime_error("File not found.");
        }
        while (getline(fpinput, line))
        {
            if (line[0] == '#')
            {
                if (prefix == "")
                    prefix = line.substr(1, 100);
                continue;
            }
            ss.clear();
            ss.str(line);
            ss >> word >> nk >> tk;
            if (seen.find(word) != seen.end())
                throw word_doubled(word, "dictionary");
            else
            {
                int temp = std::hash<std::string>{}(word);
                if (!positions.emplace(std::make_pair(temp, order.size())).second)
                    throw word_doubled(std::to_string(temp), "hash");
                order.push_back(std::make_pair(temp, nk));
                N += nk;
            }
        }
        if (fpinput.is_open())
            fpinput.close();
        if (!order.size())
            throw empty_file(file);
    };

    void book::read_from_seq(std::string file)
    {
        std::ifstream fpinput;
        std::string word;
        std::stringstream ss;
        std::vector<std::string> list;

        fpinput.open(file);
        if (fpinput.fail())
        {
            fprintf(stderr, "Sequence file  %s doesn't exist ... exiting\n", file.data());
            exit(EXIT_FAILURE);
        }
        if (isatty(0))
            std::cerr << "Reading sequence..." << std::endl;
        while (getline(fpinput, word))
        {
            if (word == "")
                break;
            if (word[0] == '#')
            {
                if (prefix == "")
                    prefix = word.substr(1, 100);
                else
                    std::cout << word << std::endl;
                continue;
            }
            list.push_back(word);
        }
        if (fpinput.is_open())
            fpinput.close();
        read_from_list(list);
    }

    void book::read_from_list(std::vector<std::string> list)
    {
        std::vector<hash_type> hash_list;
        std::hash<std::string> hasher;
        hash_list.reserve(list.size());
        if (!list.size())
            throw empty_file("sequence");
        for (auto &word : list)
        {
            hash_type temp = hasher(word);
            if (!missing_words.empty())
            {
                std::lock_guard<std::mutex> lk(missing_words_lock);
                if (missing_words.find(temp) != missing_words.end())
                    missing_words.at(temp) = word;
            }
            hash_list.push_back(temp);
        }
        read_from_list(hash_list);
    }

    void book::read_from_list(std::vector<hash_type> list)
    {
        if (!list.size())
            throw empty_file("sequence");
        for (auto &word : list)
        {
            if (positions.emplace(std::make_pair(word, order.size())).second)
                order.push_back(std::make_pair(word, 1));
            else
                order[positions[word]].second++;
        }
        N = list.size();
    }
}