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
 * @file supportlib.cpp
 * @author Giulio Tani Raffaelli (tani@cs.cas.cz)
 * @brief Functions' bodies for supportlib and globals definition.
 * @version 1.0
 * @date 2023-06-27
 * 
 * @copyright Copyright (c) 2023
 *
 */
#include "lib/supportlib.hpp"
#include <queue>
#include <algorithm>
#include <csignal>
#include <iostream>
#include <thread>
#include <unistd.h>
#include <fstream>
#include <algorithm>
#include <random>
#include <iomanip>
#include "lib/authorSplitter.hpp"
#include <climits>
#include <functional>

#include <sys/ioctl.h> //ioctl() and TIOCGWINSZ
#include <unistd.h>    // for STDOUT_FILENO
/* Definition of the list_manager methods */
#define __WORKERS_THRESHOLD__ 9999 /**< Number of workers that marks an author as completed. */
#define MAX_BUCKET_NUM 400
namespace bp = bookprob;
namespace supp
{
    bool FLAG = true;  /**< Activates the controlled shutdown when the first SIGINT arrives. */
    bool NOHUP = true; /**< Makes the signal handler ignore the first SIGHUP.*/
    std::mutex shelf_lock, queue_lock;
    task_iterator::task_iterator(dt::auth_id_t aid, std::vector<std::reference_wrapper<const supp::slice>> cake_part, size_t size) : my_auth(aid), cakelist(cake_part), _size(size)
    {
        if (!_size)
            calculate_size();
        now_slice = cakelist.begin();
        now_aut = now_slice->get().shelf.begin();
        now_book = now_aut->second.begin();
    }

    void task_iterator::calculate_size()
    {
        _size = 0;
        for (auto i : cakelist)
            for (auto j : i.get().shelf)
                _size += j.second.size();
    }

    task task_iterator::pop()
    {
        task ret;
        if (now_slice == cakelist.end())
        {
            if (returned.empty())
            {
                throw std::runtime_error("Called empty task_iterator: " + std::to_string(my_auth));
            }
            else
            {
                auto toret = returned.front();
                returned.pop();
                _size--;
                return {my_auth, toret.first, toret.second};
            }
        }
        else
        {
            ret = {my_auth, now_aut->first, *now_book};
            if (++now_book == now_aut->second.end())
            {
                if (++now_aut != now_slice->get().shelf.end())
                    now_book = now_aut->second.begin();
                else if (++now_slice != cakelist.end())
                {
                    now_aut = now_slice->get().shelf.begin();
                    now_book = now_aut->second.begin();
                }
            }
        }
        _size--;
        return ret;
    }

    void task_iterator::place_back(task job)
    {
        returned.emplace(std::make_pair(job.aut2, job.book));
        _size++;
    }

    list_manager::list_manager(const library &shelf, const std::vector<supp::slice> &cake, const std::vector<int> &slices)
    {
        create(shelf, cake, slices);
    }
    void list_manager::create(const library &shelf, const std::vector<supp::slice> &cake, const std::vector<int> &slices)
    {
        size_t n = 0;
        std::vector<std::pair<dt::auth_id_t, dt::book_id_t>> tempV, tempV2, jobs;
        std::vector<std::reference_wrapper<const supp::slice>> cake_part;
        if (slices.size() > 0)
            for (auto i : slices)
                cake_part.push_back(cake[i]);
        else
            for (auto& slice : cake)
                cake_part.push_back(slice);

        if (created)
            throw std::logic_error("Create list managers only once");
        else
        {
            for (auto my_slice : cake_part)
                for (auto &aut1 : my_slice.get().shelf)
                    n += aut1.second.size();
            taskers.clear();
            std::cout << "Scheduling tasks\r" << std::flush;
            for (auto &aut1 : shelf)
                if (aut1.first)
                    taskers.emplace(aut1.first, task_iterator(aut1.first, cake_part, n));

            Priority = std::forward_list<dt::auth_id_t>(taskers.size());
            auto autP = Priority.begin();
            for (auto autQ = taskers.begin(); autQ != taskers.end(); autQ++, autP++)
                *autP = autQ->first;
            Priority.reverse();
            actual = Priority.begin();
            first_round = true;
            std::cout << "Scheduled " << n * taskers.size() << " tasks." << std::endl;
            created = true;
        }
    };

    void list_manager::place_back(task job)
    {
        std::lock_guard<std::mutex> guard(lock);
        taskers.at(job.aut1).place_back(job);
        placed_back.emplace(job.aut1);
    }

    task list_manager::pop(dt::auth_id_t my_aut)
    {
        dt::auth_id_t new_aut;
        std::lock_guard<std::mutex> guard(lock);

        if (!my_aut || !taskers[my_aut].size())
        { // if not assigned do it
            if (first_round)
            {
                new_aut = *actual;
                clean_advance();
            }
            else
            {
                while (!Priority.empty() && !taskers[*actual].size())
                    clean_advance();
                if (Priority.empty())
                { // controllo i restituiti
                    if (!placed_back.empty())
                    {
                        for (auto i : placed_back)
                            Priority.push_front(i);
                        placed_back.clear();
                        Priority.remove_if([this](dt::auth_id_t i)
                                           { return taskers[i].size() == 0; });
                        actual = Priority.begin();
                    }
                    if (Priority.empty())
                    {
                        return {0, 0, 0};
                    }
                    else
                    {
                        new_aut = *actual;
                        clean_advance();
                    }
                }
                else
                {
                    new_aut = *actual;
                    clean_advance();
                }
            } // even the highest priority author is completed
        }
        else
            new_aut = my_aut;

        return taskers[new_aut].pop();
    };
    void list_manager::clean_advance()
    {
        if (++actual == Priority.end())
        {
            first_round = false;
            Priority.remove_if([this](dt::auth_id_t i)
                               { return taskers[i].size() == 0; });
            actual = Priority.begin();
        }
    };

    size_t list_manager::remaining()
    {
        size_t n = 0;
        std::lock_guard<std::mutex> guard(lock);
        for (auto &q : taskers)
            n += q.second.size();
        return n;
    };

    void sig_handler(int signo)
    {
        if (signo == SIGUSR1)
        {
            sleep(1e9); // The thread is paused putting it to sleep for 31 years and 9 months
        }
        if (signo == SIGINT)
        {
            fprintf(stderr, "\nReceived SIGINT\n");
            if (FLAG)
                FLAG = false;
            else
                throw user_stop("Received SIGINT.");
        }
        if (signo == SIGHUP)
        {
            if (NOHUP)
                fprintf(stderr, "\nIgnoring SIGHUP\n");
            else
                exit(SIGHUP);
        }
        if (signo == SIGCONT)
        {
        }
    }

    void load_P0(std::string fileP0, bp::glob_prob &P0)
    {
        std::ifstream inputfiles;
        P0.count = 0;
        P0.prob.clear();
        if (!std::filesystem::is_regular_file(fileP0))
            throw std::invalid_argument(fileP0 + " is missing in input folder or not a regular file.");
        if (std::filesystem::path(fileP0).extension() == ".bin")
        {
            inputfiles.open((fileP0).c_str(), std::ifstream::binary);
            if (inputfiles)
            {
                std::pair<const bookprob::hash_type, int> hash_count;
                while (!inputfiles.eof())
                {
                    inputfiles.read(reinterpret_cast<char *>(&hash_count), sizeof(hash_count));
                    P0.prob.emplace(hash_count);
                    P0.count += hash_count.second;
                }
                inputfiles.close();
            }
            else
                throw std::runtime_error("Could not open P0 input file");
        }
        else
        {
            std::stringstream ss;
            std::string line, word;
            std::unordered_set<std::string> seen;
            int molt;
            inputfiles.open((fileP0).c_str());
            if (inputfiles)
            {
                while (getline(inputfiles, line))
                {
                    ss.clear();
                    ss.str(line);
                    ss >> word >> molt;
                    if (seen.emplace(word).second)
                    {
                        bp::hash_type temp = std::hash<std::string>{}(word);
                        if (!P0.prob.emplace(temp, molt).second)
                            throw bp::word_doubled(word, "hashing");
                        P0.count += molt;
                    }
                    else
                        throw bp::word_doubled(word, "probability");
                }
                inputfiles.close();
            }
            else
                throw std::runtime_error("Could not open P0 input file");
        }
        std::cout << "Loaded P0" << std::endl;
    }

    void hash_combine(size_t & seed, size_t const& v){
        seed ^= v + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }

    void dump_P0(std::filesystem::path outputFolderPath, bp::glob_prob &P0)
    {
        std::ofstream output((outputFolderPath / "baseprobability.bin").c_str(), std::ofstream::binary);
        if (output)
        {
            for (auto p : P0.prob)
                output.write(reinterpret_cast<const char *>(&p), sizeof(p));
            output.close();
        }
        else
            throw std::runtime_error("Could not open P0 output file");
    }

    int slice::size()
    {
        int s = 0;
        for (auto v : shelf)
            s += (int)v.second.size();
        return s;
    }
    void cake_from_library(std::vector<supp::slice> *cake, const library &shelf, int slicesize)
    {
        int numbooks = 0, numslices = 0;
        std::vector<std::pair<dt::auth_id_t, dt::book_id_t>> well, authwell;
        std::vector<dt::auth_id_t> randaut;
        std::random_device rd;
        auto rng = std::default_random_engine{rd()};
        randaut.reserve(shelf.size());
        for (auto &a : shelf)
            randaut.push_back(a.first);
        std::shuffle(std::begin(randaut), std::end(randaut), rng);
        for (auto &a : randaut)
        {
            numbooks += shelf.at(a).size();
            authwell.clear();
            authwell.reserve(shelf.at(a).size());
            for (auto &b : shelf.at(a))
                authwell.push_back({a, b.first});
            std::shuffle(std::begin(authwell), std::end(authwell), rng);
            well.insert(std::end(well), std::begin(authwell), std::end(authwell));
        }
        numslices = (numbooks / slicesize);
        if (numbooks > slicesize * sqrt(double(numslices) * (numslices + 1)))
            numslices++;
        if (numslices < 0 || size_t(numslices) > std::numeric_limits<dt::slice_id_t>::max())
            throw std::invalid_argument("The number of slices is larger than the mazimum number allowed"
                                        "\nRecompile with a larger type for slice_id_t or change slice size. " +
                                        std::to_string(numslices));
        std::cout << "Cake: " << numslices << " slices of size " << slicesize << " for " << numbooks << " books." << std::endl;

        for (auto i = 0; i < numslices; i++)
            cake->push_back(supp::slice());

        for (auto i = 0; i < (int)well.size(); i++)
            if (!(cake->at(i % numslices).shelf.emplace(well[i].first, std::set<dt::book_id_t>({well[i].second})).second))
                cake->at(i % numslices).shelf[well[i].first].emplace(well[i].second); // create slicing
    }

    bp::book split_frag(std::vector<bp::hash_type> &splitBook, int A, int B, int F, int start, int end, bp::glob_prob &P0)
    {
        std::ostringstream fragname;
        std::vector<bp::hash_type> fragment;
        fragment.clear();
        std::copy(splitBook.begin() + start, splitBook.begin() + end, std::back_inserter(fragment));

        fragname.clear();
        fragname.str("");
        fragname << "A" << A << "B" << B << "F" << F;
        return bp::book(fragment, &P0, -100, -100, fragname.str());
    }

    int read_books_from_file(library &shelf_short, sequence &shelf_long, int F, std::queue<std::filesystem::path> &inputFiles, bp::glob_prob &P0)
    {
        long long int fn, tfs, res, add, totFrag = 0, B = 0;
        std::ifstream inwnt, inseq;
        std::string line, word, word2, filename;
        std::ostringstream fragname;
        std::istringstream ss;
        std::vector<std::pair<bp::hash_type, int>> newfrag;
        std::vector<std::string> tmp_splitBook;
        std::vector<bp::hash_type> splitBook;
        std::unordered_map<dt::book_id_t, std::vector<bp::hash_type>> tmp_long;
        std::unordered_map<dt::book_id_t, std::vector<std::unique_ptr<bp::book>>> tmp_short;
        std::filesystem::path inputFile;
        std::hash<std::string> hasher;
        while (!inputFiles.empty())
        {
            try
            {
                std::lock_guard<std::mutex> sl(queue_lock);
                if (!inputFiles.empty())
                {
                    inputFile = inputFiles.front();
                    inputFiles.pop();
                }
                else
                    break;
            }
            catch (const std::exception &e)
            {
                std::cerr << e.what() << '\n';
            }
            int A = std::strtol(inputFile.stem().string().substr(1).c_str(), NULL, 10);
            inseq.open(inputFile.string());
            if (inseq.is_open())
            {
                while (getline(inseq, line))
                {
                    if (line[0] != '#' || !isInteger(line.substr(2)))
                        throw std::invalid_argument("Odd lines of seq file must be \"# <book number>\", found: " + line.substr(0, 10));
                    B = std::stoi(line.substr(1));
                    if (B < 0 || size_t(B) > std::numeric_limits<dt::book_id_t>::max())
                        throw std::invalid_argument("The number of books for some author is larger than the maximum number allowed"
                                                    "\nRecompile with a larger type for booke_id_t (if impossible, split the large authors). " +
                                                    std::to_string(B));
                    if (!getline(inseq, line))
                        throw std::invalid_argument("Expected a line in seq containing the processed book.");
                    tmp_splitBook.clear();
                    splitBook.clear();
                    tmp_splitBook = asp::split(" " + line + " ", asp::ngramSize);
                    splitBook.reserve(tmp_splitBook.size());
                    for (auto &w : tmp_splitBook)
                        splitBook.push_back(hasher(w));
                    tmp_long[B] = splitBook;
                    if (splitBook.empty())
                    {
                        fragname.clear();
                        fragname.str("");
                        fragname << "A" << A << "B" << B << "F1";
                        tmp_short[B].push_back(std::make_unique<bp::book>(bp::book(std::vector<std::pair<bp::hash_type, int>>{std::make_pair(std::hash<std::string>{}(""), 0)}, &P0, -100, -100, fragname.str())));
                        continue;
                    }
                    if (F == 1)
                    {
                        fragname.clear();
                        fragname.str("");
                        fragname << "A" << A << "B" << B << "F1";
                        tmp_short[B].push_back(std::make_unique<bp::dkl_book>(bp::dkl_book(splitBook, &P0, -100, -100, fragname.str())));
                    }
                    else
                    {
                        if (F && (int)splitBook.size() > F)
                        {
                            fn = splitBook.size() / F;
                            // choses between more smaller parts or fewer bigger to avoid wasting data
                            if (splitBook.size() > F * (fn + 0.5))
                                fn++;
                            if (fn < 0 || size_t(fn) > std::numeric_limits<dt::frag_id_t>::max())
                                throw std::invalid_argument("The number of fragments in some book is larger than the maximum number allowed"
                                                            "\nRecompile with a larger type for frag_id_t or change fragment size. " +
                                                            std::to_string(fn));
                            res = splitBook.size() - F * fn;
                            if (res >= 0)
                                add = +1;
                            else
                                add = -1;
                            tfs = F + res / fn;
                            if (tfs > std::numeric_limits<dt::diff_tok_t>::max())
                                throw std::invalid_argument("The number of tokens in some fragment is larger than the maximum number allowed"
                                                            "\nRecompile with a larger type for diff_tok_t or change fragment size. " +
                                                            std::to_string(tfs));
                            res %= fn;
                        }
                        else
                        {
                            fn = 1;
                            tfs = splitBook.size();
                            res = 0;
                            add = 0;
                            if (tfs > std::numeric_limits<dt::diff_tok_t>::max())
                                throw std::invalid_argument("The number of tokens in some fragment is larger than the mazimum number allowed"
                                                            "\nRecompile with a larger type for diff_tok_t or change fragment size. " +
                                                            std::to_string(tfs));
                        }
                        tmp_short[B].reserve(fn);
                        for (auto f = 0; f < abs(res); f++)
                        {
                            auto end = ((tfs + add) * (f + 1) < (long int)splitBook.size() ? (tfs + add) * (f + 1) : splitBook.size());
                            auto start = f * (tfs + add);
                            tmp_short[B].push_back(std::make_unique<bp::book>(split_frag(splitBook, A, B, f + 1, start, end, P0)));
                        }
                        for (auto f = abs(res); f < fn; f++)
                        {
                            auto end = (tfs * (f + 1) + res < (long int)splitBook.size() ? tfs * (f + 1) + res : splitBook.size());
                            auto start = f * tfs + res;
                            tmp_short[B].push_back(std::make_unique<bp::book>(split_frag(splitBook, A, B, f + 1, start, end, P0)));
                        }
                        totFrag += fn;
                    }
                }
                inseq.close();
            }
            else
                throw std::runtime_error("Can't open " + inputFile.string());
            try
            {
                std::lock_guard<std::mutex> sl(shelf_lock);
                shelf_short[A] = std::move(tmp_short);
                shelf_long[A] = std::move(tmp_long);
            }
            catch (const std::exception &e)
            {
                std::cerr << e.what() << '\n';
            }
        }

        return totFrag;
    }

    std::string nicetime(int t)
    {
        if (t < 0 || t + 1 < 0)
            return " inf ";
        else
        {
            std::ostringstream ss;
            ss.str("");
            ss.clear();
            ss << std::setw(2) << std::setfill('0') << t / 60 << ":" << std::setw(2) << t % 60;
            return ss.str();
        }
    }

    void progress_bar(double tpass, size_t nowwri, size_t totwri, size_t remaining, size_t resCount, double eLAPsed, size_t totProbs, size_t totFra)
    {
        std::ostringstream ss, s2;
        struct winsize sizeWin;
        int hashnum, lastdgt, barlen, trem;
        bool print;
        static int calls = 0;
        try
        {
            ss.str("");
            ss.clear();
            s2.str("");
            s2.clear();
            trem = 1. * tpass / nowwri * (totProbs - totwri);
            ss << std::setw(6) << std::fixed << std::setprecision(2) << 100. * totwri / totProbs << "% ";
            s2 << "| " << std::scientific << std::setprecision(3) << 1. * totFra << " comp [" << nicetime(tpass) << "<" << nicetime(trem) << "(" << std::setfill(' ') << std::setw(int(log10(totProbs)) + 1) << remaining << "), " << std::fixed << std::setprecision(2) << std::setw(6) << resCount / eLAPsed / 1000 << " K comp/s]";
            if (isatty(STDOUT_FILENO))
            {
                ioctl(STDOUT_FILENO, TIOCGWINSZ, &sizeWin);
                print = true;
            }
            else
            {
                sizeWin.ws_col = 90;
                print = !(calls % 40);
            }
            barlen = sizeWin.ws_col - ss.tellp() - s2.tellp() - 3;
            if (barlen > 0)
            {
                ss << "|";
            }
            else
                barlen = 0;
            hashnum = 1. * barlen * totwri / totProbs;
            lastdgt = ((1. * barlen * totwri / totProbs) - hashnum) * 10;
            ss << std::string(hashnum, '#');
            if (lastdgt)
                ss << lastdgt << std::string(barlen - hashnum - 1, ' ');
            else
                ss << std::string(barlen - hashnum, ' ');
            if (print)
                std::cout << "\r" << ss.str() + s2.str() << std::flush;
            calls++;
        }
        catch (const std::exception &e)
        {
            std::cerr << "\r Progress bar error: " << e.what() << '\n';
            std::cout << "\r Progress bar error: " << e.what() << '\n';
        }
    }

    unsigned short int split_save_number(library &shelf_short, int numThreads)
    {
        unsigned short int numb_of_buck;
        int tot_book = 0;
        std::unordered_map<dt::auth_id_t, unsigned short int> author_bucket;
        for (auto &i : shelf_short)
            tot_book += i.second.size();
        numb_of_buck = sqrt(tot_book) / 100 + 1;
        numb_of_buck *= numThreads;
        numb_of_buck = numb_of_buck < shelf_short.size() ? numb_of_buck : shelf_short.size();
        numb_of_buck = numb_of_buck < MAX_BUCKET_NUM ? numb_of_buck : MAX_BUCKET_NUM;

        return numb_of_buck;
    }

    void load_main_author_parts(std::string fileName, std::map<dt::auth_id_t, dt::book_id_t> *main_author_parts, std::map<dt::auth_id_t, std::unordered_map<dt::slice_id_t, dt::book_id_t>> *slice_author_size)
    {
        std::ifstream input;
        main_author_parts->clear();
        slice_author_size->clear();

        input.open(fileName, std::ifstream::binary);
        if (input.fail() || !input.is_open())
        {
            std::cerr << "Error opening: " << fileName << std::endl;
            exit(EXIT_FAILURE);
        }
        if (input)
        {
            dt::auth_id_t N_auth, A;
            dt::book_id_t N_books, default_size;
            std::pair<dt::slice_id_t, dt::book_id_t> S_fa;
            input.read((char *)&N_auth, sizeof(N_auth));
            for (unsigned i = 0; i < N_auth; i++)
            {
                input.read((char *)&A, sizeof(A));
                slice_author_size->emplace(A, std::unordered_map<dt::slice_id_t, dt::book_id_t>());
                input.read((char *)&N_books, sizeof(dt::book_id_t));
                input.read((char *)&default_size, sizeof(dt::book_id_t));
                if (default_size)
                    main_author_parts->emplace(A, default_size);
                for (unsigned j = 0; j < N_books; j++)
                {
                    input.read((char *)&S_fa, sizeof(S_fa));
                    slice_author_size->at(A).emplace(S_fa);
                }
            }
            input.close();
        }
    }
    void catch_signals()
    {
        if (signal(SIGCONT, supp::sig_handler) == SIG_ERR)
            fprintf(stderr, "\nCan't catch SIGCONT\n");
        if (signal(SIGUSR1, supp::sig_handler) == SIG_ERR)
            fprintf(stderr, "\nCan't catch SIGUSR1\n"); // don't know why I can't effectively catch SIGSTOP, stops everything
        if (signal(SIGINT, supp::sig_handler) == SIG_ERR)
            fprintf(stderr, "\nCan't catch SIGINT\n");
        if (signal(SIGHUP, supp::sig_handler) == SIG_ERR)
            fprintf(stderr, "\nCan't catch SIGHUP\n");
    }
    std::string insert_padded_number(std::string first_part, std::string last_part, int n, int padsize)
    {
        std::stringstream padnum;
        padnum.clear();
        padnum.str("");
        padnum << first_part << std::setfill('0') << std::setw(padsize) << n << last_part;
        return padnum.str();
    }
} // namespace supp
