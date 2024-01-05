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
 * @file supportlib.hpp
 * @author Giulio Tani Raffaelli (tani@cs.cas.cz)
 * @brief This is a support library with job managing and signal handling functionalities.
 * @version 1.0
 * @date 2023-06-27
 * 
 * @copyright Copyright (c) 2023
 *
 */
#pragma once
#include "lib/bookprob.hpp" //defines namespace book
#include "lib/bookdkl.hpp"
#include "lib/datatypes.hpp"
#include <vector>
#include <set>
#include <unordered_set>
#include <mutex>
#include <queue>
#include <unordered_map>
#include <numeric>   // std::iota
#include <algorithm> // std::sort, std::stable_sort
#include <filesystem>
#include <chrono>
#include <condition_variable>
#include <forward_list>
#include <csignal>
#include <map>
#include <functional>

namespace bp = bookprob;
/**
 * @package supp
 * @brief Contains useful classes and functions.
 *
 * Although tailored for the use with the objects defined in bookprob, the functions in this namespace
 * are not really specific to books but of more general purpose.
 *
 */
namespace supp
{
    extern bool FLAG;  /**< Activates the controlled shutdown when the first SIGINT arrives. */
    extern bool NOHUP; /**< Makes the signal handler ignore the first SIGHUP.*/

    /**
     * @typedef std::unordered_map<int,std::unordered_map<int,std::vector<book>>> library
     * @brief Defines the \c library type.
     * The nested structure of authors having books having fragments is of widespread use in fragments probability computation.
     */
    typedef std::unordered_map<dt::auth_id_t, std::unordered_map<dt::book_id_t, std::vector<std::unique_ptr<bp::book>>>> library;

    typedef std::unordered_map<dt::auth_id_t, std::unordered_map<dt::book_id_t, std::vector<bp::hash_type>>> sequence;

    /**
     * @struct TASK supportlib.hpp
     * @brief Contains information on the book to compute.
     *
     */
    struct TASK
    {
        dt::auth_id_t aut1 /** Author to compare with. */, aut2 /** Author of the book to compare. */;
        dt::book_id_t book /** Number of the book to compare */;
        bool operator==(struct TASK oth) { return aut1 == oth.aut1 && aut2 == oth.aut2 && book == oth.book; };
    };

    /**
     * @typedef typedef struct TASK task
     * @brief Defines as a type the struct TASK.
     *
     */
    typedef struct TASK task;

    constexpr int ARR_SIZ()
    {
        constexpr int ret = ((2 * sizeof(dt::auth_id_t) + sizeof(dt::book_id_t)) / sizeof(unsigned short int));
        return ret;
    }
    constexpr int OFFS_AU()
    {
        constexpr int ret = sizeof(dt::auth_id_t) / sizeof(unsigned short int);
        return ret;
    }
    constexpr int OFFS_BK()
    {
        constexpr int ret = 2 * sizeof(dt::auth_id_t) / sizeof(unsigned short int);
        return ret;
    }

    class slice
    {
    public:
        std::unordered_map<dt::auth_id_t, std::set<dt::book_id_t>> shelf;
        int size();
    };

    class task_iterator
    {
    private:
        /* data */
        dt::auth_id_t my_auth;
        std::vector<std::reference_wrapper<const supp::slice>> cakelist;
        size_t _size;
        std::vector<std::reference_wrapper<const supp::slice>>::iterator now_slice;
        std::unordered_map<dt::auth_id_t, std::set<dt::book_id_t>>::const_iterator now_aut;
        std::set<dt::book_id_t>::const_iterator now_book;
        std::queue<std::pair<dt::auth_id_t, dt::book_id_t>> returned;

        void calculate_size();

    public:
        task_iterator() = default;
        task_iterator(dt::auth_id_t aid, std::vector<std::reference_wrapper<const supp::slice>> cake_part, size_t size = 0);
        task pop();
        void place_back(task job);
        size_t size() const { return _size; };
    };

    /**
     * @class list_manager
     *
     * @brief An object that manages lists of tasks.
     *
     * Assigns tasks to workers caring to not change the author they
     * are working on (if possible) or assigning them the most needing
     * author among those still having books to compute.
     *
     */
    class list_manager
    {
    private:
        std::unordered_map<dt::auth_id_t, task_iterator> taskers; /**< A \c map with the queues of tasks to do. */
        mutable std::forward_list<dt::auth_id_t> Priority;        /**< A \c vector that contains the authors sorted on priority. */
        std::mutex lock;                                          /**< A \c mutex protecting from simultaneous accesses to the queue. */
        std::forward_list<dt::auth_id_t>::iterator actual;        /**< The index ot the next author to be assigned */
        std::set<dt::auth_id_t> placed_back;
        bool first_round, created=false;

        void clean_advance();

    public:
        list_manager() = default;
        list_manager(list_manager &&) = default;
        /**
         * @brief Initialize the list_manager from a \c lybrary object.
         *
         * For every author in the \c library a task is added for every book. If an author as only one book is not possible to compare it with himself
         * and that author-book pair has no corresponding task.
         *
         * @param shelf The \c library object used to initialize.
         */
        list_manager(const library &shelf, const std::vector<supp::slice> &cake, const std::vector<int> &slices);
        void create(const library &shelf, const std::vector<supp::slice> &cake, const std::vector<int> &slices);
        /**
         * @brief Assigns a task to a worker.
         *
         * This method assigns tasks to workers. If possible assigning the same author the worker is working on at the moment to reuse caches values
         * for the author's corpus and probability. If the author has no more available tasks assigns a task from a new author selected according to
         * priority.
         *
         * @param my_aut The code of the author the worker is computing at the moment.
         * @return task The \c task containing the author to compare with and the book to compare.
         */
        task pop(dt::auth_id_t my_aut);
        void place_back(task job);
        size_t remaining();
    };

    void hash_combine(size_t & seed, size_t const& v);

    /**
     * @brief A simple signal handler.
     *
     * This serves mostly to put threads to sleep if signalled with SIGUSR1.
     *
     * @param signo The number of the signal received.
     */
    void sig_handler(int signo);

    /**
     * @class user_stop
     *
     * @brief Custom exception in case of interrupt from the user.
     *
     */
    class user_stop : public std::exception
    {
        std::string msg = "Interrupted by the user. "; /**< Base message explaining the exception. */
    public:
        user_stop(){};

        /**
         * @brief Construct a new user stop object
         *
         * @param mess Message to be appended to the def    ault one.
         */
        user_stop(std::string mess) { msg += mess; };

        /**
         * @brief Describes the exception
         *
         * @return const char* Message.
         */
        const char *what() const noexcept { return msg.c_str(); };
    };

    void catch_signals();
    /**
     * @brief Loads the P0 from a file.
     *
     * @param[in] fileP0 Path to the file containing the P0.
     * @param[out] P0  The object that will contain the P0.
     */
    void load_P0(std::string fileP0, bp::glob_prob &P0);
    void dump_P0(std::filesystem::path outputFolderPath, bp::glob_prob &P0);

    void cake_from_library(std::vector<supp::slice> *cake, const library &shelf, int slicesize);

    int read_books_from_file(library &shelf_short, sequence &shelf_long, int A, std::queue<std::filesystem::path> &inputFile, bp::glob_prob &P0);

    std::string insert_padded_number(std::string first_part, std::string last_part, int n, int padsize = 3);

    unsigned short int split_save_number(library &shelf_short, int numThreads);
    template <typename T, typename U>
    class my_hash
    {
    public:
        std::size_t operator()(const std::pair<T, U> &val) const
        {
            std::size_t seed = 0;
            hash_combine(seed, val.first);
            hash_combine(seed, val.second);
            return seed;
        }
    };
    void load_main_author_parts(std::string fileName, std::map<dt::auth_id_t, dt::book_id_t> *main_author_parts, std::map<dt::auth_id_t, std::unordered_map<dt::slice_id_t, dt::book_id_t>> *slice_author_size);

    class comp_id
    {
    public:
        dt::auth_id_t comp_aut, auth;
        dt::book_id_t book, comp_frag_num;
        comp_id() = default;
        comp_id(int c, int a, int b, int cf = -1) : comp_aut(c), auth(a), book(b), comp_frag_num(cf){};
        comp_id(const comp_id &oth) : comp_aut(oth.comp_aut), auth(oth.auth), book(oth.book), comp_frag_num(oth.comp_frag_num){};
        comp_id(std::array<unsigned short int, ARR_SIZ()> tmp)
        {
            memcpy(&comp_aut, tmp.data(), sizeof(dt::auth_id_t));
            memcpy(&auth, tmp.data() + OFFS_AU(), sizeof(dt::auth_id_t));
            memcpy(&book, tmp.data() + OFFS_BK(), sizeof(dt::book_id_t));
        }
        comp_id operator=(const comp_id &oth)
        {
            comp_aut = oth.comp_aut;
            auth = oth.auth;
            book = oth.book;
            comp_frag_num = oth.comp_frag_num;
            return *this;
        };
        ~comp_id() = default;
        bool operator==(const comp_id &oth) const
        {
            return (this->comp_aut == oth.comp_aut) && (this->auth == oth.auth) && (this->book == oth.book);
        };
        std::string print()
        {
            std::stringstream out;
            out << "A" << comp_aut << "F" << comp_frag_num << "_A" << auth << "B" << book;
            return out.str();
        };
        operator std::array<unsigned short int, ARR_SIZ()>() const
        {
            std::array<unsigned short int, ARR_SIZ()> ret;
            memcpy(ret.data(), &comp_aut, sizeof(dt::auth_id_t));
            memcpy(ret.data() + OFFS_AU(), &auth, sizeof(dt::auth_id_t));
            memcpy(ret.data() + OFFS_BK(), &book, sizeof(dt::book_id_t));
            return ret;
        };
    };
    struct comp_hash
    {
        std::size_t operator()(const comp_id &val) const
        {
            std::size_t seed = 0;
            hash_combine(seed, val.comp_aut);
            hash_combine(seed, val.auth);
            hash_combine(seed, val.book);
            return seed;
        }
    };
    class result_box
    {
    public:
        result_box() = default;
        result_box(const result_box &) = delete;
        std::unordered_map<comp_id, std::vector<unsigned char>, comp_hash> results;
        std::mutex res_lock;
    };

    class Timer
    {
    public:
        Timer() : beg_(clock_::now()) {}
        void reset() { beg_ = clock_::now(); }
        double elapsed() const
        {
            return std::chrono::duration_cast<second_>(clock_::now() - beg_).count();
        }

    private:
        typedef std::chrono::high_resolution_clock clock_;
        typedef std::chrono::duration<double, std::ratio<1>> second_;
        std::chrono::time_point<clock_> beg_;
    };

    inline bool isInteger(const std::string &s)
    {
        if (s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+')))
            return false;

        char *p;
        strtol(s.c_str(), &p, 10);

        return (*p == 0);
    }

    void progress_bar(double tpass, size_t nowwri, size_t totwri, size_t remaining, size_t resCount, double eLAPsed, size_t totProbs, size_t totFra);

    class Semaphore
    {
        std::mutex m;
        std::condition_variable cv;
        int count;

    public:
        Semaphore(int n) : count{n} {}
        void notify()
        {
            std::unique_lock<std::mutex> l(m);
            ++count;
            cv.notify_one();
        }
        void wait()
        {
            std::unique_lock<std::mutex> l(m);
            cv.wait(l, [this]
                    { return count != 0; });
            --count;
        }
    };
    class Critical_section
    {
        Semaphore &s;

    public:
        Critical_section(Semaphore &ss) : s{ss} { s.wait(); }
        ~Critical_section() { s.notify(); }
    };
} // namespace supp
