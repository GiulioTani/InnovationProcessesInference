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
 * @file bookprob.hpp
 * @author Giulio Tani Raffaelli (tani@cs.cas.cz)
 * @brief Defines classes, exceptions and functions for CP2D
 * @version 1.0
 * @date 2023-06-27
 * 
 * @copyright Copyright (c) 2023
 * 
 * @bug If betwen checking the existance of the file in the constructor and opening it in fread_from_* the file is deleted the program closes without throwing.
 */
#pragma once
#include "lib/datatypes.hpp"
#include <cstring>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <unordered_map>
#include <mutex>

/**
 * @package bookprob
 * @brief Contains classes and exceptions to work with books.
 * 
 * All the tools needed to manipulate books in order to compute their probabilities.
 * 
 */
namespace bookprob
{
    typedef size_t hash_type;
    /**
     * @brief A struct to contain the global multiplicities of words holding also their sum to avoid to recompute it if useless.
     * 
     */
    struct glob_prob
    {
        std::unordered_map<hash_type, int> prob; /**< The \c map containing the P0.*/
        int count;                               /**< The sum of multiplicities of the words defining the P0.*/
    };

    /**
     * @class not_init
     * 
     * @brief Custom exception in case of use of an uninitialised book.
     * 
     */
    class not_init : public std::exception
    {
    public:
        /**
         * @brief Construct a new not init object
         * 
         * @param wword String describing what is not initialised.
         */
        explicit not_init(const std::string &wword) : word(wword) {}

        /**
         * @brief Describes the exception.
         * 
         * @return std::string Message.
         */
        std::string what() { return "Computing the probability of a book that has not been initialised. " + word; };

    private:
        std::string word; /**< String specifying what is not initialised. */
    };

    /**
     * @class word_miss
     * @brief Custom exception in case of missing word in probability.
     * 
     * The throwing of this exception means an error in the code if auto generated multiplicities as used, otherwise
     * it means an error (a word missing) in the probability file loaded.
     * 
     */
    class word_miss : public std::exception
    {
    public:
        /**
         * @brief Construct a new word miss object
         * 
         * @param wword The word that is missing.
         */
        explicit word_miss(const std::string &wword) : word(wword) {}

        /**
         * @brief Describes the exception.
         * 
         * @return std::string Message.
         */
        std::string what() { return "Word \"" + word + "\" in dictionary missing in probability."; };

    private:
        std::string word; /**< String specifying which word is missing. */
    };
    /**
     * @class word_doubled
     * @brief Custom exception in case of doubled word.
     * 
     * The throwing of this exception usually means an error in the input files. If dictionariy files are used for to initialise
     * books or probability, the same word cannot appear twice as the multiplicity to use would then be undefined.
     * 
     */
    class word_doubled : public std::exception
    {
    public:
        /**
         * @brief Construct a new word doubled object
         * 
         * @param wword The word found doubled.
         * @param wwhere The place (dictionary or probability) where the word was found doubled.
         */
        explicit word_doubled(const std::string &wword, const std::string &wwhere = "?") : word(wword), where(wwhere) {}

        /**
         * @brief Describes the exception.
         * 
         * @return std::string Message.
         */
        std::string what() { return "Word \"" + word + "\" in \"" + where + "\" doubled."; };

    private:
        std::string word /**< String specifying which word is doubled. */, where /**< String specifying where the doubled word was found. */;
    };
    class empty_file : public std::exception
    {
    public:
        /**
         * @brief Construct a new word doubled object
         * 
         * @param fname The word found doubled.
         */
        explicit empty_file(const std::string &fname = "?") : fname(fname) {}

        /**
         * @brief Describes the exception.
         * 
         * @return std::string Message.
         */
        std::string what() { return ("Empty file: \"" + fname + "\".").c_str(); };

    private:
        std::string fname /**< String specifying where the doubled word was found. */;
    };

    /**
     * @class book bookprob.hpp
     * @brief class used to compute probabilities following the continuous PD process.
     * 
     * This class stores the dictionary with the multiplicities of the words and their order. Given these data and a probability is possible
     * to compute the probability of the sequence according to the continuous Poisson-Dirichlet process and the probability updating scheme
     * from Tria and Lalli.
     * 
     */
    class book
    {
    public:
        /**
         * @brief Construct a new book object
         * 
         * Takes care that \c alpha and \c theta are marked as uninitialised.
         * 
         */
        book() : alpha(-100), theta(-100){};

        /**
         * @brief Construct a new book object copying an existing one.
         * 
         * @param oth The book to copy.
         */
        book(const book &oth);

        /**
         * @brief Construct a new book object by moving an existing one.
         * 
         * @param oth The book to move.
         */
        book(book &&oth);

        /**
         * @brief Construct a new book object
         * 
         * @param path Path in the filesyste where to look for the book. May be a sequence file (.seq) or a dictionary file (.wnt).
         * @param inP0 A pointer to the \c map containing the global multiplicities.
         * @param Alpha The \c alpha parameter for the PDP.
         * @param Theta The \c theta parameter for the PDP.
         * @param pref The prefix, i.e. an unique name, for the book.
         * @exception invalid_argument if the input file is missing or the type is unknown (not a .seq or .wnt).
         * @sa read_from_wnt
         * @sa read_from_seq
         */
        book(const std::string &path, const glob_prob *inP0 = NULL, double Alpha = -100, double Theta = -100, std::string pref = "");

        /**
         * @brief Construct a new book object
         * 
         * @param list A list of strings that is the book as extracted from a sequence file.
         * @param inP0 A pointer to the \c map containing the global multiplicities.
         * @param Alpha The \c alpha parameter for the PDP.
         * @param Theta The \c theta parameter for the PDP.
         * @param pref The prefix, i.e. an unique name, for the book. 
         */
        book(const std::vector<std::string> &list, const glob_prob *inP0 = NULL, double Alpha = -100, double Theta = -100, std::string pref = "");
        book(const std::vector<hash_type> &list, const glob_prob *inP0 = NULL, double Alpha = -100, double Theta = -100, std::string pref = "");

        /**
         * @brief Construct a new book object
         * 
         * @param ord An ordered sequence of (word, multiplicity) pairs as extracted froma a dictionary file.
         * @param inP0 A pointer to the \c map containing the global multiplicities.
         * @param Alpha The \c alpha parameter for the PDP.
         * @param Theta The \c theta parameter for the PDP.
         * @param pref The prefix, i.e. an unique name, for the book.
         */
        book(const std::vector<std::pair<hash_type, int>> &ord, const glob_prob *inP0 = NULL, double Alpha = -100, double Theta = -100, std::string pref = "");

        /**
         * @brief Joins two books in order without modifying them.
         * 
         * The book passed is appended to this preserving the order
         * of first appearance of the words.
         * 
         * @param other The book to be appended.
         * @return book The joined book.
         */
        book join(const book &other) const;

        /**
         * @brief Appends a book to another in-place
         * 
         * @param other The book to be appended.
         */
        virtual void append(const book &other);

        /**
         * @brief Copy intialization (alpha, theta, P0) from another book.
         * 
         * @param other The book to copy from.
         */
        void copy_init(const book &other);

        /**
         * @brief Deletes the whole content of the book making it uninitialised.
         * 
         */
        void clear();

        /**
         * @brief Sets the internal reference to P0
         * 
         * @param inP0 Pointer to \c map containing P0
         */
        void add_P0(const glob_prob *inP0) { P0 = inP0; };

        /**
         * @brief Copies the \c map containing P0 and setd the internal reference.
         * 
         * @param inP0 A \c map containing P0.
         */
        void add_P0(const glob_prob &inP0);

        /**
         * @brief Sets parameters \c alpha and \c theta checking their values.
         * 
         * While setting the parameters checks if \f$\alpha \in [0,1)\f$ and \f$\theta \in [-\alpha,\infty)\f$. 
         * 
         * @param Alpha The value of the \c alpha in the PDP.
         * @param Theta The value of the \c theta in the PDP.
         * @exception invalid_argument if parameters are out of the permitted interval.
         */
        void add_par(double Alpha, double Theta);
        void add_par(std::pair<double, double> par) { add_par(par.first, par.second); }
        /**
         * @brief Set the prefix object.
         * 
         * @param pref The prefix to be used.
         */
        void set_prefix(const std::string &pref) { prefix = pref; };

        /**
         * @brief Compute the log probability of the book.
         * 
         * This function computes the probability under the hypothesis of a PD process.
         * Each word \f$ w \f$ with a multiplicity \f$ n \f$ contributes to the probability with a factor
         * \f[
         *      P_0(w)\cdot \prod_{i=1}^{n} (i-\alpha)
         * \f]
         * The probability of the sequence reads:
         * \f[
         *      P= \frac{(\theta | \alpha)_K}{(\theta)_N}\prod_{k=1}^K P_0(w_k) (1-\alpha)_{n_k-1}
         * \f]
         * Where \f$ N \f$ is the total number of words and \f$ K \f$ is the total number of differents words.
         * 
         * @return double The computed log probability.
         * @exception not_init in case any field of the book is not initialised.
         */
        std::pair<double, dt::diff_tok_t> log_prob() const;

        /**
         * @brief Compute the log conditional probability of other given the book.
         * 
         * This function computes the probability under the hypothesis of a PD process in which
         * other is the final part of book and one is interested only in its probability.
         * The probability of the sequence reads:
         * \f[
         *      P= \frac{(\theta + \alpha \cdot K | \alpha)_{K'}}{(\theta + N)_{N'}}\prod_{k=1}^K P_0(w_k) (1-\alpha)_{n'_k-n_k} \prod_{k=K}^{K+K'} P_0(w_k) (1-\alpha)_{n'_k-1}
         * \f]
         * Where \f$ N \f$ is the total number of words and \f$ K \f$ is the total number of differents words in book and the primed variables are referred to other.
         * 
         * @param other 
         * @return double The computed log probability.
         * @exception not_init in case any field of the book is not initialised.
         * @sa log_prob
         */
        std::pair<double, dt::diff_tok_t> log_prob(const book &other) const;

        /**
         * @brief Get the prefix object.
         * 
         * @return std::string The prefix.
         */
        std::string get_prefix() const { return prefix; };

        /**
         * @brief Get the dictionary object.
         * 
         * Returns a pointer to a non modifiable view of the dictionary useful e.g. to get its length.
         * If the dictionary is not initialised returns NULL.
         * 
         * @return const std::vector<std::pair<std::string,int>>* Pointer to the dictionary.
         */
        const std::vector<std::pair<hash_type, int>> *get_dictionary() const
        {
            if (order.size() == 0)
                return NULL;
            return &order;
        };

        /**
         * @brief Get the params object.
         * 
         * If \c alpha or \c theta are not initialized both their values are set to -100.
         * 
         * @return std::pair<double,double> A pair containing ( \c alpha , \c theta )
         */
        std::pair<double, double> get_params() const { return std::make_pair(alpha, theta); };

        /**
         * @brief Allows move construction when using operator=.
         * 
         * @param oth The book to be move constructed.
         * @return book& The resulting book.
         */
        book &operator=(book &&oth);

        size_t get_N(){return N;}
        virtual void log_prob_to_chararr(const book &other, void* dest) const;
        virtual void log_prob_to_chararr(void* dest) const;
        virtual size_t retsize() const;
    protected:
        double alpha /** The \c alpha parameter of the PDP*/, theta /** The \c theta parameter of the PDP*/;
        size_t N /** The length of the book.*/;
        std::vector<std::pair<hash_type, int>> order /** The order of the words in the book with their multiplicities.*/;
        const glob_prob *P0 /** Pointer to the \c struct containing the P0*/;
        std::string prefix /** The prefix associated to the book, is intended to be an unique identifier.*/;
        std::unordered_map<hash_type, int> positions /** The \c map linking the words in the book to their order in the sequence*/;
        mutable std::vector<int> tmpP0 /** The multiplicities of the words not yet appeared in the sequence. */;
        mutable int fixWeight /** The sum of the multiplicities of words missing from the book*/, nowWeight;

        /**
         * @brief Initialise the tmpP0 and fixedweight.
         * 
         * Sums in \c fixedweight the multiplicities of all the words absent in the book and in other. Creates a new \c map \c tmpP0 restricted to the words in other.
         * 
         * @param other Another book.
         * 
         * @exception word_miss if a word from the book or from other is missing in P0.
         */
        void init_tmpP0(const book &other) const;

        /**
         * @brief Returns the lop P of a word removing it from tmpP0
         * 
         * P is computed as the fraction of the occurrences of \c word over the sum of \c fixedweight and of the occurrencies of all the words in \c tmpP0.
         * Then \c word is removed from \c tmpP0.
         * 
         * @param word The word of interest.
         * @return double The associated log P.
         * @sa init_tmpP0
         */
        double logP_words() const;
    private:
        /**
         * @brief Initialise \c alpha and \c theta parameters.
         * 
         * Only if both \c alpha and \c theta have a value different from -100 (i.e. they are initialised) passes them to add_par. Otherwise, they are both
         * lef uninitialised: is not possible to have a half defined PDP.
         * @sa add_par
         * 
         * @param Alpha The \c alpha parameter from the constructor.
         * @param Theta The \c theta parameter from the constructor.
         */
        void initpar(double Alpha, double Theta);
        glob_prob myP0 /** The \c struct containing the P0 if not using an external reference*/;

        /**
         * @brief Reads the book from a dictionary file.
         * 
         * The book in a dictionary file is presented as a list of words followed by their multiplicities, one word per line.
         * The order of the words is that of first appearance in the book and no word should appear twice.
         * 
         * @param file The path to a .wnt file containing the book.
         * 
         * @exception word_doubled if using a dictionary file and the same word is found twice.
         */
        void read_from_wnt(std::string file);

        /**
         * @brief Reads the book from a sequence file.
         * 
         * The sequence file is a list of the words in the book one per line.
         * 
         * @param file The path to a .seq file containing the book.
         * @sa read_from_list
         */
        void read_from_seq(std::string file);

        /**
         * @brief Reads the book from a list of words.
         * 
         * Reads the \c vector containing the sequence of the words, either from the calling progam or from the sequence file, and builds the book.
         * 
         * @param list The vector containing the book.
         */
        void read_from_list(std::vector<std::string> list);
        void read_from_list(std::vector<hash_type> list);

        /**
         * @brief Initialise the tmpP0 and fixedweight.
         * 
         * Sums in \c fixedweight the multiplicities of all the words absent in the book. Creates a new \c map \c tmpP0 restricted to the words in the book.
         * 
         * @exception word_miss if a word from the book is missing in P0.
         */
        void init_tmpP0() const;
    };

    extern std::unordered_map<hash_type, std::string> missing_words;
    extern std::mutex missing_words_lock;
    /**
     * @typedef std::unordered_map<int,std::string> trueName
     * @brief Defines the \c trueName type.
     * Used to store the association between the codes (for authors and books) and their true names.
     */
    typedef std::unordered_map<int, std::string> trueName;
} // namespace bookprob
