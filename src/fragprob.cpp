/**
 * @file fragprob.cpp
 * @author Giulio Tani (giulio.tani@uniroma1.it)
 * @brief Computes conditional probabilities of fragments give authors
 * @version 0.1
 * @date 2020-07-12
 * @copyright Copyright (c) 2020
 *
 * Recieves from command line at least the \c -f option designing the input folder and other optional parameters as
 * \c -o a non default name for the output file, \c -p a file containing the multiplicities to be used for the probabilities,
 * \c -t the number of threads used for computation.
 *
 * The default name for the output file is ${input_folder}/tmpres.res.
 *
 * The default value of the multiplicities is taken directly from the words in the input corpus.
 *
 * The default number of threads is 36 as this has been developed for a machine with 40 cores and I had to leave some
 * of them free.
 *
 */
/* given two books and relative parameters compute the conditioned probabilities */
#include "lib/bookprob.hpp"   //defines namespace book
#include "lib/supportlib.hpp" //defines namespace supp
#include "lib/base_experiment.hpp"
#include <iostream>
#include <iomanip>
#include <stdio.h>

#define defaultNumThreads 36 /**< Default number of threads to use.*/

namespace bp = bookprob;
/**
 * @brief Computes the conditional probabilities of fragments.
 *
 * Given a corpus of books of various authors computes di conditional probability of each fragment of each book of each author given each
 * author's coprus.
 *
 * @param[in] argc Count of command line arguments.
 * @param[in] argv Command line arguments.
 * @return int Return value.
 */
int main(int argc, char **argv)
{
    supp::Timer timer_tot;
    size_t prev_res;
    be::base_experiment exp(argc, argv);
    /* Initialization */
    supp::catch_signals();

    std::cout << "Init library: " << std::fixed << std::setprecision(3) << timer_tot.elapsed() << " s" << std::endl;

    timer_tot.reset();
    std::cout << "Init  queue : " << std::fixed << std::setprecision(3) << timer_tot.elapsed() << " s" << std::endl;

    timer_tot.reset();
    prev_res = exp.load_prev_results();
    if (prev_res)
    {
        std::cout << "Init  written : " << std::fixed << std::setprecision(3) << timer_tot.elapsed() << " s. Found " << prev_res << " precomputed results." << std::endl;
    }
    else
    {
        std::cout << "No precomputed results found." << std::endl;
    }

    exp.run();

    if (supp::FLAG)
    {
        exp.order();
    }

    if (!bp::missing_words.empty())
    {
        std::string outputFile = "";
        std::cout << "Some words were missing. Dumping to " + outputFile + "_wordmiss.res" << std::endl;
        std::ofstream output_wm(outputFile + "_wordmiss.res", std::ofstream::out | std::ofstream::app);
        for (auto w : bp::missing_words)
            output_wm << w.first << " : " << w.second << "\n";
    }

    if (supp::FLAG)
    {
        std::cout << "Done!" << std::endl;
        return EXIT_SUCCESS;
    }
    throw supp::user_stop("Received SIGINT.");
    return EXIT_FAILURE;
}