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
#include "lib/bookprob.hpp"   //defines namespace book
#include "lib/supportlib.hpp" //defines namespace supp
#include "lib/datatypes.hpp"
#include <filesystem>
#include <mutex>
#include <thread>
#include <map>
#include <atomic>
#include <fstream>

namespace fs = std::filesystem;

#define defaultNumThreads 36 /**< Default number of threads to use.*/
#define nofrag 0x01
#define noauth 0x02
namespace be
{
    class base_experiment
    {
    private:
        std::atomic_int threadCount; /**< Used in main to keep track of the running threads.*/
        std::mutex par_lock, res_lock;
        std::vector<supp::result_box> results;
        std::vector<std::thread> workers;
        char flagbyte = 0;
        std::string outputFile;
        supp::library shelf_short;
        supp::sequence shelf_long;
        bp::glob_prob P0;
        std::map<dt::auth_id_t, std::unordered_map<dt::slice_id_t, dt::book_id_t>> slice_author_size;
        std::map<dt::auth_id_t, dt::book_id_t> main_parts;
        std::vector<supp::slice> cake;
        fs::path outputFolderPath;
        int totFrag;
        std::unordered_map<std::pair<dt::auth_id_t, dt::book_id_t>, std::pair<double, double>, supp::my_hash<dt::auth_id_t, dt::book_id_t>> params;
        std::set<std::pair<dt::auth_id_t, dt::slice_id_t>> known_par, computing_par;
        long int numThreads = -2000000000 /**< Number of threads to use.*/, F = -2000000000 /**< Number of tokens in fragmen*/;
        std::ofstream outpar;
        std::unordered_map<std::pair<dt::auth_id_t, dt::book_id_t>, dt::slice_id_t, supp::my_hash<dt::auth_id_t, dt::book_id_t>> which_slice;
        std::unordered_map<dt::auth_id_t, unsigned short> author_bucket;
        unsigned short number_of_buckets;
        std::unordered_set<supp::comp_id, supp::comp_hash> written;
        size_t nowwri = 0, totwri = 0, totFra = 0;
        supp::list_manager authList;
        bool resuming;

        void retrieve_number_of_buckets();
        void retrieve_results();
        void results_loader(std::string fname);
        void read_books(const fs::path &inputFolder);
        void get_P0();
        bool get_slices(const fs::path &inputFolder, unsigned int slicesize);
        bool read_slices_from_file(std::string fname);
        void dump_slices_to_file(std::string fname);
        void get_params(const fs::path &inputFolder, bool new_slices);
        void dump_fra_info();
        void dump_and_process_aut_info();
        void emplace_results(supp::result_box &results, supp::comp_id this_comp, std::vector<unsigned char> my_results);
        std::vector<std::string> output_file_names(std::string midpart);
        void results_orderer(int target);
        void worker(supp::result_box &results, supp::list_manager &authList);
        int collect_results(supp::result_box &slavres, std::vector<std::ofstream> &output_ids, std::vector<std::ofstream> &output_pro);
        void open_output_files(std::string midpart, std::vector<std::ofstream> &output_vector, bool dumpBuck=false);
        void dump_new_params();
        void lap_update(supp::Timer &timer_lap, supp::Timer &timer_tot, std::vector<std::ofstream> &output_ids, std::vector<std::ofstream> &output_pro, size_t totProbs);
        void spread_full_author_parameters();
        

    public:
        base_experiment(/* args */);
        base_experiment(int argc, char **argv);
        ~base_experiment();

        size_t load_prev_results();
        void run();
        void order();
    };

    std::pair<std::array<unsigned char, 8>, dt::diff_tok_t> trans_doub_to_chararr(std::pair<double, dt::diff_tok_t> pres);

} // namespace be
