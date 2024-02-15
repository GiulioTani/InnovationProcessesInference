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

#include "lib/base_experiment.hpp"
#include "lib/authorSplitter.hpp"
#include "lib/paramOpt.hpp" //defines namespace popt
#include <future>
#include <iostream>
#include <iomanip>
#include <regex>
#include <cmath>
#include <unistd.h>
// #include <algorithm>
#define _USE_JOINT_PARAMS_
// #define _USE_BOOK_COUNTS_
// #define _FULL_AUTHOR_PARAMETERS_
namespace be
{
    base_experiment::~base_experiment()
    {
        if (outpar.is_open())
            outpar.close();
    }

    base_experiment::base_experiment(int argc, char **argv)
    {
        char c;
        int slicesize = 10, totBooks = 0;
        std::string inputFolderName = "", line, inputTrueName, prefix, outputFolderName;
        fs::path inputFolder;
        std::vector<std::vector<int>> fragments;
        std::vector<int> slices;
        std::ostringstream fragname;
        std::ifstream inputfiles, inputseq;
        std::stringstream ss;
        P0.count = -100;
        bool new_slices, dumpP0 = false;
        resuming = false;

        while ((c = getopt(argc, argv, "f:F:o:p:t:s:d:n:S:P")) != -1 || argc == 0)
        {
            switch (c)
            {
            case 'f':
                inputFolderName = optarg;
                break;
            case 'F':
                F = std::atoi(optarg);
                break;
            case 'o':
                outputFolderName = optarg;
                break;
            case 'p':
                supp::load_P0(optarg, P0);
                break;
            case 't':
                numThreads = std::atoi(optarg);
                break;
            case 's':
                slicesize = std::atoi(optarg);
                break;
            case 'd':
                asp::authorDimension = std::atoi(optarg);
                break;
            case 'n':
                asp::ngramSize = std::atoi(optarg);
                break;
            case 'S':
                slices.push_back(std::atoi(optarg));
                break;
            case 'P':
                dumpP0 = true;
                break;
            default:
                std::cerr << "Unknown option: \"" << c << "\"\nUsage: fragprob -f [input folder name] (-o [output file name])." << std::endl;
                exit(EXIT_FAILURE);
                break;
            }
        }
        /* Check the number of threads */
        if (numThreads <= 0)
        {
            if (numThreads == -2000000000)
                numThreads = defaultNumThreads;
            else
                throw std::invalid_argument("The number of computing threads must be positive.");
        }
        /* Check the fragments size */
        if (F < 0)
        {
            if (F == -2000000000)
                F = 0;
            else
                throw std::invalid_argument("The size of fragmets must be positive or 0 (for full texts).");
        }
        if ((long unsigned)F > std::numeric_limits<dt::diff_tok_t>::max())
            throw std::invalid_argument("The size of fragmets is larger than the mazimum number of unknown tokens"
                                        "\nRecompile with a larger type for diff_tok_t or change fragment size. " +
                                        std::to_string(F));
        /* Check the size of virtual authors */
        if (asp::authorDimension < 0)
        {
            if (asp::authorDimension == -2000000000)
                asp::authorDimension = 0;
            else
                throw std::invalid_argument("The size of virtual authors must be positive.");
        }
        /* Check the size of N-grams */
        if (asp::ngramSize < 0)
        {
            if (asp::ngramSize == -2000000000)
                asp::ngramSize = 0;
            else
                throw std::invalid_argument("The size of N-grams must be positive.");
        }
        /* Check if the input folder is valid */
        if (inputFolderName == "")
            throw std::invalid_argument("You must provide an input folder name.");
        inputFolder = fs::path(inputFolderName);
        if (!fs::is_directory(inputFolder))
            throw std::invalid_argument("Input folder must be a directory.");
        /* Check if the output folder is valid */
        if (outputFolderName == "")
            throw std::invalid_argument("You must provide an output folder name.");
        outputFolderPath = fs::path(outputFolderName);
        if (!fs::is_directory(outputFolderPath))
            throw std::invalid_argument("Output folder must be a directory.");
        /* Check output file or creates a default output file if needed*/
        outputFile = (outputFolderPath / std::string("tmpres")).string();
        std::cout << "Saving results in: " << outputFile << std::endl;

        read_books(inputFolder);

        get_P0();

        if (dumpP0)
            supp::dump_P0(outputFolderPath, P0);
        /* Check the slice size */
        for (auto &A : shelf_long)
            totBooks += A.second.size();
        if (slicesize <= 0 || slicesize > totBooks / sqrt(2))
        {
            throw std::invalid_argument("The slice size must be in [1, #books/sqrt(2)] (found " + std::to_string(slicesize) + " with " + std::to_string(totBooks) + " books).");
        }
        new_slices = get_slices(inputFolder, slicesize);

        get_params(inputFolder, new_slices);
#ifdef _FULL_AUTHOR_PARAMETERS_
        spread_full_author_parameters();
#endif

        flagbyte += (!F) * nofrag + (!asp::authorDimension) * noauth;

        dump_fra_info();

        dump_and_process_aut_info();

        authList.create(shelf_short, cake, slices);
        std::cout<<"Total fragments: "<< totFrag << std::endl;
    }

    void base_experiment::read_books(const fs::path &inputFolder)
    {
        /* Reads the books */
        fs::directory_iterator dit{inputFolder};
        std::queue<fs::path> file_names;
        std::vector<std::future<int>> tasks;
        for (const auto &fname : dit)
            if (fs::path(fname).extension() == ".seq")
                file_names.push(fname.path());
        if (file_names.size() > std::numeric_limits<dt::auth_id_t>::max())
            throw std::invalid_argument("The number of authors is larger than the mazimum number allowed"
                                        "\nRecompile with a larger type for auth_id_t. " +
                                        std::to_string(file_names.size()));

        totFrag = 0;
        for (auto i = 0; i < numThreads; i++)
        {
            tasks.emplace_back(std::async(std::launch::async, supp::read_books_from_file, std::ref(shelf_short), std::ref(shelf_long), std::ref(F), std::ref(file_names), std::ref(P0)));
        }
        for (auto &&task : tasks)
            totFrag += task.get();
    }

    void base_experiment::get_P0()
    {
#ifdef _USE_BOOK_COUNTS_
        std::unordered_set<bp::hash_type> tmp;
#endif // _USE_BOOK_COUNTS_
        if (P0.prob.size() == 0)
        {
            P0.count = 0;
            for (auto &auth : shelf_short)
                for (auto &volume : auth.second)
#ifdef _USE_BOOK_COUNTS_
                {
                    tmp.clear();
                    for (auto &fragm : volume.second)
                        for (auto p : *(fragm.second.get_dictionary()))
                            tmp.insert(p.first);
                    for (auto p : tmp)
                        if (!P0.prob.emplace(p, 1).second)
                            P0.prob[p] += 1;
                    P0.count += tmp.size();
                }
#else
                    for (auto &fragm : volume.second)
                    {
                        for (auto p : *(fragm->get_dictionary()))
                        {
                            if (!P0.prob.emplace(p).second)
                                P0.prob[p.first] += p.second;
                            P0.count += p.second;
                        }
                    }
#endif // _USE_BOOK_COUNTS_
            std::cout << "Computed P0" << std::endl;
        }
        P0.prob.emplace(std::make_pair(std::hash<std::string>{}(""), 0));
    }

    bool base_experiment::get_slices(const fs::path &inputFolder, unsigned int slicesize)
    {
        bool new_slices;

        /* Reads the slicing if provided */
        try
        {
            new_slices = read_slices_from_file((inputFolder / "slices.bin").string());
        }
        catch (bp::word_doubled &e)
        {
            std::cerr << e.what() << std::endl;
            throw e;
        }
        catch (bp::word_miss &e)
        {
            std::cerr << e.what() << std::endl;
            throw e;
        }

        if (new_slices)
            supp::cake_from_library(&cake, shelf_short, slicesize);

        dump_slices_to_file((outputFolderPath / "slices.bin").string());

        for (dt::slice_id_t i = 0; i < cake.size(); i++)
            for (auto &a : cake[i].shelf)
                for (auto &b : a.second)
                    if (a.first > 0)
                        which_slice[{a.first, b}] = i;
                    else
                        which_slice[{a.first, b}] = -1;

        for (auto &auth : shelf_long)
            for (auto &book : auth.second)
                if (which_slice.find({auth.first, book.first}) == which_slice.end())
                    throw bp::word_miss("A book has no place in any slice: " + std::to_string(auth.first) + " " + std::to_string(book.first));

        return new_slices;
    }

    bool base_experiment::read_slices_from_file(std::string fname)
    {
        std::ifstream inputfiles;
        inputfiles.open(fname.c_str(), std::ifstream::binary);
        if (inputfiles.is_open())
        {
            dt::slice_id_t N_slices;
            dt::auth_id_t A, N_auth;
            dt::book_id_t N_books, B;
            inputfiles.read((char *)&N_slices, sizeof(N_slices));
            for (dt::slice_id_t i = 0; i < N_slices; i++)
            {
                cake.push_back(supp::slice());
                inputfiles.read((char *)&N_auth, sizeof(N_auth));
                for (dt::auth_id_t j = 0; j < N_auth; j++)
                {
                    inputfiles.read((char *)&A, sizeof(A));
                    cake[i].shelf.emplace(A, std::set<dt::book_id_t>());
                    inputfiles.read((char *)&N_books, sizeof(N_books));
                    for (dt::book_id_t k = 0; k < N_books; k++)
                    {
                        inputfiles.read((char *)&B, sizeof(B));
                        if ((shelf_long.find(A) == shelf_long.end()) || (shelf_long[A].find(B) == shelf_long[A].end()))
                            throw bp::word_miss("Adding in slice a book that does not exist: " + std::to_string(i) + " " + std::to_string(A) + " " + std::to_string(B));
                        if (!(cake[i].shelf.at(A).emplace(B).second))
                            throw bp::word_doubled("Same book twice in slice: " + std::to_string(i) + " " + std::to_string(A) + " " + std::to_string(B));
                    }
                }
            }
        }
        if (inputfiles.is_open())
            inputfiles.close();
        return cake.size() == 0;
    }

    void base_experiment::dump_slices_to_file(std::string fname)
    {
        std::ofstream out_support;
        out_support.open(fname.c_str(), std::ofstream::out | std::ofstream::binary);
        if (out_support.is_open())
        {
            dt::slice_id_t N_slices;
            dt::auth_id_t A, N_auth;
            dt::book_id_t N_books;
            N_slices = (dt::slice_id_t)cake.size();
            out_support.write((char *)&N_slices, sizeof(N_slices));
            for (dt::slice_id_t i = 0; i < N_slices; i++)
            {
                N_auth = cake[i].shelf.size();
                out_support.write((char *)&N_auth, sizeof(N_auth));
                for (auto &auth : cake[i].shelf)
                {
                    A = auth.first;
                    out_support.write((char *)&A, sizeof(A));
                    N_books = (dt::book_id_t)auth.second.size();
                    out_support.write((char *)&N_books, sizeof(N_books));
                    for (auto B : auth.second)
                        out_support.write((char *)&B, sizeof(B));
                }
            }
            out_support.close();
        }
    }

    void base_experiment::get_params(const fs::path &inputFolder, bool new_slices)
    {
        /* Reads parameters for all authors from params.par file if provided*/
        std::string inputTrueName = "params.bin";
        std::pair<dt::auth_id_t, dt::book_id_t> par_id;
        std::pair<double, double> par_val;
        std::ifstream inputfiles;
        inputfiles.open((inputFolder / inputTrueName).c_str(), std::ofstream::binary);
        if (inputfiles.is_open())
        {
            while (!inputfiles.eof())
            {
                inputfiles.read((char *)&par_id.first, sizeof(par_id.first));
                inputfiles.read((char *)&par_id.second, sizeof(par_id.second));
                inputfiles.read((char *)&par_val, sizeof(par_val));
#ifdef _FULL_AUTHOR_PARAMETERS_
                if (par_id.second==-1)
#endif
                    if (!(params.emplace(par_id, par_val).second) && !inputfiles.eof())
                        throw bp::word_doubled("Same parameters twice: " + std::to_string(par_id.first) + " " + std::to_string(par_id.second));
            }
        }
        if (inputfiles.is_open())
        {
            inputfiles.close();
            if (new_slices)
                std::remove((inputFolder / inputTrueName).c_str()); // removes outdated params file
        }

        // write valid params
        outpar.open((outputFolderPath / inputTrueName).c_str(), std::ofstream::out | std::ofstream::binary);
        if (outpar.fail() || !outpar.is_open())
        {
            std::cerr << "Error opening: " << (outputFolderPath / inputTrueName).c_str() << std::endl;
            exit(EXIT_FAILURE);
        }
        for (auto &add : params)
        {
            known_par.emplace(add.first);
            outpar.write((char *)&add.first.first, sizeof(add.first.first));
            outpar.write((char *)&add.first.second, sizeof(add.first.second));
            outpar.write((char *)&add.second, sizeof(add.second));
        }
    }

    void base_experiment::dump_fra_info()
    {
        std::ofstream out_support;
        if (!fs::exists(outputFile + "_fra.res"))
        {
            out_support.open(outputFile + "_fra.res", std::ofstream::out | std::ofstream::binary);
            if (out_support.fail() || !out_support.is_open())
            {
                std::cerr << "Error opening: " << outputFile << std::endl;
                exit(EXIT_FAILURE);
            }
            if (out_support.write((char *)&flagbyte, sizeof(char)).good())
                out_support.flush();
            for (auto &auth : shelf_short)
            {
                for (auto &book : auth.second)
                {
                    out_support.write((char *)(&auth.first), sizeof(auth.first));
                    out_support.write((char *)(&book.first), sizeof(book.first));
                    if (flagbyte & nofrag)
                    {
                        dt::diff_tok_t fsiz = book.second[0]->get_N();
                        out_support.write((char *)(&fsiz), sizeof(fsiz));
                    }
                    else
                    {
                        dt::frag_id_t fnum = book.second.size();
                        out_support.write((char *)(&fnum), sizeof(fnum));
                        for (dt::frag_id_t i = 0; i < fnum; i++)
                        {
                            dt::diff_tok_t fsiz = book.second[i]->get_N();
                            out_support.write((char *)(&fsiz), sizeof(fsiz));
                        }
                    }
                }
            }
            out_support.close();
        }
    }

    void base_experiment::dump_and_process_aut_info()
    {
        std::ofstream out_support;
        main_parts.clear();
        if (!(flagbyte & noauth))
        {
            if (fs::exists(outputFile + "_aut.res"))
                supp::load_main_author_parts(outputFile + "_aut.res", &main_parts, &slice_author_size);
            else
            {
                std::set<dt::auth_id_t> sortaut;
                for (auto &S_auth : shelf_long)
                    sortaut.emplace(S_auth.first);
                out_support.open(outputFile + "_aut.res", std::ofstream::out | std::ofstream::binary);
                if (out_support.fail() || !out_support.is_open())
                {
                    std::cerr << "Error opening: " << outputFile << std::endl;
                    exit(EXIT_FAILURE);
                }
                dt::auth_id_t N_auth = (dt::auth_id_t)shelf_long.size();
                if (shelf_long.find(0) != shelf_long.end())
                    N_auth--;
                out_support.write((char *)&N_auth, sizeof(N_auth));
                for (auto A : sortaut)
                {
                    if (!A)
                        continue;
                    out_support.write((char *)&A, sizeof(A));
                    dt::book_id_t N_slices = 0;
                    for (dt::slice_id_t i = 0; i < cake.size(); i++)
                        if (cake[i].shelf.find(A) != cake[i].shelf.end())
                            N_slices++;
                    out_support.write((char *)&N_slices, sizeof(N_slices));

                    dt::book_id_t default_size = asp::parts_num(shelf_long.at(A).size(), 0, asp::get_fullLength(shelf_long, A, std::set<dt::book_id_t>()));
                    main_parts[A] = default_size;
                    out_support.write((char *)&default_size, sizeof(default_size));

                    for (dt::slice_id_t i = 0; i < cake.size(); i++)
                    {
                        if (cake[i].shelf.find(A) != cake[i].shelf.end())
                        {
                            std::set<dt::book_id_t> excl = cake[i].shelf.at(A);
                            dt::book_id_t sli_siz = asp::parts_num(shelf_long.at(A).size(), excl.size(), asp::get_fullLength(shelf_long, A, excl));
                            slice_author_size[A][i] = sli_siz;
                            out_support.write((char *)&i, sizeof(i));
                            out_support.write((char *)&sli_siz, sizeof(sli_siz));
                        }
                    }
                }
                out_support.close();
            }
        }
    }

    void base_experiment::retrieve_number_of_buckets()
    {
        std::ifstream input;
        fs::directory_iterator dit{outputFolderPath};
        std::set<std::string> ids_file_names, pro_file_names;
        for (const auto &fname : dit)
        {
            if (std::regex_search(fname.path().c_str(), std::regex("_uo_ids[\\d]{3}")))
                ids_file_names.emplace(fname.path().string());
            if (std::regex_search(fname.path().c_str(), std::regex("_uo_pro[\\d]{3}")))
                pro_file_names.emplace(fname.path().string());
        }
        std::cout << "Found " << ids_file_names.size() << " ids and " << pro_file_names.size() << "pro files." << std::endl;
        number_of_buckets = 0;
        for (auto &fname : ids_file_names)
        {
            if (!number_of_buckets)
            {
                input.open(fname, std::ifstream::binary);
                if (input)
                {
                    input.read((char *)&number_of_buckets, sizeof(number_of_buckets));
                    input.close();
                }
                else
                {
                    if (fs::is_regular_file(fname))
                        std::remove(fname.c_str()); // remove ill formed ids
                }
            }
            if (!fs::is_regular_file(std::regex_replace(fname, std::regex("ids"), "pro")) && fs::is_regular_file(fname))
                std::remove(fname.c_str()); // remove non matched files
        }
        for (auto &fname : ids_file_names)
        {
            if (!fs::is_regular_file(std::regex_replace(fname, std::regex("pro"), "ids")) && fs::is_regular_file(fname))
                std::remove(fname.c_str()); // remove non matched files
        }
        std::cout << "Number of buckets: " << number_of_buckets << std::endl;
    }

    void base_experiment::results_loader(std::string fname)
    {
        std::ifstream input;
        std::unordered_set<supp::comp_id, supp::comp_hash> tmp_written;
        const unsigned int recorLen = shelf_short.begin()->second.begin()->second.front()->retsize();
        std::string pro_file = std::regex_replace(fname, std::regex("ids"), "pro");
        size_t offNow = 0, tmp_totfra = 0;
        input.open(fname, std::ifstream::binary);
        if (input)
        {
            input.seekg(sizeof(number_of_buckets)); // to skip the bytes storing the number of buckets
            std::array<unsigned short, supp::ARR_SIZ()> tmp;
            while (!input.eof())
            {
                input.read((char *)&tmp, sizeof(unsigned short int) * supp::ARR_SIZ());
                if (!input.eof())
                {
                    supp::comp_id tmp_id(tmp);
                    dt::frag_id_t no_of_fragm = 0;
                    dt::book_id_t no_of_auth_pieces = 0;
                    try
                    {
                        no_of_fragm = (flagbyte & nofrag) ? 1 : shelf_short.at(tmp_id.auth).at(tmp_id.book).size();
                    }
                    catch (std::out_of_range &e)
                    {
                        std::cerr << std::this_thread::get_id() << " no_of_fragm " << e.what() << std::endl;
                    }
                    try
                    {
                        no_of_auth_pieces = (flagbyte & noauth) ? 1 : (((tmp_id.auth == 0) || (cake[which_slice[{tmp_id.auth, tmp_id.book}]].shelf.find(tmp_id.comp_aut) == cake[which_slice[{tmp_id.auth, tmp_id.book}]].shelf.end())) ? main_parts.at(tmp_id.comp_aut) : slice_author_size.at(tmp_id.comp_aut).at(which_slice[{tmp_id.auth, tmp_id.book}]));
                    }
                    catch (std::out_of_range &e)
                    {
                        std::cerr << std::this_thread::get_id() << " no_of_auth_pieces " << e.what() << std::endl;
                    }
                    offNow += no_of_fragm * no_of_auth_pieces;
                    tmp_written.emplace(tmp_id);
                    tmp_totfra += shelf_short.at(tmp_id.auth).at(tmp_id.book).size();
                }
            }
            input.close();
        }
        else
            throw std::runtime_error("Missing file " + fname);

        if (fs::file_size(pro_file) != offNow * recorLen)
        {
            if (fs::is_regular_file(fname))
                std::remove(fname.c_str());
            if (fs::is_regular_file(pro_file))
                std::remove(pro_file.c_str());
        }
        else
        {
            std::lock_guard<std::mutex> gurs(res_lock);
            totwri += tmp_written.size();
            for (auto &&id : tmp_written)
                written.emplace(id);
            totFra += tmp_totfra;
        }
    }

    void base_experiment::retrieve_results()
    {
        fs::directory_iterator dit{outputFolderPath};
        std::set<std::string> ids_file_names;
        for (const auto &fname : dit)
        {
            if (std::regex_search(fname.path().c_str(), std::regex("_uo_ids[\\d]{3}")))
                ids_file_names.emplace(fname.path().string());
        }
        supp::Semaphore task_limiter(numThreads);
        std::vector<std::future<void>> tasks;

        std::cout << "\nLoading" << std::endl;
        for (auto &fname : ids_file_names)
        {
            tasks.emplace_back(std::async(std::launch::async, [fname, this, &task_limiter]
                                          {
                supp::Critical_section __(task_limiter);
                results_loader(fname); }));
        }
        for (auto &&task : tasks)
            task.wait();
        resuming = true;
    }

    size_t base_experiment::load_prev_results()
    {
        std::ifstream input;
        retrieve_number_of_buckets();
        if (!number_of_buckets) // I failed retrieving no more files left
        {
            number_of_buckets = supp::split_save_number(shelf_short, numThreads);
        }
        else
        {
            retrieve_results();
        }
        return totwri;
    }

    std::vector<std::string> base_experiment::output_file_names(std::string midpart)
    {
        std::vector<std::string> result;
        for (auto i = 0; i < number_of_buckets; i++)
        {
            result.push_back(supp::insert_padded_number(outputFile + midpart, ".res", i));
        }
        return result;
    }

    void base_experiment::emplace_results(supp::result_box &results, supp::comp_id this_comp, std::vector<unsigned char> my_results)
    {
        std::lock_guard<std::mutex> guard(results.res_lock);
        results.results.emplace(this_comp, my_results);
    }

    std::pair<std::array<unsigned char, 8>, dt::diff_tok_t> trans_doub_to_chararr(std::pair<double, dt::diff_tok_t> pres)
    {
        return {*reinterpret_cast<std::array<unsigned char, 8> *>(&pres.first), pres.second};
    }

    int base_experiment::collect_results(supp::result_box &slavres, std::vector<std::ofstream> &output_ids, std::vector<std::ofstream> &output_pro)
    {
        std::lock_guard<std::mutex> guard(slavres.res_lock);
        supp::my_hash<dt::auth_id_t, dt::book_id_t> hasher;
        int tmpres = 0;
        for (auto &res : slavres.results)
        {
            output_ids[hasher({res.first.auth, res.first.book}) % number_of_buckets].write((char *)((std::array<unsigned short, supp::ARR_SIZ()>)res.first).data(), sizeof(unsigned short int) * supp::ARR_SIZ());
            output_pro[hasher({res.first.auth, res.first.book}) % number_of_buckets].write((char *)res.second.data(), sizeof(unsigned char) * res.second.size());
            dt::frag_id_t fnum = shelf_short.at(res.first.auth).at(res.first.book).size();
            totFra += fnum * res.first.comp_frag_num;
            tmpres += fnum * res.first.comp_frag_num;
        }
        totwri += (int)slavres.results.size();
        nowwri += (int)slavres.results.size();
        slavres.results.clear();

        return tmpres;
    }

    void base_experiment::open_output_files(std::string midpart, std::vector<std::ofstream> &output_vector, bool dumpBuck)
    {
        for (auto oname : output_file_names(midpart))
        {
            output_vector.push_back(std::ofstream(oname, std::ofstream::out | std::ofstream::app | std::ofstream::binary));
            if (output_vector.back().fail() || !output_vector.back().is_open())
            {
                std::cerr << "Error opening: " << oname << std::endl;
                exit(EXIT_FAILURE);
            }
            if (dumpBuck && !resuming)
                output_vector.back().write((char *)&number_of_buckets, sizeof(number_of_buckets));
        }
    }

    void base_experiment::dump_new_params()
    {
        for (auto &par : params)
            if (known_par.find(par.first) == known_par.end())
            {
                known_par.emplace(par.first);
                outpar.write((char *)&par.first.first, sizeof(par.first.first));
                outpar.write((char *)&par.first.second, sizeof(par.first.second));
                outpar.write((char *)&par.second, sizeof(par.second));
            }
    }

    void base_experiment::lap_update(supp::Timer &timer_lap, supp::Timer &timer_tot, std::vector<std::ofstream> &output_ids, std::vector<std::ofstream> &output_pro, size_t totProbs)
    {
        size_t resCount = 0;
        for (auto &slavres : results)
        {
            resCount += collect_results(slavres, output_ids, output_pro);
        }
#ifndef _FULL_AUTHOR_PARAMETERS_
        dump_new_params();
#endif

        supp::progress_bar(timer_tot.elapsed(), nowwri, totwri, authList.remaining(), resCount, timer_lap.elapsed(), totProbs, totFra);
        timer_lap.reset();
    }

    void base_experiment::run()
    {
        supp::Timer timer_lap, timer_tot;
        std::vector<std::ofstream> output_ids, output_pro;
        size_t totProbs = authList.remaining();

        open_output_files("_uo_ids", output_ids, true);
        open_output_files("_uo_pro", output_pro);

        results = std::vector<supp::result_box>(numThreads);
        threadCount = 0;
        for (auto i = 0; i < numThreads; i++)
            workers.push_back(std::thread(&base_experiment::worker, this, std::ref(results[i]), std::ref(authList)));
        std::this_thread::sleep_for(std::chrono::milliseconds(250));

        timer_lap.reset();
        timer_tot.reset();
        while (threadCount)
        {
            std::this_thread::sleep_for(std::chrono::milliseconds(250));
            lap_update(timer_lap, timer_tot, output_ids, output_pro, totProbs);
            if (!supp::FLAG)
                break;
        }
        if (supp::FLAG)
        {
            for (auto &t : workers)
                if (t.joinable())
                    t.join();
            lap_update(timer_lap, timer_tot, output_ids, output_pro, totProbs);
        }

        for (auto i = 0; i < number_of_buckets; i++)
        {
            output_ids[i].close();
            output_pro[i].close();
        }
        outpar.close();
    }

    void base_experiment::worker(supp::result_box &results, supp::list_manager &authList)
    { //
        std::pair<double, double> par;
        supp::task job;
        supp::comp_id this_comp;
        size_t res_len = shelf_short.begin()->second.begin()->second.front()->retsize();
        unsigned char *pos_now;
        std::vector<unsigned char> my_results;
        std::map<dt::book_id_t, std::unique_ptr<bp::book>> newAuthors;
        std::string tmpAuthor;
        std::vector<std::string> tmpSeq;
        dt::auth_id_t my_aut = 0;
        dt::slice_id_t my_slice = -1, new_slice, par_slice;
        threadCount++;

        job = authList.pop(my_aut);
        while (job.aut1 > 0 && job.book > 0)
        {
            this_comp = supp::comp_id(job.aut1, job.aut2, job.book);
            if (written.find(this_comp) != written.end())
            {
                job = authList.pop(my_aut);
                continue;
            }
            new_slice = which_slice.at({job.aut2, job.book});
            bool author1InOneOfSlices = (my_slice != dt::slice_id_t(-1) && (cake[my_slice].shelf.find(job.aut1) != cake[my_slice].shelf.end())) || (new_slice != dt::slice_id_t(-1) && cake[new_slice].shelf.find(job.aut1) != cake[new_slice].shelf.end()); //(my_slice == dt::slice_id_t(-1) || (cake[my_slice].shelf.find(job.aut1) != cake[my_slice].shelf.end())) || (new_slice == dt::slice_id_t(-1) || cake[new_slice].shelf.find(job.aut1) != cake[new_slice].shelf.end());//
#ifdef _MANUALLY_SETTABLE_PARAMS
            par_lock.lock();
            bool paramsDefinedForSomeReason = params.find({job.aut1, new_slice}) != params.end() || params.find({my_aut, my_slice}) != params.end();
            par_lock.unlock();
#else
            bool paramsDefinedForSomeReason = false;
#endif //_MANUALLY_SETTABLE_PARAMS
            if (job.aut1 != my_aut || (my_slice != new_slice && (author1InOneOfSlices || paramsDefinedForSomeReason)))
            { // if changing author or changing slice (and one of the previous or the latter had books from the author) or comparing with books of the same author I have to recompute
                my_aut = job.aut1;
                my_slice = new_slice;
                if (my_slice != dt::slice_id_t(-1))
                {
#ifdef _MANUALLY_SETTABLE_PARAMS
                    if (cake[my_slice].shelf.find(my_aut) != cake[my_slice].shelf.end() || params.find({my_aut, new_slice}) != params.end())
#else
                    if (cake[my_slice].shelf.find(my_aut) != cake[my_slice].shelf.end())
#endif                                        //_MANUALLY_SETTABLE_PARAMS
                        par_slice = my_slice; // I may have defined parameters for this combination of author and slice even if no books of the author are present, should not treat it as -1
                    else
                        par_slice = -1;
                }
                else
                    par_slice = -1;
                bool authorNotInThisSlice = par_slice != dt::slice_id_t(-1) && cake[par_slice].shelf.find(my_aut) == cake[par_slice].shelf.end(); // needed if params defined for some reason
                try
                {
                    // as I may need the full author (possibly without some book) to compute the parameters I take this as a starting point also for the splitting
                    std::set<dt::book_id_t> bookExcl;
                    for (auto &volume : shelf_long.at(my_aut))
                    {
                        bool authorInThisSliceButThisBookNot = par_slice != dt::slice_id_t(-1) && (cake[par_slice].shelf.find(my_aut) != cake[par_slice].shelf.end() && cake[par_slice].shelf.at(my_aut).find(volume.first) == cake[par_slice].shelf.at(my_aut).end());
                        if (!(par_slice == dt::slice_id_t(-1) || authorNotInThisSlice || authorInThisSliceButThisBookNot))
                            // for every fragment of the book if its author is not in my_slice or this book of the author is not
                            bookExcl.emplace(volume.first);
                    }
                    if (bookExcl.size() == shelf_long.at(my_aut).size())
                    { // Author has one book, just empty results
                        my_results = std::vector<unsigned char>(shelf_short.at(job.aut2).at(job.book).size() * res_len);
                        for (auto dest = my_results.data(); dest < my_results.data() + my_results.size(); dest += res_len)
                            shelf_short.at(job.aut2).at(job.book).front()->log_prob_to_chararr(dest);
                        this_comp.comp_frag_num = 1;
                        emplace_results(results, this_comp, my_results);
                        newAuthors.clear();
                        job = authList.pop(my_aut);
                        continue;
                    }
                    newAuthors = asp::authorSplitter(shelf_long, my_aut, bookExcl, F);
                    if (!newAuthors.size())
                    { // Author is empty, just empty results
                        my_results = std::vector<unsigned char>(shelf_short.at(job.aut2).at(job.book).size() * res_len);
                        for (auto dest = my_results.data(); dest < my_results.data() + my_results.size(); dest += res_len)
                            shelf_short.at(job.aut2).at(job.book).front()->log_prob_to_chararr(dest);
                        this_comp.comp_frag_num = 1;
                        emplace_results(results, this_comp, my_results);
                        job = authList.pop(my_aut);
                        continue;
                    }
#ifdef _USE_JOINT_PARAMS_
                    par_lock.lock();
                    if (params.find({my_aut, par_slice}) == params.end())
                    {
                        if (computing_par.find({my_aut, par_slice}) != computing_par.end())
                        {
                            my_aut = 0;
                            authList.place_back(job);
                            par_lock.unlock();
                            job = authList.pop(my_aut);
                            continue;
                        }
                        computing_par.emplace(std::make_pair(my_aut, par_slice));
                        par_lock.unlock();
                        bp::book tmpSelf;
                        for (auto &au : newAuthors)
                            tmpSelf.append(*au.second);
                        par = popt::param_opt(tmpSelf);
                        par_lock.lock();
                        params.emplace(std::make_pair(std::make_pair(my_aut, par_slice), par));
                    }
                    for (auto &au : newAuthors)
                    {
                        au.second->add_par(params.at({my_aut, par_slice}));
                        au.second->add_P0(&P0);
                    }
                    par_lock.unlock();
#else
                    bool goodPars = true;
                    for (auto &au : newAuthors)
                    {
                        par_lock.lock();
                        if (params.find({au.first, par_slice}) == params.end())
                        {
                            if (computing_par.find({au.first, par_slice}) != computing_par.end())
                            {
                                my_aut = 0;
                                authList.place_back(job);
                                par_lock.unlock();
                                job = authList.pop(my_aut);
                                goodPars = false;
                                break;
                            }
                            computing_par.emplace(std::make_pair(au.first, par_slice));
                            par_lock.unlock();
                            par = popt::param_opt(au.second);
                            par_lock.lock();
                            params.emplace(std::make_pair(std::make_pair(au.first, par_slice), par));
                        }
                        au.second.add_par(params.at({au.first, par_slice}));
                        au.second.add_P0(&P0);

                        par_lock.unlock();
                    }
                    if (!goodPars)
                        continue;
#endif // _USE_JOINT_PARAMS_
                }
                catch (std::runtime_error &e)
                {
                    std::cerr << "In making corpus. "
                              << " " << my_aut << " " << e.what() << std::endl;
                }
            }
            else
            {
                my_slice = new_slice;
            }
            if (!newAuthors.size())
            { // Author is empty, just empty results
                my_results = std::vector<unsigned char>(shelf_short.at(job.aut2).at(job.book).size() * res_len);
                for (auto dest = my_results.data(); dest < my_results.data() + my_results.size(); dest += res_len)
                    shelf_short.at(job.aut2).at(job.book).front()->log_prob_to_chararr(dest);
                this_comp.comp_frag_num = 1;
                emplace_results(results, this_comp, my_results);
                job = authList.pop(my_aut);
                continue;
            }
            this_comp.comp_frag_num = newAuthors.size();
            // Compare all the fragments of the author
            my_results = std::vector<unsigned char>(shelf_short.at(job.aut2).at(job.book).size() * this_comp.comp_frag_num * res_len);
            pos_now = my_results.data();
            for (auto &fragm1 : newAuthors)
            {
                for (auto &fragm2 : shelf_short.at(job.aut2).at(job.book))
                {
                    try
                    {
                        fragm1.second->log_prob_to_chararr(*fragm2, pos_now);
                    }
                    catch (bp::word_miss &e)
                    {
                        std::cerr << "Word miss in compare. " << job.aut2 << " " << job.book << " " << e.what() << std::endl;
                        fragm1.second->log_prob_to_chararr(pos_now);
                    }
                    catch (bp::not_init &e)
                    {
                        std::cerr << "Not init in compare. " << job.aut2 << " " << job.book << " " << e.what() << std::endl;
                        fragm1.second->log_prob_to_chararr(pos_now);
                    }
                    pos_now += res_len;
                }
            }
            emplace_results(results, this_comp, my_results);
            job = authList.pop(my_aut);
        }
        threadCount--;
    }

    void base_experiment::order()
    {
        std::ofstream out_ids, out_pro;
        std::ifstream input;
        supp::Semaphore task_limiter(numThreads);
        std::vector<std::future<void>> tasks;

        std::cout << "\nOrdering" << std::endl;
        for (int i = 0; i < number_of_buckets; i++)
        {
            tasks.emplace_back(std::async(std::launch::async, [i, this, &task_limiter]
                                          {
                supp::Critical_section __(task_limiter);
                results_orderer(i); }));
        }
        for (auto &&task : tasks)
            task.wait();

        if (number_of_buckets <= numThreads)
        {
            std::cout << "Joining" << std::endl;
            out_ids.open(outputFile + "_ids.res", std::ofstream::out | std::ofstream::binary);
            out_pro.open(outputFile + "_pro.res", std::ofstream::out | std::ofstream::binary);
            if (out_ids.fail() || !out_ids.is_open() || out_pro.fail() || !out_pro.is_open())
            {
                std::cerr << "Error opening: " << outputFile << "_???" << std::endl;
                exit(EXIT_FAILURE);
            }
            for (int i = 0; i < number_of_buckets; i++)
            {
                if (fs::is_regular_file(supp::insert_padded_number(outputFile + "_ids", ".res", i)) && fs::file_size(supp::insert_padded_number(outputFile + "_ids", ".res", i)))
                {
                    input.open(supp::insert_padded_number(outputFile + "_ids", ".res", i), std::ifstream::binary);
                    auto qwe = input.rdbuf();
                    out_ids << qwe << std::flush;
                    input.close();
                    input.open(supp::insert_padded_number(outputFile + "_pro", ".res", i), std::ifstream::binary);
                    out_pro << input.rdbuf();
                    input.close();
                }
                std::remove(supp::insert_padded_number(outputFile + "_ids", ".res", i).c_str());
                std::remove(supp::insert_padded_number(outputFile + "_pro", ".res", i).c_str());
            }
            out_ids.close();
            out_pro.close();
        }
    }

    void base_experiment::results_orderer(int target)
    {
        std::unordered_map<std::pair<dt::auth_id_t, dt::book_id_t>, std::map<dt::auth_id_t, std::pair<unsigned int, unsigned int>>, supp::my_hash<dt::auth_id_t, dt::book_id_t>> pieces;
        const unsigned int recorLen = shelf_short.begin()->second.begin()->second.front()->retsize();
        std::ifstream input;
        std::ofstream out_ids, out_pro;
        unsigned int offNow = 0;
        try
        {
            input.open(supp::insert_padded_number(outputFile + "_uo_ids", ".res", target), std::ofstream::binary);
            if (input)
            {
                input.seekg(sizeof(number_of_buckets)); // to skip the bytes storing the number of buckets
                std::array<unsigned short, supp::ARR_SIZ()> tmp;
                while (!input.eof())
                {
                    input.read((char *)&tmp, sizeof(unsigned short int) * supp::ARR_SIZ());
                    if (!input.eof())
                    {
                        supp::comp_id tmp_id(tmp);
                        dt::frag_id_t no_of_fragm = 0;
                        dt::book_id_t no_of_auth_pieces = 0;
                        try
                        {
                            no_of_fragm = (flagbyte & nofrag) ? 1 : shelf_short.at(tmp_id.auth).at(tmp_id.book).size();
                        }
                        catch (std::out_of_range &e)
                        {
                            std::ostringstream oss;
                            oss << std::this_thread::get_id() << " no_of_fragm " << e.what() << std::endl;
                            if (shelf_short.find(tmp_id.auth) == shelf_short.end())
                                oss << "Missing author.\n";
                            if (shelf_short.at(tmp_id.auth).find(tmp_id.book) == shelf_short.at(tmp_id.auth).end())
                                oss << "Missing book.\n";
                            oss << "Comp\nA: " << tmp_id.auth << " B: " << tmp_id.book << " CA: " << tmp_id.comp_aut << "\n";
                            oss << "Target: " << target << " Pieces: " << pieces.size() << " Off: " << offNow;
                            std::cerr << oss.str() << std::endl;
                            throw e;
                        }
                        try
                        {
                            no_of_auth_pieces = (flagbyte & noauth) ? 1 : (((tmp_id.auth == 0) || (cake[which_slice[{tmp_id.auth, tmp_id.book}]].shelf.find(tmp_id.comp_aut) == cake[which_slice[{tmp_id.auth, tmp_id.book}]].shelf.end())) ? main_parts.at(tmp_id.comp_aut) : slice_author_size.at(tmp_id.comp_aut).at(which_slice[{tmp_id.auth, tmp_id.book}]));
                        }
                        catch (std::out_of_range &e)
                        {
                            std::cerr << std::this_thread::get_id() << " no_of_auth_pieces " << e.what() << std::endl;
                        }
                        pieces[{tmp_id.auth, tmp_id.book}][tmp_id.comp_aut] = {offNow, no_of_fragm * no_of_auth_pieces};
                        offNow += no_of_fragm * no_of_auth_pieces;
                    }
                }
                input.close();
            }
            else
                throw std::runtime_error("Missing file " + supp::insert_padded_number(outputFile + "_uo_ids", ".res", target));
            input.open(supp::insert_padded_number(outputFile + "_uo_pro", ".res", target), std::ifstream::binary);
            out_ids.open(supp::insert_padded_number(outputFile + "_ids", ".res", target), std::ofstream::out | std::ofstream::app | std::ofstream::binary);
            out_pro.open(supp::insert_padded_number(outputFile + "_pro", ".res", target), std::ofstream::out | std::ofstream::app | std::ofstream::binary);
            if (input && out_ids && out_pro)
            {
                for (auto &ab : pieces)
                {
                    unsigned int tot_this = 0;
                    for (auto &ca : ab.second)
                    {
                        int bytes_to_read = ca.second.second * recorLen;
                        input.seekg(ca.second.first * recorLen);
                        char *memblock = new char[bytes_to_read];
                        input.read(memblock, bytes_to_read);
                        out_pro.write(memblock, bytes_to_read);
                        delete[] memblock;
                        tot_this += ca.second.second;
                    }
                    std::array<char, 24> buffer;
                    unsigned short used = sizeof(ab.first.first) + sizeof(ab.first.second);
                    memcpy(buffer.data(), &ab.first.first, sizeof(ab.first.first));
                    memcpy(buffer.data() + sizeof(ab.first.first), &ab.first.second, sizeof(ab.first.second));
                    memcpy(buffer.data() + used, &tot_this, sizeof(tot_this));
                    used += sizeof(tot_this);
                    out_ids.write((char *)buffer.data(), used);
                }
                input.close();
                out_ids.close();
                out_pro.close();
            }
            else
                throw std::runtime_error("Missing some file");
            std::remove(supp::insert_padded_number(outputFile + "_uo_ids", ".res", target).c_str());
            std::remove(supp::insert_padded_number(outputFile + "_uo_pro", ".res", target).c_str());
        }
        catch (std::exception &e)
        {
            std::cerr << std::this_thread::get_id() << " " << e.what() << std::endl;
            throw e;
        }
    }

    void base_experiment::spread_full_author_parameters()
    {
        std::vector<dt::auth_id_t> missing_parameters;
        supp::Semaphore task_limiter(numThreads);
        std::vector<std::future<std::pair<dt::auth_id_t, std::pair<double, double>>>> tasks;

        for (auto aut : shelf_long)
            if (params.find({aut.first, -1}) == params.end())
                missing_parameters.push_back(aut.first);

        for (auto aut : missing_parameters)
        {
            tasks.emplace_back(std::async(std::launch::async, [this, &task_limiter, aut]() -> std::pair<dt::auth_id_t, std::pair<double, double>>
                                          {       
                                            supp::Critical_section __(task_limiter);
                                            bp::book tmpSelf;
                                            for (auto &au : shelf_short.at(aut)) for (auto &bk : au.second)
                                                tmpSelf.append(*bk);
                                            return std::make_pair(aut, popt::param_opt(tmpSelf)); }));
        }
        for (auto &&task : tasks)
        {
            std::pair<dt::auth_id_t, std::pair<double, double>> tmp;
            tmp = task.get();
            params.emplace(std::make_pair(std::make_pair(tmp.first, -1), tmp.second));
        }
        dump_new_params();
        for (auto aut : shelf_long)
            for (dt::slice_id_t i = 0; i < cake.size(); i++)
                if (!params.emplace(std::make_pair(std::make_pair(aut.first, i), params.at({aut.first, -1}))).second)
                    throw bp::word_doubled("Same parameters twice: " + std::to_string(aut.first) + " " + std::to_string(i));
    }
}