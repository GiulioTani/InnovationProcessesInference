#include "lib/attributor_Ctypes.hpp"
#include <map>
#include <vector>
#include <istream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <limits>
#include <cmath>
#include <cstring>
#include <sstream>
#include <iomanip>

#define INF std::numeric_limits<double>::infinity()

/**
 * @brief Loads the probabilities and number of new tokens per fragment to allow attribution.
 * 
 * The input data are the 'pro' files produced by base_experiment.
 * Format of the loaded data (output):
 * - 1 double: number of authors (n);
 * - 2*n doubles: pairs (author id, number of author parts)
 * - for every author part:
 * -  for every book fragment:
 * -   2 doubles: probability per token of the fragment (with delta =1), fraction of new tokens.
 * 
 * Return codes:
 * - 0, no error.
 * - 1, failed opening file.
 * - 2, file ended before expected.
 * - 3, other I/O error.
 * - 4, uncaught exception.
 * 
 * @param data Structure containing the information for data loading.
 * @param results Pointer to the memory that will contain the data for attribution with attribute().
 * @return int Return code for error notification to python (positive numbers) or the number of authors (negative numbers).
 */
int load(struct dataLoad data, double **results)
{
    try
    {
        std::ifstream inputFile;
        if (data.bucket < 1000)/*results in multiple buckets, recover bucket file name.*/
        {
            std::stringstream padnum;
            padnum.clear();
            padnum.str("");
            padnum << std::setfill('0') << std::setw(3) << data.bucket;
            inputFile.open(data.fname + padnum.str() + ".res", std::ifstream::binary);
        }
        else /*results in a single file.*/
            inputFile.open(data.fname + std::string(".res"), std::ifstream::binary);
        std::vector<double> pro = std::vector<double>(data.ncids * data.fn, 0);
        std::vector<double> wc = std::vector<double>(data.ncids * data.fn, 0);
        const int recorLen = sizeof(double) + (data.F == 1 ? sizeof(double) : sizeof(dt::diff_tok_t)); //records may contain the probability of the fragment and the number of new tokens or, when using single token fragments, the fraction of new tokens.
        int infcheck, oldcid = -1, sgood = 0, length = data.N * recorLen;
        bool orderauth = true;
        double *_results;
        char *buffer;
        if (!inputFile.is_open())
            return 1;
        inputFile.seekg(data.offs * recorLen);
        buffer = new char[length];
        inputFile.read(buffer, length);
        if (inputFile.eof())
            return 2;
        if (!inputFile.good())
            return 3;
        if (data.F == 1)
            for (auto i = 0; i < data.N; i++)
            {
                pro[i] = *reinterpret_cast<double *>(buffer + i * recorLen);
                wc[i] = *reinterpret_cast<double *>(buffer + i * recorLen + sizeof(pro.front()));
            }
        else
            for (auto i = 0; i < data.N; i++)
            {
                pro[i] = *reinterpret_cast<double *>(buffer + i * recorLen);
                wc[i] = *reinterpret_cast<dt::diff_tok_t *>(buffer + i * recorLen + sizeof(pro.front()));
            }
        delete[] buffer;
        for (auto i = 0; i < data.ncids; i++)
            sgood += data.goods[i];
        _results = (double *)malloc(sizeof(double) * (sgood * (data.fn + 1) + 1) * 2);
        *results = _results;
        for (auto i = 0; i < data.ncids; i++)
        {
            if (data.goods[i])
            {
                if ((int)data.cids[i] < oldcid)
                {
                    orderauth = false;
                    break;
                }
                oldcid = data.cids[i];
            }
        }
        if (orderauth)
        {
            int cpos = 0, num_auts, aut_parts = 0;
            _results[0] = -1;
            for (auto i = 0; i < data.ncids; i++)
            {
                if (data.goods[i])
                {
                    if (data.cids[i] == _results[2 * cpos])
                    {
                        _results[2 * cpos + 1]++;
                    }
                    else
                    {
                        cpos++;
                        _results[2 * cpos] = data.cids[i];
                        _results[2 * cpos + 1] = 1;
                    }
                }
            }
            _results[0] = num_auts = cpos;
            cpos = 1;
            for (auto i = 0; i < data.ncids; i++)
            {
                if (data.goods[i])
                {
                    cpos += _results[2 * cpos] != data.cids[i];
                    for (infcheck = 0; pro[i * data.fn + infcheck] == -INF && infcheck < data.fn; infcheck++)
                        ;
                    if (infcheck == data.fn && pro[i * data.fn + infcheck - 1] == -INF)
                    {
                        _results[2 * cpos + 1]--;
                        continue;
                    }
                    if (data.F == 1)
                        for (auto j = 0; j < data.fn; j++)
                        {
                            _results[2 * (num_auts + 1 + aut_parts * data.fn + j)] = pro[i * data.fn + j];
                            _results[2 * (num_auts + 1 + aut_parts * data.fn + j) + 1] = wc[i * data.fn + j];
                        }
                    else
                        for (auto j = 0; j < data.fn; j++)
                        {
                            _results[2 * (num_auts + 1 + aut_parts * data.fn + j)] = pro[i * data.fn + j] / data.lengths[j];
                            _results[2 * (num_auts + 1 + aut_parts * data.fn + j) + 1] = wc[i * data.fn + j] / data.lengths[j];
                        }
                    aut_parts++;
                }
            }
            _results[1] = aut_parts;
        }
        else
        {
            int cpos = 0, num_auts, aut_parts = 0;
            std::map<int, std::pair<std::vector<double>, std::vector<double>>> baseres;
            for (auto i = 0; i < data.ncids; i++)
            {
                if (data.goods[i])
                {
                    for (infcheck = 0; pro[i * data.fn + infcheck] == -INF && infcheck < data.fn; infcheck++)
                        ;
                    if (infcheck == data.fn && pro[i * data.fn + infcheck - 1] == -INF)
                        continue;
                    baseres.emplace_hint(baseres.end(), std::make_pair(data.cids[i], std::pair<std::vector<double>, std::vector<double>>()));
                    std::copy(pro.begin() + (i * data.fn), pro.begin() + ((i + 1) * data.fn), std::back_inserter(baseres.at(data.cids[i]).first));
                    std::copy(wc.begin() + (i * data.fn), wc.begin() + ((i + 1) * data.fn), std::back_inserter(baseres.at(data.cids[i]).second));
                }
            }
            _results[0] = num_auts = baseres.size();
            for (auto &cid : baseres)
            {
                _results[2 * ++cpos] = cid.first;
                _results[2 * cpos + 1] = cid.second.first.size();
                for (auto i = 0; i < (int)cid.second.first.size(); i++)
                {
                    _results[2 * (num_auts + 1 + aut_parts + i)] = cid.second.first[i];
                    _results[2 * (num_auts + 1 + aut_parts + i) + 1] = cid.second.second[i];
                }
                aut_parts += (int)cid.second.first.size() / data.fn;
            }
            _results[1] = aut_parts;
        }
        return 0;
    }
    catch (std::exception &e)
    {
        std::cerr << e.what() << std::endl;
        return 4;
    }
}

/**
 * @brief Attributes books to authors.
 * Given the results loaded with load() and a value of delta assign the books to the authors in 2 to 6 different ways.
 * 
 * Format of the output in procssed:
 * It contains the results for every author in sequence, ordered by descending FNN probability. For every author it contains:
 * - dt::auth_id_t, the id number of the author.
 * - 3*double, the probability according to FNN (or Maximum Likelihood), TOP (the best of the author slices), and WP (weighted average of the author slices).
 * - 3*int, the number of fragments assigned to the author according to MR, TMR (TOP majority rule), and WMR (Weighted profile majority rule).
 * 
 * Return codes:
 * - 4, uncaught exception.
 * - -n, minus the number of authors to help parsing the output.
 * 
 * @param results Pointer to the results loaded with load(), see for the structure.
 * @param fn Number of fragments of the book.
 * @param delta Value of delta for the attribution.
 * @param processed Address for the results of the attribution.
 * @return int Return code for error notification to python (positive numbers) or the number of authors (negative numbers).
 */
int attribute(double *results, long fn, double delta, char **processed)
{
    try
    {
        int num_auts = (int)results[0], aut_parts = (int)results[1], rpos = 0;
        std::vector<std::pair<int, double>> MRM_mean(fn, {-1, -INF}), MRM_max(fn, {-1, -INF}), MRM_wp(fn, {-1, -INF});
        std::vector<double> perares_d(num_auts * 3), partial(aut_parts * fn);
        std::vector<int> perares_i(num_auts * 3);
        std::vector<dt::auth_id_t> authOrder;
        std::map<int, int> authplaces;
        double *nowpart;
        char *_processed;
        authOrder.reserve(num_auts);
        if (!aut_parts)
        { // empty book all -inf
            double dumb = -INF;
            short int n = 2, p = 0;
            int d = 1;
            _processed = (char *)malloc((size_t)(38 * 2));
            *processed = _processed;
            *_processed = *reinterpret_cast<char *>(&n);
            for (short int i = 0; i < n; i++)
            {
                memcpy(_processed + p, &i, sizeof(i));
                p += sizeof(i);
                for (auto j = 0; j < 3; j++)
                {
                    memcpy(_processed + p, &dumb, sizeof(dumb));
                    p += sizeof(dumb);
                }
                for (auto j = 0; j < 3; j++)
                {
                    memcpy(_processed + p, &d, sizeof(d));
                    p += sizeof(d);
                }
            }
            return -n;
        }
        for (auto i = 0; i < aut_parts * fn; i++)
            partial[i] = results[2 * (num_auts + 1 + i)] + results[2 * (num_auts + 1 + i) + 1] * delta;

        for (auto cpos = 1; cpos <= num_auts; cpos++)
        {
            nowpart = partial.data() + rpos;
            int thisCidLen = results[2 * cpos + 1];
            if (!thisCidLen)
                continue;
            std::vector<double> prosd(thisCidLen), ML_tmax(fn), ML_twp(fn), ML_tmean(fn);
            for (auto fnow = 0; fnow < fn; fnow++)
            {
                int infcount = 0;
                for (auto auth = 0; auth < thisCidLen; auth += 1)
                {
                    prosd[auth] = nowpart[fn * auth + fnow];
                    infcount += (prosd[auth] == -INF);
                }

                if (thisCidLen == 1)
                {
                    ML_tmean[fnow] = infcount ? -INF : prosd[0];
                    if (prosd[0] > MRM_mean[fnow].second)
                        MRM_mean[fnow] = {results[2 * cpos], prosd[0]};
                    if (prosd[0] > MRM_max[fnow].second)
                        MRM_max[fnow] = {results[2 * cpos], prosd[0]};
                    if (prosd[0] > MRM_wp[fnow].second)
                        MRM_wp[fnow] = {results[2 * cpos], prosd[0]};
                }
                else
                {
                    std::stable_sort(prosd.begin(), prosd.end());
                    ML_tmax[fnow] = prosd.back();
                    if (ML_tmax[fnow] > MRM_max[fnow].second)
                        MRM_max[fnow] = {results[2 * cpos], ML_tmax[fnow]};

                    double wp = 0, w = 0, j = 1;
                    for (auto i = (int)prosd.size(); i-- > infcount; j++)
                    {
                        ML_tmean[fnow] += prosd[i];
                        wp += prosd[i] / j;
                        w += 1 / j;
                    }
                    ML_tmean[fnow] /= thisCidLen - infcount;
                    if (ML_tmean[fnow] > MRM_mean[fnow].second)
                        MRM_mean[fnow] = {results[2 * cpos], ML_tmean[fnow]};

                    ML_twp[fnow] = wp / w;
                    if (ML_twp[fnow] > MRM_wp[fnow].second)
                        MRM_wp[fnow] = {results[2 * cpos], ML_twp[fnow]};
                }
            }
            double mlmax = 0, mlwp = 0, mlmean = 0;
            for (auto p : ML_tmean)
                mlmean += p;
            perares_d[(cpos - 1) * 3] = mlmean / fn;
            if (thisCidLen == 1)
            {
                perares_d[(cpos - 1) * 3 + 1] = perares_d[(cpos - 1) * 3 + 2] = perares_d[(cpos - 1) * 3];
            }
            else
            {
                for (auto p : ML_tmax)
                    mlmax += p;
                perares_d[(cpos - 1) * 3 + 1] = mlmax / fn;
                for (auto p : ML_twp)
                    mlwp += p;
                perares_d[(cpos - 1) * 3 + 2] = mlwp / fn;
            }
            rpos += fn * thisCidLen;
        }
        for (auto i = 0; i < num_auts; i++)
        {
            if (results[2 * (i + 1) + 1] > 0)
            {
                authOrder.push_back(results[2 * (i + 1)]);
                authplaces[results[2 * (i + 1)]] = i;
            }
        }
        std::stable_sort(authOrder.begin(), authOrder.end(), [&perares_d, &authplaces](size_t i1, size_t i2)
                         { return perares_d[authplaces.at(i1) * 3] < perares_d[authplaces.at(i2) * 3]; });

        for (auto aid : MRM_mean)
            perares_i[authplaces.at(aid.first) * 3]++;
        for (auto aid : MRM_max)
            perares_i[authplaces.at(aid.first) * 3 + 1]++;
        for (auto aid : MRM_wp)
            perares_i[authplaces.at(aid.first) * 3 + 2]++;
        int recordSize = sizeof(authOrder[0]) + 3 * sizeof(perares_d[0]) + 3 * sizeof(perares_i[0]);
        _processed = (char *)malloc(recordSize * authOrder.size());
        *processed = _processed;
        int p = 0, bs;
        for (auto i = 0; i < (int)authOrder.size(); i++)
        {
            memcpy(_processed + p, authOrder.data() + i, sizeof(authOrder[0]));
            p += sizeof(authOrder[0]);
            bs = authplaces.at(authOrder[i]) * 3;
            for (auto j = 0; j < 3; j++)
            {
                memcpy(_processed + p, perares_d.data() + bs + j, sizeof(perares_d[0]));
                p += sizeof(perares_d[0]);
            }
            for (auto j = 0; j < 3; j++)
            {
                memcpy(_processed + p, perares_i.data() + bs + j, sizeof(perares_i[0]));
                p += sizeof(perares_i[0]);
            }
        }
        return -(int)authOrder.size();
    }
    catch (std::exception &e)
    {
        std::cerr << e.what() << std::endl;
        return 4;
    }
}
/**
 * @brief Simple interface to free() for cleaning the memory.
 * 
 * @param data Pointer to the memory to free.
 */
void clean(void *data)
{
    free(data);
}