#include "lib/paramOpt.hpp"
#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <map>

#define LEN 1000000
int main(int argc, char **argv)
{
    std::map<std::string, double> args;
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        size_t pos = arg.find('=');
        if (pos != std::string::npos)
        {
            std::string key = arg.substr(0, pos);
            double val = std::stod(arg.substr(pos + 1));
            args[key] = val;
        }
    }
    double alpha = args.count("alpha") ? args["alpha"] : 0.4;
    double theta = args.count("theta") ? args["theta"] : 50;
    std::vector<double> Dt;
    Dt.reserve(LEN);
    int D = 0;
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0., 1.0);

    std::cout << "Start simulation with:" << alpha << ", " << theta << std::endl;
    std::vector<int> counts(LEN, 0);
    std::vector<double> cumsum(LEN);
    std::ptrdiff_t last = 0;
    double n = 0;
    for (auto lim : std::vector<double>{10000, 100000, 1000000})
    {
        for (; n < lim; n++)
        {
            auto r = dis(gen) - (theta + D * alpha) / (theta + n);
            if (r < 0)
            {
                counts[D]++;
                last = D++;
            }
            else
            {
                for (auto ptr = cumsum.end() - last - 1; ptr < cumsum.end(); ptr++)
                    *ptr = *(counts.begin() + (cumsum.end() - ptr - 1)) - alpha + *(ptr - 1);
                auto ind = std::upper_bound(cumsum.end() - D, cumsum.end(), r * (theta + n));
                auto last = cumsum.end() - ind - 1;
                counts[last]++;
            }
            Dt.push_back(D);
        }
        std::cout << alpha << " " << theta << " " << lim << " " << D;
        auto res1 = popt::param_opt_fit(Dt, alpha - 0.05, theta - 2);
        std::cout << " " << res1.first << " " << res1.second;
        std::vector<int> mults = std::vector<int>(counts.begin(), counts.begin() + D);
        auto res2 = popt::param_opt_ML(mults, alpha - 0.05, theta - 2);
        std::cout << " " << res2.first << " " << res2.second;
        auto res3 = popt::param_opt_fit2D(Dt, alpha - 0.05, theta - 2);
        std::cout << " " << res3.first << " " << res3.second << std::endl;
    }
    return 0;
}