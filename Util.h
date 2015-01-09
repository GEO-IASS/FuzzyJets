#ifndef UTIL_H
#define UTIL_H

#include <TStyle.h>

#include <string>
#include <vector>

namespace util {
    extern const std::vector<Color_t> colors;
    extern const std::vector<int> NPVs;
    extern const std::vector<std::string> processes;
    extern const std::vector<int> EJWs;
    extern const std::vector<int> EJOs;
    extern const std::vector<int> PPs;
    extern const std::vector<int> seed_pT_cuts;

    std::vector<Color_t> get_colors(unsigned int n);

    std::string format_float(const float x, const int digits);

    std::string time_generated(std::string process, int NPV, int EJW, int EJO,
                               int PP, int seed_pT_cut);

    std::string filename_base(std::string process, int NPV, int EJW, int EJO,
                              int PP, int seed_pT_cut);

    std::string filename(std::string process, int NPV, int EJW, int EJO,
                         int PP, int seed_pT_cut);

    std::string filename_w(std::string process, int NPV, int EJW, int EJO,
                           int PP, int seed_pT_cut);

    std::string eventname(std::string process, int NPV, int EJW, int EJO,
                          int PP, int seed_pT_cut);
}

#endif
