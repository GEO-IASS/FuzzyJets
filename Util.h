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
    extern const std::vector<int> TBSs;
    extern const std::vector<int> seed_pT_cuts;
    extern const std::vector<int> seed_noises;

    std::vector<Color_t> get_colors(unsigned int n);

    std::string format_float(const float x, const int digits);

    std::string time_generated(std::string process, int NPV, int EJW, int EJO,
                               int PP, int TBS, int seed_pT_cut, int seed_noise);

    std::string filename_base(std::string process, int NPV, int EJW, int EJO,
                              int PP, int TBS, int seed_pT_cut, int seed_noise);

    std::string filename(std::string process, int NPV, int EJW, int EJO,
                         int PP, int TBS, int seed_pT_cut, int seed_noise);

    std::string filename_w(std::string process, int NPV, int EJW, int EJO,
                           int PP, int TBS, int seed_pT_cut, int seed_noise);

    std::string eventname(std::string process, int NPV, int EJW, int EJO,
                          int PP, int TBS, int seed_pT_cut, int seed_noise);

    std::string fancy_var_name(std::string var_name);

    std::string sanitize(std::string s);

    std::string latexed_pT_range(float low, float high);

    class ParameterSet {
    public:
        std::string _process;
        int _NPV;
        int _EJW;
        int _EJO;
        int _PP;
        int _TBS;
        int _seed_pT_cut;
        int _seed_noise;

        ParameterSet(std::string process, int NPV, int EJW, int EJO,
                     int PP, int TBS, int seed_pT_cut, int seed_noise) :
            _process(process), _NPV(NPV), _EJW(EJW), _EJO(EJO),
            _PP(PP), _TBS(TBS), _seed_pT_cut(seed_pT_cut), _seed_noise(seed_noise) {}
        ParameterSet();

        std::string Print();
        std::string PrintExceptDistinguished(std::string distinguished);
        std::string EventName();
    };

    std::string fancy_distinguished_param_name(util::ParameterSet ps, std::string distinguished_var_type);
}


#endif
