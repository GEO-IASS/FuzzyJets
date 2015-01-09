#include <string>
#include <sstream>

#include <cassert>
#include <TStyle.h>

#include "Util.h"

namespace util {
    const std::vector<Color_t> colors = {kBlack, kBlue, kRed, kGreen+1, kMagenta+1, kOrange+1, kYellow, kGray, kTeal};
    const std::vector<int> NPVs {0, 20, 40};
    const std::vector<std::string> processes {"qcd"}; // run one at a time
    const std::vector<int> EJWs {0, 10, 20, 40, 100};
    const std::vector<int> EJOs {-2, 0, 2};
    const std::vector<int> PPs {0, 1};
    const std::vector<int> seed_pT_cuts {5, 10, 15, 20};


    std::vector<Color_t> get_colors(unsigned int n) {
        assert(n <= util::colors.size());
        std::vector<Color_t>::const_iterator s = util::colors.begin();
        std::vector<Color_t>::const_iterator f = util::colors.begin() + n;
        std::vector<Color_t> out(s, f);
        return out;
    }

    std::string format_float(const float x, const int digits) {
        std::stringstream ss;
        ss.str(std::string());
        ss << std::fixed;
        ss.precision(digits);
        ss << x;
        return ss.str();
    }

    std::string time_generated(std::string process,
                               __attribute__((unused)) int NPV,
                               __attribute__((unused)) int EJW,
                               __attribute__((unused)) int EJO,
                               int PP,
                               __attribute__((unused)) int seed_pT_cut) {
        if (process == "qcd") {
            if (PP == 0) {
                return "2015_01_05_23h47m01s";
            } else {
                return "2015_01_06_14h43m15s";
            }
        } else if(process == "zprime") {
            if (PP == 0) {
                return "2015_01_07_16h38m28s";
            } else {
                return "2015_01_07_23h57m27s";
            }
        }
        return "";
    }

    std::string dataset_name (std::string process,
                              __attribute__((unused)) int NPV,
                              __attribute__((unused)) int EJW,
                              __attribute__((unused)) int EJO,
                              int PP,
                              __attribute__((unused)) int seed_pT_cut) {
        std::stringstream ss;
        ss.str(std::string());
        ss << "2500evts_" << process << ((PP == 0) ? "_ndm_" : "_ydm_");
        if (process == "qcd") {
            // nothing
        } else {
            ss << "ntbs_";
        }
        ss << "offset_test";
        return ss.str();
    }

    std::string filename_base(std::string process, int NPV, int EJW, int EJO,
                                  int PP, int seed_pT_cut) {
        std::string base = "/u/at/chstan/nfs/summer_2014/ForConrad/files/";
        std::stringstream ss;
        ss.str(std::string());
        ss << base
           << util::dataset_name(process, NPV, EJW, EJO, PP, seed_pT_cut)
           << "/"
           << util::time_generated(process, NPV, EJW, EJO, PP, seed_pT_cut)
           << "/10s_" << NPV << "mu_0lw_0rec_"
           << seed_pT_cut << "cut_"
           << EJW << "ejw_"
           << EJO << "off_"
           << PP << "PP_0TBS";
        return ss.str();
    }

    std::string filename(std::string process, int NPV, int EJW, int EJO,
                         int PP, int seed_pT_cut) {
        return util::filename_base(process, NPV, EJW,
                                   EJO, PP, seed_pT_cut) + ".root";
    }

    std::string filename_w(std::string process, int NPV, int EJW, int EJO,
                           int PP, int seed_pT_cut) {
        return util::filename_base(process, NPV, EJW,
                                   EJO, PP, seed_pT_cut) + "_w.root";
    }

    std::string eventname(std::string process, int NPV, int EJW, int EJO,
                          int PP, int seed_pT_cut) {
        std::stringstream ss;
        ss.str(std::string());
        ss << process << "_" << NPV << "_" << EJW << "_" << EJO << "_"
           << PP << "_" << seed_pT_cut;
        return ss.str();
    }
}
