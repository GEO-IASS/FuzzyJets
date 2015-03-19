#include <string>
#include <sstream>
#include <map>

#include <cassert>
#include <TStyle.h>

#include "Util.h"

#include <boost/algorithm/string.hpp>

namespace util {
    const std::vector<Color_t> colors = {kBlue, kRed, kBlack, kGreen+1, kMagenta+1, kOrange+1, kYellow, kGray, kTeal};
    const std::vector<int> NPVs {0};
    const std::vector<std::string> processes {"qcd", "zprime", "wprime"};
    const std::vector<int> EJWs {0};
    const std::vector<int> EJOs {0};
    const std::vector<int> PPs {0};
    const std::vector<int> TBSs {0};
    const std::vector<int> seed_pT_cuts {5};
    const std::vector<int> seed_noises {0, 20, 100};

    std::map<std::string, std::string> fancy_names =
        {
            {"mGMMc_r", "Leading Learned #sigma"},
            {"mGMMc_r*mGMMc_r", "Leading Learned #sigma^{2}"},
            {"mGMMc_m", "Fuzzy Mass [GeV]"},
            {"mGMMc_pt", "Fuzzy p_{T} [GeV]"},
            {"antikt_m", "Anti-k_{t} Mass [GeV]"},
            {"antikt_pt", "Anti-k_{t} p_{T} [GeV]"},
            {"deltatop_mGMMc", "#Delta R to truth top"},
            {"mGMMc_r_avg", "Average Learned #sigma"},
            {"mGMMc_dr", "#DeltaR Between Fuzzy Jets"},
            {"antikt_dr", "#DeltaR Between anti-k_{t} (R=1.0) Jets"},
            {"antikt_area_trimmed_three",
             "Area of Leading anti-k_{t} Trimmed (R=0.3) Jet"},
            {"mGMMc_r_second", "Sub-leading Learned #sigma"},
            {"mGMMc_r_second/mGMMc_r",
             "Relative Size of Second Fuzzy Jet (#sigma_{2}/#sigma_{1})"},

            {"antikt_m/antikt_pt", "Leading anti-k_{t} Jet Mass/Jet p_{T}"},
            {"mGMMc_m/mGMMc_pt", "Leading Fuzzy Jet Mass/Jet p_{T}"},
            {"qcd", "QCD"},
            {"zprime", "Z' #rightarrow t#bar{t}"},
            {"wprime", "W' #rightarrow qq'"},
            {"pt", "p_{T}"},
            {"antikt_nsubjettiness:0", "Anti-k_{t} #tau_{1}"},
            {"antikt_nsubjettiness:1", "Anti-k_{t} #tau_{2}"},
            {"antikt_nsubjettiness:2", "Anti-k_{t} #tau_{3}"}
        };

    std::string fancy_var_name(std::string var_name) {
        auto it = fancy_names.find(var_name);
        if (it == fancy_names.end())
            return var_name;
        return it->second;
    }

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

    std::string time_generated(__attribute__((unused)) std::string process,
                               __attribute__((unused)) int NPV,
                               __attribute__((unused)) int EJW,
                               __attribute__((unused)) int EJO,
                               __attribute__((unused)) int PP,
                               __attribute__((unused)) int TBS,
                               __attribute__((unused)) int seed_pT_cut,
                               __attribute__((unused)) int seed_noise) {
        if (process == "wprime") {
            // only have the 5 cut for this sample
            return "2015_02_25_00h25m23s";
        }
        if (process == "qcd") {
            return "2015_02_25_00h25m43s";
        }
        if (process == "zprime") {
            return "2015_02_25_00h26m07s";
        }
        return "";
    }

    std::string dataset_name (__attribute__((unused)) std::string process,
                              __attribute__((unused)) int NPV,
                              __attribute__((unused)) int EJW,
                              __attribute__((unused)) int EJO,
                              __attribute__((unused)) int PP,
                              __attribute__((unused)) int TBS,
                              __attribute__((unused)) int seed_pT_cut,
                              __attribute__((unused)) int seed_noise) {
        return "100k_" + process + "_noise";
    }

    std::string filename_base(std::string process, int NPV, int EJW, int EJO,
                              int PP, int TBS, int seed_pT_cut, int seed_noise) {
        std::string base = "/u/at/chstan/nfs/summer_2014/ForConrad/files/";
        std::stringstream ss;
        ss.str(std::string());
        ss << base
           << util::dataset_name(process, NPV, EJW, EJO, PP, TBS, seed_pT_cut, seed_noise)
           << "/"
           << util::time_generated(process, NPV, EJW, EJO, PP, TBS, seed_pT_cut, seed_noise)
           << "/10s_" << NPV << "mu_0lw_0rec_"
           << seed_pT_cut << "cut_"
           << EJW << "ejw_"
           << EJO << "off_"
           << PP << "PP_"
           << TBS << "TBS_"
           << seed_noise << "SN";
        return ss.str();
    }

    std::string sanitize(std::string s) {
        std::string t = s;
        boost::replace_all(t, "/", "_div_");
        return t;
    }

    std::string filename(std::string process, int NPV, int EJW, int EJO,
                         int PP, int TBS, int seed_pT_cut, int seed_noise) {
        return util::filename_base(process, NPV, EJW,
                                   EJO, PP, TBS, seed_pT_cut, seed_noise) + ".root";
    }

    std::string filename_w(std::string process, int NPV, int EJW, int EJO,
                           int PP, int TBS, int seed_pT_cut, int seed_noise) {
        return util::filename_base(process, NPV, EJW,
                                   EJO, PP, TBS, seed_pT_cut, seed_noise) + "_w.root";
    }

    std::string eventname(std::string process, int NPV, int EJW, int EJO,
                          int PP, int TBS, int seed_pT_cut, int seed_noise) {
        std::stringstream ss;
        ss.str(std::string());
        ss << process << "_" << NPV << "_" << EJW << "_" << EJO << "_"
           << PP << "_" << TBS << "_" << seed_pT_cut << "_" << seed_noise;
        return ss.str();
    }

    ParameterSet::ParameterSet() {};

    std::string ParameterSet::EventName() {
        return util::eventname(_process, _NPV, _EJW, _EJO, _PP, _TBS, _seed_pT_cut, _seed_noise);
    }

    std::string latexed_pT_range(float low, float high) {
        unsigned int i_low = static_cast<unsigned int>(low);
        unsigned int i_high = static_cast<unsigned int>(high);
        std::stringstream ss;
        ss.str(std::string());
        ss << i_low << " GeV < "
           << "p_{T}^{Leading anti-k_{t} R=1.0} < "
           << i_high << " GeV";
        return ss.str();
    }

    std::string ParameterSet::Print() {
        std::stringstream ss;
        ss.str(std::string());
        ss << _process << "_"
           << _NPV << "mu_"
           << _EJW << "ejw_"
           << _EJO << "off_"
           << _PP << "PP_"
           << _TBS << "TBS_"
           << _seed_pT_cut << "cut_"
           << _seed_noise << "noise";
        return ss.str();
    }
    std::string ParameterSet::PrintExceptDistinguished(std::string distinguished) {
        std::stringstream ss;
        ss.str(std::string());
        if (distinguished != "process") {
            ss << _process;
            if (distinguished != "NPV") {
                ss << "_"; // print the _
            }
        }
        if (distinguished != "NPV") {
            ss << _NPV << "mu"; // note no _
        }
        if (distinguished != "EJW") {
            ss << "_" << _EJW << "ejw";
        }
        if (distinguished != "EJO") {
            ss << "_" << _EJO << "off";
        }
        if (distinguished != "PP") {
            ss << "_" << _PP << "PP";
        }
        if (distinguished != "TBS") {
            ss << "_" << _TBS << "TBS";
        }
        if (distinguished != "cut") {
            ss << "_" << _seed_pT_cut << "cut";
        }
        if (distinguished != "noise") {
            ss << "_" << _seed_noise << "noise";
        }

        return ss.str();
    }

    std::string fancy_distinguished_param_name(util::ParameterSet ps, std::string distinguished_var_type) {
        std::stringstream ss;
        ss.str(std::string());
        if (distinguished_var_type == "process")
            return util::fancy_var_name(ps._process);
        if (distinguished_var_type == "NPV") {
            ss << "#mu = " << ps._NPV;
            return ss.str();
        }
        if (distinguished_var_type == "EJW") {
            ss << "Event Jet Strength  = " << ps._EJW;
            return ss.str();
        }
        if (distinguished_var_type == "EJO") {
            ss << "#rho Offset = " << ps._EJO;
            return ss.str();
        }
        if (distinguished_var_type == "PP") {
            return ps._PP == 0 ? "No Post Processing" : "Final Merger";
        }
        if (distinguished_var_type == "TBS") {
            return ps._TBS == 0 ? "Tower Subtraction Off" : "Tower Subtraction On";
        }
        if (distinguished_var_type == "cut") {
            ss << "Seed "<< util::fancy_var_name("pt")
               << " cut = "
               << ps._seed_pT_cut << " GeV";
            return ss.str();
        }
        if (distinguished_var_type == "noise") {
            ss << "Seed Noise = " << format_float(0.01 * ps._seed_noise, 3);
            return ss.str();
        }
        assert(0);
        return "";
    }
}
