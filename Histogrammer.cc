#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <tuple>
#include <utility>
#include <unordered_map>

#include <TColor.h>

#include "AnalyzeFuzzyTools.h"
#include "AtlasStyle.h"

#include "boost/program_options.hpp"

namespace po = boost::program_options;

typedef std::tuple<int, int, int> map_key_t;

struct key_hash : public std::unary_function<map_key_t, std::size_t> {
    std::size_t operator()(const map_key_t& key) const {
        return std::get<0>(key) ^ std::get<1>(key) ^ std::get<2>(key);
    }
};

struct key_equal : public std::binary_function<map_key_t, map_key_t, bool> {
    bool operator()(const map_key_t& key1, const map_key_t& key2) const {
        return (std::get<0>(key1) == std::get<0>(key2) &&
                std::get<1>(key1) == std::get<1>(key2) &&
                std::get<2>(key1) == std::get<2>(key2));
    }
};

typedef std::unordered_map<map_key_t, std::string, key_hash, key_equal> file_map_t;

TString buildFileLocation(int size, int learn, int n_pileup_vertices, std::string prefix) {
  return TString::Format("%s%ds_%dmu_%d.root", prefix.c_str(), size, n_pileup_vertices, learn);
}

file_map_t
constructFileMap(std::vector<int> const& sizes, std::vector<int> const& learns,
                 std::vector<int> const& pileup_vertices, std::string const& prefix) {
  file_map_t m;
  for (unsigned int size_iter = 0; size_iter < sizes.size(); size_iter++) {
    for (unsigned int learn_iter = 0; learn_iter < learns.size(); learn_iter++) {
      for (unsigned int pileup_iter = 0; pileup_iter < pileup_vertices.size(); pileup_iter++) {
        int current_size = sizes[size_iter];
        int current_learn = learns[learn_iter];
        int current_n_pileup = pileup_vertices[pileup_iter];
        TString file_location = buildFileLocation(current_size, current_learn, 
                                                  current_n_pileup, prefix);
        map_key_t new_key = std::make_tuple(current_size, current_learn, current_n_pileup);

        m[new_key] = file_location;                                           
      }
    }
  }
  return m;
}

int main(int argc, char *argv[]) {
    std::cout << "Called as: ";
    for (int i = 0; i < argc; i++) {
        std::cout << argv[i] << " ";
    }
    std::cout << std::endl;

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "Produces help message");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help") > 0) {
        std::cout << desc << std::endl;
        return 1;
    }

    SetAtlasStyle();

    static const int sizes_arr[] = {7, 8, 9, 10};
    static const int NPVs_arr[] = {0, 10, 20, 30};
    static const int learns_arr[] = {0, 1};
    static const std::string algs_arr[] = {"mGMM", "mGMMs", "mUMM", "mTGMM", "mTGMMs"};

    static const Int_t colors_arr[] = {kMagenta + 3, kViolet + 9, kTeal - 5, kGray + 3};
    static const Style_t styles_arr[] = {kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown};

    std::vector<int> sizes(sizes_arr, sizes_arr+sizeof(sizes_arr) / sizeof(sizes_arr[0]));
    std::vector<int> NPVs(NPVs_arr, NPVs_arr+sizeof(NPVs_arr) / sizeof(NPVs_arr[0]));
    std::vector<int> learns(learns_arr, learns_arr+sizeof(learns_arr) / sizeof(learns_arr[0]));
    std::vector<std::string> algs(algs_arr, algs_arr+sizeof(algs_arr) / sizeof(algs_arr[0]));

    std::string file_prefix = "/u/at/chstan/nfs/summer_2014/ForConrad/files/2014_08_04_19h02m14s/";

    file_map_t file_m = constructFileMap(sizes, learns, NPVs, file_prefix);

    for (auto iter = file_m.begin(); iter != file_m.end(); iter++) {
      std::cout << ":" << iter->second << std::endl;
    }

    std::string out_dir = "/u/at/chstan/nfs/summer_2014/ForConrad/results/plots/";
    
    CanvasHelper c_dec_test("Mass [GeV]", "", "Test Title", out_dir, 800, 800);
    

    // KEYS ARE BUILT BY (SIZE, LEARN, PILEUP)
    HistHelper hh_one(file_m[std::make_tuple(10, 0, 0)], "CA_m", "CA Mass", 510, 0, 400, 50,
                      StyleTypes::NONE, kBlue, kFullCircle, "");
    HistHelper hh_two(file_m[std::make_tuple(10, 0, 0)], "antikt_m", "Anti-kt Mass", 510, 0, 400, 50,
                      StyleTypes::DASHED, kBlack, kFullCircle, "");
    
    std::vector<HistHelper> v_hist_decs_test;
    v_hist_decs_test.push_back(hh_one);
    v_hist_decs_test.push_back(hh_two);
    
    prettyHist<float>(v_hist_decs_test, c_dec_test);


    // DO PILEUP COMPARISONS: MASS RESOLUTION
    std::vector<HistHelper> v_hist_decs;

    std::stringstream ss;

    int n_bins = 20;
    
    // god I'd love a for each loop guys please update the system gcc!
    for (unsigned int alg_iter = 0; alg_iter < algs.size(); alg_iter++) {
        for (unsigned int size_iter = 0; size_iter < sizes.size(); size_iter++) {
            for(unsigned int learn_iter = 0; learn_iter < learns.size(); learn_iter++) {
                v_hist_decs.clear();
                int size = sizes[size_iter];
                int learn = learns[learn_iter];
                std::string alg = algs[alg_iter];
                ss.str(std::string());
                ss << "Mass resolution with pileup " << alg << " jets: sz " << size << " lw " << learn;
                std::string canvas_title = ss.str();
                CanvasHelper c_dec("Mass [GeV]", "", canvas_title, out_dir, 800, 800);
                c_dec.diff_scale = false;
                for(unsigned int npv_iter = 0; npv_iter < NPVs.size(); npv_iter++) {
                    int npv = NPVs[npv_iter];
                    ss.str(std::string());
                    ss << alg << " Mass - PuV " << npv;
                    std::string title = ss.str();
                    ss.str(std::string());
                    ss << alg << "_m";
                    std::string branch = ss.str();
                    HistHelper hist_helper_temp(file_m[std::make_tuple(size, learn, npv)],
                                                branch, title, 510, 0, 400, n_bins, 
                                                StyleTypes::NONE, colors_arr[npv_iter],
                                                styles_arr[npv_iter], "p");
                    v_hist_decs.push_back(hist_helper_temp);
                }
                HistHelper hist_helper_temp(file_m[std::make_tuple(size, learn, 0)],
                                            "antikt_m", "Anti-kt Mass - PuV 0", 
                                            510, 0, 400, n_bins, StyleTypes::STRIPED, kBlue,
                                            kFullCircle, "");
                v_hist_decs.push_back(hist_helper_temp);
                prettyHist<float>(v_hist_decs, c_dec);
            }
        }
    }

    // PT RESOLUTION
    for (unsigned int alg_iter = 0; alg_iter < algs.size(); alg_iter++) {
        for (unsigned int size_iter = 0; size_iter < sizes.size(); size_iter++) {
            for(unsigned int learn_iter = 0; learn_iter < learns.size(); learn_iter++) {
                v_hist_decs.clear();
                int size = sizes[size_iter];
                int learn = learns[learn_iter];
                std::string alg = algs[alg_iter];
                ss.str(std::string());
                ss << "pT resolution with pileup " << alg << " jets: sz " << size << " lw " << learn;
                std::string canvas_title = ss.str();
                CanvasHelper c_dec("pT [GeV]", "", canvas_title, out_dir, 800, 800);
                c_dec.diff_scale = false;
                for(unsigned int npv_iter = 0; npv_iter < NPVs.size(); npv_iter++) {
                    int npv = NPVs[npv_iter];
                    ss.str(std::string());
                    ss << alg << " pT - PuV " << npv;
                    std::string title = ss.str();
                    ss.str(std::string());
                    ss << alg << "_pt";
                    std::string branch = ss.str();
                    HistHelper hist_helper_temp(file_m[std::make_tuple(size, learn, npv)],
                                                branch, title, 510, 0, 800, n_bins, 
                                                StyleTypes::NONE, colors_arr[npv_iter],
                                                styles_arr[npv_iter], "p");
                    v_hist_decs.push_back(hist_helper_temp);
                }
                HistHelper hist_helper_temp(file_m[std::make_tuple(size, learn, 0)],
                                            "antikt_pt", "Anti-kt pT - PuV 0", 
                                            510, 0, 800, n_bins, StyleTypes::STRIPED, kBlue,
                                            kFullCircle, "");
                v_hist_decs.push_back(hist_helper_temp);
                prettyHist<float>(v_hist_decs, c_dec);
            }
        }
    }

    // HISTOGRAMS FOR SIZE COMPARISON
    for (unsigned int alg_iter = 0; alg_iter < algs.size(); alg_iter++) {
        for (unsigned int npv_iter = 0; npv_iter < NPVs.size(); npv_iter++) {
            for(unsigned int learn_iter = 0; learn_iter < learns.size(); learn_iter++) {
                v_hist_decs.clear();
                int learn = learns[learn_iter];
                int npv = NPVs[npv_iter];
                std::string alg = algs[alg_iter];
                ss.str(std::string());
                ss << "Mass resolution by size " << alg << " lw " << learn << " npu " << npv;
                std::string canvas_title = ss.str();
                CanvasHelper c_dec("Mass [GeV]", "", canvas_title, out_dir, 800, 800);
                c_dec.diff_scale = false;
                for(unsigned int size_iter = 0; size_iter < sizes.size(); size_iter++) {
                    int size = sizes[size_iter];
                    ss.str(std::string());
                    ss << alg << " Mass - Size " << size;
                    std::string title = ss.str();
                    ss.str(std::string());
                    ss << alg << "_m";
                    std::string branch = ss.str();
                    HistHelper hist_helper_temp(file_m[std::make_tuple(size, learn, npv)],
                                                branch, title, 510, 0, 400, n_bins, 
                                                StyleTypes::NONE, colors_arr[size_iter],
                                                styles_arr[size_iter], "p");
                    v_hist_decs.push_back(hist_helper_temp);
                }
                HistHelper hist_helper_temp(file_m[std::make_tuple(10, learn, npv)],
                                            "antikt_m_trimmed_three", "Anti-kt trimmed Mass", 
                                            510, 0, 400, n_bins, StyleTypes::STRIPED, kBlue,
                                            kFullCircle, "");
                v_hist_decs.push_back(hist_helper_temp);
                prettyHist<float>(v_hist_decs, c_dec);
            }
        }
    }

    // HISTOGRAMS FOR WEIGHT LEARNING COMPARISON
    for (unsigned int alg_iter = 0; alg_iter < algs.size(); alg_iter++) {
        for (unsigned int npv_iter = 0; npv_iter < NPVs.size(); npv_iter++) {
            for(unsigned int size_iter = 0; size_iter < sizes.size(); size_iter++) {
                v_hist_decs.clear();
                int size = sizes[size_iter];
                int npv = NPVs[npv_iter];
                std::string alg = algs[alg_iter];
                ss.str(std::string());
                ss << "Mass resolution weight learning " << alg << " sz " << size << " npu " << npv;
                std::string canvas_title = ss.str();
                CanvasHelper c_dec("Mass [GeV]", "", canvas_title, out_dir, 800, 800);
                c_dec.diff_scale = false;
                for(unsigned int learn_iter = 0; learn_iter < learns.size(); learn_iter++) {
                    int learn = learns[learn_iter];
                    ss.str(std::string());
                    ss << alg << " Mass - Learning " << learn;
                    std::string title = ss.str();
                    ss.str(std::string());
                    ss << alg << "_m";
                    std::string branch = ss.str();
                    HistHelper hist_helper_temp(file_m[std::make_tuple(size, learn, npv)],
                                                branch, title, 510, 0, 400, n_bins, 
                                                StyleTypes::NONE, colors_arr[learn_iter],
                                                styles_arr[learn_iter], "p");
                    v_hist_decs.push_back(hist_helper_temp);
                }
                HistHelper hist_helper_temp(file_m[std::make_tuple(size, 0, npv)],
                                            "antikt_m_trimmed_three", "Anti-kt trimmed Mass", 
                                            510, 0, 400, n_bins, StyleTypes::STRIPED, kBlue,
                                            kFullCircle, "");
                v_hist_decs.push_back(hist_helper_temp);
                prettyHist<float>(v_hist_decs, c_dec);
            }
        }
    }

    // =============================SOFT MASS=======================================
    
    for (unsigned int alg_iter = 0; alg_iter < algs.size(); alg_iter++) {
        for (unsigned int size_iter = 0; size_iter < sizes.size(); size_iter++) {
            for(unsigned int learn_iter = 0; learn_iter < learns.size(); learn_iter++) {
                v_hist_decs.clear();
                int size = sizes[size_iter];
                int learn = learns[learn_iter];
                std::string alg = algs[alg_iter];
                ss.str(std::string());
                ss << "Soft mass resolution with pileup " << alg << " jets: sz " << size << " lw " << learn;
                std::string canvas_title = ss.str();
                CanvasHelper c_dec("Mass [GeV]", "", canvas_title, out_dir, 800, 800);
                c_dec.diff_scale = false;
                for(unsigned int npv_iter = 0; npv_iter < NPVs.size(); npv_iter++) {
                    int npv = NPVs[npv_iter];
                    ss.str(std::string());
                    ss << alg << " Soft Mass - PuV " << npv;
                    std::string title = ss.str();
                    ss.str(std::string());
                    ss << alg << "_m_soft";
                    std::string branch = ss.str();
                    HistHelper hist_helper_temp(file_m[std::make_tuple(size, learn, npv)],
                                                branch, title, 510, 0, 400, n_bins, 
                                                StyleTypes::NONE, colors_arr[npv_iter],
                                                styles_arr[npv_iter], "p");
                    v_hist_decs.push_back(hist_helper_temp);
                }
                HistHelper hist_helper_temp(file_m[std::make_tuple(size, learn, 0)],
                                            "antikt_m", "Anti-kt Mass - PuV 0", 
                                            510, 0, 400, n_bins, StyleTypes::STRIPED, kBlue,
                                            kFullCircle, "");
                v_hist_decs.push_back(hist_helper_temp);
                prettyHist<float>(v_hist_decs, c_dec);
            }
        }
    }

    // HISTOGRAMS FOR SIZE COMPARISON
    for (unsigned int alg_iter = 0; alg_iter < algs.size(); alg_iter++) {
        for (unsigned int npv_iter = 0; npv_iter < NPVs.size(); npv_iter++) {
            for(unsigned int learn_iter = 0; learn_iter < learns.size(); learn_iter++) {
                v_hist_decs.clear();
                int learn = learns[learn_iter];
                int npv = NPVs[npv_iter];
                std::string alg = algs[alg_iter];
                ss.str(std::string());
                ss << "Soft mass resolution by size " << alg << " lw " << learn << " npu " << npv;
                std::string canvas_title = ss.str();
                CanvasHelper c_dec("Mass [GeV]", "", canvas_title, out_dir, 800, 800);
                c_dec.diff_scale = false;
                for(unsigned int size_iter = 0; size_iter < sizes.size(); size_iter++) {
                    int size = sizes[size_iter];
                    ss.str(std::string());
                    ss << alg << " Soft Mass - Size " << size;
                    std::string title = ss.str();
                    ss.str(std::string());
                    ss << alg << "_m_soft";
                    std::string branch = ss.str();
                    HistHelper hist_helper_temp(file_m[std::make_tuple(size, learn, npv)],
                                                branch, title, 510, 0, 400, n_bins, 
                                                StyleTypes::NONE, colors_arr[size_iter],
                                                styles_arr[size_iter], "p");
                    v_hist_decs.push_back(hist_helper_temp);
                }
                HistHelper hist_helper_temp(file_m[std::make_tuple(10, learn, npv)],
                                            "antikt_m_trimmed_three", "Anti-kt trimmed Mass", 
                                            510, 0, 400, n_bins, StyleTypes::STRIPED, kBlue,
                                            kFullCircle, "");
                v_hist_decs.push_back(hist_helper_temp);
                prettyHist<float>(v_hist_decs, c_dec);
            }
        }
    }

    // HISTOGRAMS FOR WEIGHT LEARNING COMPARISON
    for (unsigned int alg_iter = 0; alg_iter < algs.size(); alg_iter++) {
        for (unsigned int npv_iter = 0; npv_iter < NPVs.size(); npv_iter++) {
            for(unsigned int size_iter = 0; size_iter < sizes.size(); size_iter++) {
                v_hist_decs.clear();
                int size = sizes[size_iter];
                int npv = NPVs[npv_iter];
                std::string alg = algs[alg_iter];
                ss.str(std::string());
                ss << "Soft mass resolution weight learning " << alg << " sz " << size << " npu " << npv;
                std::string canvas_title = ss.str();
                CanvasHelper c_dec("Mass [GeV]", "", canvas_title, out_dir, 800, 800);
                c_dec.diff_scale = false;
                for(unsigned int learn_iter = 0; learn_iter < learns.size(); learn_iter++) {
                    int learn = learns[learn_iter];
                    ss.str(std::string());
                    ss << alg << " Soft Mass - Learning " << learn;
                    std::string title = ss.str();
                    ss.str(std::string());
                    ss << alg << "_m_soft";
                    std::string branch = ss.str();
                    HistHelper hist_helper_temp(file_m[std::make_tuple(size, learn, npv)],
                                                branch, title, 510, 0, 400, n_bins, 
                                                StyleTypes::NONE, colors_arr[learn_iter],
                                                styles_arr[learn_iter], "p");
                    v_hist_decs.push_back(hist_helper_temp);
                }
                HistHelper hist_helper_temp(file_m[std::make_tuple(size, 0, npv)],
                                            "antikt_m_trimmed_three", "Anti-kt trimmed Mass", 
                                            510, 0, 400, n_bins, StyleTypes::STRIPED, kBlue,
                                            kFullCircle, "");
                v_hist_decs.push_back(hist_helper_temp);
                prettyHist<float>(v_hist_decs, c_dec);
            }
        }
    }

    /// ==========================PILEUP STUDIES=============================
    for (unsigned int alg_iter = 0; alg_iter < algs.size(); alg_iter++) {
        for (unsigned int size_iter = 0; size_iter < sizes.size(); size_iter++) {
            for(unsigned int learn_iter = 0; learn_iter < learns.size(); learn_iter++) {
                v_hist_decs.clear();
                int size = sizes[size_iter];
                int learn = learns[learn_iter];
                std::string alg = algs[alg_iter];
                ss.str(std::string());
                ss << "Unclustered pileup " << alg << " jets: sz " << size << " lw " << learn;
                std::string canvas_title = ss.str();
                CanvasHelper c_dec("Fraction of Clustered Particles", "", canvas_title, out_dir, 800, 800);
                c_dec.diff_scale = false;
                for(unsigned int npv_iter = 1; npv_iter < NPVs.size(); npv_iter++) {

                    int npv = NPVs[npv_iter];
                    ss.str(std::string());
                    ss << alg << " - PuV " << npv;
                    std::string title = ss.str();
                    ss.str(std::string());
                    ss << alg << "_ucpu";
                    std::string branch = ss.str();
                    HistHelper hist_helper_temp(file_m[std::make_tuple(size, learn, npv)],
                                                branch, title, 510, 0, 1, n_bins, 
                                                StyleTypes::NONE, colors_arr[npv_iter],
                                                styles_arr[npv_iter], "p");
                    v_hist_decs.push_back(hist_helper_temp);
                }
                prettyHist<float>(v_hist_decs, c_dec);
            }
        }
    }
    
}
