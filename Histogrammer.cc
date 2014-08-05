#include <iostream>
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

    std::vector<int> sizes(sizes_arr, sizes_arr+sizeof(sizes_arr) / sizeof(sizes_arr[0]));
    std::vector<int> NPVs(NPVs_arr, NPVs_arr+sizeof(NPVs_arr) / sizeof(NPVs_arr[0]));
    std::vector<int> learns(learns_arr, learns_arr+sizeof(learns_arr) / sizeof(learns_arr[0]));
    std::string file_prefix = "/u/at/chstan/nfs/summer_2014/ForConrad/files/2014_08_04_19h02m14s/";

    file_map_t file_m = constructFileMap(sizes, learns, NPVs, file_prefix);

    for (auto iter = file_m.begin(); iter != file_m.end(); iter++) {
      std::cout << ":" << iter->second << std::endl;
    }

    std::string out_dir = "/u/at/chstan/nfs/summer_2014/ForConrad/results/plots/";
    
    //CanvasHelper c_dec_test("Mass [GeV]", "", "Test Title", out_dir, 800, 800);
    //
    //HistHelper hh_one(zero_loc, "CA_m", "CA Mass", 510, 0, 400, 50,
    //                  StyleTypes::NONE, kBlue);
    //HistHelper hh_two(zero_loc, "antikt_m", "Anti-kt Mass", 510, 0, 400, 50,
    //                  StyleTypes::DASHED, kBlack);
    //
    //std::vector<HistHelper> v_hist_decs_test;
    //v_hist_decs_test.push_back(hh_one);
    //v_hist_decs_test.push_back(hh_two);
    //
    //prettyHist<float>(v_hist_decs_test, c_dec_test);
}
