#include <iostream>
#include <vector>
#include <string>

#include <TColor.h>

#include "AnalyzeFuzzyTools.h"
#include "AtlasStyle.h"

#include "boost/program_options.hpp"

namespace po = boost::program_options;

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

    std::string zero_loc = "/u/at/chstan/nfs/summer_2014/ForConrad/results/with_pu/wno_s10/0.root";
    std::string ten_loc = "/u/at/chstan/nfs/summer_2014/ForConrad/results/with_pu/wno_s10/10.root";
    std::string twenty_loc = "/u/at/chstan/nfs/summer_2014/ForConrad/results/with_pu/wno_s10/20.root";
    std::string thirty_loc = "/u/at/chstan/nfs/summer_2014/ForConrad/results/with_pu/wno_s10/30.root";

    std::string out_dir = "/u/at/chstan/nfs/summer_2014/ForConrad/results/plots/";

    CanvasHelper c_dec_test("Mass [GeV]", "", "Test Title", out_dir, 800, 800);

    HistHelper hh_one(zero_loc, "CA_m", "CA Mass", 510, 0, 400, 50,
                      StyleTypes::NONE, kBlue);
    HistHelper hh_two(zero_loc, "antikt_m", "Anti-kt Mass", 510, 0, 400, 50,
                      StyleTypes::DASHED, kBlack);

    std::vector<HistHelper> v_hist_decs_test;
    v_hist_decs_test.push_back(hh_one);
    v_hist_decs_test.push_back(hh_two);

    prettyHist<float>(v_hist_decs_test, c_dec_test);
}
