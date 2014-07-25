#include <iostream>
#include <vector>
#include <string>

#include <TColor.h>

#include "AnalyzeFuzzyTools.h"
#include "AtlasStyle.h"

#include "boost/program_options.hpp"

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::map;

namespace po = boost::program_options;

int main(int argc, char *argv[]) {
    cout << "Called as: ";
    for (int i = 0; i < argc; i++) {
        cout << argv[i] << " ";
    }
    cout << endl;

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "Produces help message");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help") > 0) {
        cout << desc << "\n";
        return 1;
    }

    SetAtlasStyle();

    string zero_loc = "/u/at/chstan/nfs/summer_2014/ForConrad/results/with_pu/wno_s10/0.root";
    string ten_loc = "/u/at/chstan/nfs/summer_2014/ForConrad/results/with_pu/wno_s10/10.root";
    string twenty_loc = "/u/at/chstan/nfs/summer_2014/ForConrad/results/with_pu/wno_s10/20.root";
    string thirty_loc = "/u/at/chstan/nfs/summer_2014/ForConrad/results/with_pu/wno_s10/30.root";

    string outdir = "/u/at/chstan/nfs/summer_2014/ForConrad/results/plots/";

    CanvasHelper cdec_test("Mass [GeV]", "", "Test Title", outdir, 800, 800);

    HistHelper hh_one(zero_loc, "CA_m", "CA Mass", 510, 0, 400, 50,
                      StyleTypes::NONE, kBlue);
    HistHelper hh_two(zero_loc, "antikt_m", "Anti-kt Mass", 510, 0, 400, 50,
                      StyleTypes::DASHED, kBlack);

    vector<HistHelper> histdecs_test;
    histdecs_test.push_back(hh_one);
    histdecs_test.push_back(hh_two);

    prettyHist<float>(histdecs_test, cdec_test);
}
