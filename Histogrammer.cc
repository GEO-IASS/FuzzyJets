#include <iostream>
#include <vector>
#include <string>

#include <TColor.h>

#include "AnalyzeFuzzyTools.h"

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

    CanvasHelper cdec_default("Mass [GeV]", "", "Test Title", 600, 600);
}
