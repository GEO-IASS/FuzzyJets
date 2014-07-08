#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iterator>
#include <algorithm>


#include "TString.h"
#include "TSystem.h"
#include "TError.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TDatabasePDG.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"

#include "Pythia8/Pythia.h"

#include "FuzzyTools.h"
#include "FuzzyAnalysis.h"

#include "boost/program_options.hpp"

using std::cout;
using std::endl;
using std::string;
using std::map;
using namespace std;
namespace po = boost::program_options;


int main(int argc, char* argv[]){
    // argument parsing  ------------------------
    cout << "Called as: ";
    for(int ii=0; ii<argc; ++ii){
        cout << argv[ii] << " ";
    }
    cout << endl;

    // agruments
    int nEvents = 0;
    int fDebug  = 0;
    string outName = "FuzzyJets.root";
    string pythiaConfigName = "configs/default.pythia";
    string directory = "results/tmp/";

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("NEvents", po::value<int>(&nEvents)->default_value(1000) ,    "Number of Events ")
        ("Debug",   po::value<int>(&fDebug) ->default_value(0) ,     "Debug flag")
        ("OutFile", po::value<string>(&outName)->default_value("test.root"), "output file name")
        ("PythiaConfig", po::value<string>(&pythiaConfigName)->default_value("configs/default.pythia"), "Pythia configuration file location")
        ("Directory", po::value<string>(&directory)->default_value("results/tmp/"), "Directory in which to place results (.root files etc.)");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")>0){
        cout << desc << "\n";
        return 1;
    }

    // Configure and initialize pythia
    Pythia8::Pythia* pythia8 = new Pythia8::Pythia();

    vector<string> pythiaCommands;

    ifstream pythiaConfigFile(pythiaConfigName.c_str());
    string line;
    while (getline(pythiaConfigFile, line, '\n')) {
        pythiaCommands.push_back(line);
    }

    string command;
    for(unsigned int i = 0; i < pythiaCommands.size(); i++) {
        command = pythiaCommands[i];
        if (command.size() == 0) continue;
        if (command[0] == '#') {
            // its a comment
            cout << command << endl;
            continue;
        }
        pythia8->readString(command);
    }

    if (fDebug == 0) {
        pythia8->readString("Next:numberShowEvent = 0");
    }

    //this has to be the last line!
    pythia8->init(2212 /* p */, 2212 /* p */, 14000. /* GeV */);

    // FuzzyAnalysis
    FuzzyAnalysis * analysis = new FuzzyAnalysis();
    analysis->SetOutName(outName);
    analysis->SetPrefix(directory);
    analysis->Begin();
    analysis->Debug(fDebug);

    // Event loop
    cout << "running on " << nEvents << " events " << endl;
    for (Int_t iev = 0; iev < nEvents; iev++) {
        if (iev%100==0) cout << iev << " " << nEvents << endl;
        analysis->AnalyzeEvent(iev, pythia8);
    }

    analysis->End();

    // that was it
    delete pythia8;
    delete analysis;

    return 0;
}
