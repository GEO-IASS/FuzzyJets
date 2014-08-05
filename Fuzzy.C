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
#include <assert.h>
#include <time.h>

#include "boost/program_options.hpp"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"

#include "Pythia8/Pythia.h"

#include "FuzzyTools.h"
#include "FuzzyAnalysis.h"
#include "ROOTConf.h"

#ifdef WITHROOT
#include "TString.h"
#include "TSystem.h"
#include "TError.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#endif

namespace po = boost::program_options;

void eventLoop (int n_events, Pythia8::Pythia *pythia8,
                Pythia8::Pythia *pythia_MB, int NPV,
                FuzzyAnalysis *analysis) {
    clock_t start = clock();
    clock_t now;
    double secs;
    for (int event_iter = 0; event_iter < n_events; event_iter++) {
        analysis->AnalyzeEvent(event_iter, pythia8, pythia_MB, NPV);
        now = clock();
        secs = ((double) (now - start) / (event_iter + 1)) / CLOCKS_PER_SEC;
        cout << "Seconds per event (sample: " << (event_iter + 1)
             << " events): " << secs << endl;
    }
}

int main(int argc, char* argv[]){
    // argument parsing  ------------------------
    cout << "Called as: ";
    for(int ii=0; ii<argc; ++ii){
        cout << argv[ii] << " ";
    }
    cout << endl;

    // agruments
    int n_events = 0;
    int f_debug  = 0;
    int NPV = -1;
    double size = -1;
    bool learn_weights = false;
    bool is_batch = false;
    string out_name = "FuzzyJets.root";
    string pythia_config_name = "configs/default.pythia";
    string directory = "results/tmp/";

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("NEvents", po::value<int>(&n_events)->default_value(1000) ,    "Number of Events ")
        ("Debug",   po::value<int>(&f_debug) ->default_value(0) ,     "Debug flag")
        ("OutFile", po::value<string>(&out_name)->default_value("test.root"), "Output file name")
        ("NPV",     po::value<int>(&NPV)->default_value(-1), "Number of primary vertices (pile-up)")
        ("Size",    po::value<double>(&size)->default_value(-1), "Internal size variable for Fuzzy Clustering")
        ("LearnWeights", po::value<bool>(&learn_weights)->default_value(false), "Whether to learn cluster weights")
        ("PythiaConfig", po::value<string>(&pythia_config_name)->default_value("configs/default.pythia"), "Pythia configuration file location")
        ("Batch",   po::value<bool>(&is_batch)->default_value(false), "Is this running on the batch?")
        ("Directory", po::value<string>(&directory)->default_value("results/tmp/"), "Directory in which to place results (.root files etc.)");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")>0){
        cout << desc << "\n";
        return 1;
    } else {
        cout << "LearnWeights: " << learn_weights << endl;
        cout << "Size: " << size << endl;
        cout << "Directory: " << directory << endl;
        cout << "NPV: " << NPV << endl;
    }
    if (is_batch) {
        assert(NPV != -1);
    }
    assert(size > 0);

    // Configure and initialize pythia
    Pythia8::Pythia* pythia8 = new Pythia8::Pythia();
    Pythia8::Pythia* pythia_MB = new Pythia8::Pythia();

    vector<string> pythia_commands;

    ifstream pythia_config_file(pythia_config_name.c_str());
    string line;
    while (getline(pythia_config_file, line, '\n')) {
        pythia_commands.push_back(line);
    }

    string command;
    for(unsigned int i = 0; i < pythia_commands.size(); i++) {
        command = pythia_commands[i];
        if (command.size() == 0) continue;
        if (command[0] == '#') {
            // its a comment
            cout << command << endl;
            continue;
        }
        pythia8->readString(command);
    }

    int seed = 0;
    stringstream ss;
    ss.clear();
    ss.str("");
    ss << "Random:seed = " << seed;
    //pythia_MB->readstring("Random:setSeed = on");
    //pythia_MB->readstring(ss.str());
    pythia_MB->readString("SoftQCD:nonDiffractive = on");
    pythia_MB->readString("HardQCD:all = off");
    pythia_MB->readString("PhaseSpace:pTHatMin = .1");
    pythia_MB->readString("PhaseSpace:pTHatMax = 20000");


    if (f_debug == 0) {
        pythia8->readString("Next:numberShowEvent = 0");
        pythia_MB->readString("Next:numberShowEvent = 0");
    }



    //this has to be the last line!
    pythia8->init(2212 /* p */, 2212 /* p */, 14000. /* GeV */);
    pythia_MB->init(2212, 2212, 14000.);

    // FuzzyAnalysis
    FuzzyAnalysis *analysis = new FuzzyAnalysis();
    if(NPV == -1) {
        // compute for NPV = 0, 10, 20, 30
        int NPVs [4] = {0, 10, 20, 30};
        for (unsigned int npv_iter = 0; npv_iter < 4; npv_iter++) {
            int current_npv = NPVs[npv_iter];
            ss.clear();
            ss.str("");
            ss << current_npv << ".root";

            analysis->SetOutName(ss.str());
            analysis->SetPrefix(directory);
            analysis->SetSize(size);
            analysis->SetBatched(false);
            analysis->SetLearnWeights(learn_weights);
            analysis->Begin();
            analysis->Debug(f_debug);

            eventLoop(n_events, pythia8, pythia_MB, current_npv, analysis);

            analysis->End();
        }
    } else {
        ss.clear();
        ss.str("");
        ss << NPV << ".root";
        if(is_batch) {
            analysis->SetOutName(out_name);
            analysis->SetPrefix("");
            analysis->SetBatched(true);
        } else {
            analysis->SetOutName(ss.str());
            analysis->SetPrefix(directory);
            analysis->SetBatched(false);
        }
        analysis->SetShouldPrint(false);
        analysis->SetSize(size);
        analysis->SetLearnWeights(learn_weights);
        analysis->Begin();
        analysis->Debug(f_debug);

        cout << "running on " << n_events << " events " << endl;
        eventLoop(n_events, pythia8, pythia_MB, NPV, analysis);

        analysis->End();
    }

    // that was it
    delete pythia8;
    delete analysis;

    return 0;
}
