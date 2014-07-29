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
    int NPV = -1;
    double size = -1;
    bool learnWeights = false;
    bool isBatch = false;
    string outName = "FuzzyJets.root";
    string pythiaConfigName = "configs/default.pythia";
    string directory = "results/tmp/";

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("NEvents", po::value<int>(&nEvents)->default_value(1000) ,    "Number of Events ")
        ("Debug",   po::value<int>(&fDebug) ->default_value(0) ,     "Debug flag")
        ("OutFile", po::value<string>(&outName)->default_value("test.root"), "Output file name")
        ("NPV",     po::value<int>(&NPV)->default_value(-1), "Number of primary vertices (pile-up)")
        ("Size",    po::value<double>(&size)->default_value(-1), "Internal size variable for Fuzzy Clustering")
        ("LearnWeights", po::value<bool>(&learnWeights)->default_value(false), "Whether to learn cluster weights")
        ("PythiaConfig", po::value<string>(&pythiaConfigName)->default_value("configs/default.pythia"), "Pythia configuration file location")
        ("Batch",   po::value<bool>(&isBatch)->default_value(false), "Is this running on the batch?")
        ("Directory", po::value<string>(&directory)->default_value("results/tmp/"), "Directory in which to place results (.root files etc.)");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")>0){
        cout << desc << "\n";
        return 1;
    } else {
        cout << "LearnWeights: " << learnWeights << endl;
        cout << "Size: " << size << endl;
        cout << "Directory: " << directory << endl;
        cout << "NPV: " << NPV << endl;
    }
    if (isBatch) {
        assert(NPV != -1);
    }
    assert(size > 0);

    // Configure and initialize pythia
    Pythia8::Pythia* pythia8 = new Pythia8::Pythia();
    Pythia8::Pythia* pythia_MB = new Pythia8::Pythia();

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


    if (fDebug == 0) {
        pythia8->readString("Next:numberShowEvent = 0");
        pythia_MB->readString("Next:numberShowEvent = 0");
    }



    //this has to be the last line!
    pythia8->init(2212 /* p */, 2212 /* p */, 14000. /* GeV */);
    pythia_MB->init(2212, 2212, 14000.);

    // FuzzyAnalysis
    FuzzyAnalysis * analysis = new FuzzyAnalysis();
    if(NPV == -1) {
        // compute for NPV = 0, 10, 20, 30
        int NPVs [4] = {0, 10, 20, 30};
        for (unsigned int npviter = 0; npviter < 4; npviter++) {
            int currentNPV = NPVs[npviter];
            ss.clear();
            ss.str("");
            ss << currentNPV << ".root";
            analysis->SetOutName(ss.str());
            analysis->SetPrefix(directory);
            analysis->SetSize(size);
            analysis->SetBatched(false);
            analysis->SetLearnWeights(learnWeights);
            analysis->Begin();
            analysis->Debug(fDebug);
            for (int iev = 0; iev < nEvents; iev++) {
                analysis->AnalyzeEvent(iev, pythia8, pythia_MB, NPV);
            }
            analysis->End();
        }
    } else {
        ss.clear();
        ss.str("");
        ss << NPV << ".root";
        if(isBatch) {
            analysis->SetOutName(outName);
            analysis->SetPrefix("");
            analysis->SetBatched(true);
        } else {
            analysis->SetOutName(ss.str());
            analysis->SetPrefix(directory);
            analysis->SetBatched(false);
        }
        analysis->SetSize(size);
        analysis->SetLearnWeights(learnWeights);
        analysis->Begin();
        analysis->Debug(fDebug);

        // Event loop
        cout << "running on " << nEvents << " events " << endl;
        for (Int_t iev = 0; iev < nEvents; iev++) {
            if (iev%100==0) cout << iev << " " << nEvents << endl;
            analysis->AnalyzeEvent(iev, pythia8, pythia_MB, NPV);
        }

        analysis->End();
    }

    // that was it
    delete pythia8;
    delete analysis;

    return 0;
}
