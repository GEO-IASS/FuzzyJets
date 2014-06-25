#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>

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

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("NEvents", po::value<int>(&nEvents)->default_value(1000) ,    "Number of Events ")
        ("Debug",   po::value<int>(&fDebug) ->default_value(0) ,     "Debug flag")
        ("OutFile", po::value<string>(&outName)->default_value("test.root"), "output file name")
        ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")>0){
        cout << desc << "\n";
        return 1;
    }



    //------

    // Configure and initialize pythia
    Pythia8::Pythia* pythia8 = new Pythia8::Pythia();
    //pythia8->init(2212,2212,8000.);

    //Signal: Z' -> ttbar -> all hadronic


    pythia8->readString("NewGaugeBoson:ffbar2gmZZprime= on");
    pythia8->readString("32:m0=1500");
    pythia8->readString("Zprime:gmZmode=3");
    pythia8->readString("32:onMode = off");
    pythia8->readString("32:onIfAny = 6");
    pythia8->readString("24:onMode = off");
    pythia8->readString("24:onIfAny = 1 2 3 4");

    if (fDebug == 0) {

        pythia8->readString("Next:numberShowEvent = 0");

    }



    /*
      pythia8->readString("NewGaugeBoson:ffbar2Wprime = on");
      pythia8->readString("Wprime:coup2WZ=1");
      pythia8->readString("34:m0=600");
      pythia8->readString("34:onMode = off");
      pythia8->readString("34:onIfAny = 23 24");
      pythia8->readString("24:onMode = off");
      pythia8->readString("24:onIfAny = 1 2 3 4");
      pythia8->readString("23:onMode = off");
      pythia8->readString("23:onIfAny = 11");
    */

    //this has to be the last line!
    pythia8->init(2212 /* p */, 2212 /* p */, 14000. /* GeV */);


    //Background: QCD
    //pythia8->readString("HardQCD:all = on ");
    //pythia8->readString("PhaseSpace:pTHatMax = 700");
    //pythia8->readString("PhaseSpace:pTHatMin = 300");

    //pythia8->init(2212,2212,8000.);

    // FuzzyAnalysis
    FuzzyAnalysis * analysis = new FuzzyAnalysis();
    analysis->SetOutName(outName);
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
