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
#include "TROOT.h"
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
    float rho = 0;
    float rho_sum = 0;
    analysis->SetNEvents(n_events);
    for (int event_iter = 0; event_iter < n_events; event_iter++) {
        analysis->AnalyzeEvent(event_iter, pythia8, pythia_MB, NPV);
        if (analysis->QueryFloatValue("rho", &rho)) {
            rho_sum += rho;
        } else {
            assert(0 && "Could not find rho in analysis float map.");
        }
        now = clock();
        secs = ((double) (now - start) / (event_iter + 1)) / CLOCKS_PER_SEC;
        if (!analysis->IsBatched() || ((event_iter + 1) % 500) == 0) {
        cout << "Seconds per event (sample: " << (event_iter + 1)
             << " events): " << secs << endl;
        }
        //cout << "Average rho: " << rho_sum / (event_iter + 1) << endl;
    }
}

int main(int argc, char* argv[]){
    // Make sure we can write vectors
    #ifdef WITHROOT
    gROOT->ProcessLine("#include <vector>");
    #endif

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
    double ej_strength = -1;
    double rho_offset = 0;
    float pT_min = 5;
    float alpha = 1;
    float seed_location_noise;
    bool learn_weights = false;
    bool is_batch = false;
    bool do_recombination = true;
    bool do_mod_pi_test = false;
    bool do_jet_multiplicity_study = false;
    bool do_wandering_study = false;
    bool do_sigma_dependence_on_mu_study = false;
    bool do_sigma_dependence_on_mu_complete_study = false;
    bool do_sigma_mass_study = false;
    bool do_noise_study = false;
    bool do_final_recombination = false;
    bool do_tower_subtraction = false;
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
        ("SeedNoise",    po::value<float>(&seed_location_noise)->default_value(0), "Gaussian noise used to perturb seeds.")
        ("LearnWeights", po::value<bool>(&learn_weights)->default_value(false), "Whether to learn cluster weights")
        ("FinalRecombination", po::value<bool>(&do_final_recombination)->default_value(false), "Whether to perform one final recombination step.")
        ("TowerSubtraction", po::value<bool>(&do_tower_subtraction)->default_value(false), "Whether to perform tower based pT subtraction to suppress pileup.")
        ("PythiaConfig", po::value<string>(&pythia_config_name)->default_value("configs/default.pythia"), "Pythia configuration file location")
        ("pTMin", po::value<float>(&pT_min)->default_value(5), "Minimum pT for standard jets. (Including those used as seeds")
        ("Alpha", po::value<float>(&alpha)->default_value(1), "Fuzzy Clustering parameter controlling the pT power scaling exponent.")
        ("EventJetStrength", po::value<double>(&ej_strength)->default_value(-1), "Scalar used to determine how strong the likelihood for the event jet is.")
        ("RhoOffset", po::value<double>(&rho_offset)->default_value(0), "Scalar offset used to adjust event jet strength (together with EventJetStrength they give a pair of affine variables).")
        ("Batch",   po::value<bool>(&is_batch)->default_value(false), "Is this running on the batch?")
        ("Recombination", po::value<bool>(&do_recombination)->default_value(false), "Should we use fixed clusters or all particles with recombination?")
        ("ModPiTest", po::value<bool>(&do_mod_pi_test)->default_value(false), "Whether to run the mod pi testing suite.")
        ("JetMultiplicityStudy", po::value<bool>(&do_jet_multiplicity_study)->default_value(false), "Run the run jet multiplicity study.")
        ("NoiseStudy", po::value<bool>(&do_noise_study)->default_value(false), "Run the seed location noise study.")
        ("WanderingStudy", po::value<bool>(&do_wandering_study)->default_value(false), "Run jet wandering study.")
        ("SigmaDependenceOnMuStudy", po::value<bool>(&do_sigma_dependence_on_mu_study)->default_value(false), "Run sigma dependence study for first event.")
        ("SigmaDependenceOnMuCompleteStudy", po::value<bool>(&do_sigma_dependence_on_mu_complete_study)->default_value(false), "Run sigma dependence study for many events.")
        ("SigmaMassStudy", po::value<bool>(&do_sigma_mass_study)->default_value(false), "Mass inside 1 sigma study.")
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
    analysis->SetOutName(out_name);
    analysis->SetSeedLocationNoise(seed_location_noise);

    analysis->SetToolAlpha(alpha);

    analysis->SetDoTowerSubtraction(do_tower_subtraction);
    if (do_final_recombination) {
        analysis->SetPostProcessingMethod(FuzzyTools::ONE_DISTANCE_MERGER);
    } else {
        analysis->SetPostProcessingMethod(FuzzyTools::NO_POST_PROCESS);
    }
    if (ej_strength > 0) {
        analysis->SetEventJetType(FuzzyTools::FLAT);
        analysis->SetEventJetStrength(ej_strength);
        analysis->SetEventJetRhoOffset(rho_offset);
    } else {
        analysis->SetEventJetType(FuzzyTools::NONE);
        analysis->SetEventJetStrength(0);
        analysis->SetEventJetRhoOffset(0);
    }

    if (do_sigma_dependence_on_mu_complete_study) {
        analysis->IsNonStandardStudy();
        if (is_batch) {
            analysis->SetPrefix("");
        } else {
            analysis->SetPrefix(directory);
        }
        analysis->Begin();
        analysis->SigmaDependenceOnMuCompleteStudy(pythia8, pythia_MB);
        analysis->End();
        // shutdown
        delete analysis;
        delete pythia8;
        delete pythia_MB;
        return 0;
    }
    if (do_wandering_study) {
        analysis->SetRecombination(do_recombination);
        analysis->SetShouldPrint(false);
        analysis->SetSize(size);
        analysis->SetLearnWeights(learn_weights);

        analysis->Begin();

        analysis->IsNonStandardStudy();
        analysis->SetBatched(true);
        if (is_batch) {
            analysis->SetPrefix("");
        } else {
            analysis->SetPrefix(directory);
        }
        analysis->Begin();
        analysis->PtCutWanderingStudy(pythia8, pythia_MB, NPV, n_events);
        analysis->End();
        // shutdown
        delete analysis;
        delete pythia8;
        delete pythia_MB;
        return 0;
    }
    if (do_noise_study) {
        analysis->SetRecombination(do_recombination);
        analysis->SetShouldPrint(false);
        analysis->SetSize(size);
        analysis->SetLearnWeights(learn_weights);
        analysis->SetPtMin(pT_min);

        analysis->Begin();

        analysis->IsNonStandardStudy();
        analysis->SetBatched(true);
        if (is_batch) {
            analysis->SetPrefix("");
        } else {
            analysis->SetPrefix(directory);
        }
        analysis->Begin();
        analysis->NoiseStudy(pythia8, pythia_MB, NPV, n_events);
        analysis->End();
        // shutdown
        delete analysis;
        delete pythia8;
        delete pythia_MB;
        return 0;
    }
    if (do_sigma_dependence_on_mu_study) {
        analysis->IsNonStandardStudy();
        if (is_batch) {
            analysis->SetPrefix("");
        } else {
            analysis->SetPrefix(directory);
        }
        analysis->Begin();
        analysis->SigmaDependenceOnMuStudy(pythia8, pythia_MB);
        analysis->End();
        // shutdown
        delete analysis;
        delete pythia8;
        delete pythia_MB;
        return 0;
    }

    if (do_sigma_mass_study) {
        analysis->IsNonStandardStudy();
        if (is_batch) {
            analysis->SetPrefix("");
        } else {
            analysis->SetPrefix(directory);
        }
        analysis->Begin();
        analysis->SigmaMassStudy(pythia8, pythia_MB);
        analysis->End();
        // shutdown
        delete analysis;
        delete pythia8;
        delete pythia_MB;
        return 0;
    }

    analysis->SetToolAlpha(alpha);

    analysis->SetDoTowerSubtraction(do_tower_subtraction);
    if (do_final_recombination) {
        analysis->SetPostProcessingMethod(FuzzyTools::ONE_DISTANCE_MERGER);
    } else {
        analysis->SetPostProcessingMethod(FuzzyTools::NO_POST_PROCESS);
    }
    if (ej_strength > 0) {
        analysis->SetEventJetType(FuzzyTools::FLAT);
        analysis->SetEventJetStrength(ej_strength);
        analysis->SetEventJetRhoOffset(rho_offset);
    } else {
        analysis->SetEventJetType(FuzzyTools::NONE);
        analysis->SetEventJetStrength(0);
        analysis->SetEventJetRhoOffset(0);
    }


    if (do_jet_multiplicity_study) {
        analysis->SetPrefix(directory);
        analysis->SetSize(size);
        analysis->SetBatched(true);
        analysis->Begin();
        analysis->JetMultiplicityStudy(pythia8, pythia_MB);
        analysis->End();
        return 0;
    }
    if(NPV == -1) {
        // compute for NPV = 0, 10, 20, 30
        int NPVs [4] = {0, 10, 20, 30};
        for (unsigned int npv_iter = 0; npv_iter < 4; npv_iter++) {
            int current_npv = NPVs[npv_iter];
            ss.clear();
            ss.str("");
            ss << current_npv << ".root";
            analysis->SetPtMin(pT_min);
            analysis->SetRecombination(do_recombination);
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
        analysis->SetPtMin(pT_min);
        analysis->SetRecombination(do_recombination);
        analysis->SetShouldPrint(false);
        analysis->SetSize(size);
        analysis->SetLearnWeights(learn_weights);
        analysis->Begin();
        analysis->Debug(f_debug);

        cout << "running on " << n_events << " events " << endl;
        if (!do_mod_pi_test)
            eventLoop(n_events, pythia8, pythia_MB, NPV, analysis);

        if (do_mod_pi_test && !is_batch) {
            // don't do this on the batch, and rewrite necessary variables
            analysis->SetRecombination(false);
            analysis->SetPrefix(directory);
            analysis->SetSize(size);
            analysis->SetBatched(false);
            analysis->SetLearnWeights(true);
            analysis->PiFixStudy();
        }
        analysis->End();
    }

    // that was it
    delete pythia8;
    delete pythia_MB;
    delete analysis;

    return 0;
}
