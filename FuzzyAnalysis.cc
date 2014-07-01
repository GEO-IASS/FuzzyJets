#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <set>


#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TMath.h"
#include "TVector3.h"
#include "TMatrix.h"

#include "FuzzyAnalysis.h"
#include "FuzzyTools.h"

#include "myFastJetBase.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/Selector.hh"

#include "Pythia8/Pythia.h"

using namespace std;

// Constructor
FuzzyAnalysis::FuzzyAnalysis(){
    fDebug = false;
    if(fDebug) cout << "FuzzyAnalysis::FuzzyAnalysis Start " << endl;
    ftest = 0;
    fOutName = "test.root";
    tool = new FuzzyTools();

    tool->SetClusteringMode(FuzzyTools::RECOMBINATION);
    tool->SetKernelType(FuzzyTools::GAUSSIAN);
    // jet def
    m_jet_def                = new fastjet::JetDefinition(fastjet::antikt_algorithm, 0.4);
    m_jet_def_largeR_antikt  = new fastjet::JetDefinition(fastjet::antikt_algorithm, 1.0);
    m_jet_def_largeR_ca      = new fastjet::JetDefinition(fastjet::cambridge_algorithm, 1.0);

    if(fDebug) cout << "FuzzyAnalysis::FuzzyAnalysis End " << endl;
}

// Destructor
FuzzyAnalysis::~FuzzyAnalysis(){
    delete tool;
    delete m_jet_def;
    delete m_jet_def_largeR_antikt;
    delete m_jet_def_largeR_ca;
}

// Begin method
void FuzzyAnalysis::Begin(){
    // Declare TTree
    tF = new TFile(fOutName.c_str(), "RECREATE");
    tT = new TTree("EventTree", "Event Tree for Fuzzy");

    DeclareBranches();
    ResetBranches();


    return;
}

// End
void FuzzyAnalysis::End(){

    tT->Write();
    tF->Close();
    return;
}

// Analyze
void FuzzyAnalysis::AnalyzeEvent(int ievt, Pythia8::Pythia* pythia8){
    if(fDebug) cout << "FuzzyAnalysis::AnalyzeEvent Begin " << endl;

    // generate a new event
    if (!pythia8->next()) return;
    if(fDebug) cout << "FuzzyAnalysis::AnalyzeEvent Event Number " << ievt << endl;

    // reset branches
    ResetBranches();

    // new event-----------------------
    fTEventNumber = ievt;
    std::vector <fastjet::PseudoJet>           particlesForJets;
    vector<fastjet::PseudoJet> tops;
    fastjet::PseudoJet dl;
    tops.push_back(dl);
    fastjet::PseudoJet dl2;
    tops.push_back(dl);


    // Particle loop -----------------------------------------------------------
    // The Pythia event listing contains a lot more than we want to process,
    // we prune out certain particles (muons / neutrinos) and only add final
    // state particles

    double px, py, pz, e;
    for (unsigned int ip=0; ip < (unsigned) pythia8->event.size(); ++ip){
        px = pythia8->event[ip].px();
        py = pythia8->event[ip].py();
        pz = pythia8->event[ip].pz();
        e  = pythia8->event[ip].e();
        fastjet::PseudoJet p(px, py, pz, e);
        p.set_user_info(new MyUserInfo(pythia8->event[ip].id(),ip,pythia8->event[ip].charge()));

        // In reality we should be more careful about finding tops,
        // but this will do for now. In the future consider refactoring
        // and tracing to find a top quark with no daughters
        if (pythia8->event[ip].id()  ==6) tops[0]=p;
        if (pythia8->event[ip].id()  ==-6) tops[1]=p;

        // prune uninteresting particles
        if (!pythia8->event[ip].isFinal() )      continue; // only final state
        if (fabs(pythia8->event[ip].id())  ==12) continue; // prune nu-e
        if (fabs(pythia8->event[ip].id())  ==13) continue; // ...   mu
        if (fabs(pythia8->event[ip].id())  ==14) continue; // ...   nu-mu
        if (fabs(pythia8->event[ip].id())  ==16) continue; // ...   nu-tau
        if (pythia8->event[ip].pT()       < 0.5) continue; // ...   low pT

        particlesForJets.push_back(p);

    } // end particle loop -----------------------------------------------

    // large-R jets: C/A --------------------
    fastjet::ClusterSequence csLargeR_ca(particlesForJets, *m_jet_def_largeR_ca);
    vector<fastjet::PseudoJet> myJetsLargeR_ca = fastjet::sorted_by_pt(csLargeR_ca.inclusive_jets(1.0)); //was 10
    for (unsigned int ij = 0; ij < myJetsLargeR_ca.size(); ij++){
        if (ij > 0) break;
        fastjet::PseudoJet tj = myJetsLargeR_ca[ij];
        //std::cout << "CA " << tj.pt() << " " << tj.m() << std::endl;
        fTCA_m = tj.m();
        fTCA_pt = tj.pt();
    }

    // Fuzzy Jets: mGMM --------------------
    vector<fastjet::PseudoJet> parts = particlesForJets;
    tool->SetSeeds(parts);

    vector<vector<double> > Weights;
    vector<TMatrix> mGMMjetsparams;
    vector<fastjet::PseudoJet> mGMMjets = tool->ClusterFuzzyGaussian(parts, &Weights, &mGMMjetsparams);

    int learnedClusterCount = Weights[0].size();
    if (learnedClusterCount == 0) {
        cout << "At event " << ievt << " returned no clusters. Continuing." << endl;
        return;
    }

    int leadmGMM=-1;
    double maxpt = -1;
    for (unsigned int i= 0; i<mGMMjets.size(); i++){
        double holdpt = tool->MLpT(parts,Weights,i,Weights[0].size(),0);
        if (holdpt>maxpt){
            maxpt = holdpt;
            leadmGMM = i;
        }
    }


    std::cout << mGMMjets.size() << " " << myJetsLargeR_ca.size() << " " << Weights.size() << std::endl;
    tool->EventDisplay(parts,myJetsLargeR_ca,tops,mGMMjets,Weights,leadmGMM, mGMMjetsparams,TString::Format("%i",ievt));
    tool->Qjetmass(parts, Weights, leadmGMM, TString::Format("%i",ievt));

    fTmGMM_m = tool->MLpT(parts,Weights,leadmGMM,Weights[0].size(),1);
    fTmGMM_pt = tool->MLpT(parts,Weights,leadmGMM,Weights[0].size(),0);
    //std::cout << "mGMM " << tool->MLpT(particlesForJets,Weights,leadmGMM,newk,0) << " " << tool->MLpT(particlesForJets,Weights,leadmGMM,newk,1) << std::endl;
    int mytop=0;
    if (tops[1].pt()> tops[0].pt()){
        mytop=1;
    }
    fTtoppt = tops[0].pt();
    fTdeltatop = tops[0].delta_R(mGMMjets[leadmGMM]);
    if (tops[1].delta_R(mGMMjets[leadmGMM]) < fTdeltatop){
        fTdeltatop = tops[1].delta_R(mGMMjets[leadmGMM]);
        fTtoppt = tops[1].pt();
    }

    tT->Fill();

    if(fDebug) cout << "FuzzyAnalysis::AnalyzeEvent End " << endl;
    return;
}



// declate branches
void FuzzyAnalysis::DeclareBranches(){

    // Event Properties
    tT->Branch("EventNumber",               &fTEventNumber,            "EventNumber/I");
    tT->Branch("CA_m,",                  &fTCA_m,                 "CA_m/F");
    tT->Branch("CA_pt,",                  &fTCA_pt,                 "CA_pt/F");
    tT->Branch("mGMM_m,",                  &fTmGMM_m,                 "mGMM_m/F");
    tT->Branch("mGMM_pt,",                  &fTmGMM_pt,                 "mGMM_pt/F");
    tT->Branch("deltatop", &fTdeltatop, "deltatop/F");
    tT->Branch("toppt", &fTtoppt, "toppt/F");

    tT->GetListOfBranches()->ls();

    return;
}


// resets vars
void FuzzyAnalysis::ResetBranches(){
    // reset branches
    fTEventNumber                 = -999;
    fTCA_m = 0;
    fTCA_pt = 0;
    fTmGMM_m = 0;
    fTmGMM_pt = 0;
    fTdeltatop = 0;
    fTtoppt = 0.;

}
