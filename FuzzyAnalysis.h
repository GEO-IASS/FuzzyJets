#ifndef  FuzzyAnalysis_H
#define  FuzzyAnalysis_H

#include <vector>
#include <math.h>
#include <string>

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/Selector.hh"

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"

#include "FuzzyTools.h"
#include "myFastJetBase.h"
#include "Pythia8/Pythia.h"

using namespace std;

class FuzzyAnalysis{
 private:
    int  ftest;
    int  fDebug;
    string fOutName;

    TFile *tF;
    TTree *tT;
    FuzzyTools *tool;

    // Tree Vars ---------------------------------------
    int              fTEventNumber;
    float fTCA_m;
    float fTCA_pt;
    float fTmGMM_m;
    float fTmGMM_pt;
    float fTdeltatop;
    float fTtoppt;

    fastjet::JetDefinition     *m_jet_def;
    fastjet::JetDefinition     *m_jet_def_largeR_antikt;
    fastjet::JetDefinition     *m_jet_def_largeR_ca;

 public:
    FuzzyAnalysis ();
    ~FuzzyAnalysis ();

    void Begin();
    void AnalyzeEvent(int iEvt, Pythia8::Pythia *pythia8);
    void End();
    void DeclareBranches();
    void ResetBranches();
    void Debug(int debug){
        fDebug = debug;
    }
    void SetOutName(string outname){
        fOutName = outname;
    }

};

#endif
