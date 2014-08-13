#ifndef  FuzzyAnalysis_H
#define  FuzzyAnalysis_H

#include <vector>
#include <string>
#include <map>
#include <math.h>

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/Selector.hh"
#include "FuzzyTools.h"
#include "myFastJetBase.h"

#include "Pythia8/Pythia.h"

#include "ROOTConf.h"

#ifdef WITHROOT
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"
#endif

class FuzzyAnalysis{
 private:
    int  f_test;
    int  f_debug;
    double  f_size;
    bool f_learn_weights;
    bool batched;
    bool should_print;
    bool do_recombination;
    string f_out_name;
    string directory_prefix;

    std::map<string, float*> tree_vars;
    std::map<string, std::vector<float> > map_weight_vecs;

    #ifdef WITHROOT
    TFile *t_f;
    TTree *t_t;
    #endif
    
    FuzzyTools *tool;

    // Tree Vars ---------------------------------------
    int   fTEventNumber;
    int   fTNPV;
    float fTtoppt;

    float fTCA_pt;
    float fTCA_m;
    float fTCA_pufrac;
    float fTCA_m_pu;

    float fTantikt_m;
    float fTantikt_pt;
    float fTantikt_m_trimmed_two;
    float fTantikt_pt_trimmed_two;
    float fTantikt_m_trimmed_three;
    float fTantikt_pt_trimmed_three;
    float fTantikt_pufrac_trimmed_two;
    float fTantikt_pufrac_trimmed_three;
    float fTantikt_m_pu;
    float fTantikt_m_pu_trimmed_two;
    float fTantikt_m_pu_trimmed_three;

    float fTmUMM_m;
    float fTmUMM_pt;
    float fTdeltatop_mUMM;
    float fTmUMM_m_mean;
    float fTmUMM_m_var;
    float fTmUMM_m_skew;
    float fTmUMM_pt_mean;
    float fTmUMM_pt_var;
    float fTmUMM_pt_skew;
    float fTmUMM_ptl;
    float fTmUMM_ml;
    float fTmUMM_m_soft;
    float fTmUMM_pt_soft;
    float fTmUMM_ucpu;
    float fTmUMM_clpu;
    float fTmUMM_pufrac_soft;
    float fTmUMM_pufrac_hard;
    float fTmUMM_m_pu_soft;
    float fTmUMM_m_pu_hard;

    float fTmGMMs_m;
    float fTmGMMs_pt;
    float fTdeltatop_mGMMs;
    float fTmGMMs_m_mean;
    float fTmGMMs_m_var;
    float fTmGMMs_m_skew;
    float fTmGMMs_pt_mean;
    float fTmGMMs_pt_var;
    float fTmGMMs_pt_skew;
    float fTmGMMs_ptl;
    float fTmGMMs_ml;
    float fTmGMMs_m_soft;
    float fTmGMMs_pt_soft;
    float fTmGMMs_ucpu;
    float fTmGMMs_clpu;
    float fTmGMMs_pufrac_soft;
    float fTmGMMs_pufrac_hard;
    float fTmGMMs_m_pu_soft;
    float fTmGMMs_m_pu_hard;

    float fTmGMMc_m;
    float fTmGMMc_pt;
    float fTdeltatop_mGMMc;
    float fTmGMMc_m_mean;
    float fTmGMMc_m_var;
    float fTmGMMc_m_skew;
    float fTmGMMc_pt_mean;
    float fTmGMMc_pt_var;
    float fTmGMMc_pt_skew;
    float fTmGMMc_ptl;
    float fTmGMMc_ml;
    float fTmGMMc_m_soft;
    float fTmGMMc_pt_soft;
    float fTmGMMc_ucpu;
    float fTmGMMc_clpu;
    float fTmGMMc_pufrac_soft;
    float fTmGMMc_pufrac_hard;
    float fTmGMMc_m_pu_soft;
    float fTmGMMc_m_pu_hard;
    
    float fTmGMMc_r;
    float fTmGMMc_r_avg;
    float fTmGMMc_r_weighted_avg;

    float fTmTGMMs_m;
    float fTmTGMMs_pt;
    float fTdeltatop_mTGMMs;
    float fTmTGMMs_m_mean;
    float fTmTGMMs_m_var;
    float fTmTGMMs_m_skew;
    float fTmTGMMs_pt_mean;
    float fTmTGMMs_pt_var;
    float fTmTGMMs_pt_skew;
    float fTmTGMMs_ptl;
    float fTmTGMMs_ml;
    float fTmTGMMs_m_soft;
    float fTmTGMMs_pt_soft;
    float fTmTGMMs_ucpu;
    float fTmTGMMs_clpu;
    float fTmTGMMs_pufrac_soft;
    float fTmTGMMs_pufrac_hard;
    float fTmTGMMs_m_pu_soft;
    float fTmTGMMs_m_pu_hard;

    float fTmGMM_m;
    float fTmGMM_pt;
    float fTdeltatop_mGMM;
    float fTmGMM_m_mean;
    float fTmGMM_m_var;
    float fTmGMM_m_skew;
    float fTmGMM_pt_mean;
    float fTmGMM_pt_var;
    float fTmGMM_pt_skew;
    float fTmGMM_ptl;
    float fTmGMM_ml;
    float fTmGMM_m_soft;
    float fTmGMM_pt_soft;
    float fTmGMM_ucpu;
    float fTmGMM_clpu;
    float fTmGMM_pufrac_soft;
    float fTmGMM_pufrac_hard;
    float fTmGMM_m_pu_soft;
    float fTmGMM_m_pu_hard;

    float fTmTGMM_m;
    float fTmTGMM_pt;
    float fTdeltatop_mTGMM;
    float fTmTGMM_m_mean;
    float fTmTGMM_m_var;
    float fTmTGMM_m_skew;
    float fTmTGMM_pt_mean;
    float fTmTGMM_pt_var;
    float fTmTGMM_pt_skew;
    float fTmTGMM_ptl;
    float fTmTGMM_ml;
    float fTmTGMM_m_soft;
    float fTmTGMM_pt_soft;
    float fTmTGMM_ucpu;
    float fTmTGMM_clpu;
    float fTmTGMM_pufrac_soft;
    float fTmTGMM_pufrac_hard;
    float fTmTGMM_m_pu_soft;
    float fTmTGMM_m_pu_hard;

    fastjet::JetDefinition     *m_jet_def;
    fastjet::JetDefinition     *m_jet_def_large_r_antikt;
    fastjet::JetDefinition     *m_jet_def_large_r_ca;

 public:
    FuzzyAnalysis ();
    ~FuzzyAnalysis ();

    void Begin();
    void SetupHistosMap();
    void WriteHistosMap();
    void AnalyzeEvent(int i_evt, Pythia8::Pythia *pythia8, Pythia8::Pythia *pythia_MB, int NPV);
    void End();
    void DeclareBranches();
    void ResetBranches();
    void Debug(int debug){
        f_debug = debug;
    }
    void SetOutName(string outname){
        f_out_name = outname;
    }
    void SetRecombination(bool r) {
        do_recombination = r;
    }
    void SetPrefix(string prefix) {
        directory_prefix = prefix;
    }
    void SetSize(double s) {
        f_size = s;
    }
    void SetLearnWeights(bool w) {
        f_learn_weights = w;
    }
    void SetBatched(bool b) {
        batched = b;
    }
    void SetShouldPrint(bool s) {
        should_print = s;
    }
};

#endif
