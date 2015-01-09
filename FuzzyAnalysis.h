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

typedef std::vector<fastjet::PseudoJet> vecPseudoJet;

class FuzzyAnalysis{
 private:
    int  f_test;
    int  f_debug;
    unsigned int _n_events;
    double  f_size;
    float pT_min;
    bool f_learn_weights;
    bool batched;
    bool should_print;
    bool do_recombination;
    string f_out_name;
    string directory_prefix;

    bool do_protobuf_serialize;

    std::map<string, float*> tree_vars;
    std::map<string, std::vector<float> > map_weight_vecs;

    std::vector<float> mGMM_weight_vec;
    std::vector<float> mGMMs_weight_vec;
    std::vector<float> mGMMc_weight_vec;
    std::vector<float> mTGMM_weight_vec;
    std::vector<float> mTGMMs_weight_vec;
    std::vector<float> mUMM_weight_vec;

    std::vector<float> mGMM_distance_vec;
    std::vector<float> mGMMs_distance_vec;
    std::vector<float> mGMMc_distance_vec;
    std::vector<float> mTGMM_distance_vec;
    std::vector<float> mTGMMs_distance_vec;
    std::vector<float> mUMM_distance_vec;

    std::vector<bool> is_pileup_vec;
    std::vector<bool> is_lead_antikt_constituent_vec;
    std::vector<bool> is_lead_antikt_trimmed_two_constituent_vec;

    #ifdef WITHROOT
    TFile *t_f;
    TTree *t_t;
    #endif

    FuzzyTools *tool;

    double _event_jet_strength;
    double _event_jet_rho_offset;
    FuzzyTools::EventJet _event_jet_type;
    FuzzyTools::PostProcess _post_processing_method;

    bool _do_tower_subtraction;

    // Tree Vars ---------------------------------------
    int   fTEventNumber;
    int   fTNPV;
    float fTtoppt;
    float fTrho;
    float fTevent_jet_strength;
    float fTevent_jet_rho_offset;

    float fTCA_pt;
    float fTCA_m;
    float fTCA_pufrac;
    float fTCA_m_pu;
    float fTCA_dr;
    float fTCA_area;
    std::vector<float> CA_nsubjettiness;

    float fTantikt_m;
    float fTantikt_m_pu;
    float fTantikt_pt;
    float fTantikt_pufrac;
    float fTantikt_dr;
    float fTantikt_area;
    std::vector<float> antikt_nsubjettiness;

    float fTantikt_m_trimmed_two;
    float fTantikt_m_pu_trimmed_two;
    float fTantikt_pt_trimmed_two;
    float fTantikt_pufrac_trimmed_two;
    float fTantikt_area_trimmed_two;

    float fTantikt_m_trimmed_three;
    float fTantikt_m_pu_trimmed_three;
    float fTantikt_pt_trimmed_three;
    float fTantikt_pufrac_trimmed_three;
    float fTantikt_area_trimmed_three;

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
    float fTmUMM_dr;
    unsigned int fTmUMM_iter;
    std::vector<float> fTmUMM_etas;
    std::vector<float> fTmUMM_phis;

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
    float fTmGMMs_dr;
    unsigned int fTmGMMs_iter;
    std::vector<float> fTmGMMs_etas;
    std::vector<float> fTmGMMs_phis;

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
    float fTmGMMc_dr;
    unsigned int fTmGMMc_iter;
    std::vector<float> fTmGMMc_etas;
    std::vector<float> fTmGMMc_phis;

    float fTmGMMc_r;
    float fTmGMMc_r_avg;
    float fTmGMMc_r_weighted_avg;
    float fTmGMMc_r_second;
    float fTmGMMc_r_third;

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
    float fTmTGMMs_dr;
    unsigned int fTmTGMMs_iter;
    std::vector<float> fTmTGMMs_etas;
    std::vector<float> fTmTGMMs_phis;

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
    float fTmGMM_dr;
    unsigned int fTmGMM_iter;
    std::vector<float> fTmGMM_etas;
    std::vector<float> fTmGMM_phis;

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
    float fTmTGMM_dr;
    unsigned int fTmTGMM_iter;
    std::vector<float> fTmTGMM_etas;
    std::vector<float> fTmTGMM_phis;

    fastjet::JetDefinition     *m_jet_def;
    fastjet::JetDefinition     *m_jet_def_large_r_antikt;
    fastjet::JetDefinition     *m_jet_def_large_r_ca;
    fastjet::JetDefinition     *m_jet_def_vlarge_r_antikt;
    fastjet::JetDefinition     *m_jet_def_vlarge_r_ca;



 public:
    FuzzyAnalysis ();
    ~FuzzyAnalysis ();

    void Begin();
    void SetupHistosMap();
    void WriteHistosMap();
    void PiFixStudy();
    vecPseudoJet DiscretizeParticles(vecPseudoJet const& arg_particles,
                                                    int n_phi, int n_eta,
                                                    float eta_cutoff);
    vecPseudoJet ChargedParticles(vecPseudoJet const& arg_particles, double eta_max);
    vecPseudoJet UnchargedParticles(vecPseudoJet const& arg_particles, double eta_max);
    void GroomingStudy(vecPseudoJet leading_particles, int event_iter);
    void SubstructureStudy(vecPseudoJet ca_jets, vecPseudoJet antikt_jets, int event_iter);
    void AnalyzeEventNew(int i_evt, Pythia8::Pythia *pythia8, Pythia8::Pythia *pythia_MB, int NPV);
    void AnalyzeEvent(int i_evt, Pythia8::Pythia *pythia8, Pythia8::Pythia *pythia_MB, int NPV);
    void End();
    void DeclareBranches();
    void ResetBranches();
    void Debug(int debug){
        f_debug = debug;
    }

    bool QueryFloatValue(std::string name, float *f_out) {
        std::map<string, float*>::iterator it = tree_vars.find(name);
        if (it == tree_vars.end()) {
            // no such value in tree
            return false;
        } else {
            // found it! write to output parameter
            *f_out = *(it->second);
            return true;
        }
    }

    void SetNEvents(unsigned int n) {
        _n_events = n;
        if (tool) {
            tool->SetNEvents(n);
        }
    }

    void SetDoTowerSubtraction(bool d) {
        _do_tower_subtraction = d;
    }

    void SetEventJetRhoOffset(double off) {
        _event_jet_rho_offset = off;
    }

    void SetEventJetType(FuzzyTools::EventJet ejt)  {
        _event_jet_type = ejt;
    }

    void SetPostProcessingMethod(FuzzyTools::PostProcess pj) {
        _post_processing_method = pj;
    }

    void SetEventJetStrength(float s) {
        _event_jet_strength = s;
    }

    void SetPtMin(float p) {
        pT_min = p;
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
    bool IsBatched() {
        return batched;
    }

    void SetProtobufSerialize(bool p) {
        do_protobuf_serialize = p;
    }

    void SetShouldPrint(bool s) {
        should_print = s;
    }
};

#endif
