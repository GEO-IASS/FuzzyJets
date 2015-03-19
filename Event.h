#ifndef EVENT_H
#define EVENT_H

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

#include <vector>
#include <string>
#include <assert.h>

#include "Histogram.h"

#include <boost/algorithm/string.hpp>

class EventBuffer {
public:
    TTree          *_tree;

    // ouchies
    std::map<std::string, void *> _raw_member_locations;

    // Declaration of leaf types
    Int_t           EventNumber;
    Int_t           NPV;
    Int_t           n_jet_seeds;

    std::vector<float>   *mGMM_weight_vec;
    std::vector<float>   *mGMMs_weight_vec;
    std::vector<float>   *mGMMc_weight_vec;
    std::vector<float>   *mTGMM_weight_vec;
    std::vector<float>   *mTGMMs_weight_vec;
    std::vector<float>   *mUMM_weight_vec;
    std::vector<float>   *mGMM_distance_vec;
    std::vector<float>   *mGMMs_distance_vec;
    std::vector<float>   *mGMMc_distance_vec;
    std::vector<float>   *mTGMM_distance_vec;
    std::vector<float>   *mTGMMs_distance_vec;
    std::vector<float>   *mUMM_distance_vec;
    std::vector<bool>    *is_pileup_vec;
    std::vector<bool>    *is_lead_antikt_constituent_vec;
    std::vector<bool>    *is_lead_antikt_trimmed_two_constituent_vec;

    std::vector<float>   *CA_nsubjettiness;
    std::vector<float>   *antikt_nsubjettiness;

    std::vector<float>   *mGMM_etas;
    std::vector<float>   *mGMM_phis;
    std::vector<float>   *mGMMc_etas;
    std::vector<float>   *mGMMc_phis;
    std::vector<float>   *mGMMs_etas;
    std::vector<float>   *mGMMs_phis;
    std::vector<float>   *mTGMM_etas;
    std::vector<float>   *mTGMM_phis;
    std::vector<float>   *mTGMMs_etas;
    std::vector<float>   *mTGMMs_phis;
    std::vector<float>   *mUMM_etas;
    std::vector<float>   *mUMM_phis;

    Float_t         CA_area;
    Float_t         antikt_area;
    Float_t         antikt_area_trimmed_two;
    Float_t         antikt_area_trimmed_three;

    Float_t         CA_dr;
    Float_t         CA_m;
    Float_t         CA_m_pu;
    Float_t         CA_pt;
    Float_t         CA_pufrac;
    Float_t         antikt_dr;
    Float_t         antikt_m;
    Float_t         antikt_m_pu;
    Float_t         antikt_m_pu_trimmed_three;
    Float_t         antikt_m_pu_trimmed_two;
    Float_t         antikt_m_trimmed_three;
    Float_t         antikt_m_trimmed_two;
    Float_t         antikt_pt;
    Float_t         antikt_pt_trimmed_three;
    Float_t         antikt_pt_trimmed_two;
    Float_t         antikt_pufrac;
    Float_t         antikt_pufrac_trimmed_three;
    Float_t         antikt_pufrac_trimmed_two;
    Float_t         deltatop_mGMM;
    Float_t         deltatop_mGMMc;
    Float_t         deltatop_mGMMs;
    Float_t         deltatop_mTGMM;
    Float_t         deltatop_mTGMMs;
    Float_t         deltatop_mUMM;
    Float_t         mGMM_clpu;
    Float_t         mGMM_dr;
    Float_t         mGMM_m;
    Float_t         mGMM_m_mean;
    Float_t         mGMM_m_pu_hard;
    Float_t         mGMM_m_pu_soft;
    Float_t         mGMM_m_skew;
    Float_t         mGMM_m_soft;
    Float_t         mGMM_m_var;
    Float_t         mGMM_ml;
    Float_t         mGMM_pt;
    Float_t         mGMM_pt_mean;
    Float_t         mGMM_pt_skew;
    Float_t         mGMM_pt_soft;
    Float_t         mGMM_pt_var;
    Float_t         mGMM_ptl;
    Float_t         mGMM_pufrac_hard;
    Float_t         mGMM_pufrac_soft;
    Float_t         mGMM_ucpu;
    Float_t         mGMMc_clpu;
    Float_t         mGMMc_dr;
    Float_t         mGMMc_m;
    Float_t         mGMMc_m_mean;
    Float_t         mGMMc_m_pu_hard;
    Float_t         mGMMc_m_pu_soft;
    Float_t         mGMMc_m_skew;
    Float_t         mGMMc_m_soft;
    Float_t         mGMMc_m_var;
    Float_t         mGMMc_ml;
    Float_t         mGMMc_pt;
    Float_t         mGMMc_pt_mean;
    Float_t         mGMMc_pt_skew;
    Float_t         mGMMc_pt_soft;
    Float_t         mGMMc_pt_var;
    Float_t         mGMMc_ptl;
    Float_t         mGMMc_pufrac_hard;
    Float_t         mGMMc_pufrac_soft;
    Float_t         mGMMc_r;
    Float_t         mGMMc_r_avg;
    Float_t         mGMMc_r_second;
    Float_t         mGMMc_r_third;
    Float_t         mGMMc_r_weighted_avg;
    Float_t         mGMMc_ucpu;
    Float_t         mGMMs_clpu;
    Float_t         mGMMs_dr;
    Float_t         mGMMs_m;
    Float_t         mGMMs_m_mean;
    Float_t         mGMMs_m_pu_hard;
    Float_t         mGMMs_m_pu_soft;
    Float_t         mGMMs_m_skew;
    Float_t         mGMMs_m_soft;
    Float_t         mGMMs_m_var;
    Float_t         mGMMs_ml;
    Float_t         mGMMs_pt;
    Float_t         mGMMs_pt_mean;
    Float_t         mGMMs_pt_skew;
    Float_t         mGMMs_pt_soft;
    Float_t         mGMMs_pt_var;
    Float_t         mGMMs_ptl;
    Float_t         mGMMs_pufrac_hard;
    Float_t         mGMMs_pufrac_soft;
    Float_t         mGMMs_ucpu;
    Float_t         mTGMM_clpu;
    Float_t         mTGMM_dr;
    Float_t         mTGMM_m;
    Float_t         mTGMM_m_mean;
    Float_t         mTGMM_m_pu_hard;
    Float_t         mTGMM_m_pu_soft;
    Float_t         mTGMM_m_skew;
    Float_t         mTGMM_m_soft;
    Float_t         mTGMM_m_var;
    Float_t         mTGMM_ml;
    Float_t         mTGMM_pt;
    Float_t         mTGMM_pt_mean;
    Float_t         mTGMM_pt_skew;
    Float_t         mTGMM_pt_soft;
    Float_t         mTGMM_pt_var;
    Float_t         mTGMM_ptl;
    Float_t         mTGMM_pufrac_hard;
    Float_t         mTGMM_pufrac_soft;
    Float_t         mTGMM_ucpu;
    Float_t         mTGMMs_clpu;
    Float_t         mTGMMs_dr;
    Float_t         mTGMMs_m;
    Float_t         mTGMMs_m_mean;
    Float_t         mTGMMs_m_pu_hard;
    Float_t         mTGMMs_m_pu_soft;
    Float_t         mTGMMs_m_skew;
    Float_t         mTGMMs_m_soft;
    Float_t         mTGMMs_m_var;
    Float_t         mTGMMs_ml;
    Float_t         mTGMMs_pt;
    Float_t         mTGMMs_pt_mean;
    Float_t         mTGMMs_pt_skew;
    Float_t         mTGMMs_pt_soft;
    Float_t         mTGMMs_pt_var;
    Float_t         mTGMMs_ptl;
    Float_t         mTGMMs_pufrac_hard;
    Float_t         mTGMMs_pufrac_soft;
    Float_t         mTGMMs_ucpu;
    Float_t         mUMM_clpu;
    Float_t         mUMM_dr;
    Float_t         mUMM_m;
    Float_t         mUMM_m_mean;
    Float_t         mUMM_m_pu_hard;
    Float_t         mUMM_m_pu_soft;
    Float_t         mUMM_m_skew;
    Float_t         mUMM_m_soft;
    Float_t         mUMM_m_var;
    Float_t         mUMM_ml;
    Float_t         mUMM_pt;
    Float_t         mUMM_pt_mean;
    Float_t         mUMM_pt_skew;
    Float_t         mUMM_pt_soft;
    Float_t         mUMM_pt_var;
    Float_t         mUMM_ptl;
    Float_t         mUMM_pufrac_hard;
    Float_t         mUMM_pufrac_soft;
    Float_t         mUMM_ucpu;
    Float_t         toppt;

    // List of branches
    TBranch        *b_EventNumber;   //!
    TBranch        *b_NPV;   //!
    TBranch        *b_n_jet_seeds; //!

    TBranch        *b_mGMM_weight_vec;   //!
    TBranch        *b_mGMMs_weight_vec;   //!
    TBranch        *b_mGMMc_weight_vec;   //!
    TBranch        *b_mTGMM_weight_vec;   //!
    TBranch        *b_mTGMMs_weight_vec;   //!
    TBranch        *b_mUMM_weight_vec;   //!
    TBranch        *b_mGMM_distance_vec;   //!
    TBranch        *b_mGMMs_distance_vec;   //!
    TBranch        *b_mGMMc_distance_vec;   //!
    TBranch        *b_mTGMM_distance_vec;   //!
    TBranch        *b_mTGMMs_distance_vec;   //!
    TBranch        *b_mUMM_distance_vec;   //!
    TBranch        *b_is_pileup_vec;   //!
    TBranch        *b_is_lead_antikt_constituent_vec;   //!
    TBranch        *b_is_lead_antikt_trimmed_two_constituent_vec;   //!

    TBranch        *b_mGMM_etas;
    TBranch        *b_mGMM_phis;
    TBranch        *b_mGMMc_etas;
    TBranch        *b_mGMMc_phis;
    TBranch        *b_mGMMs_etas;
    TBranch        *b_mGMMs_phis;
    TBranch        *b_mTGMM_etas;
    TBranch        *b_mTGMM_phis;
    TBranch        *b_mTGMMs_etas;
    TBranch        *b_mTGMMs_phis;
    TBranch        *b_mUMM_etas;
    TBranch        *b_mUMM_phis;

    TBranch        *b_CA_nsubjettiness;
    TBranch        *b_antikt_nsubjettiness;

    TBranch        *b_CA_area;
    TBranch        *b_antikt_area;
    TBranch        *b_antikt_area_trimmed_two;
    TBranch        *b_antikt_area_trimmed_three;

    TBranch        *b_CA_dr;   //!
    TBranch        *b_CA_m;   //!
    TBranch        *b_CA_m_pu;   //!
    TBranch        *b_CA_pt;   //!
    TBranch        *b_CA_pufrac;   //!
    TBranch        *b_antikt_dr;   //!
    TBranch        *b_antikt_m;   //!
    TBranch        *b_antikt_m_pu;   //!
    TBranch        *b_antikt_m_pu_trimmed_three;   //!
    TBranch        *b_antikt_m_pu_trimmed_two;   //!
    TBranch        *b_antikt_m_trimmed_three;   //!
    TBranch        *b_antikt_m_trimmed_two;   //!
    TBranch        *b_antikt_pt;   //!
    TBranch        *b_antikt_pt_trimmed_three;   //!
    TBranch        *b_antikt_pt_trimmed_two;   //!
    TBranch        *b_antikt_pufrac;   //!
    TBranch        *b_antikt_pufrac_trimmed_three;   //!
    TBranch        *b_antikt_pufrac_trimmed_two;   //!
    TBranch        *b_deltatop_mGMM;   //!
    TBranch        *b_deltatop_mGMMc;   //!
    TBranch        *b_deltatop_mGMMs;   //!
    TBranch        *b_deltatop_mTGMM;   //!
    TBranch        *b_deltatop_mTGMMs;   //!
    TBranch        *b_deltatop_mUMM;   //!
    TBranch        *b_mGMM_clpu;   //!
    TBranch        *b_mGMM_dr;   //!
    TBranch        *b_mGMM_m;   //!
    TBranch        *b_mGMM_m_mean;   //!
    TBranch        *b_mGMM_m_pu_hard;   //!
    TBranch        *b_mGMM_m_pu_soft;   //!
    TBranch        *b_mGMM_m_skew;   //!
    TBranch        *b_mGMM_m_soft;   //!
    TBranch        *b_mGMM_m_var;   //!
    TBranch        *b_mGMM_ml;   //!
    TBranch        *b_mGMM_pt;   //!
    TBranch        *b_mGMM_pt_mean;   //!
    TBranch        *b_mGMM_pt_skew;   //!
    TBranch        *b_mGMM_pt_soft;   //!
    TBranch        *b_mGMM_pt_var;   //!
    TBranch        *b_mGMM_ptl;   //!
    TBranch        *b_mGMM_pufrac_hard;   //!
    TBranch        *b_mGMM_pufrac_soft;   //!
    TBranch        *b_mGMM_ucpu;   //!
    TBranch        *b_mGMMc_clpu;   //!
    TBranch        *b_mGMMc_dr;   //!
    TBranch        *b_mGMMc_m;   //!
    TBranch        *b_mGMMc_m_mean;   //!
    TBranch        *b_mGMMc_m_pu_hard;   //!
    TBranch        *b_mGMMc_m_pu_soft;   //!
    TBranch        *b_mGMMc_m_skew;   //!
    TBranch        *b_mGMMc_m_soft;   //!
    TBranch        *b_mGMMc_m_var;   //!
    TBranch        *b_mGMMc_ml;   //!
    TBranch        *b_mGMMc_pt;   //!
    TBranch        *b_mGMMc_pt_mean;   //!
    TBranch        *b_mGMMc_pt_skew;   //!
    TBranch        *b_mGMMc_pt_soft;   //!
    TBranch        *b_mGMMc_pt_var;   //!
    TBranch        *b_mGMMc_ptl;   //!
    TBranch        *b_mGMMc_pufrac_hard;   //!
    TBranch        *b_mGMMc_pufrac_soft;   //!
    TBranch        *b_mGMMc_r;   //!
    TBranch        *b_mGMMc_r_avg;   //!
    TBranch        *b_mGMMc_r_second;   //!
    TBranch        *b_mGMMc_r_third;   //!
    TBranch        *b_mGMMc_r_weighted_avg;   //!
    TBranch        *b_mGMMc_ucpu;   //!
    TBranch        *b_mGMMs_clpu;   //!
    TBranch        *b_mGMMs_dr;   //!
    TBranch        *b_mGMMs_m;   //!
    TBranch        *b_mGMMs_m_mean;   //!
    TBranch        *b_mGMMs_m_pu_hard;   //!
    TBranch        *b_mGMMs_m_pu_soft;   //!
    TBranch        *b_mGMMs_m_skew;   //!
    TBranch        *b_mGMMs_m_soft;   //!
    TBranch        *b_mGMMs_m_var;   //!
    TBranch        *b_mGMMs_ml;   //!
    TBranch        *b_mGMMs_pt;   //!
    TBranch        *b_mGMMs_pt_mean;   //!
    TBranch        *b_mGMMs_pt_skew;   //!
    TBranch        *b_mGMMs_pt_soft;   //!
    TBranch        *b_mGMMs_pt_var;   //!
    TBranch        *b_mGMMs_ptl;   //!
    TBranch        *b_mGMMs_pufrac_hard;   //!
    TBranch        *b_mGMMs_pufrac_soft;   //!
    TBranch        *b_mGMMs_ucpu;   //!
    TBranch        *b_mTGMM_clpu;   //!
    TBranch        *b_mTGMM_dr;   //!
    TBranch        *b_mTGMM_m;   //!
    TBranch        *b_mTGMM_m_mean;   //!
    TBranch        *b_mTGMM_m_pu_hard;   //!
    TBranch        *b_mTGMM_m_pu_soft;   //!
    TBranch        *b_mTGMM_m_skew;   //!
    TBranch        *b_mTGMM_m_soft;   //!
    TBranch        *b_mTGMM_m_var;   //!
    TBranch        *b_mTGMM_ml;   //!
    TBranch        *b_mTGMM_pt;   //!
    TBranch        *b_mTGMM_pt_mean;   //!
    TBranch        *b_mTGMM_pt_skew;   //!
    TBranch        *b_mTGMM_pt_soft;   //!
    TBranch        *b_mTGMM_pt_var;   //!
    TBranch        *b_mTGMM_ptl;   //!
    TBranch        *b_mTGMM_pufrac_hard;   //!
    TBranch        *b_mTGMM_pufrac_soft;   //!
    TBranch        *b_mTGMM_ucpu;   //!
    TBranch        *b_mTGMMs_clpu;   //!
    TBranch        *b_mTGMMs_dr;   //!
    TBranch        *b_mTGMMs_m;   //!
    TBranch        *b_mTGMMs_m_mean;   //!
    TBranch        *b_mTGMMs_m_pu_hard;   //!
    TBranch        *b_mTGMMs_m_pu_soft;   //!
    TBranch        *b_mTGMMs_m_skew;   //!
    TBranch        *b_mTGMMs_m_soft;   //!
    TBranch        *b_mTGMMs_m_var;   //!
    TBranch        *b_mTGMMs_ml;   //!
    TBranch        *b_mTGMMs_pt;   //!
    TBranch        *b_mTGMMs_pt_mean;   //!
    TBranch        *b_mTGMMs_pt_skew;   //!
    TBranch        *b_mTGMMs_pt_soft;   //!
    TBranch        *b_mTGMMs_pt_var;   //!
    TBranch        *b_mTGMMs_ptl;   //!
    TBranch        *b_mTGMMs_pufrac_hard;   //!
    TBranch        *b_mTGMMs_pufrac_soft;   //!
    TBranch        *b_mTGMMs_ucpu;   //!
    TBranch        *b_mUMM_clpu;   //!
    TBranch        *b_mUMM_dr;   //!
    TBranch        *b_mUMM_m;   //!
    TBranch        *b_mUMM_m_mean;   //!
    TBranch        *b_mUMM_m_pu_hard;   //!
    TBranch        *b_mUMM_m_pu_soft;   //!
    TBranch        *b_mUMM_m_skew;   //!
    TBranch        *b_mUMM_m_soft;   //!
    TBranch        *b_mUMM_m_var;   //!
    TBranch        *b_mUMM_ml;   //!
    TBranch        *b_mUMM_pt;   //!
    TBranch        *b_mUMM_pt_mean;   //!
    TBranch        *b_mUMM_pt_skew;   //!
    TBranch        *b_mUMM_pt_soft;   //!
    TBranch        *b_mUMM_pt_var;   //!
    TBranch        *b_mUMM_ptl;   //!
    TBranch        *b_mUMM_pufrac_hard;   //!
    TBranch        *b_mUMM_pufrac_soft;   //!
    TBranch        *b_mUMM_ucpu;   //!
    TBranch        *b_toppt;   //!

    EventBuffer() {
        _raw_member_locations["CA_area"] = &CA_area;
        _raw_member_locations["antikt_area"] = &antikt_area;
        _raw_member_locations["antikt_area_trimmed_two"] = &antikt_area_trimmed_two;
        _raw_member_locations["antikt_area_trimmed_three"] = &antikt_area_trimmed_three;

        _raw_member_locations["CA_nsubjettiness"] = &CA_nsubjettiness;
        _raw_member_locations["antikt_nsubjettiness"] = &antikt_nsubjettiness;

        _raw_member_locations["mGMM_etas"] = &mGMM_etas;
        _raw_member_locations["mGMM_phis"] = &mGMM_phis;
        _raw_member_locations["mGMMc_etas"] = &mGMMc_etas;
        _raw_member_locations["mGMMc_phis"] = &mGMMc_phis;
        _raw_member_locations["mGMMs_etas"] = &mGMMs_etas;
        _raw_member_locations["mGMMs_phis"] = &mGMMs_phis;
        _raw_member_locations["mTGMM_etas"] = &mTGMM_etas;
        _raw_member_locations["mTGMM_phis"] = &mTGMM_phis;
        _raw_member_locations["mTGMMs_etas"] = &mTGMMs_etas;
        _raw_member_locations["mTGMMs_phis"] = &mTGMMs_phis;
        _raw_member_locations["mUMM_etas"] = &mUMM_etas;
        _raw_member_locations["mUMM_phis"] = &mUMM_phis;

        _raw_member_locations["EventNumber"] = &EventNumber;
        _raw_member_locations["NPV"] = &NPV;
        _raw_member_locations["n_jet_seeds"] = &n_jet_seeds;
        _raw_member_locations["mGMM_weight_vec"] = &mGMM_weight_vec;
        _raw_member_locations["mGMMs_weight_vec"] = &mGMMs_weight_vec;
        _raw_member_locations["mGMMc_weight_vec"] = &mGMMc_weight_vec;
        _raw_member_locations["mTGMM_weight_vec"] = &mTGMM_weight_vec;
        _raw_member_locations["mTGMMs_weight_vec"] = &mTGMMs_weight_vec;
        _raw_member_locations["mUMM_weight_vec"] = &mUMM_weight_vec;
        _raw_member_locations["mGMM_distance_vec"] = &mGMM_distance_vec;
        _raw_member_locations["mGMMs_distance_vec"] = &mGMMs_distance_vec;
        _raw_member_locations["mGMMc_distance_vec"] = &mGMMc_distance_vec;
        _raw_member_locations["mTGMM_distance_vec"] = &mTGMM_distance_vec;
        _raw_member_locations["mTGMMs_distance_vec"] = &mTGMMs_distance_vec;
        _raw_member_locations["mUMM_distance_vec"] = &mUMM_distance_vec;
        _raw_member_locations["is_pileup_vec"] = &is_pileup_vec;
        _raw_member_locations["is_lead_antikt_constituent_vec"] = &is_lead_antikt_constituent_vec;
        _raw_member_locations["is_lead_antikt_trimmed_two_constituent_vec"] = &is_lead_antikt_trimmed_two_constituent_vec;
        _raw_member_locations["CA_dr"] = &CA_dr;
        _raw_member_locations["CA_m"] = &CA_m;
        _raw_member_locations["CA_m_pu"] = &CA_m_pu;
        _raw_member_locations["CA_pt"] = &CA_pt;
        _raw_member_locations["CA_pufrac"] = &CA_pufrac;
        _raw_member_locations["antikt_dr"] = &antikt_dr;
        _raw_member_locations["antikt_m"] = &antikt_m;
        _raw_member_locations["antikt_m_pu"] = &antikt_m_pu;
        _raw_member_locations["antikt_m_pu_trimmed_three"] = &antikt_m_pu_trimmed_three;
        _raw_member_locations["antikt_m_pu_trimmed_two"] = &antikt_m_pu_trimmed_two;
        _raw_member_locations["antikt_m_trimmed_three"] = &antikt_m_trimmed_three;
        _raw_member_locations["antikt_m_trimmed_two"] = &antikt_m_trimmed_two;
        _raw_member_locations["antikt_pt"] = &antikt_pt;
        _raw_member_locations["antikt_pt_trimmed_three"] = &antikt_pt_trimmed_three;
        _raw_member_locations["antikt_pt_trimmed_two"] = &antikt_pt_trimmed_two;
        _raw_member_locations["antikt_pufrac"] = &antikt_pufrac;
        _raw_member_locations["antikt_pufrac_trimmed_three"] = &antikt_pufrac_trimmed_three;
        _raw_member_locations["antikt_pufrac_trimmed_two"] = &antikt_pufrac_trimmed_two;
        _raw_member_locations["deltatop_mGMM"] = &deltatop_mGMM;
        _raw_member_locations["deltatop_mGMMc"] = &deltatop_mGMMc;
        _raw_member_locations["deltatop_mGMMs"] = &deltatop_mGMMs;
        _raw_member_locations["deltatop_mTGMM"] = &deltatop_mTGMM;
        _raw_member_locations["deltatop_mTGMMs"] = &deltatop_mTGMMs;
        _raw_member_locations["deltatop_mUMM"] = &deltatop_mUMM;
        _raw_member_locations["mGMM_clpu"] = &mGMM_clpu;
        _raw_member_locations["mGMM_dr"] = &mGMM_dr;
        _raw_member_locations["mGMM_m"] = &mGMM_m;
        _raw_member_locations["mGMM_m_mean"] = &mGMM_m_mean;
        _raw_member_locations["mGMM_m_pu_hard"] = &mGMM_m_pu_hard;
        _raw_member_locations["mGMM_m_pu_soft"] = &mGMM_m_pu_soft;
        _raw_member_locations["mGMM_m_skew"] = &mGMM_m_skew;
        _raw_member_locations["mGMM_m_soft"] = &mGMM_m_soft;
        _raw_member_locations["mGMM_m_var"] = &mGMM_m_var;
        _raw_member_locations["mGMM_ml"] = &mGMM_ml;
        _raw_member_locations["mGMM_pt"] = &mGMM_pt;
        _raw_member_locations["mGMM_pt_mean"] = &mGMM_pt_mean;
        _raw_member_locations["mGMM_pt_skew"] = &mGMM_pt_skew;
        _raw_member_locations["mGMM_pt_soft"] = &mGMM_pt_soft;
        _raw_member_locations["mGMM_pt_var"] = &mGMM_pt_var;
        _raw_member_locations["mGMM_ptl"] = &mGMM_ptl;
        _raw_member_locations["mGMM_pufrac_hard"] = &mGMM_pufrac_hard;
        _raw_member_locations["mGMM_pufrac_soft"] = &mGMM_pufrac_soft;
        _raw_member_locations["mGMM_ucpu"] = &mGMM_ucpu;
        _raw_member_locations["mGMMc_clpu"] = &mGMMc_clpu;
        _raw_member_locations["mGMMc_dr"] = &mGMMc_dr;
        _raw_member_locations["mGMMc_m"] = &mGMMc_m;
        _raw_member_locations["mGMMc_m_mean"] = &mGMMc_m_mean;
        _raw_member_locations["mGMMc_m_pu_hard"] = &mGMMc_m_pu_hard;
        _raw_member_locations["mGMMc_m_pu_soft"] = &mGMMc_m_pu_soft;
        _raw_member_locations["mGMMc_m_skew"] = &mGMMc_m_skew;
        _raw_member_locations["mGMMc_m_soft"] = &mGMMc_m_soft;
        _raw_member_locations["mGMMc_m_var"] = &mGMMc_m_var;
        _raw_member_locations["mGMMc_ml"] = &mGMMc_ml;
        _raw_member_locations["mGMMc_pt"] = &mGMMc_pt;
        _raw_member_locations["mGMMc_pt_mean"] = &mGMMc_pt_mean;
        _raw_member_locations["mGMMc_pt_skew"] = &mGMMc_pt_skew;
        _raw_member_locations["mGMMc_pt_soft"] = &mGMMc_pt_soft;
        _raw_member_locations["mGMMc_pt_var"] = &mGMMc_pt_var;
        _raw_member_locations["mGMMc_ptl"] = &mGMMc_ptl;
        _raw_member_locations["mGMMc_pufrac_hard"] = &mGMMc_pufrac_hard;
        _raw_member_locations["mGMMc_pufrac_soft"] = &mGMMc_pufrac_soft;
        _raw_member_locations["mGMMc_r"] = &mGMMc_r;
        _raw_member_locations["mGMMc_r_avg"] = &mGMMc_r_avg;
        _raw_member_locations["mGMMc_r_second"] = &mGMMc_r_second;
        _raw_member_locations["mGMMc_r_third"] = &mGMMc_r_third;
        _raw_member_locations["mGMMc_r_weighted_avg"] = &mGMMc_r_weighted_avg;
        _raw_member_locations["mGMMc_ucpu"] = &mGMMc_ucpu;
        _raw_member_locations["mGMMs_clpu"] = &mGMMs_clpu;
        _raw_member_locations["mGMMs_dr"] = &mGMMs_dr;
        _raw_member_locations["mGMMs_m"] = &mGMMs_m;
        _raw_member_locations["mGMMs_m_mean"] = &mGMMs_m_mean;
        _raw_member_locations["mGMMs_m_pu_hard"] = &mGMMs_m_pu_hard;
        _raw_member_locations["mGMMs_m_pu_soft"] = &mGMMs_m_pu_soft;
        _raw_member_locations["mGMMs_m_skew"] = &mGMMs_m_skew;
        _raw_member_locations["mGMMs_m_soft"] = &mGMMs_m_soft;
        _raw_member_locations["mGMMs_m_var"] = &mGMMs_m_var;
        _raw_member_locations["mGMMs_ml"] = &mGMMs_ml;
        _raw_member_locations["mGMMs_pt"] = &mGMMs_pt;
        _raw_member_locations["mGMMs_pt_mean"] = &mGMMs_pt_mean;
        _raw_member_locations["mGMMs_pt_skew"] = &mGMMs_pt_skew;
        _raw_member_locations["mGMMs_pt_soft"] = &mGMMs_pt_soft;
        _raw_member_locations["mGMMs_pt_var"] = &mGMMs_pt_var;
        _raw_member_locations["mGMMs_ptl"] = &mGMMs_ptl;
        _raw_member_locations["mGMMs_pufrac_hard"] = &mGMMs_pufrac_hard;
        _raw_member_locations["mGMMs_pufrac_soft"] = &mGMMs_pufrac_soft;
        _raw_member_locations["mGMMs_ucpu"] = &mGMMs_ucpu;
        _raw_member_locations["mTGMM_clpu"] = &mTGMM_clpu;
        _raw_member_locations["mTGMM_dr"] = &mTGMM_dr;
        _raw_member_locations["mTGMM_m"] = &mTGMM_m;
        _raw_member_locations["mTGMM_m_mean"] = &mTGMM_m_mean;
        _raw_member_locations["mTGMM_m_pu_hard"] = &mTGMM_m_pu_hard;
        _raw_member_locations["mTGMM_m_pu_soft"] = &mTGMM_m_pu_soft;
        _raw_member_locations["mTGMM_m_skew"] = &mTGMM_m_skew;
        _raw_member_locations["mTGMM_m_soft"] = &mTGMM_m_soft;
        _raw_member_locations["mTGMM_m_var"] = &mTGMM_m_var;
        _raw_member_locations["mTGMM_ml"] = &mTGMM_ml;
        _raw_member_locations["mTGMM_pt"] = &mTGMM_pt;
        _raw_member_locations["mTGMM_pt_mean"] = &mTGMM_pt_mean;
        _raw_member_locations["mTGMM_pt_skew"] = &mTGMM_pt_skew;
        _raw_member_locations["mTGMM_pt_soft"] = &mTGMM_pt_soft;
        _raw_member_locations["mTGMM_pt_var"] = &mTGMM_pt_var;
        _raw_member_locations["mTGMM_ptl"] = &mTGMM_ptl;
        _raw_member_locations["mTGMM_pufrac_hard"] = &mTGMM_pufrac_hard;
        _raw_member_locations["mTGMM_pufrac_soft"] = &mTGMM_pufrac_soft;
        _raw_member_locations["mTGMM_ucpu"] = &mTGMM_ucpu;
        _raw_member_locations["mTGMMs_clpu"] = &mTGMMs_clpu;
        _raw_member_locations["mTGMMs_dr"] = &mTGMMs_dr;
        _raw_member_locations["mTGMMs_m"] = &mTGMMs_m;
        _raw_member_locations["mTGMMs_m_mean"] = &mTGMMs_m_mean;
        _raw_member_locations["mTGMMs_m_pu_hard"] = &mTGMMs_m_pu_hard;
        _raw_member_locations["mTGMMs_m_pu_soft"] = &mTGMMs_m_pu_soft;
        _raw_member_locations["mTGMMs_m_skew"] = &mTGMMs_m_skew;
        _raw_member_locations["mTGMMs_m_soft"] = &mTGMMs_m_soft;
        _raw_member_locations["mTGMMs_m_var"] = &mTGMMs_m_var;
        _raw_member_locations["mTGMMs_ml"] = &mTGMMs_ml;
        _raw_member_locations["mTGMMs_pt"] = &mTGMMs_pt;
        _raw_member_locations["mTGMMs_pt_mean"] = &mTGMMs_pt_mean;
        _raw_member_locations["mTGMMs_pt_skew"] = &mTGMMs_pt_skew;
        _raw_member_locations["mTGMMs_pt_soft"] = &mTGMMs_pt_soft;
        _raw_member_locations["mTGMMs_pt_var"] = &mTGMMs_pt_var;
        _raw_member_locations["mTGMMs_ptl"] = &mTGMMs_ptl;
        _raw_member_locations["mTGMMs_pufrac_hard"] = &mTGMMs_pufrac_hard;
        _raw_member_locations["mTGMMs_pufrac_soft"] = &mTGMMs_pufrac_soft;
        _raw_member_locations["mTGMMs_ucpu"] = &mTGMMs_ucpu;
        _raw_member_locations["mUMM_clpu"] = &mUMM_clpu;
        _raw_member_locations["mUMM_dr"] = &mUMM_dr;
        _raw_member_locations["mUMM_m"] = &mUMM_m;
        _raw_member_locations["mUMM_m_mean"] = &mUMM_m_mean;
        _raw_member_locations["mUMM_m_pu_hard"] = &mUMM_m_pu_hard;
        _raw_member_locations["mUMM_m_pu_soft"] = &mUMM_m_pu_soft;
        _raw_member_locations["mUMM_m_skew"] = &mUMM_m_skew;
        _raw_member_locations["mUMM_m_soft"] = &mUMM_m_soft;
        _raw_member_locations["mUMM_m_var"] = &mUMM_m_var;
        _raw_member_locations["mUMM_ml"] = &mUMM_ml;
        _raw_member_locations["mUMM_pt"] = &mUMM_pt;
        _raw_member_locations["mUMM_pt_mean"] = &mUMM_pt_mean;
        _raw_member_locations["mUMM_pt_skew"] = &mUMM_pt_skew;
        _raw_member_locations["mUMM_pt_soft"] = &mUMM_pt_soft;
        _raw_member_locations["mUMM_pt_var"] = &mUMM_pt_var;
        _raw_member_locations["mUMM_ptl"] = &mUMM_ptl;
        _raw_member_locations["mUMM_pufrac_hard"] = &mUMM_pufrac_hard;
        _raw_member_locations["mUMM_pufrac_soft"] = &mUMM_pufrac_soft;
        _raw_member_locations["mUMM_ucpu"] = &mUMM_ucpu;
        _raw_member_locations["toppt"] = &toppt;
    }

    ~EventBuffer() { }

    template <typename T>
    const T Get(std::string name) const {
        // VERY DANGEROUS! BE CAREFUL USING THIS
        return *(T *) _raw_member_locations.find(name)->second;
    }

    template <typename T>
    const T GetComposite(std::string name) const {
        // VERY DANGEROUS! BE CAREFUL USING THIS
        std::vector<std::string> components;
        boost::split(components, name, boost::is_any_of("/"));
        if (components.size() == 1) {
            std::vector<std::string> components_mul;
            boost::split(components_mul, name, boost::is_any_of("*"));
            if (components_mul.size() == 1) {
                std::vector<std::string> components_2;
                boost::split(components_2, name, boost::is_any_of(":"));
                if (components_2.size() == 1) {
                    return *(T *) _raw_member_locations.find(name)->second;
                } else {
                    int p_idx = std::stoi(components_2.at(1));
                    std::vector<T> *v = Get<std::vector<T> *>(components_2.at(0));
                    return v->at(p_idx);
                }
            } else {
                T a = GetComposite<T>(components_mul.at(0));
                T b = GetComposite<T>(components_mul.at(1));
                return a * b;
            }
        } else {
            // assume there are only two for now
            T a = GetComposite<T>(components.at(0));
            T b = GetComposite<T>(components.at(1));
            return a/b;
        }

    }

    void SetTree(TTree *tree) { _tree = tree; }
    void    Init();
    Int_t   LoadEvent(Long64_t entry, Int_t getall = 0) {
        return _tree ? _tree->GetTree()->GetEntry(entry, getall) : 0;
    }
};

class EventManager {
protected:
    std::vector<UpdatesOnEvent*> _to_update;

    std::vector<EventBuffer*> _event_buffers;
    std::vector<std::string> _event_buffer_locations;
    std::vector<std::string> _event_buffer_labels;
    std::map<std::string, EventBuffer*> _event_buffer_map;

    std::vector<TH1F *> _pT_spectra;
    std::map<std::string, TH1F *> _pT_spectrum_map;

    unsigned int _n_pT_bins;

    float _pT_low;
    float _pT_high;

    std::string _reweight_rel;

public:
    unsigned int _n_events;
    unsigned int _event_iter;

    bool _do_reweighting;

    EventManager() {
        _n_events = 0;

        _do_reweighting = false;

        _reweight_rel = "";

        _n_pT_bins = 100;
        _pT_low = 0;
        _pT_high = 500;
    }

    ~EventManager() {
        for (unsigned int iter = 0; iter < _to_update.size(); iter++) {
            delete _to_update[iter];
        }
        for (unsigned int iter = 0; iter < _pT_spectra.size(); iter++) {
            delete _pT_spectra[iter];
        }
        for (unsigned int iter = 0; iter < _event_buffers.size(); iter++) {
            delete _event_buffers[iter];
        }
    }

    float Reweight(std::string event_name) const {
        if (!_do_reweighting) return 1;
        if (event_name == _reweight_rel) return 1;
        auto it = _event_buffer_map.find(event_name);
        if (it == _event_buffer_map.end()) {
            assert(0 && "No event found with supplied name.");
        }

        TH1F *pT_spectrum = _pT_spectrum_map.find(event_name)->second;
        double pT = it->second->antikt_pt;
        Int_t bin_number = pT_spectrum->FindBin(pT);
        return _pT_spectrum_map.find(_reweight_rel)->second->GetBinContent(bin_number) /
            pT_spectrum->GetBinContent(bin_number);
    }

    void PreparePtReweighting(std::string relative) {
        assert(_n_events);

        _reweight_rel = relative;
        _do_reweighting = true;

        for(unsigned int iter = 0; iter < _event_buffers.size(); iter++) {
            _event_buffers.at(iter)->_tree->SetBranchStatus("*", 0);
            _event_buffers.at(iter)->_tree->SetBranchStatus("antikt_pt", 1);
        }

        for (unsigned int temp_event_iter = 0; temp_event_iter < _n_events; temp_event_iter++) {
            for(unsigned int iter = 0; iter < _event_buffers.size(); iter++) {
                _event_buffers.at(iter)->LoadEvent(temp_event_iter);
                _pT_spectra.at(iter)->Fill(_event_buffers.at(iter)->antikt_pt);
            }
        }

        for(unsigned int iter = 0; iter < _event_buffers.size(); iter++) {
            _event_buffers.at(iter)->_tree->SetBranchStatus("*", 1);
        }
    }

    EventBuffer const& GetEventByName(std::string name) const {
        auto it = _event_buffer_map.find(name);
        if (it != _event_buffer_map.end()) {
            return *(it->second);
        }
        assert(0 && "No event found with supplied name.");

        // shuts up gcc
        return *(_event_buffers.at(0));
    }

    void SetEventCount(unsigned int n_events) {
        _n_events = n_events;
        _event_iter = 0;
    }

    bool NextEvent() {
        if (_n_events == _event_iter) return false;
        if (!(_event_iter % 1000)) std::cout << "Event " << _event_iter << std::endl;
        LoadEvent(_event_iter);
        _event_iter++;
        return true;
    }

    void UpdateUpdaters() {
        for (unsigned int iter = 0; iter < _to_update.size(); iter++) {
            UpdatesOnEvent *next_to_update = _to_update[iter];
            next_to_update->Update(this);
        }
    }

    void InsertUpdater(UpdatesOnEvent* to_insert) {
        _to_update.push_back(to_insert);
    }

    void StartUpdaters() {
        for (unsigned int iter = 0; iter < _to_update.size(); iter++) {
            UpdatesOnEvent *next_to_start = _to_update[iter];
            next_to_start->Start(this);
        }
    }

    void FinishUpdaters() {
        for (unsigned int iter = 0; iter < _to_update.size(); iter++) {
            UpdatesOnEvent *next_to_finish = _to_update[iter];
            next_to_finish->Finish(this);
        }
    }

    void Init() {
        for (unsigned int iter = 0; iter < _event_buffers.size(); iter++) {
            _event_buffers.at(iter)->Init();
        }
    }

    int LoadEvent(Long64_t entry) {
        int read_bytes = 0;
        for (unsigned int iter = 0; iter < _event_buffers.size(); iter++) {
            read_bytes += _event_buffers.at(iter)->LoadEvent(entry);
        }
        return read_bytes;
    }

    EventBuffer const& operator[] (const std::string name) const {
        return GetEventByName(name);
    }

    void InstallEvent(std::string name, std::string loation);
};

EventManager& operator<< (EventManager &manager, UpdatesOnEvent* updater);
#endif
