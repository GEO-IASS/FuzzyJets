#include <sstream>

#include "Event.h"

#include <TH2.h>
#include <TStyle.h>

void EventBuffer::Init()
{
   mGMM_weight_vec = 0;
   mGMMs_weight_vec = 0;
   mGMMc_weight_vec = 0;
   mTGMM_weight_vec = 0;
   mTGMMs_weight_vec = 0;
   mUMM_weight_vec = 0;
   mGMM_distance_vec = 0;
   mGMMs_distance_vec = 0;
   mGMMc_distance_vec = 0;
   mTGMM_distance_vec = 0;
   mTGMMs_distance_vec = 0;
   mUMM_distance_vec = 0;
   is_pileup_vec = 0;
   is_lead_antikt_constituent_vec = 0;
   is_lead_antikt_trimmed_two_constituent_vec = 0;

   CA_nsubjettiness = 0;
   antikt_nsubjettiness = 0;

   mGMM_etas = 0;
   mGMM_phis = 0;
   mGMMc_etas = 0;
   mGMMc_phis = 0;
   mGMMs_etas = 0;
   mGMMs_phis = 0;
   mTGMM_etas = 0;
   mTGMM_phis = 0;
   mTGMMs_etas = 0;
   mTGMMs_phis = 0;
   mUMM_etas = 0;
   mUMM_phis = 0;

   _tree->SetBranchAddress("antikt_nsubjettiness", &antikt_nsubjettiness, &b_antikt_nsubjettiness);
   _tree->SetBranchAddress("CA_nsubjettiness", &CA_nsubjettiness, &b_CA_nsubjettiness);

   _tree->SetBranchAddress("CA_area", &CA_area, &b_CA_area);
   _tree->SetBranchAddress("antikt_area", &antikt_area, &b_antikt_area);
   _tree->SetBranchAddress("antikt_area_trimmed_two", &antikt_area_trimmed_two, &b_antikt_area_trimmed_two);
   _tree->SetBranchAddress("antikt_area_trimmed_three", &antikt_area_trimmed_three, &b_antikt_area_trimmed_three);

   _tree->SetBranchAddress("mGMM_etas", &mGMM_etas, &b_mGMM_etas);
   _tree->SetBranchAddress("mGMM_phis", &mGMM_phis, &b_mGMM_phis);
   _tree->SetBranchAddress("mGMMc_etas", &mGMMc_etas, &b_mGMMc_etas);
   _tree->SetBranchAddress("mGMMc_phis", &mGMMc_phis, &b_mGMMc_phis);
   _tree->SetBranchAddress("mGMMs_etas", &mGMMs_etas, &b_mGMMs_etas);
   _tree->SetBranchAddress("mGMMs_phis", &mGMMs_phis, &b_mGMMs_phis);
   _tree->SetBranchAddress("mTGMM_etas", &mTGMM_etas, &b_mTGMM_etas);
   _tree->SetBranchAddress("mTGMM_phis", &mTGMM_phis, &b_mTGMM_phis);
   _tree->SetBranchAddress("mTGMMs_etas", &mTGMMs_etas, &b_mTGMMs_etas);
   _tree->SetBranchAddress("mTGMMs_phis", &mTGMMs_phis, &b_mTGMMs_phis);
   _tree->SetBranchAddress("mUMM_etas", &mUMM_etas, &b_mUMM_etas);
   _tree->SetBranchAddress("mUMM_phis", &mUMM_phis, &b_mUMM_phis);

   _tree->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   _tree->SetBranchAddress("NPV", &NPV, &b_NPV);
   _tree->SetBranchAddress("mGMM_weight_vec", &mGMM_weight_vec, &b_mGMM_weight_vec);
   _tree->SetBranchAddress("mGMMs_weight_vec", &mGMMs_weight_vec, &b_mGMMs_weight_vec);
   _tree->SetBranchAddress("mGMMc_weight_vec", &mGMMc_weight_vec, &b_mGMMc_weight_vec);
   _tree->SetBranchAddress("mTGMM_weight_vec", &mTGMM_weight_vec, &b_mTGMM_weight_vec);
   _tree->SetBranchAddress("mTGMMs_weight_vec", &mTGMMs_weight_vec, &b_mTGMMs_weight_vec);
   _tree->SetBranchAddress("mUMM_weight_vec", &mUMM_weight_vec, &b_mUMM_weight_vec);
   _tree->SetBranchAddress("mGMM_distance_vec", &mGMM_distance_vec, &b_mGMM_distance_vec);
   _tree->SetBranchAddress("mGMMs_distance_vec", &mGMMs_distance_vec, &b_mGMMs_distance_vec);
   _tree->SetBranchAddress("mGMMc_distance_vec", &mGMMc_distance_vec, &b_mGMMc_distance_vec);
   _tree->SetBranchAddress("mTGMM_distance_vec", &mTGMM_distance_vec, &b_mTGMM_distance_vec);
   _tree->SetBranchAddress("mTGMMs_distance_vec", &mTGMMs_distance_vec, &b_mTGMMs_distance_vec);
   _tree->SetBranchAddress("mUMM_distance_vec", &mUMM_distance_vec, &b_mUMM_distance_vec);
   _tree->SetBranchAddress("is_pileup_vec", &is_pileup_vec, &b_is_pileup_vec);
   _tree->SetBranchAddress("is_lead_antikt_constituent_vec", &is_lead_antikt_constituent_vec, &b_is_lead_antikt_constituent_vec);
   _tree->SetBranchAddress("is_lead_antikt_trimmed_two_constituent_vec", &is_lead_antikt_trimmed_two_constituent_vec, &b_is_lead_antikt_trimmed_two_constituent_vec);
   _tree->SetBranchAddress("CA_dr", &CA_dr, &b_CA_dr);
   _tree->SetBranchAddress("CA_m", &CA_m, &b_CA_m);
   _tree->SetBranchAddress("CA_m_pu", &CA_m_pu, &b_CA_m_pu);
   _tree->SetBranchAddress("CA_pt", &CA_pt, &b_CA_pt);
   _tree->SetBranchAddress("CA_pufrac", &CA_pufrac, &b_CA_pufrac);
   _tree->SetBranchAddress("antikt_dr", &antikt_dr, &b_antikt_dr);
   _tree->SetBranchAddress("antikt_m", &antikt_m, &b_antikt_m);
   _tree->SetBranchAddress("antikt_m_pu", &antikt_m_pu, &b_antikt_m_pu);
   _tree->SetBranchAddress("antikt_m_pu_trimmed_three", &antikt_m_pu_trimmed_three, &b_antikt_m_pu_trimmed_three);
   _tree->SetBranchAddress("antikt_m_pu_trimmed_two", &antikt_m_pu_trimmed_two, &b_antikt_m_pu_trimmed_two);
   _tree->SetBranchAddress("antikt_m_trimmed_three", &antikt_m_trimmed_three, &b_antikt_m_trimmed_three);
   _tree->SetBranchAddress("antikt_m_trimmed_two", &antikt_m_trimmed_two, &b_antikt_m_trimmed_two);
   _tree->SetBranchAddress("antikt_pt", &antikt_pt, &b_antikt_pt);
   _tree->SetBranchAddress("antikt_pt_trimmed_three", &antikt_pt_trimmed_three, &b_antikt_pt_trimmed_three);
   _tree->SetBranchAddress("antikt_pt_trimmed_two", &antikt_pt_trimmed_two, &b_antikt_pt_trimmed_two);
   _tree->SetBranchAddress("antikt_pufrac", &antikt_pufrac, &b_antikt_pufrac);
   _tree->SetBranchAddress("antikt_pufrac_trimmed_three", &antikt_pufrac_trimmed_three, &b_antikt_pufrac_trimmed_three);
   _tree->SetBranchAddress("antikt_pufrac_trimmed_two", &antikt_pufrac_trimmed_two, &b_antikt_pufrac_trimmed_two);
   _tree->SetBranchAddress("deltatop_mGMM", &deltatop_mGMM, &b_deltatop_mGMM);
   _tree->SetBranchAddress("deltatop_mGMMc", &deltatop_mGMMc, &b_deltatop_mGMMc);
   _tree->SetBranchAddress("deltatop_mGMMs", &deltatop_mGMMs, &b_deltatop_mGMMs);
   _tree->SetBranchAddress("deltatop_mTGMM", &deltatop_mTGMM, &b_deltatop_mTGMM);
   _tree->SetBranchAddress("deltatop_mTGMMs", &deltatop_mTGMMs, &b_deltatop_mTGMMs);
   _tree->SetBranchAddress("deltatop_mUMM", &deltatop_mUMM, &b_deltatop_mUMM);
   _tree->SetBranchAddress("mGMM_clpu", &mGMM_clpu, &b_mGMM_clpu);
   _tree->SetBranchAddress("mGMM_dr", &mGMM_dr, &b_mGMM_dr);
   _tree->SetBranchAddress("mGMM_m", &mGMM_m, &b_mGMM_m);
   _tree->SetBranchAddress("mGMM_m_mean", &mGMM_m_mean, &b_mGMM_m_mean);
   _tree->SetBranchAddress("mGMM_m_pu_hard", &mGMM_m_pu_hard, &b_mGMM_m_pu_hard);
   _tree->SetBranchAddress("mGMM_m_pu_soft", &mGMM_m_pu_soft, &b_mGMM_m_pu_soft);
   _tree->SetBranchAddress("mGMM_m_skew", &mGMM_m_skew, &b_mGMM_m_skew);
   _tree->SetBranchAddress("mGMM_m_soft", &mGMM_m_soft, &b_mGMM_m_soft);
   _tree->SetBranchAddress("mGMM_m_var", &mGMM_m_var, &b_mGMM_m_var);
   _tree->SetBranchAddress("mGMM_ml", &mGMM_ml, &b_mGMM_ml);
   _tree->SetBranchAddress("mGMM_pt", &mGMM_pt, &b_mGMM_pt);
   _tree->SetBranchAddress("mGMM_pt_mean", &mGMM_pt_mean, &b_mGMM_pt_mean);
   _tree->SetBranchAddress("mGMM_pt_skew", &mGMM_pt_skew, &b_mGMM_pt_skew);
   _tree->SetBranchAddress("mGMM_pt_soft", &mGMM_pt_soft, &b_mGMM_pt_soft);
   _tree->SetBranchAddress("mGMM_pt_var", &mGMM_pt_var, &b_mGMM_pt_var);
   _tree->SetBranchAddress("mGMM_ptl", &mGMM_ptl, &b_mGMM_ptl);
   _tree->SetBranchAddress("mGMM_pufrac_hard", &mGMM_pufrac_hard, &b_mGMM_pufrac_hard);
   _tree->SetBranchAddress("mGMM_pufrac_soft", &mGMM_pufrac_soft, &b_mGMM_pufrac_soft);
   _tree->SetBranchAddress("mGMM_ucpu", &mGMM_ucpu, &b_mGMM_ucpu);
   _tree->SetBranchAddress("mGMMc_clpu", &mGMMc_clpu, &b_mGMMc_clpu);
   _tree->SetBranchAddress("mGMMc_dr", &mGMMc_dr, &b_mGMMc_dr);
   _tree->SetBranchAddress("mGMMc_m", &mGMMc_m, &b_mGMMc_m);
   _tree->SetBranchAddress("mGMMc_m_mean", &mGMMc_m_mean, &b_mGMMc_m_mean);
   _tree->SetBranchAddress("mGMMc_m_pu_hard", &mGMMc_m_pu_hard, &b_mGMMc_m_pu_hard);
   _tree->SetBranchAddress("mGMMc_m_pu_soft", &mGMMc_m_pu_soft, &b_mGMMc_m_pu_soft);
   _tree->SetBranchAddress("mGMMc_m_skew", &mGMMc_m_skew, &b_mGMMc_m_skew);
   _tree->SetBranchAddress("mGMMc_m_soft", &mGMMc_m_soft, &b_mGMMc_m_soft);
   _tree->SetBranchAddress("mGMMc_m_var", &mGMMc_m_var, &b_mGMMc_m_var);
   _tree->SetBranchAddress("mGMMc_ml", &mGMMc_ml, &b_mGMMc_ml);
   _tree->SetBranchAddress("mGMMc_pt", &mGMMc_pt, &b_mGMMc_pt);
   _tree->SetBranchAddress("mGMMc_pt_mean", &mGMMc_pt_mean, &b_mGMMc_pt_mean);
   _tree->SetBranchAddress("mGMMc_pt_skew", &mGMMc_pt_skew, &b_mGMMc_pt_skew);
   _tree->SetBranchAddress("mGMMc_pt_soft", &mGMMc_pt_soft, &b_mGMMc_pt_soft);
   _tree->SetBranchAddress("mGMMc_pt_var", &mGMMc_pt_var, &b_mGMMc_pt_var);
   _tree->SetBranchAddress("mGMMc_ptl", &mGMMc_ptl, &b_mGMMc_ptl);
   _tree->SetBranchAddress("mGMMc_pufrac_hard", &mGMMc_pufrac_hard, &b_mGMMc_pufrac_hard);
   _tree->SetBranchAddress("mGMMc_pufrac_soft", &mGMMc_pufrac_soft, &b_mGMMc_pufrac_soft);
   _tree->SetBranchAddress("mGMMc_r", &mGMMc_r, &b_mGMMc_r);
   _tree->SetBranchAddress("mGMMc_r_avg", &mGMMc_r_avg, &b_mGMMc_r_avg);
   _tree->SetBranchAddress("mGMMc_r_second", &mGMMc_r_second, &b_mGMMc_r_second);
   _tree->SetBranchAddress("mGMMc_r_third", &mGMMc_r_third, &b_mGMMc_r_third);
   _tree->SetBranchAddress("mGMMc_r_weighted_avg", &mGMMc_r_weighted_avg, &b_mGMMc_r_weighted_avg);
   _tree->SetBranchAddress("mGMMc_ucpu", &mGMMc_ucpu, &b_mGMMc_ucpu);
   _tree->SetBranchAddress("mGMMs_clpu", &mGMMs_clpu, &b_mGMMs_clpu);
   _tree->SetBranchAddress("mGMMs_dr", &mGMMs_dr, &b_mGMMs_dr);
   _tree->SetBranchAddress("mGMMs_m", &mGMMs_m, &b_mGMMs_m);
   _tree->SetBranchAddress("mGMMs_m_mean", &mGMMs_m_mean, &b_mGMMs_m_mean);
   _tree->SetBranchAddress("mGMMs_m_pu_hard", &mGMMs_m_pu_hard, &b_mGMMs_m_pu_hard);
   _tree->SetBranchAddress("mGMMs_m_pu_soft", &mGMMs_m_pu_soft, &b_mGMMs_m_pu_soft);
   _tree->SetBranchAddress("mGMMs_m_skew", &mGMMs_m_skew, &b_mGMMs_m_skew);
   _tree->SetBranchAddress("mGMMs_m_soft", &mGMMs_m_soft, &b_mGMMs_m_soft);
   _tree->SetBranchAddress("mGMMs_m_var", &mGMMs_m_var, &b_mGMMs_m_var);
   _tree->SetBranchAddress("mGMMs_ml", &mGMMs_ml, &b_mGMMs_ml);
   _tree->SetBranchAddress("mGMMs_pt", &mGMMs_pt, &b_mGMMs_pt);
   _tree->SetBranchAddress("mGMMs_pt_mean", &mGMMs_pt_mean, &b_mGMMs_pt_mean);
   _tree->SetBranchAddress("mGMMs_pt_skew", &mGMMs_pt_skew, &b_mGMMs_pt_skew);
   _tree->SetBranchAddress("mGMMs_pt_soft", &mGMMs_pt_soft, &b_mGMMs_pt_soft);
   _tree->SetBranchAddress("mGMMs_pt_var", &mGMMs_pt_var, &b_mGMMs_pt_var);
   _tree->SetBranchAddress("mGMMs_ptl", &mGMMs_ptl, &b_mGMMs_ptl);
   _tree->SetBranchAddress("mGMMs_pufrac_hard", &mGMMs_pufrac_hard, &b_mGMMs_pufrac_hard);
   _tree->SetBranchAddress("mGMMs_pufrac_soft", &mGMMs_pufrac_soft, &b_mGMMs_pufrac_soft);
   _tree->SetBranchAddress("mGMMs_ucpu", &mGMMs_ucpu, &b_mGMMs_ucpu);
   _tree->SetBranchAddress("mTGMM_clpu", &mTGMM_clpu, &b_mTGMM_clpu);
   _tree->SetBranchAddress("mTGMM_dr", &mTGMM_dr, &b_mTGMM_dr);
   _tree->SetBranchAddress("mTGMM_m", &mTGMM_m, &b_mTGMM_m);
   _tree->SetBranchAddress("mTGMM_m_mean", &mTGMM_m_mean, &b_mTGMM_m_mean);
   _tree->SetBranchAddress("mTGMM_m_pu_hard", &mTGMM_m_pu_hard, &b_mTGMM_m_pu_hard);
   _tree->SetBranchAddress("mTGMM_m_pu_soft", &mTGMM_m_pu_soft, &b_mTGMM_m_pu_soft);
   _tree->SetBranchAddress("mTGMM_m_skew", &mTGMM_m_skew, &b_mTGMM_m_skew);
   _tree->SetBranchAddress("mTGMM_m_soft", &mTGMM_m_soft, &b_mTGMM_m_soft);
   _tree->SetBranchAddress("mTGMM_m_var", &mTGMM_m_var, &b_mTGMM_m_var);
   _tree->SetBranchAddress("mTGMM_ml", &mTGMM_ml, &b_mTGMM_ml);
   _tree->SetBranchAddress("mTGMM_pt", &mTGMM_pt, &b_mTGMM_pt);
   _tree->SetBranchAddress("mTGMM_pt_mean", &mTGMM_pt_mean, &b_mTGMM_pt_mean);
   _tree->SetBranchAddress("mTGMM_pt_skew", &mTGMM_pt_skew, &b_mTGMM_pt_skew);
   _tree->SetBranchAddress("mTGMM_pt_soft", &mTGMM_pt_soft, &b_mTGMM_pt_soft);
   _tree->SetBranchAddress("mTGMM_pt_var", &mTGMM_pt_var, &b_mTGMM_pt_var);
   _tree->SetBranchAddress("mTGMM_ptl", &mTGMM_ptl, &b_mTGMM_ptl);
   _tree->SetBranchAddress("mTGMM_pufrac_hard", &mTGMM_pufrac_hard, &b_mTGMM_pufrac_hard);
   _tree->SetBranchAddress("mTGMM_pufrac_soft", &mTGMM_pufrac_soft, &b_mTGMM_pufrac_soft);
   _tree->SetBranchAddress("mTGMM_ucpu", &mTGMM_ucpu, &b_mTGMM_ucpu);
   _tree->SetBranchAddress("mTGMMs_clpu", &mTGMMs_clpu, &b_mTGMMs_clpu);
   _tree->SetBranchAddress("mTGMMs_dr", &mTGMMs_dr, &b_mTGMMs_dr);
   _tree->SetBranchAddress("mTGMMs_m", &mTGMMs_m, &b_mTGMMs_m);
   _tree->SetBranchAddress("mTGMMs_m_mean", &mTGMMs_m_mean, &b_mTGMMs_m_mean);
   _tree->SetBranchAddress("mTGMMs_m_pu_hard", &mTGMMs_m_pu_hard, &b_mTGMMs_m_pu_hard);
   _tree->SetBranchAddress("mTGMMs_m_pu_soft", &mTGMMs_m_pu_soft, &b_mTGMMs_m_pu_soft);
   _tree->SetBranchAddress("mTGMMs_m_skew", &mTGMMs_m_skew, &b_mTGMMs_m_skew);
   _tree->SetBranchAddress("mTGMMs_m_soft", &mTGMMs_m_soft, &b_mTGMMs_m_soft);
   _tree->SetBranchAddress("mTGMMs_m_var", &mTGMMs_m_var, &b_mTGMMs_m_var);
   _tree->SetBranchAddress("mTGMMs_ml", &mTGMMs_ml, &b_mTGMMs_ml);
   _tree->SetBranchAddress("mTGMMs_pt", &mTGMMs_pt, &b_mTGMMs_pt);
   _tree->SetBranchAddress("mTGMMs_pt_mean", &mTGMMs_pt_mean, &b_mTGMMs_pt_mean);
   _tree->SetBranchAddress("mTGMMs_pt_skew", &mTGMMs_pt_skew, &b_mTGMMs_pt_skew);
   _tree->SetBranchAddress("mTGMMs_pt_soft", &mTGMMs_pt_soft, &b_mTGMMs_pt_soft);
   _tree->SetBranchAddress("mTGMMs_pt_var", &mTGMMs_pt_var, &b_mTGMMs_pt_var);
   _tree->SetBranchAddress("mTGMMs_ptl", &mTGMMs_ptl, &b_mTGMMs_ptl);
   _tree->SetBranchAddress("mTGMMs_pufrac_hard", &mTGMMs_pufrac_hard, &b_mTGMMs_pufrac_hard);
   _tree->SetBranchAddress("mTGMMs_pufrac_soft", &mTGMMs_pufrac_soft, &b_mTGMMs_pufrac_soft);
   _tree->SetBranchAddress("mTGMMs_ucpu", &mTGMMs_ucpu, &b_mTGMMs_ucpu);
   _tree->SetBranchAddress("mUMM_clpu", &mUMM_clpu, &b_mUMM_clpu);
   _tree->SetBranchAddress("mUMM_dr", &mUMM_dr, &b_mUMM_dr);
   _tree->SetBranchAddress("mUMM_m", &mUMM_m, &b_mUMM_m);
   _tree->SetBranchAddress("mUMM_m_mean", &mUMM_m_mean, &b_mUMM_m_mean);
   _tree->SetBranchAddress("mUMM_m_pu_hard", &mUMM_m_pu_hard, &b_mUMM_m_pu_hard);
   _tree->SetBranchAddress("mUMM_m_pu_soft", &mUMM_m_pu_soft, &b_mUMM_m_pu_soft);
   _tree->SetBranchAddress("mUMM_m_skew", &mUMM_m_skew, &b_mUMM_m_skew);
   _tree->SetBranchAddress("mUMM_m_soft", &mUMM_m_soft, &b_mUMM_m_soft);
   _tree->SetBranchAddress("mUMM_m_var", &mUMM_m_var, &b_mUMM_m_var);
   _tree->SetBranchAddress("mUMM_ml", &mUMM_ml, &b_mUMM_ml);
   _tree->SetBranchAddress("mUMM_pt", &mUMM_pt, &b_mUMM_pt);
   _tree->SetBranchAddress("mUMM_pt_mean", &mUMM_pt_mean, &b_mUMM_pt_mean);
   _tree->SetBranchAddress("mUMM_pt_skew", &mUMM_pt_skew, &b_mUMM_pt_skew);
   _tree->SetBranchAddress("mUMM_pt_soft", &mUMM_pt_soft, &b_mUMM_pt_soft);
   _tree->SetBranchAddress("mUMM_pt_var", &mUMM_pt_var, &b_mUMM_pt_var);
   _tree->SetBranchAddress("mUMM_ptl", &mUMM_ptl, &b_mUMM_ptl);
   _tree->SetBranchAddress("mUMM_pufrac_hard", &mUMM_pufrac_hard, &b_mUMM_pufrac_hard);
   _tree->SetBranchAddress("mUMM_pufrac_soft", &mUMM_pufrac_soft, &b_mUMM_pufrac_soft);
   _tree->SetBranchAddress("mUMM_ucpu", &mUMM_ucpu, &b_mUMM_ucpu);
   _tree->SetBranchAddress("toppt", &toppt, &b_toppt);
}

void EventManager::InstallEvent(std::string name, std::string location) {
    std::stringstream ss;
    ss.str(std::string());
    ss << "_pT_spectrum_" << name;
    if (_pT_spectrum_map.find(ss.str()) != _pT_spectrum_map.end())
        return;
    TH1F *spectrum = new TH1F(ss.str().c_str(), "", _n_pT_bins, _pT_low, _pT_high);
    _pT_spectra.push_back(spectrum);
    _pT_spectrum_map[name] = spectrum;
    TFile *f = TFile::Open(location.c_str(), "READ");
    TTree *t = (TTree *) f->Get("EventTree");

    EventBuffer *buffer = new EventBuffer();
    buffer->SetTree(t);
    _event_buffers.push_back(buffer);
    _event_buffer_map[name] = buffer;
    _event_buffer_labels.push_back(name);
    _event_buffer_locations.push_back(location);
}

EventManager& operator<< (EventManager &manager, UpdatesOnEvent* updater) {
    manager.InsertUpdater(updater);
    return manager;
}
