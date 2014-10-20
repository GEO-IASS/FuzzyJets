#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <tuple>
#include <utility>
#include <unordered_map>

#include <TROOT.h>
#include <TColor.h>
#include <TLine.h>

#include "AnalyzeFuzzyTools.h"
#include "AtlasStyle.h"
#include "Event.h"

#include "boost/program_options.hpp"

namespace po = boost::program_options;

typedef std::tuple<int, int, int, int> map_key_t;

int mass_low = 0;
int mass_high = 200;

struct key_hash : public std::unary_function<map_key_t, std::size_t> {
    std::size_t operator()(const map_key_t& key) const {
        return std::get<0>(key) ^ std::get<1>(key) ^ std::get<2>(key) ^ std::get<3>(key);
    }
};

struct key_equal : public std::binary_function<map_key_t, map_key_t, bool> {
    bool operator()(const map_key_t& key1, const map_key_t& key2) const {
        return (std::get<0>(key1) == std::get<0>(key2) &&
                std::get<1>(key1) == std::get<1>(key2) &&
                std::get<2>(key1) == std::get<2>(key2) &&
                std::get<3>(key1) == std::get<3>(key2));
    }
};

typedef std::unordered_map<map_key_t, std::string, key_hash, key_equal> file_map_t;

TString buildFileLocation(int size, int learn, int n_pileup_vertices, int do_rec, std::string prefix) {
    return TString::Format("%s%ds_%dmu_%dlw_%drec.root", prefix.c_str(), size, n_pileup_vertices, learn, do_rec);
}

file_map_t
constructFileMap(std::vector<int> const& sizes, std::vector<int> const& learns,
                 std::vector<int> const& do_recs, std::vector<int> const& pileup_vertices,
                 std::string const& prefix) {
  file_map_t m;
  for (unsigned int size_iter = 0; size_iter < sizes.size(); size_iter++) {
      for (unsigned int learn_iter = 0; learn_iter < learns.size(); learn_iter++) {
          for (unsigned int pileup_iter = 0; pileup_iter < pileup_vertices.size(); pileup_iter++) {
              for (unsigned int rec_iter = 0; rec_iter < do_recs.size(); rec_iter++) {
                  int current_do_rec = do_recs[rec_iter];
                  int current_size = sizes[size_iter];
                  int current_learn = learns[learn_iter];
                  int current_n_pileup = pileup_vertices[pileup_iter];
                  TString file_location = buildFileLocation(current_size, current_learn,
                                                            current_n_pileup, current_do_rec,
                                                            prefix);
                  map_key_t new_key = std::make_tuple(current_size, current_learn, current_n_pileup, current_do_rec);
                  m[new_key] = file_location;
              }
          }
      }
  }
  return m;
}

void doSanityPlot(std::string branch, std::string old_loc, std::string new_loc, std::string title,
                  std::string out) {
    std::vector<HistHelper> v_hist_decs;
    CanvasHelper c_dec("X", "", title, out, 800, 800);
    c_dec.diff_scale = false;
    HistHelper h_one(old_loc, branch, "OLD", 510, mass_low, mass_high, 25, StyleTypes::STRIPED, kMagenta + 3, 20, "");
    HistHelper h_two(new_loc, branch, "NEW", 510, mass_low, mass_high, 25, StyleTypes::STRIPED, kViolet + 9, 20, "");
    v_hist_decs.push_back(h_one);
    v_hist_decs.push_back(h_two);
    prettyHist<float>(v_hist_decs, c_dec);
}

void sanityTests() {
    std::string sanity_out = "/u/at/chstan/nfs/summer_2014/ForConrad/results/plots/sanity/";
    std::string old_data_pref = "/u/at/chstan/nfs/summer_2014/ForConrad/results/with_pu/";
    std::string new_data_pref = "/u/at/chstan/nfs/summer_2014/ForConrad/files/2014_08_04_19h02m14s/";
    static const std::string sets[] = {"wno_s10/", "wno_s8/", "wyes_s10/", "wyes_s8/"};
    static const int sizes[] = {10, 8};
    static const int learns[] = {0, 1};
    static const int NPVs[] = {0, 10, 20, 30};
    static const std::string sporadic[] =
        { "CA_m", "CA_pt", "toppt", "antikt_m", "antikt_pt" };
    static const std::string algs[] = { "mGMM", "mGMMs", "mUMM", "mTGMM", "mTGMMs" };
    static const std::string vars[] = { "_m", "_pt", "_m_soft", "_pt_soft" };
    for (unsigned int size_iter = 0; size_iter < 2; size_iter++) {
        for (unsigned int learn_iter = 0; learn_iter < 2; learn_iter++) {
            unsigned int set_idx = size_iter + 2*learn_iter;
            for (unsigned int npv_iter = 0; npv_iter < 4; npv_iter++) {
                int npv = NPVs[npv_iter];
                int size = sizes[size_iter];
                int learn = learns[learn_iter];
                std::stringstream ss;
                ss.str("");
                ss << old_data_pref << sets[set_idx] << npv << ".root";
                std::string old_loc = ss.str();
                ss.str("");
                ss << new_data_pref << size << "s_" << npv << "mu_" << learn << ".root";
                std::string new_loc = ss.str();

                for (unsigned int sporadic_iter = 0; sporadic_iter < 5; sporadic_iter++) {
                    ss.str("");
                    ss << "SANITY CHECK MU " << npv << " SZ "
                       << size << " LEARN " << learn << " BRANCH "
                       << sporadic[sporadic_iter];
                    std::string title = ss.str();
                    doSanityPlot(sporadic[sporadic_iter], old_loc, new_loc, title, sanity_out);
                }
                for (unsigned int alg_iter = 0; alg_iter < 5; alg_iter++) {
                    for (unsigned int branch_iter = 0; branch_iter < 4; branch_iter++) {
                        ss.str("");
                        ss << algs[alg_iter] << vars[branch_iter];
                        std::string branch_name = ss.str();
                        ss.str("");
                        ss << "SANITY CHECK MU " << npv << " SZ "
                           << size << " LEARN " << learn << " BRANCH "
                           << branch_name;
                        std::string title = ss.str();
                        doSanityPlot(branch_name, old_loc, new_loc, title, sanity_out);
                    }
                }
            }
        }
    }
}

void doFancyWeightHisto(std::string alg, int do_rec, int npv, int learn, int size,
                        std::string out_dir, file_map_t file_m, int num) {
    std::string file_name = file_m[std::make_tuple(size, learn, npv, do_rec)];

    std::vector<float> antikt_trimmed_pufracs =
        loadSingleBranch<float>(file_name, "antikt_pufrac_trimmed_three");
    float antikt_pufrac = 0;
    for (unsigned int antikt_iter = 0; antikt_iter < antikt_trimmed_pufracs.size(); antikt_iter++) {
        antikt_pufrac += antikt_trimmed_pufracs[antikt_iter];
    }
    antikt_pufrac /= antikt_trimmed_pufracs.size();

    std::stringstream ss;
    ss << "Hard Scatter Weights " << alg;
    std::string hs_hist_name = ss.str();
    ss.str(std::string());
    ss << "Pileup Weights " << alg;
    std::string pu_hist_name = ss.str();

    TFile *f = new TFile(file_name.c_str());
    TH1F *pu_hist = (TH1F *)f->Get(pu_hist_name.c_str());
    TH1F *hs_hist = (TH1F *)f->Get(hs_hist_name.c_str());
    TH1F *pu_hist_int = new TH1F("pu_int", "", 48, 0.04, 1);
    TH1F *hs_hist_int = new TH1F("hs_int", "", 48, 0.04, 1);
    TH1F *total_hist_int = new TH1F("total_int", "", 48, 0.04, 1);
    TH1F *pu_frac_hist = new TH1F("pu_frac", "", 48, 0.04, 1);
    Float_t pu_sum = 0;
    Float_t hs_sum = 0;
    for (Int_t iter = 48; iter >= 1; iter--) {
        pu_sum += pu_hist->GetBinContent(iter);
        hs_sum += hs_hist->GetBinContent(iter);
        pu_hist_int->SetBinContent(iter, pu_sum);
        hs_hist_int->SetBinContent(iter, hs_sum);
    }
    total_hist_int->Add(pu_hist_int, hs_hist_int, 1., 1.);
    pu_frac_hist->Divide(pu_hist_int, total_hist_int, 1., 1., "b");

    if(!(hs_hist->Integral(-1, hs_hist->GetNbinsX()+1) > 0) ||
       !(pu_hist->Integral(-1, pu_hist->GetNbinsX()+1) > 0)) {
        std::cout << "PU OR HS HIST HAS NO CONTENT." << std::endl;
        return;
    }

    pu_hist->Scale(1./pu_hist->Integral(-1, pu_hist->GetNbinsX()+1));
    hs_hist->Scale(1./hs_hist->Integral(-1, hs_hist->GetNbinsX()+1));

    TCanvas *c = new TCanvas("", "Particle Composition by Weight Comp. Antikt", 800, 800);
    TLine *antikt_pufrac_l = new TLine(c->GetUxmin(), antikt_pufrac, c->GetUxmax(), antikt_pufrac);
    antikt_pufrac_l->SetLineColor(kBlack);
    antikt_pufrac_l->SetLineStyle(2);

    pu_frac_hist->SetLineColor(kBlue);
    pu_frac_hist->SetFillStyle(3004);
    pu_frac_hist->SetFillColor(kBlue);
    pu_frac_hist->GetYaxis()->SetRangeUser(0, 1);

    pu_frac_hist->Draw("");
    antikt_pufrac_l->Draw("same");

    TLegend *leggaa = new TLegend(0.6, 0.7, 0.9, 0.8);
    leggaa->SetTextFont(42);
    leggaa->AddEntry(pu_frac_hist, "Pu part. frac w/ w > bin", "l");
    leggaa->AddEntry(antikt_pufrac_l, "Antikt pu part. frac", "l");

    leggaa->SetFillStyle(0);
    leggaa->SetFillColor(0);
    leggaa->SetBorderSize(0);
    leggaa->Draw();
    //DrawAtlasLabel("Particle Comp");

    ss.str(std::string());
    ss << out_dir << num << "_Pufrac_antikt_" << alg << "_sz_" << size << "_lw_" << learn
       << "_mu_" << npv << "_rec_" << do_rec << ".pdf";
    std::string out_file_name = ss.str();
    c->Print(out_file_name.c_str());

    delete leggaa;
    delete c;

    // also do the cumulative version of the old plot
    c = new TCanvas("", "PU and HS Integral by Weight", 800, 800);
    pu_hist->SetLineColor(kBlue);
    pu_hist->SetFillStyle(3004);
    pu_hist->SetFillColor(kBlue);
    pu_hist->SetMinimum(0);

    hs_hist->SetLineColor(kRed);
    hs_hist->SetFillStyle(3004);
    hs_hist->SetFillColor(kRed);
    hs_hist->SetMinimum(0);

    hs_hist->Draw("");
    pu_hist->Draw("same");

    leggaa = new TLegend(0.6, 0.7, 0.9, 0.8);
    leggaa->SetTextFont(42);
    leggaa->AddEntry(pu_hist, "PU Integral", "l");
    leggaa->AddEntry(hs_hist, "HS Integral", "l");

    leggaa->SetFillStyle(0);
    leggaa->SetFillColor(0);
    leggaa->SetBorderSize(0);
    leggaa->Draw();
    //DrawAtlasLabel(c_dec.title);

    ss.str(std::string());
    ss << out_dir << (num+1) << "_Cumul_norm_pufrac_" << alg << "_sz_" << size << "_lw_" << learn
       << "_mu_" << npv << "_rec_" << do_rec << ".pdf";
    out_file_name = ss.str();
    c->Print(out_file_name.c_str());

    delete leggaa;
    delete c;

    delete pu_hist;
    delete hs_hist;
    delete pu_hist_int;
    delete hs_hist_int;
    delete total_hist_int;
    delete pu_frac_hist;
    delete antikt_pufrac_l;
}

int fancyWeightHistos() {
    SetAtlasStyle();
    SetupATLASStyle();

    static const int sizes_arr[] = {7, 8, 9, 10};
    static const int NPVs_arr[] = {5};
    static const int learns_arr[] = {0, 1};
    static const int do_recs_arr[] = {0};
    static const std::string algs_arr[] = {"mGMM", "mGMMs", "mGMMc", "mUMM", "mTGMM", "mTGMMs"};

    //static const Int_t colors_arr[] = {kMagenta + 3, kViolet + 9, kTeal - 5, kGray + 3};
    //static const Style_t styles_arr[] = {kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown};

    std::vector<int> sizes(sizes_arr, sizes_arr+sizeof(sizes_arr) / sizeof(sizes_arr[0]));
    std::vector<int> NPVs(NPVs_arr, NPVs_arr+sizeof(NPVs_arr) / sizeof(NPVs_arr[0]));
    std::vector<int> learns(learns_arr, learns_arr+sizeof(learns_arr) / sizeof(learns_arr[0]));
    std::vector<int> do_recs(do_recs_arr, do_recs_arr+sizeof(do_recs_arr) / sizeof(do_recs_arr[0]));
    std::vector<std::string> algs(algs_arr, algs_arr+sizeof(algs_arr) / sizeof(algs_arr[0]));
    int num = 0;

    std::string file_prefix = "/u/at/chstan/nfs/summer_2014/ForConrad/files/20kevts_wprime_mu0_and_mu5_norec/2014_08_13_14h24m44s/";

    file_map_t file_m = constructFileMap(sizes, learns, do_recs, NPVs, file_prefix);

    std::string out_dir = "/u/at/chstan/nfs/summer_2014/ForConrad/results/plots/20kevts_wprime_mu0_and_mu5_norec/";

    for (unsigned int alg_iter = 0; alg_iter < algs.size(); alg_iter++) {
        for (unsigned int do_rec_iter = 0; do_rec_iter < do_recs.size(); do_rec_iter++) {
            for (unsigned int npv_iter = 0; npv_iter < NPVs.size(); npv_iter++) {
                for (unsigned int learn_iter = 0; learn_iter < learns.size(); learn_iter++) {
                    for (unsigned int size_iter = 0; size_iter < sizes.size(); size_iter++) {
                        std::string alg = algs[alg_iter];
                        int do_rec = do_recs[do_rec_iter];
                        int npv = NPVs[npv_iter];
                        int learn = learns[learn_iter];
                        int size = sizes[size_iter];
                        doFancyWeightHisto(alg, do_rec, npv, learn, size, out_dir, file_m, num);
                        num += 2;
                    }
                }
            }
        }
    }
    return num;
}

void pTOverlay() {
    SetAtlasStyle();
    SetupATLASStyle();

    static const std::string datasets_arr[] = {
        "/u/at/chstan/nfs/summer_2014/ForConrad/files/20kevts_wprime_mu0_and_mu5_norec/2014_08_13_14h24m44s/",
        "/u/at/chstan/nfs/summer_2014/ForConrad/files/20kevts_qcd_mu0_and_mu5_norec/2014_08_13_14h28m03s/",
        "/u/at/chstan/nfs/summer_2014/ForConrad/files/500kevts_mu0_bothrec/2014_08_11_19h46m58s/"};

    static const std::string dataset_names_arr[] = {"W' m 600 GeV", "QCD: pT 300-700 GeV", "Z' m 1500 GeV"};

    std::vector<std::string> datasets(datasets_arr, datasets_arr + sizeof(datasets_arr) / sizeof(datasets_arr[0]));

    std::string out_dir = "/u/at/chstan/nfs/summer_2014/ForConrad/results/plots/ptcomp/";

    static const Int_t colors_arr[] = {kMagenta + 3, kViolet + 9, kTeal - 5, kGray + 3};

    CanvasHelper c_dec_trimmed("pT [GeV]", "", "pT Comparison, W', Z', QCD", out_dir, 800, 800);
    c_dec_trimmed.diff_scale = true;

    std::vector<HistHelper> v_hist_decs;
    for (unsigned int dataset_iter = 0; dataset_iter < datasets.size(); dataset_iter++) {
        std::string dataset = datasets[dataset_iter];
        std::string dataset_name = dataset_names_arr[dataset_iter];

        std::stringstream ss;
        ss.str(std::string());
        ss << dataset_name << ": anti-kt trimmed pT";
        std::string title = ss.str();
        ss.str(std::string());
        ss << dataset << "10s_0mu_0lw_0rec.root";
        std::string file_location = ss.str();
        HistHelper hist_helper_temp(file_location,
                                    "antikt_pt_trimmed_three", title, 510, 0, 800, 100,
                                    StyleTypes::NONE, colors_arr[dataset_iter], kFullCircle, "");
        v_hist_decs.push_back(hist_helper_temp);
    }
    prettyHist<float>(v_hist_decs, c_dec_trimmed);

    v_hist_decs.clear();

    CanvasHelper c_dec_untrimmed("pT [GeV]", "", "Untrimmed pT Comparison, W', Z', QCD", out_dir, 800, 800);
    c_dec_untrimmed.diff_scale = true;

    for (unsigned int dataset_iter = 0; dataset_iter < datasets.size(); dataset_iter++) {
        std::string dataset = datasets[dataset_iter];
        std::string dataset_name = dataset_names_arr[dataset_iter];

        std::stringstream ss;
        ss.str(std::string());
        ss << dataset_name << ": anti-kt untrimmed pT";
        std::string title = ss.str();
        ss.str(std::string());
        ss << dataset << "10s_0mu_0lw_0rec.root";
        std::string file_location = ss.str();
        HistHelper hist_helper_temp(file_location,
                                    "antikt_pt", title, 510, 0, 800, 100,
                                    StyleTypes::NONE, colors_arr[dataset_iter], kFullCircle, "");
        v_hist_decs.push_back(hist_helper_temp);
    }
    prettyHist<float>(v_hist_decs, c_dec_untrimmed);
}

void EventTest() {
    std::string qcd_location_5 = "/u/at/chstan/nfs/summer_2014/ForConrad/files/150kevts_qcd_area/2014_10_14_11h23m44s/10s_0mu_0lw_0rec_5cut.root";
    std::string qcd_location_15 = "/u/at/chstan/nfs/summer_2014/ForConrad/files/150kevts_qcd_area/2014_10_14_11h23m44s/10s_0mu_0lw_0rec_15cut.root";
    std::string qcd_location_25 = "/u/at/chstan/nfs/summer_2014/ForConrad/files/150kevts_qcd_area/2014_10_14_11h23m44s/10s_0mu_0lw_0rec_25cut.root";
    std::string qcd_location_50 = "/u/at/chstan/nfs/summer_2014/ForConrad/files/150kevts_qcd_area/2014_10_14_11h23m44s/10s_0mu_0lw_0rec_50cut.root";
    std::string wprime_location_5 = "/u/at/chstan/nfs/summer_2014/ForConrad/files/150kevts_wprime_area/2014_10_15_00h40m54s/10s_0mu_0lw_0rec_5cut.root";
    std::string wprime_location_15 = "/u/at/chstan/nfs/summer_2014/ForConrad/files/150kevts_wprime_area/2014_10_15_00h40m54s/10s_0mu_0lw_0rec_15cut.root";
    std::string wprime_location_25 = "/u/at/chstan/nfs/summer_2014/ForConrad/files/150kevts_wprime_area/2014_10_15_00h40m54s/10s_0mu_0lw_0rec_25cut.root";
    std::string wprime_location_50 = "/u/at/chstan/nfs/summer_2014/ForConrad/files/150kevts_wprime_area/2014_10_15_00h40m54s/10s_0mu_0lw_0rec_50cut.root";
    std::string zprime_location_5 = "/u/at/chstan/nfs/summer_2014/ForConrad/files/150kevts_zprime_area/2014_10_13_22h17m12s/10s_0mu_0lw_0rec_5cut.root";
    std::string zprime_location_15 = "/u/at/chstan/nfs/summer_2014/ForConrad/files/150kevts_zprime_area/2014_10_13_22h17m12s/10s_0mu_0lw_0rec_15cut.root";
    std::string zprime_location_25 = "/u/at/chstan/nfs/summer_2014/ForConrad/files/150kevts_zprime_area/2014_10_13_22h17m12s/10s_0mu_0lw_0rec_25cut.root";
    std::string zprime_location_50 = "/u/at/chstan/nfs/summer_2014/ForConrad/files/150kevts_zprime_area/2014_10_13_22h17m12s/10s_0mu_0lw_0rec_50cut.root";

    static const std::vector<std::string> all_algs {"mGMM", "mGMMs", "mGMMc", "mTGMM", "mTGMMs", "mUMM", "CA", "antikt"};
    static const std::vector<std::string> cs_algs {"CA", "antikt"};
    static const std::vector<std::string> cs_algs_all
    {"CA", "antikt", "antikt_trimmed_two", "antikt_trimmed_three"};

    static const std::vector<std::string> algs {"mGMM", "mGMMs", "mGMMc", "mTGMM", "mTGMMs", "mUMM"};
    static const std::vector<std::string> event_labels {"qcd", "wprime", "zprime"};

    static const std::vector<float> pT_bins {0, 250, 300, 350, 400, 450, 500, 550, 600, 700, 1000};

    EventManager manager;
    manager.InstallEvent("qcd_5", qcd_location_5);
    manager.InstallEvent("qcd_15", qcd_location_15);
    manager.InstallEvent("qcd_25", qcd_location_25);
    //manager.InstallEvent("qcd_50", qcd_location_50);

    manager.InstallEvent("wprime_5", wprime_location_5);
    manager.InstallEvent("wprime_15", wprime_location_15);
    manager.InstallEvent("wprime_25", wprime_location_25);
    //manager.InstallEvent("wprime_50", wprime_location_50);

    manager.InstallEvent("zprime_5", zprime_location_5);
    manager.InstallEvent("zprime_15", zprime_location_15);
    manager.InstallEvent("zprime_25", zprime_location_25);
    //manager.InstallEvent("zprime_50", zprime_location_50);

    // INSTALL HISTOGRAMS
    // SOME TESTS FIRST
    manager << new FuzzyAntiktPtCorrelation()
            << new EfficiencyGenTest()

    // RADIUS PLOTS
            << new RadiusComparisonHistogram()
            << new AverageRadiusComparisonHistogram()


    // POSTER STUFF
            << new SigmaJetSizeCorrelationPoster()
            << new SigmaEfficiencyPosterPlot()
            << new SkewEfficiencyPlot()
            << new SkewHistogram();

    // Tagging efficiencies using sigma
    for (unsigned int event_label_iter = 0; event_label_iter < event_labels.size(); event_label_iter++) {
        for (unsigned int other_event_label_iter = 0; other_event_label_iter < event_labels.size(); other_event_label_iter++) {
            for (unsigned int pT_bin_iter = 0; pT_bin_iter < pT_bins.size() - 1; pT_bin_iter++) {
                if (event_label_iter == other_event_label_iter) continue;
                std::string event_label = event_labels[event_label_iter];
                std::string other_event_label = event_labels[other_event_label_iter];
                float cut_low = pT_bins[pT_bin_iter];
                float cut_high = pT_bins[pT_bin_iter + 1];
                manager << new SigmaEfficiencyPlot(event_label + "_5",
                                                   other_event_label + "_5",
                                                   cut_low, cut_high);
                std::vector<std::string> postfixes = {"_5", "_15", "_25"};

                for (unsigned int postfix_iter = 0; postfix_iter < postfixes.size();
                     postfix_iter++) {
                    std::string cpf = postfixes.at(postfix_iter);

                    manager << new SigmaNSubjettinessEfficiencyPlot(event_label + cpf,
                                                                    other_event_label + cpf,
                                                                    cut_low, cut_high);
                    manager << new SigmaImprovementEfficiencyTau32(event_label + cpf,
                                                                   other_event_label + cpf,
                                                                   cut_low, cut_high);
                    manager << new SigmaImprovementEfficiencyTau21(event_label + cpf,
                                                                   other_event_label + cpf,
                                                                   cut_low, cut_high);
                    manager << new SigmaImprovementEfficiencyMultiTau(event_label + cpf,
                                                                      other_event_label + cpf,
                                                                      cut_low, cut_high);
                    manager << new AreaEfficiency(event_label + cpf,
                                                  other_event_label + cpf,
                                                  cut_low, cut_high);
                }
                for (unsigned int alg_iter = 0; alg_iter < algs.size(); alg_iter++) {
                    std::string alg = algs[alg_iter];
                    manager << new FuzzyJetMassEfficiencyPlot(event_label + "_5",
                                                              other_event_label + "_5",
                                                              alg, cut_low, cut_high);
                }
            }
        }
    }

    // DeltaR plots
    for (unsigned int alg_iter = 0; alg_iter < algs.size(); alg_iter++) {
        std::string alg = algs[alg_iter];
        manager << new DeltaRHistogram(alg);
    }

    // Weight distance correlation plots + jet multiplicity
    for (unsigned int event_label_iter = 0; event_label_iter < event_labels.size(); event_label_iter++) {
        std::string event_label = event_labels[event_label_iter];
        manager << new JetMultiplicityPtCut(event_label); // deliberately no cut postfix, see ::Update
        manager << new RadiusPtSeedHistogram(event_label); // deliberately no cut postfix
        manager << new AverageRadiusPtSeedHistogram(event_label); // deliberately no cut postfix

        for (unsigned int alg_iter = 0; alg_iter < algs.size(); alg_iter++) {
            std::string alg = algs[alg_iter];
            manager << new WeightDistanceCorrelation(event_label + "_5", alg);
        }
    }

    // Sigma Area Correlation
    for (unsigned int event_label_iter = 0; event_label_iter < event_labels.size(); event_label_iter++) {
        for (unsigned int cs_alg_iter = 0; cs_alg_iter < cs_algs_all.size(); cs_alg_iter++) {
            std::string event_label = event_labels.at(event_label_iter);
            std::string cs_alg = cs_algs_all.at(cs_alg_iter);
            manager << new SigmaAreaCorrelation(event_label + "_5", cs_alg);
        }
    }

    // Mass and pT correlation plots
    for (unsigned int event_label_iter = 0; event_label_iter < event_labels.size(); event_label_iter++) {
        for (unsigned int alg_iter = 0; alg_iter < all_algs.size(); alg_iter++) {
            for (unsigned int other_alg_iter = 0; other_alg_iter < all_algs.size(); other_alg_iter++) {
                std::string event_label = event_labels[event_label_iter];
                std::string alg = all_algs[alg_iter];
                std::string other_alg = all_algs[other_alg_iter];
                manager << new MassCorrelation(event_label + "_5", alg, other_alg);
                manager << new PtCorrelation(event_label + "_5", alg, other_alg);
            }
        }
    }

    for (unsigned int event_label_iter = 0; event_label_iter < event_labels.size(); event_label_iter++) {
        for (unsigned int alg_iter = 0; alg_iter < all_algs.size(); alg_iter++) {
                std::string event_label = event_labels[event_label_iter];
                std::string alg = all_algs[alg_iter];
                manager << new SigmaJetSizeCorrelation(event_label + "_5", alg);
        }
    }

    // Start running the analysis
    manager.SetEventCount(5000);
    manager.Init();
    manager.PreparePtReweighting("qcd_5");

    manager.StartUpdaters();
    while (manager.NextEvent()) {
        manager.UpdateUpdaters();
    }
    manager.FinishUpdaters();
}

int main(int argc, char *argv[]) {
    std::cout << "Called as: ";
    for (int i = 0; i < argc; i++) {
        std::cout << argv[i] << " ";
    }
    std::cout << std::endl;

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "Produces help message");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help") > 0) {
        std::cout << desc << std::endl;
        return 1;
    }



    SetAtlasStyle();

    gROOT->ProcessLine("#include <vector>");
    EventTest();
    return 0;

    pTOverlay();

    int current_hist = fancyWeightHistos();

    static const int sizes_arr[] = {7, 8, 9, 10};
    static const int NPVs_arr[] = {0, 5};
    static const int learns_arr[] = {0, 1};
    static const int do_recs_arr[] = {0};
    static const std::string algs_arr[] = {"mGMM", "mGMMs", "mGMMc", "mUMM", "mTGMM", "mTGMMs"};

    static const Int_t colors_arr[] = {kMagenta + 3, kViolet + 9, kTeal - 5, kGray + 3};
    static const Style_t styles_arr[] = {kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown};

    std::vector<int> sizes(sizes_arr, sizes_arr+sizeof(sizes_arr) / sizeof(sizes_arr[0]));
    std::vector<int> NPVs(NPVs_arr, NPVs_arr+sizeof(NPVs_arr) / sizeof(NPVs_arr[0]));
    std::vector<int> learns(learns_arr, learns_arr+sizeof(learns_arr) / sizeof(learns_arr[0]));
    std::vector<int> do_recs(do_recs_arr, do_recs_arr+sizeof(do_recs_arr) / sizeof(do_recs_arr[0]));
    std::vector<std::string> algs(algs_arr, algs_arr+sizeof(algs_arr) / sizeof(algs_arr[0]));

    std::string file_prefix = "/u/at/chstan/nfs/summer_2014/ForConrad/files/20kevts_wprime_mu0_and_mu5_norec/2014_08_13_14h24m44s/";

    file_map_t file_m = constructFileMap(sizes, learns, do_recs, NPVs, file_prefix);

    for (auto iter = file_m.begin(); iter != file_m.end(); iter++) {
      std::cout << ":" << iter->second << std::endl;
    }

    std::string out_dir = "/u/at/chstan/nfs/summer_2014/ForConrad/results/plots/20kevts_wprime_mu0_and_mu5_norec/";

    // KEYS ARE BUILT BY (SIZE, LEARN, PILEUP)
    // DO PILEUP COMPARISONS: MASS RESOLUTION
    std::vector<HistHelper> v_hist_decs;

    std::stringstream ss;

    int n_bins = 100;

    std::string draw_opt = "p";
    std::string draw_opt_ref = "";
    StyleTypes::HistOptions options = StyleTypes::NONE;
    StyleTypes::HistOptions options_ref = StyleTypes::STRIPED;

    // god I'd love a for each loop guys please update the system gcc!
    for (unsigned int do_rec_iter = 0; do_rec_iter < do_recs.size(); do_rec_iter++) {
        for (unsigned int alg_iter = 0; alg_iter < algs.size(); alg_iter++) {
            for (unsigned int size_iter = 0; size_iter < sizes.size(); size_iter++) {
                for(unsigned int learn_iter = 0; learn_iter < learns.size(); learn_iter++) {
                    v_hist_decs.clear();

                    int do_rec = do_recs[do_rec_iter];
                    int size = sizes[size_iter];
                    int learn = learns[learn_iter];

                    std::string alg = algs[alg_iter];
                    ss.str(std::string());
                    ss << current_hist << " Mass resolution with pileup " << alg << " jets sz " << size << " lw " << learn << " rec " << do_rec;
                    current_hist++;
                    std::string canvas_title = ss.str();
                    CanvasHelper c_dec("Mass [GeV]", "", canvas_title, out_dir, 800, 800);
                    c_dec.diff_scale = false;
                    for(unsigned int npv_iter = 0; npv_iter < NPVs.size(); npv_iter++) {
                        int npv = NPVs[npv_iter];
                        ss.str(std::string());
                        ss << alg << " Mass - PuV " << npv;
                        std::string title = ss.str();
                        ss.str(std::string());
                        ss << alg << "_m";
                        std::string branch = ss.str();
                        HistHelper hist_helper_temp(file_m[std::make_tuple(size, learn, npv, do_rec)],
                                                    branch, title, 510, mass_low, mass_high, n_bins,
                                                    options, colors_arr[npv_iter],
                                                    styles_arr[npv_iter], draw_opt);
                        v_hist_decs.push_back(hist_helper_temp);
                    }
                    HistHelper hist_helper_temp(file_m[std::make_tuple(size, learn, 0, do_rec)],
                                                "antikt_m_trimmed_three", "Anti-kt Trimmed Mass - PuV 0",
                                                510, mass_low, mass_high, n_bins, options_ref, kBlue,
                                                kFullCircle, draw_opt_ref);
                    v_hist_decs.push_back(hist_helper_temp);
                    prettyHist<float>(v_hist_decs, c_dec);
                }
            }
        }
    }

    // HISTOGRAMS FOR SIZE COMPARISON
    for (unsigned int do_rec_iter = 0; do_rec_iter < do_recs.size(); do_rec_iter++) {
        for (unsigned int alg_iter = 0; alg_iter < algs.size(); alg_iter++) {
            for (unsigned int npv_iter = 0; npv_iter < NPVs.size(); npv_iter++) {
                for(unsigned int learn_iter = 0; learn_iter < learns.size(); learn_iter++) {
                    v_hist_decs.clear();

                    int do_rec = do_recs[do_rec_iter];
                    int learn = learns[learn_iter];
                    int npv = NPVs[npv_iter];
                    std::string alg = algs[alg_iter];
                    ss.str(std::string());
                    ss << current_hist << " Mass resolution by size " << alg << " lw " << learn << " rec " << do_rec << " npu " << npv;
                    current_hist++;
                    std::string canvas_title = ss.str();
                    CanvasHelper c_dec("Mass [GeV]", "", canvas_title, out_dir, 800, 800);
                    c_dec.diff_scale = false;
                    for(unsigned int size_iter = 0; size_iter < sizes.size(); size_iter++) {
                        int size = sizes[size_iter];
                        ss.str(std::string());
                        ss << alg << " Mass - Size " << size;
                        std::string title = ss.str();
                        ss.str(std::string());
                        ss << alg << "_m";
                        std::string branch = ss.str();
                        HistHelper hist_helper_temp(file_m[std::make_tuple(size, learn, npv, do_rec)],
                                                    branch, title, 510, mass_low, mass_high, n_bins,
                                                    options, colors_arr[size_iter],
                                                    styles_arr[size_iter], draw_opt);
                        v_hist_decs.push_back(hist_helper_temp);
                    }
                    HistHelper hist_helper_temp(file_m[std::make_tuple(10, learn, npv, do_rec)],
                                                "antikt_m_trimmed_three", "Anti-kt trimmed Mass",
                                                510, mass_low, mass_high, n_bins, options_ref, kBlue,
                                                kFullCircle, draw_opt_ref);
                    v_hist_decs.push_back(hist_helper_temp);
                    prettyHist<float>(v_hist_decs, c_dec);
                }
            }
        }
    }

    // HISTOGRAMS FOR WEIGHT LEARNING COMPARISON
    for (unsigned int do_rec_iter = 0; do_rec_iter < do_recs.size(); do_rec_iter++) {
        for (unsigned int alg_iter = 0; alg_iter < algs.size(); alg_iter++) {
            for (unsigned int npv_iter = 0; npv_iter < NPVs.size(); npv_iter++) {
                for(unsigned int size_iter = 0; size_iter < sizes.size(); size_iter++) {
                    v_hist_decs.clear();

                    int do_rec = do_recs[do_rec_iter];
                    int size = sizes[size_iter];
                    int npv = NPVs[npv_iter];
                    std::string alg = algs[alg_iter];
                    ss.str(std::string());
                    ss << current_hist << " Mass resolution weight learning " << alg << " sz " << size << " rec " << do_rec << " npu " << npv;
                    current_hist++;
                    std::string canvas_title = ss.str();
                    CanvasHelper c_dec("Mass [GeV]", "", canvas_title, out_dir, 800, 800);
                    c_dec.diff_scale = false;
                    for(unsigned int learn_iter = 0; learn_iter < learns.size(); learn_iter++) {
                        int learn = learns[learn_iter];
                        ss.str(std::string());
                        ss << alg << " Mass - Learning " << learn;
                        std::string title = ss.str();
                        ss.str(std::string());
                        ss << alg << "_m";
                        std::string branch = ss.str();
                        HistHelper hist_helper_temp(file_m[std::make_tuple(size, learn, npv, do_rec)],
                                                    branch, title, 510, mass_low, mass_high, n_bins,
                                                    options, colors_arr[learn_iter],
                                                    styles_arr[learn_iter], draw_opt);
                        v_hist_decs.push_back(hist_helper_temp);
                    }
                    HistHelper hist_helper_temp(file_m[std::make_tuple(size, 0, npv, do_rec)],
                                                "antikt_m_trimmed_three", "Anti-kt trimmed Mass",
                                                510, mass_low, mass_high, n_bins, options_ref, kBlue,
                                                kFullCircle, draw_opt_ref);
                    v_hist_decs.push_back(hist_helper_temp);
                    prettyHist<float>(v_hist_decs, c_dec);
                }
            }
        }
    }

    // HISTOGRAMS FOR RECOMBINATION SCHEME
    for (unsigned int learn_iter = 0; learn_iter < learns.size(); learn_iter++) {
        for (unsigned int alg_iter = 0; alg_iter < algs.size(); alg_iter++) {
            for (unsigned int npv_iter = 0; npv_iter < NPVs.size(); npv_iter++) {
                for(unsigned int size_iter = 0; size_iter < sizes.size(); size_iter++) {
                    v_hist_decs.clear();

                    int learn = learns[learn_iter];
                    int size = sizes[size_iter];
                    int npv = NPVs[npv_iter];
                    std::string alg = algs[alg_iter];
                    ss.str(std::string());
                    ss << current_hist << " Mass resolution recombination scheme " << alg << " sz " << size << " lw " << learn << " npu " << npv;
                    current_hist++;
                    std::string canvas_title = ss.str();
                    CanvasHelper c_dec("Mass [GeV]", "", canvas_title, out_dir, 800, 800);
                    c_dec.diff_scale = false;
                    for(unsigned int do_rec_iter = 0; do_rec_iter < do_recs.size(); do_rec_iter++) {
                        int do_rec = do_recs[do_rec_iter];
                        ss.str(std::string());
                        std::string scheme = do_rec_iter == 1 ? "rec all particles" : "fixed C/A locations";
                        ss << alg << " Mass - " << scheme;
                        std::string title = ss.str();
                        ss.str(std::string());
                        ss << alg << "_m";
                        std::string branch = ss.str();
                        HistHelper hist_helper_temp(file_m[std::make_tuple(size, learn, npv, do_rec)],
                                                    branch, title, 510, mass_low, mass_high, n_bins,
                                                    options, colors_arr[do_rec_iter],
                                                    styles_arr[do_rec_iter], draw_opt);
                        v_hist_decs.push_back(hist_helper_temp);
                    }
                    HistHelper hist_helper_temp(file_m[std::make_tuple(size, learn, npv, 0)],
                                                "antikt_m_trimmed_three", "Anti-kt trimmed Mass",
                                                510, mass_low, mass_high, n_bins, options_ref, kBlue,
                                                kFullCircle, draw_opt_ref);
                    v_hist_decs.push_back(hist_helper_temp);
                    prettyHist<float>(v_hist_decs, c_dec);
                }
            }
        }
    }

    /// ==========================SOFT MASS==================================
    for (unsigned int do_rec_iter = 0; do_rec_iter < do_recs.size(); do_rec_iter++) {
        for (unsigned int alg_iter = 0; alg_iter < algs.size(); alg_iter++) {
            for (unsigned int size_iter = 0; size_iter < sizes.size(); size_iter++) {
                for(unsigned int learn_iter = 0; learn_iter < learns.size(); learn_iter++) {
                    v_hist_decs.clear();

                    int do_rec = do_recs[do_rec_iter];
                    int size = sizes[size_iter];
                    int learn = learns[learn_iter];

                    std::string alg = algs[alg_iter];
                    ss.str(std::string());
                    ss << current_hist << " Soft mass resolution with pileup " << alg << " jets sz " << size << " lw " << learn << " rec " << do_rec;
                    current_hist++;
                    std::string canvas_title = ss.str();
                    CanvasHelper c_dec("Mass [GeV]", "", canvas_title, out_dir, 800, 800);
                    c_dec.diff_scale = false;
                    for(unsigned int npv_iter = 0; npv_iter < NPVs.size(); npv_iter++) {
                        int npv = NPVs[npv_iter];
                        ss.str(std::string());
                        ss << alg << " Soft Mass - PuV " << npv;
                        std::string title = ss.str();
                        ss.str(std::string());
                        ss << alg << "_m_soft";
                        std::string branch = ss.str();
                        HistHelper hist_helper_temp(file_m[std::make_tuple(size, learn, npv, do_rec)],
                                                    branch, title, 510, mass_low, mass_high, n_bins,
                                                    options, colors_arr[npv_iter],
                                                    styles_arr[npv_iter], draw_opt);
                        v_hist_decs.push_back(hist_helper_temp);
                    }
                    HistHelper hist_helper_temp(file_m[std::make_tuple(size, learn, 0, do_rec)],
                                                "antikt_m_trimmed_three", "Anti-kt Trimmed Mass - PuV 0",
                                                510, mass_low, mass_high, n_bins, options_ref, kBlue,
                                                kFullCircle, draw_opt_ref);
                    v_hist_decs.push_back(hist_helper_temp);
                    prettyHist<float>(v_hist_decs, c_dec);
                }
            }
        }
    }

    // HISTOGRAMS FOR SIZE COMPARISON
    for (unsigned int do_rec_iter = 0; do_rec_iter < do_recs.size(); do_rec_iter++) {
        for (unsigned int alg_iter = 0; alg_iter < algs.size(); alg_iter++) {
            for (unsigned int npv_iter = 0; npv_iter < NPVs.size(); npv_iter++) {
                for(unsigned int learn_iter = 0; learn_iter < learns.size(); learn_iter++) {
                    v_hist_decs.clear();

                    int do_rec = do_recs[do_rec_iter];
                    int learn = learns[learn_iter];
                    int npv = NPVs[npv_iter];
                    std::string alg = algs[alg_iter];
                    ss.str(std::string());
                    ss << current_hist << " Soft mass resolution by size " << alg << " lw " << learn << " rec " << do_rec << " npu " << npv;
                    current_hist++;
                    std::string canvas_title = ss.str();
                    CanvasHelper c_dec("Mass [GeV]", "", canvas_title, out_dir, 800, 800);
                    c_dec.diff_scale = false;
                    for(unsigned int size_iter = 0; size_iter < sizes.size(); size_iter++) {
                        int size = sizes[size_iter];
                        ss.str(std::string());
                        ss << alg << " Soft Mass - Size " << size;
                        std::string title = ss.str();
                        ss.str(std::string());
                        ss << alg << "_m_soft";
                        std::string branch = ss.str();
                        HistHelper hist_helper_temp(file_m[std::make_tuple(size, learn, npv, do_rec)],
                                                    branch, title, 510, mass_low, mass_high, n_bins,
                                                    options, colors_arr[size_iter],
                                                    styles_arr[size_iter], draw_opt);
                        v_hist_decs.push_back(hist_helper_temp);
                    }
                    HistHelper hist_helper_temp(file_m[std::make_tuple(10, learn, npv, do_rec)],
                                                "antikt_m_trimmed_three", "Anti-kt trimmed Mass",
                                                510, mass_low, mass_high, n_bins, options_ref, kBlue,
                                                kFullCircle, draw_opt_ref);
                    v_hist_decs.push_back(hist_helper_temp);
                    prettyHist<float>(v_hist_decs, c_dec);
                }
            }
        }
    }

    // HISTOGRAMS FOR WEIGHT LEARNING COMPARISON
    for (unsigned int do_rec_iter = 0; do_rec_iter < do_recs.size(); do_rec_iter++) {
        for (unsigned int alg_iter = 0; alg_iter < algs.size(); alg_iter++) {
            for (unsigned int npv_iter = 0; npv_iter < NPVs.size(); npv_iter++) {
                for(unsigned int size_iter = 0; size_iter < sizes.size(); size_iter++) {
                    v_hist_decs.clear();

                    int do_rec = do_recs[do_rec_iter];
                    int size = sizes[size_iter];
                    int npv = NPVs[npv_iter];
                    std::string alg = algs[alg_iter];
                    ss.str(std::string());
                    ss << current_hist << " Soft mass resolution weight learning " << alg << " sz " << size << " rec " << do_rec << " npu " << npv;
                    current_hist++;
                    std::string canvas_title = ss.str();
                    CanvasHelper c_dec("Mass [GeV]", "", canvas_title, out_dir, 800, 800);
                    c_dec.diff_scale = false;
                    for(unsigned int learn_iter = 0; learn_iter < learns.size(); learn_iter++) {
                        int learn = learns[learn_iter];
                        ss.str(std::string());
                        ss << alg << " Soft Mass - Learning " << learn;
                        std::string title = ss.str();
                        ss.str(std::string());
                        ss << alg << "_m_soft";
                        std::string branch = ss.str();
                        HistHelper hist_helper_temp(file_m[std::make_tuple(size, learn, npv, do_rec)],
                                                    branch, title, 510, mass_low, mass_high, n_bins,
                                                    options, colors_arr[learn_iter],
                                                    styles_arr[learn_iter], draw_opt);
                        v_hist_decs.push_back(hist_helper_temp);
                    }
                    HistHelper hist_helper_temp(file_m[std::make_tuple(size, 0, npv, do_rec)],
                                                "antikt_m_trimmed_three", "Anti-kt trimmed Mass",
                                                510, mass_low, mass_high, n_bins, options_ref, kBlue,
                                                kFullCircle, draw_opt_ref);
                    v_hist_decs.push_back(hist_helper_temp);
                    prettyHist<float>(v_hist_decs, c_dec);
                }
            }
        }
    }

    // HISTOGRAMS FOR RECOMBINATION SCHEME
    for (unsigned int learn_iter = 0; learn_iter < learns.size(); learn_iter++) {
        for (unsigned int alg_iter = 0; alg_iter < algs.size(); alg_iter++) {
            for (unsigned int npv_iter = 0; npv_iter < NPVs.size(); npv_iter++) {
                for(unsigned int size_iter = 0; size_iter < sizes.size(); size_iter++) {
                    v_hist_decs.clear();

                    int learn = learns[learn_iter];
                    int size = sizes[size_iter];
                    int npv = NPVs[npv_iter];
                    std::string alg = algs[alg_iter];
                    ss.str(std::string());
                    ss << current_hist << " Soft mass resolution recombination scheme " << alg << " sz " << size << " lw " << learn << " npu " << npv;
                    current_hist++;
                    std::string canvas_title = ss.str();
                    CanvasHelper c_dec("Mass [GeV]", "", canvas_title, out_dir, 800, 800);
                    c_dec.diff_scale = false;
                    for(unsigned int do_rec_iter = 0; do_rec_iter < do_recs.size(); do_rec_iter++) {
                        int do_rec = do_recs[do_rec_iter];
                        ss.str(std::string());
                        std::string scheme = do_rec_iter == 1 ? "rec all particles" : "fixed C/A locations";
                        ss << alg << " Soft Mass - " << scheme;
                        std::string title = ss.str();
                        ss.str(std::string());
                        ss << alg << "_m_soft";
                        std::string branch = ss.str();
                        HistHelper hist_helper_temp(file_m[std::make_tuple(size, learn, npv, do_rec)],
                                                    branch, title, 510, mass_low, mass_high, n_bins,
                                                    options, colors_arr[do_rec_iter],
                                                    styles_arr[do_rec_iter], draw_opt);
                        v_hist_decs.push_back(hist_helper_temp);
                    }
                    HistHelper hist_helper_temp(file_m[std::make_tuple(size, learn, npv, 0)],
                                                "antikt_m_trimmed_three", "Anti-kt trimmed Mass",
                                                510, mass_low, mass_high, n_bins, options_ref, kBlue,
                                                kFullCircle, draw_opt_ref);
                    v_hist_decs.push_back(hist_helper_temp);
                    prettyHist<float>(v_hist_decs, c_dec);
                }
            }
        }
    }

    /// ======================RADIUS LEARNING STUDIES========================
    static const std::string branches_arr[] =
        {"mGMMc_r", "mGMMc_r_avg", "mGMMc_r_weighted_avg"};
    static const std::string radius_titles_arr[] =
        {"Lead jet radius", "Average jet radius", "Weighted average jet radius"};
    for (unsigned int do_rec_iter = 0; do_rec_iter < do_recs.size(); do_rec_iter++) {
        for (unsigned int size_iter = 0; size_iter <  sizes.size(); size_iter++) {
            for (unsigned int learn_iter = 0; learn_iter < learns.size(); learn_iter++) {
                for (unsigned int npv_iter = 0; npv_iter < NPVs.size(); npv_iter++) {
                    for (unsigned int branch_iter = 0; branch_iter < 3; branch_iter++) {
                        v_hist_decs.clear();

                        int do_rec = do_recs[do_rec_iter];
                        int size = sizes[size_iter];
                        int learn = learns[learn_iter];
                        int npv = NPVs[npv_iter];
                        std::string branch = branches_arr[branch_iter];
                        std::string partial_title = radius_titles_arr[branch_iter];

                        ss.str(std::string());
                        ss << current_hist << " " << partial_title << " mGMMc jets sz "
                           << size << " lw " << learn << " rec " << do_rec << " pu " << npv;
                        current_hist++;
                        std::string canvas_title = ss.str();

                        CanvasHelper c_dec(partial_title, "", canvas_title, out_dir, 800, 800);
                        c_dec.diff_scale = false;
                        HistHelper hist_helper_temp(file_m[std::make_tuple(size, learn, npv, do_rec)],
                                                    branch, partial_title, 510,
                                                    0, 1.5, n_bins,
                                                    options_ref, kBlue, kFullCircle, draw_opt_ref);
                        v_hist_decs.push_back(hist_helper_temp);
                        prettyHist<float>(v_hist_decs, c_dec);
                    }
                }
            }
        }
    }

    /// ========================FUZZINESS STUDIES============================
    static const std::string branch_postfix_arr[] =
        {"_m_mean", "_m_var", "_m_skew", "_pt_mean", "_pt_var", "_pt_skew"};
    static const std::string fuzzy_titles_arr[] =
        {"Mass Mean ", "Mass Variance ", "Mass Skew ", "pT Mean ", "pT Variance ", "pT Skew"};
    static const double low_edges[] =
        {0, -20, -20, 0, -20, 20};
    static const double high_edges[] =
        {400, 20, 20, 800, 20, 20};

    for (unsigned int do_rec_iter = 0; do_rec_iter < do_recs.size(); do_rec_iter++) {
        for (unsigned int alg_iter = 0; alg_iter < algs.size(); alg_iter++) {
            for (unsigned int size_iter = 0; size_iter < sizes.size(); size_iter++) {
                for(unsigned int learn_iter = 0; learn_iter < learns.size(); learn_iter++) {
                    for (unsigned int npv_iter = 0; npv_iter < NPVs.size(); npv_iter++) {
                        for (unsigned int branch_iter = 0; branch_iter < 6; branch_iter++) {
                            v_hist_decs.clear();

                            int do_rec = do_recs[do_rec_iter];
                            int size = sizes[size_iter];
                            int learn = learns[learn_iter];
                            int npv = NPVs[npv_iter];
                            std::string alg = algs[alg_iter];
                            double low_edge = low_edges[branch_iter];
                            double high_edge = high_edges[branch_iter];

                            ss.str(std::string());
                            ss << current_hist << " " << fuzzy_titles_arr[branch_iter] << alg << " jets sz "
                               << size << " lw " << learn << " rec " << do_rec << " pu " << npv;
                            current_hist++;
                            std::string canvas_title = ss.str();
                            CanvasHelper c_dec(fuzzy_titles_arr[branch_iter], "", canvas_title, out_dir, 800, 800);
                            c_dec.diff_scale = false;
                            ss.str(std::string());
                            ss << alg << branch_postfix_arr[branch_iter];
                            std::string branch = ss.str();
                            HistHelper hist_helper_temp(file_m[std::make_tuple(size, learn, npv, do_rec)],
                                                        branch, fuzzy_titles_arr[branch_iter], 510,
                                                        low_edge, high_edge, n_bins,
                                                        options_ref, kBlue, kFullCircle, draw_opt_ref);
                            v_hist_decs.push_back(hist_helper_temp);
                            prettyHist<float>(v_hist_decs, c_dec);
                        }
                    }
                }
            }
        }
    }

    /// ==========================PILEUP STUDIES=============================
    for (unsigned int do_rec_iter = 0; do_rec_iter < do_recs.size(); do_rec_iter++) {
        for (unsigned int alg_iter = 0; alg_iter < algs.size(); alg_iter++) {
            for (unsigned int size_iter = 0; size_iter < sizes.size(); size_iter++) {
                for(unsigned int learn_iter = 0; learn_iter < learns.size(); learn_iter++) {
                    v_hist_decs.clear();

                    int do_rec = do_recs[do_rec_iter];
                    int size = sizes[size_iter];
                    int learn = learns[learn_iter];
                    std::string alg = algs[alg_iter];
                    ss.str(std::string());
                    ss << current_hist << " Unclustered pileup " << alg << " jets sz " << size << " lw " << learn << " rec " << do_rec;
                    std::string canvas_title = ss.str();
                    current_hist++;
                    CanvasHelper c_dec("Fraction of Clustered Particles", "", canvas_title, out_dir, 800, 800);
                    c_dec.diff_scale = false;
                    for(unsigned int npv_iter = 1; npv_iter < NPVs.size(); npv_iter++) {

                        int npv = NPVs[npv_iter];
                        ss.str(std::string());
                        ss << alg << " - PuV " << npv;
                        std::string title = ss.str();
                        ss.str(std::string());
                        ss << alg << "_ucpu";
                        std::string branch = ss.str();
                        HistHelper hist_helper_temp(file_m[std::make_tuple(size, learn, npv, do_rec)],
                                                    branch, title, 510, 0, 1, n_bins,
                                                    options, colors_arr[npv_iter],
                                                    styles_arr[npv_iter], draw_opt_ref);
                        v_hist_decs.push_back(hist_helper_temp);
                    }
                    prettyHist<float>(v_hist_decs, c_dec);
                }
            }
        }
    }

    return 0;
}
