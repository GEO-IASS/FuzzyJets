#include <iostream.h>
#include <string.h>
#include <cstddef>
#include <sstream>
#include <map>
#include <vector>

#include <TStyle.h>
#include <THStack.h>
#include <TH1F.h>
#include <TH1D.h>

#include "TColor.h"
#include "AtlasUtils.cc"
#include "AtlasStyle.h"

//double HistMaximum(std::vector<TH1D *> const& hists) {
//    double max = -1;
//
//    const size_t n_hists = hists.size();
//    for (size_t hist_iter = 0; hist_iter < n_hists; hist_iter++) {
//        double n_max = hists[hist_iter]->GetMaximum();
//        if (n_max > max) {
//            max = n_max;
//        }
//    }
//    return max;
//}

size_t RoundDoubleUp(double d, size_t q) {
    size_t i = (size_t) ceil(d);
    size_t r = i % q;
    size_t o = i - r;
    if (r) o += q;
    return o;
}

std::string fileify(std::string text) {
    for(std::string::iterator it = text.begin(); it != text.end(); ++it) {
        if(isspace(*it)) {
            *it = '_';
        }
    }
    return text;
}

void SetupATLASStyle() {
    gStyle->SetOptStat(0);
    gROOT->Reset();
    AtlasStyle();
    gROOT->SetStyle("ATLAS");
    gROOT->ForceStyle();
    gStyle->SetPadLeftMargin(0.16);
}

void DrawAtlasLabel(std::string title, double x, double y) {
    ATLAS_LABEL(x, y, kBlack);
    myText(x, y-0.08, kBlack, "#sqrt{s} = 8 TeV");
    myText(x, y-0.14, kBlack, title.c_str(), 0.4);
}

void prettify(TH1F *h, Color_t like_color, Int_t smashes) {
    h->SetLineColor(like_color);
    TColor *col = gROOT->GetColor(smashes);
    TColor *col_ref = gROOT->GetColor(like_color);
    col->SetRGB(col_ref->GetRed(), col_ref->GetGreen(), col_ref->GetBlue());
    col->SetAlpha(0.28);

    h->SetLineStyle(1);
    h->SetFillColor(smashes);
    h->SetFillStyle(1001);
    h->SetLineWidth(2);

    float norm = h->Integral(-1, h->GetNbinsX() + 1);
    h->Scale(1./norm);
}

std::string gen_x_label(Int_t npv) {
    std::stringstream ss;
    ss.str(std::string());

    //ss << "Coefficient of Variation: c_{v, <#sigma_{" << npv << "}> - #sigma_{0}}";
    ss << "#frac{#sigma_{" << npv << "} - #sigma_{0}}{#sigma_{0}}";

    return ss.str();
}

std::string gen_x_mass_label(Int_t npv) {
    return "Mass [GeV]";
}

std::string titley_bit() {
    return "norm_dev";
}

void NoiseStudy() {
    std::string base = "/u/at/chstan/nfs/summer_2014/ForConrad/";
    std::string out_directory = base + "results/plots/Simple/";
    base = base + "files/100k_zprime_noise/2015_03_10_18h48m07s/";
    TFile *f_ej0 = TFile::Open((base + "10s_0mu_0lw_0rec_0ejw_0off_0PP_0TBS_5cut.root").c_str());
    TFile *f_ej1 = TFile::Open((base + "10s_0mu_0lw_0rec_1ejw_0off_0PP_0TBS_5cut.root").c_str());
    TFile *f_ej10 = TFile::Open((base + "10s_0mu_0lw_0rec_10ejw_0off_0PP_0TBS_5cut.root").c_str());

    std::cout << f_ej0 << std::endl;
    std::cout << f_ej1 << std::endl;
    std::cout << f_ej10 << std::endl;

    size_t n_bins = 40;

    TH1F *hs_ndev_0 = new TH1F("hs_ndev_0", "", n_bins, -0.1, 0.1);
    TH1F *hs_ndev_1 = new TH1F("hs_ndev_1", "", n_bins, -0.1, 0.1);
    TH1F *hs_ndev_10 = new TH1F("hs_ndev_10", "", n_bins, -0.1, 0.1);

    float s_ndev_0;
    float s_ndev_1;
    float s_ndev_10;

    TTree *t_ej0 = (TTree *) f_ej0->Get("EventTree");
    TTree *t_ej1 = (TTree *) f_ej1->Get("EventTree");
    TTree *t_ej10 = (TTree *) f_ej10->Get("EventTree");

    t_ej0->SetBranchAddress("delta_sigma", &s_ndev_0);
    t_ej1->SetBranchAddress("delta_sigma", &s_ndev_1);
    t_ej10->SetBranchAddress("delta_sigma", &s_ndev_10);

    Long64_t nentries = t_ej0->GetEntries();
    std::cout << nentries << std::endl;
    Int_t nbytes = 0;
    for (Long64_t i = 0; i < nentries; i++) {
        nbytes += t_ej0->GetEntry(i);
        nbytes += t_ej1->GetEntry(i);
        nbytes += t_ej10->GetEntry(i);

        hs_ndev_0->Fill(s_ndev_0);
        hs_ndev_1->Fill(s_ndev_1);
        hs_ndev_10->Fill(s_ndev_10);
    }

    TCanvas canvas("canv", "", 600, 600);
    canvas.cd();

    {
        THStack hist_stack("temporary_stack", "");
        hist_stack.Add(hs_ndev_0);
        hist_stack.Add(hs_ndev_1);
        hist_stack.Add(hs_ndev_10);

        prettify(hs_ndev_0, kBlue, 3);
        prettify(hs_ndev_1, kRed, 4);
        prettify(hs_ndev_10, kBlack, 5);

        hist_stack.Draw("nostack");

        hist_stack.GetHistogram()->
            GetXaxis()->SetTitle("#Delta#sigma");
        hist_stack.GetHistogram()->GetYaxis()->SetTitle("Normalized to Unity");
        hist_stack.GetHistogram()->GetXaxis()->SetNdivisions(505);

        hist_stack.SetMaximum(1.4 * hist_stack.GetMaximum("nostack"));
        //hist_stack.SetMaximum(0.25);
        hist_stack.Draw("nostack");

        TLegend legend(0.65, 0.75, 0.95, 0.92);
        legend.AddEntry(hs_ndev_0, "No Event Jet", "f");
        legend.AddEntry(hs_ndev_1, "#gamma = 0.001", "f");
        legend.AddEntry(hs_ndev_10, "#gamma = 0.01", "f");

        legend.SetTextFont(42);
        legend.SetFillStyle(0);
        legend.SetFillColor(0);
        legend.SetBorderSize(0);
        legend.SetTextSize(1.3 * legend.GetTextSize());

        legend.Draw();

        //DrawAtlasLabel("");
        canvas.Update();

        canvas.Print((out_directory + "noise_sigma_" + titley_bit() + ".pdf").c_str(), "pdf");
        canvas.Clear();
    }

    {
        THStack hist_stack("temporary_stack", "");
        hist_stack.Add(hs_ndev_0);
        //hist_stack.Add(hs_ndev_1);
        //hist_stack.Add(hs_ndev_10);

        prettify(hs_ndev_0, kBlue, 3);
        //prettify(hs_ndev_1, kRed, 4);
        //prettify(hs_ndev_10, kBlack, 5);

        hist_stack.Draw("nostack");

        hist_stack.GetHistogram()->
            GetXaxis()->SetTitle("#Delta#sigma");
        hist_stack.GetHistogram()->GetYaxis()->SetTitle("Normalized to Unity");
        hist_stack.GetHistogram()->GetXaxis()->SetNdivisions(505);

        hist_stack.SetMaximum(1.4 * hist_stack.GetMaximum("nostack"));
        //hist_stack.SetMaximum(0.25);
        hist_stack.Draw("nostack");

//TLegend legend(0.65, 0.75, 0.95, 0.92);
//legend.AddEntry(hs_ndev_0, "No Event Jet", "f");
//legend.AddEntry(hs_ndev_1, "#gamma = 0.001", "f");
//legend.AddEntry(hs_ndev_10, "#gamma = 0.01", "f");
//
//legend.SetTextFont(42);
//legend.SetFillStyle(0);
//legend.SetFillColor(0);
//legend.SetBorderSize(0);
//legend.SetTextSize(1.3 * legend.GetTextSize());
//
//legend.Draw();

        //DrawAtlasLabel("");
        canvas.Update();

        canvas.Print((out_directory + "noise_sigma_simple_" + titley_bit() + ".pdf").c_str(), "pdf");
        canvas.Clear();
    }

    delete f_ej0;
    delete f_ej1;
    delete f_ej10;
}

void WanderingStudy() {
    std::string base = "/u/at/chstan/nfs/summer_2014/ForConrad/";
    std::string out_directory = base + "results/plots/Simple/";
    base = base + "files/30k_zprime_wandering/2015_03_10_21h15m56s/";
    TFile *f_ej0 = TFile::Open((base + "10s_0mu_0lw_0rec_0ejw_0off_0PP_0TBS_0SN.root").c_str());
    TFile *f_ej1 = TFile::Open((base + "10s_0mu_0lw_0rec_1ejw_0off_0PP_0TBS_0SN.root").c_str());
    TFile *f_ej10 = TFile::Open((base + "10s_0mu_0lw_0rec_10ejw_0off_0PP_0TBS_0SN.root").c_str());

    std::cout << f_ej0 << std::endl;
    std::cout << f_ej1 << std::endl;
    std::cout << f_ej10 << std::endl;

    size_t n_bins = 40;

    TH1F *hs_ndev_0_10 = new TH1F("hs_ndev_0_10", "", n_bins, -0.1, 0.1);
    TH1F *hs_ndev_1_10 = new TH1F("hs_ndev_1_10", "", n_bins, -0.1, 0.1);
    TH1F *hs_ndev_10_10 = new TH1F("hs_ndev_10_10", "", n_bins, -0.1, 0.1);

    TH1F *hs_phi_dev_0_10 = new TH1F("hs_phi_0_10", "", n_bins, -0.1, 0.1);
    TH1F *hs_phi_dev_1_10 = new TH1F("hs_phi_1_10", "", n_bins, -0.1, 0.1);
    TH1F *hs_phi_dev_10_10 = new TH1F("hs_phi_10_10", "", n_bins, -0.1, 0.1);

    TH1F *hs_eta_dev_0_10 = new TH1F("hs_eta_dev_0_10", "", n_bins, -0.1, 0.1);
    TH1F *hs_eta_dev_1_10 = new TH1F("hs_eta_dev_1_10", "", n_bins, -0.1, 0.1);
    TH1F *hs_eta_dev_10_10 = new TH1F("hs_eta_dev_10_10", "", n_bins, -0.1, 0.1);


    TH1F *hs_ndev_0_20 = new TH1F("hs_ndev_0_20", "", n_bins, -0.1, 0.1);
    TH1F *hs_ndev_1_20 = new TH1F("hs_ndev_1_20", "", n_bins, -0.1, 0.1);
    TH1F *hs_ndev_10_20 = new TH1F("hs_ndev_10_20", "", n_bins, -0.1, 0.1);

    TH1F *hs_phi_dev_0_20 = new TH1F("hs_phi_0_20", "", n_bins, -0.1, 0.1);
    TH1F *hs_phi_dev_1_20 = new TH1F("hs_phi_1_20", "", n_bins, -0.1, 0.1);
    TH1F *hs_phi_dev_10_20 = new TH1F("hs_phi_10_20", "", n_bins, -0.1, 0.1);

    TH1F *hs_eta_dev_0_20 = new TH1F("hs_eta_dev_0_20", "", n_bins, -0.1, 0.1);
    TH1F *hs_eta_dev_1_20 = new TH1F("hs_eta_dev_1_20", "", n_bins, -0.1, 0.1);
    TH1F *hs_eta_dev_10_20 = new TH1F("hs_eta_dev_10_20", "", n_bins, -0.1, 0.1);


    TH1F *hs_ndev_0_30 = new TH1F("hs_ndev_0_30", "", n_bins, -0.1, 0.1);
    TH1F *hs_ndev_1_30 = new TH1F("hs_ndev_1_30", "", n_bins, -0.1, 0.1);
    TH1F *hs_ndev_10_30 = new TH1F("hs_ndev_10_30", "", n_bins, -0.1, 0.1);

    TH1F *hs_phi_dev_0_30 = new TH1F("hs_phi_0_30", "", n_bins, -0.1, 0.1);
    TH1F *hs_phi_dev_1_30 = new TH1F("hs_phi_1_30", "", n_bins, -0.1, 0.1);
    TH1F *hs_phi_dev_10_30 = new TH1F("hs_phi_10_30", "", n_bins, -0.1, 0.1);

    TH1F *hs_eta_dev_0_30 = new TH1F("hs_eta_dev_0_30", "", n_bins, -0.1, 0.1);
    TH1F *hs_eta_dev_1_30 = new TH1F("hs_eta_dev_1_30", "", n_bins, -0.1, 0.1);
    TH1F *hs_eta_dev_10_30 = new TH1F("hs_eta_dev_10_30", "", n_bins, -0.1, 0.1);

    float s_ndev_0_10;
    float s_ndev_1_10;
    float s_ndev_10_10;

    float s_phi_dev_0_10;
    float s_phi_dev_1_10;
    float s_phi_dev_10_10;

    float s_eta_dev_0_10;
    float s_eta_dev_1_10;
    float s_eta_dev_10_10;

    float s_ndev_0_20;
    float s_ndev_1_20;
    float s_ndev_10_20;

    float s_phi_dev_0_20;
    float s_phi_dev_1_20;
    float s_phi_dev_10_20;

    float s_eta_dev_0_20;
    float s_eta_dev_1_20;
    float s_eta_dev_10_20;

    float s_ndev_0_30;
    float s_ndev_1_30;
    float s_ndev_10_30;

    float s_phi_dev_0_30;
    float s_phi_dev_1_30;
    float s_phi_dev_10_30;

    float s_eta_dev_0_30;
    float s_eta_dev_1_30;
    float s_eta_dev_10_30;

    TTree *t_ej0 = (TTree *) f_ej0->Get("EventTree");
    TTree *t_ej1 = (TTree *) f_ej1->Get("EventTree");
    TTree *t_ej10 = (TTree *) f_ej10->Get("EventTree");

    t_ej0->SetBranchAddress("sigma_norm_deviation_10", &s_ndev_0_10);
    t_ej1->SetBranchAddress("sigma_norm_deviation_10", &s_ndev_1_10);
    t_ej10->SetBranchAddress("sigma_norm_deviation_10", &s_ndev_10_10);

    t_ej0->SetBranchAddress("phi_deviation_10", &s_phi_dev_0_10);
    t_ej1->SetBranchAddress("phi_deviation_10", &s_phi_dev_1_10);
    t_ej10->SetBranchAddress("phi_deviation_10", &s_phi_dev_10_10);

    t_ej0->SetBranchAddress("eta_deviation_10", &s_eta_dev_0_10);
    t_ej1->SetBranchAddress("eta_deviation_10", &s_eta_dev_1_10);
    t_ej10->SetBranchAddress("eta_deviation_10", &s_eta_dev_10_10);

    t_ej0->SetBranchAddress("sigma_norm_deviation_20", &s_ndev_0_20);
    t_ej1->SetBranchAddress("sigma_norm_deviation_20", &s_ndev_1_20);
    t_ej10->SetBranchAddress("sigma_norm_deviation_20", &s_ndev_10_20);

    t_ej0->SetBranchAddress("phi_deviation_20", &s_phi_dev_0_20);
    t_ej1->SetBranchAddress("phi_deviation_20", &s_phi_dev_1_20);
    t_ej10->SetBranchAddress("phi_deviation_20", &s_phi_dev_10_20);

    t_ej0->SetBranchAddress("eta_deviation_20", &s_eta_dev_0_20);
    t_ej1->SetBranchAddress("eta_deviation_20", &s_eta_dev_1_20);
    t_ej10->SetBranchAddress("eta_deviation_20", &s_eta_dev_10_20);

    t_ej0->SetBranchAddress("sigma_norm_deviation_30", &s_ndev_0_30);
    t_ej1->SetBranchAddress("sigma_norm_deviation_30", &s_ndev_1_30);
    t_ej10->SetBranchAddress("sigma_norm_deviation_30", &s_ndev_10_30);

    t_ej0->SetBranchAddress("phi_deviation_30", &s_phi_dev_0_30);
    t_ej1->SetBranchAddress("phi_deviation_30", &s_phi_dev_1_30);
    t_ej10->SetBranchAddress("phi_deviation_30", &s_phi_dev_10_30);

    t_ej0->SetBranchAddress("eta_deviation_30", &s_eta_dev_0_30);
    t_ej1->SetBranchAddress("eta_deviation_30", &s_eta_dev_1_30);
    t_ej10->SetBranchAddress("eta_deviation_30", &s_eta_dev_10_30);

    Long64_t nentries = t_ej0->GetEntries();
    std::cout << nentries << std::endl;
    Int_t nbytes = 0;
    for (Long64_t i = 0; i < nentries; i++) {
        nbytes += t_ej0->GetEntry(i);
        nbytes += t_ej1->GetEntry(i);
        nbytes += t_ej10->GetEntry(i);

        hs_ndev_0_10->Fill(s_ndev_0_10);
        hs_ndev_0_20->Fill(s_ndev_0_20);
        hs_ndev_0_30->Fill(s_ndev_0_30);
        hs_ndev_1_10->Fill(s_ndev_1_10);
        hs_ndev_1_20->Fill(s_ndev_1_20);
        hs_ndev_1_30->Fill(s_ndev_1_30);
        hs_ndev_10_10->Fill(s_ndev_10_10);
        hs_ndev_10_20->Fill(s_ndev_10_20);
        hs_ndev_10_30->Fill(s_ndev_10_30);

        hs_phi_dev_0_10->Fill(s_phi_dev_0_10);
        hs_phi_dev_0_20->Fill(s_phi_dev_0_20);
        hs_phi_dev_0_30->Fill(s_phi_dev_0_30);
        hs_phi_dev_1_10->Fill(s_phi_dev_1_10);
        hs_phi_dev_1_20->Fill(s_phi_dev_1_20);
        hs_phi_dev_1_30->Fill(s_phi_dev_1_30);
        hs_phi_dev_10_10->Fill(s_phi_dev_10_10);
        hs_phi_dev_10_20->Fill(s_phi_dev_10_20);
        hs_phi_dev_10_30->Fill(s_phi_dev_10_30);

        hs_eta_dev_0_10->Fill(s_eta_dev_0_10);
        hs_eta_dev_0_20->Fill(s_eta_dev_0_20);
        hs_eta_dev_0_30->Fill(s_eta_dev_0_30);
        hs_eta_dev_1_10->Fill(s_eta_dev_1_10);
        hs_eta_dev_1_20->Fill(s_eta_dev_1_20);
        hs_eta_dev_1_30->Fill(s_eta_dev_1_30);
        hs_eta_dev_10_10->Fill(s_eta_dev_10_10);
        hs_eta_dev_10_20->Fill(s_eta_dev_10_20);
        hs_eta_dev_10_30->Fill(s_eta_dev_10_30);
    }

    TCanvas canvas("canv", "", 600, 600);
    canvas.cd();
    {
        THStack hist_stack("temporary_stack", "");
        hist_stack.Add(hs_ndev_0_10);
        hist_stack.Add(hs_ndev_1_10);
        hist_stack.Add(hs_ndev_10_10);

        prettify(hs_ndev_0_10, kBlue, 3);
        prettify(hs_ndev_1_10, kRed, 4);
        prettify(hs_ndev_10_10, kBlack, 5);

        hist_stack.Draw("nostack");

        hist_stack.GetHistogram()->
            GetXaxis()->SetTitle("#Delta#sigma");
        hist_stack.GetHistogram()->GetYaxis()->SetTitle("Normalized to Unity");
        hist_stack.GetHistogram()->GetXaxis()->SetNdivisions(505);

        hist_stack.SetMaximum(1.4 * hist_stack.GetMaximum("nostack"));
        //hist_stack.SetMaximum(0.25);
        hist_stack.Draw("nostack");

        TLegend legend(0.65, 0.75, 0.95, 0.92);
        legend.AddEntry(hs_ndev_0_10, "No Event Jet", "f");
        legend.AddEntry(hs_ndev_1_10, "#gamma = 0.001", "f");
        legend.AddEntry(hs_ndev_10_10, "#gamma = 0.01", "f");

        legend.SetTextFont(42);
        legend.SetFillStyle(0);
        legend.SetFillColor(0);
        legend.SetBorderSize(0);
        legend.SetTextSize(1.3 * legend.GetTextSize());

        legend.Draw();

        //DrawAtlasLabel("");
        canvas.Update();

        canvas.Print((out_directory + "wandering_pT_10_sigma_" + titley_bit() + ".pdf").c_str(), "pdf");
        canvas.Clear();
    }

    {
        THStack hist_stack("temporary_stack", "");
        hist_stack.Add(hs_ndev_0_20);
        hist_stack.Add(hs_ndev_1_20);
        hist_stack.Add(hs_ndev_10_20);

        prettify(hs_ndev_0_20, kBlue, 3);
        prettify(hs_ndev_1_20, kRed, 4);
        prettify(hs_ndev_10_20, kBlack, 5);

        hist_stack.Draw("nostack");

        hist_stack.GetHistogram()->
            GetXaxis()->SetTitle("#Delta#sigma");
        hist_stack.GetHistogram()->GetYaxis()->SetTitle("Normalized to Unity");
        hist_stack.GetHistogram()->GetXaxis()->SetNdivisions(505);

        hist_stack.SetMaximum(1.4 * hist_stack.GetMaximum("nostack"));
        //hist_stack.SetMaximum(0.25);
        hist_stack.Draw("nostack");

        TLegend legend(0.65, 0.75, 0.95, 0.92);
        legend.AddEntry(hs_ndev_0_20, "No Event Jet", "f");
        legend.AddEntry(hs_ndev_1_20, "#gamma = 0.001", "f");
        legend.AddEntry(hs_ndev_10_20, "#gamma = 0.01", "f");

        legend.SetTextFont(42);
        legend.SetFillStyle(0);
        legend.SetFillColor(0);
        legend.SetBorderSize(0);
        legend.SetTextSize(1.3 * legend.GetTextSize());

        legend.Draw();

        //DrawAtlasLabel("");
        canvas.Update();

        canvas.Print((out_directory + "wandering_pT_10_sigma_" + titley_bit() + ".pdf").c_str(), "pdf");
        canvas.Clear();
    }

    {
        THStack hist_stack("temporary_stack", "");
        hist_stack.Add(hs_ndev_0_30);
        hist_stack.Add(hs_ndev_1_30);
        hist_stack.Add(hs_ndev_10_30);

        prettify(hs_ndev_0_30, kBlue, 3);
        prettify(hs_ndev_1_30, kRed, 4);
        prettify(hs_ndev_10_30, kBlack, 5);

        hist_stack.Draw("nostack");

        hist_stack.GetHistogram()->
            GetXaxis()->SetTitle("#Delta#sigma");
        hist_stack.GetHistogram()->GetYaxis()->SetTitle("Normalized to Unity");
        hist_stack.GetHistogram()->GetXaxis()->SetNdivisions(505);

        hist_stack.SetMaximum(1.4 * hist_stack.GetMaximum("nostack"));
        //hist_stack.SetMaximum(0.25);
        hist_stack.Draw("nostack");

        TLegend legend(0.65, 0.75, 0.95, 0.92);
        legend.AddEntry(hs_ndev_0_30, "No Event Jet", "f");
        legend.AddEntry(hs_ndev_1_30, "#gamma = 0.001", "f");
        legend.AddEntry(hs_ndev_10_30, "#gamma = 0.01", "f");

        legend.SetTextFont(42);
        legend.SetFillStyle(0);
        legend.SetFillColor(0);
        legend.SetBorderSize(0);
        legend.SetTextSize(1.3 * legend.GetTextSize());

        legend.Draw();

        //DrawAtlasLabel("");
        canvas.Update();

        canvas.Print((out_directory + "wandering_pT_10_sigma_" + titley_bit() + ".pdf").c_str(), "pdf");
        canvas.Clear();
    }



    delete f_ej0;
    delete f_ej1;
    delete f_ej10;
}

void PileupStudy() {
    std::string base = "/u/at/chstan/nfs/summer_2014/ForConrad/";
    std::string file_loc_corr = base + "files/20k_zprime_pileup_effect_corrected/2015_03_15_17h04m17s/";
    std::string file_loc_uncorr = base + "files/20k_zprime_pileup_effect_uncorrected/2015_03_15_15h54m31s/";
    std::string out_directory = base + "results/plots/Simple/";

    std::map<size_t, TH1F *> s1m;
    std::map<size_t, TH1F *> mt1m;
    std::map<size_t, TH1F *> m1m;
    std::map<size_t, TH1F *> ns1m;
    std::map<size_t, float> bs1m;
    std::map<size_t, float> bm1m;
    std::map<size_t, float> bmt1m;

    std::map<size_t, TH1F *> s1um;
    std::map<size_t, TH1F *> mt1um;
    std::map<size_t, TH1F *> m1um;
    std::map<size_t, TH1F *> ns1um;
    std::map<size_t, float> bs1um;
    std::map<size_t, float> bm1um;
    std::map<size_t, float> bmt1um;

    TFile *f_corr = TFile::Open((file_loc_corr + "fin.root").c_str());
    TFile *f_uncorr = TFile::Open((file_loc_uncorr + "fin.root").c_str());
    TTree *t_corr = (TTree *) f_corr->Get("EventTree");
    TTree *t_uncorr = (TTree *) f_uncorr->Get("EventTree");
    size_t n_bins = 40;

    const size_t NPVs[10] = {0, 2, 5, 10, 15, 20, 25, 30, 35, 40};
    size_t npv_iter_max = sizeof(NPVs)/sizeof(NPVs[0]);

    for (size_t npv_iter = 0; npv_iter < npv_iter_max; npv_iter++) {
        std::cout << NPVs[npv_iter] << " ";
    }

    std::stringstream ss;
    for (size_t npv_iter = 0; npv_iter < npv_iter_max; npv_iter++) {
        std::cout << NPVs[npv_iter] << " ";
    }
    std::cout << std::endl;

    for (size_t npv_iter = 0; npv_iter < npv_iter_max; npv_iter++) {
        size_t NPV = NPVs[npv_iter];

        ss.str(std::string());
        ss << "hs1_" << NPV;
        s1m[NPV] = new TH1F(ss.str().c_str(), ss.str().c_str(), n_bins, 0, 1.5);

        ss.str(std::string());
        ss << "hns1_" << NPV;
        ns1m[NPV] = new TH1F(ss.str().c_str(), ss.str().c_str(), n_bins, -1.0, 1.5);

        ss.str(std::string());
        ss << "hmt1_" << NPV;
        mt1m[NPV] = new TH1F(ss.str().c_str(), ss.str().c_str(), n_bins, 0, 300);

        ss.str(std::string());
        ss << "hm1_" << NPV;
        m1m[NPV] = new TH1F(ss.str().c_str(), ss.str().c_str(), n_bins, 0, 300);

        ss.str(std::string());
        ss << "s1_" << NPV;
        bs1m[NPV] = 0;

        ss.str(std::string());
        ss << "m1_" << NPV;
        bm1m[NPV] = 0;

        ss.str(std::string());
        ss << "mt1_" << NPV;
        bmt1m[NPV] = 0;

        // uncorrected
        ss.str(std::string());
        ss << "hs1u_" << NPV;
        s1um[NPV] = new TH1F(ss.str().c_str(), ss.str().c_str(), n_bins, 0, 1.5);

        ss.str(std::string());
        ss << "hns1u_" << NPV;
        ns1um[NPV] = new TH1F(ss.str().c_str(), ss.str().c_str(), n_bins, -1.0, 1.5);

        ss.str(std::string());
        ss << "hmt1u_" << NPV;
        mt1um[NPV] = new TH1F(ss.str().c_str(), ss.str().c_str(), n_bins, 0, 300);

        ss.str(std::string());
        ss << "hm1u_" << NPV;
        m1um[NPV] = new TH1F(ss.str().c_str(), ss.str().c_str(), n_bins, 0, 300);

        ss.str(std::string());
        ss << "s1u_" << NPV;
        bs1um[NPV] = 0;

        ss.str(std::string());
        ss << "m1u_" << NPV;
        bm1um[NPV] = 0;

        ss.str(std::string());
        ss << "mt1u_" << NPV;
        bmt1um[NPV] = 0;
    }

    for (size_t npv_iter = 0; npv_iter < npv_iter_max; npv_iter++) {
        size_t NPV = NPVs[npv_iter];
        std::cout << NPV << std::endl;
        ss.str(std::string());
        ss << "sigma1_" << NPV;
        t_corr->SetBranchAddress(ss.str().c_str(), &bs1m[NPV]);

        ss.str(std::string());
        ss << "mass1_" << NPV;
        t_corr->SetBranchAddress(ss.str().c_str(), &bm1m[NPV]);

        ss.str(std::string());
        ss << "mass_trimmed1_" << NPV;
        t_corr->SetBranchAddress(ss.str().c_str(), &bmt1m[NPV]);

        // uncorrected
        ss.str(std::string());
        ss << "sigma1_" << NPV;
        t_uncorr->SetBranchAddress(ss.str().c_str(), &bs1um[NPV]);

        ss.str(std::string());
        ss << "mass1_" << NPV;
        t_uncorr->SetBranchAddress(ss.str().c_str(), &bm1um[NPV]);

        ss.str(std::string());
        ss << "mass_trimmed1_" << NPV;
        t_uncorr->SetBranchAddress(ss.str().c_str(), &bmt1um[NPV]);
    }

    Long64_t nentries = t_corr->GetEntries();
    std::cout << nentries << std::endl;
    Int_t nbytes = 0;
    for (Long64_t i=0; i<nentries;i++) {
        nbytes += t_corr->GetEntry(i);
        nbytes += t_uncorr->GetEntry(i);

        for (size_t npv_iter = 0; npv_iter < npv_iter_max; npv_iter++) {
            size_t NPV = NPVs[npv_iter];

            s1m[NPV]->Fill(bs1m[NPV]);
            m1m[NPV]->Fill(bm1m[NPV]);
            mt1m[NPV]->Fill(bmt1m[NPV]);
            s1um[NPV]->Fill(bs1um[NPV]);
            m1um[NPV]->Fill(bm1um[NPV]);
            mt1um[NPV]->Fill(bmt1um[NPV]);

            ns1m[NPV]->Fill((bs1m[NPV] - bs1m[0])/bs1m[0]);
            ns1um[NPV]->Fill((bs1um[NPV] - bs1um[0])/bs1um[0]);
        }
    }

    TCanvas canvas("canv", "", 600, 600);
    canvas.cd();

    std::vector<float> uncorrected_w;
    std::vector<float> corrected_w;
    std::vector<float> xs;

    for (size_t npv_iter = 0; npv_iter < npv_iter_max; npv_iter++) {
        size_t NPV = NPVs[npv_iter];
        if (NPV != 0) {
            xs.push_back(static_cast<float>(NPV));
            uncorrected_w.push_back(s1um[NPV]->GetRMS() / s1um[NPV]->GetMean());
            corrected_w.push_back(s1m[NPV]->GetRMS() / s1m[NPV]->GetMean());
        }
        THStack hist_stack("temporary_stack", "");

        hist_stack.Add(ns1m[NPV]);
        hist_stack.Add(ns1um[NPV]);

        prettify(ns1m[NPV], kBlue, 3);
        prettify(ns1um[NPV], kRed, 4);

        hist_stack.Draw("nostack");
        hist_stack.GetHistogram()->GetXaxis()->SetTitle(gen_x_label(NPV).c_str());
        hist_stack.GetHistogram()->GetYaxis()->SetTitle("Normalized to Unity");
        hist_stack.GetHistogram()->GetXaxis()->SetNdivisions(505);

        hist_stack.SetMaximum(0.35);
        hist_stack.Draw("nostack");

        TLegend legend(0.65, 0.75, 0.95, 0.92);
        ss.str(std::string());
        ss << "Corrected at NPV = " << NPV;
        legend.AddEntry(ns1m[NPV], ss.str().c_str(), "f");
        ss.str(std::string());
        ss << "Uncorrected at NPV = " << NPV;
        legend.AddEntry(ns1um[NPV], ss.str().c_str(), "f");

        legend.SetTextFont(42);
        legend.SetFillStyle(0);
        legend.SetFillColor(0);
        legend.SetBorderSize(0);
        legend.SetTextSize(1.3 * legend.GetTextSize());

        legend.Draw();

        //DrawAtlasLabel("");
        canvas.Update();

        ss.str(std::string());
        ss << out_directory << "ns1_" << NPV << "_comp_" << titley_bit() << ".pdf";
        canvas.Print(ss.str().c_str(), "pdf");
        canvas.Clear();
    }

    {
        TMultiGraph mg;
        TGraph *g_corrected = new TGraph(xs.size(), &xs.at(0), &corrected_w.at(0));
        g_corrected->SetMarkerColor(kBlue);
        g_corrected->SetMarkerStyle(22);
        TGraph *g_uncorrected = new TGraph(xs.size(), &xs.at(0), &uncorrected_w.at(0));
        g_uncorrected->SetMarkerColor(kRed);
        g_uncorrected->SetMarkerStyle(23);

        mg.Add(g_corrected, "p");
        mg.Add(g_uncorrected, "p");

        mg.Draw("ap");
        mg.GetXaxis()->SetTitle("Truth Primary Vertex Count");
        mg.GetYaxis()->SetTitle("#sigma_{#sigma} / #mu_{#sigma}");

        TLegend legend(0.2, 0.25, 0.45, 0.35);
        legend.SetTextFont(42);
        legend.SetFillStyle(0);
        legend.SetFillColor(0);
        legend.SetBorderSize(0);
        legend.SetTextSize(1.3 * legend.GetTextSize());
        legend.AddEntry(g_corrected, "Corrected", "p");
        legend.AddEntry(g_uncorrected, "Uncorrected", "p");

        canvas.Update();
        legend.Draw();

        ss.str(std::string());
        ss << out_directory << "ns1_wom" << ".pdf";
        canvas.Print(ss.str().c_str(), "pdf");
        canvas.Clear();
    }

    // MASS
    //{
    //    THStack hist_stack("temporary_stack", "");
    //    hist_stack.Add(hmt1_80);
    //    hist_stack.Add(hm1_80);
    //
    //    prettify(hmt1_80, kBlue, 3);
    //    prettify(hm1_80, kRed, 4);
    //
    //    hist_stack.Draw("nostack");
    //    hist_stack.GetHistogram()->GetXaxis()->SetTitle(gen_x_mass_label(80).c_str());
    //    hist_stack.GetHistogram()->GetYaxis()->SetTitle("Normalized to Unity");
    //    hist_stack.GetHistogram()->GetXaxis()->SetNdivisions(505);
    //
    //    hist_stack.SetMaximum(1.4 * hist_stack.GetMaximum("nostack"));
    //
    //    hist_stack.Draw("nostack");
    //
    //    TLegend legend(0.65, 0.75, 0.95, 0.92);
    //    legend.AddEntry(hmt1_80, "R=0.3 Trimmed Anti-k_{t} Mass", "f");
    //    legend.AddEntry(hm1_80, "Untrimmed Anti-k_{t} Mass", "f");
    //
    //    legend.SetTextFont(42);
    //    legend.SetFillStyle(0);
    //    legend.SetFillColor(0);
    //    legend.SetBorderSize(0);
    //    legend.SetTextSize(1.3 * legend.GetTextSize());
    //
    //    legend.Draw();
    //
    //    //DrawAtlasLabel("");
    //    canvas.Update();
    //
    //    canvas.Print((out_directory+"trimming_80_comp_" + titley_bit() + ".pdf").c_str(), "pdf");
    //    canvas.Clear();
    //}

    for(size_t npv_iter = 0; npv_iter < npv_iter_max; npv_iter++) {
        size_t NPV = NPVs[npv_iter];
        std::cout << NPV << ": " << s1m[NPV]->GetMean() << std::endl;
        std::cout << NPV << ": " << s1um[NPV]->GetMean() << std::endl << std::endl;

    }

    delete f_corr;
    delete f_uncorr;

}


void Simple() {
    SetupATLASStyle();

    bool do_noise_study = false;
    bool do_wandering_study = false;
    bool do_pileup_study = true;

    if (do_noise_study) {
        NoiseStudy();
    }
    if (do_wandering_study) {
        WanderingStudy();
    }
    if (do_pileup_study) {
        PileupStudy();
    }
    return;



}
