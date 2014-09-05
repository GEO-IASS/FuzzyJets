#include <iostream>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <assert.h>

#include <TColor.h>
#include <THStack.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TLegend.h>

#include "Histogram.h"
#include "Event.h"
#include "AtlasUtils.h"
#include "AnalyzeFuzzyTools.h"

void gen_string(char *s, const unsigned int len) {
    static const char alphabet[] = 
        "0123456789"
        "abcdefghijklmnopqrstuvwxyz"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    for (unsigned int iter = 0; iter < len; iter++) {
        s[iter] = alphabet[rand() % (sizeof(alphabet) - 1)];
    }
    s[len] = 0;
}

void StackedEfficiencyHistogramBase::Finish(__attribute__((unused)) EventManager const* event_manager) {
    SetupATLASStyle();

    // scale histograms
    TMultiGraph *multi = new TMultiGraph();

    std::vector<float> random_x = {0, 1};
    std::vector<float> random_y = {1, 0};
    TGraph *random = new TGraph(2, &random_x.at(0), &random_y.at(0));
    random->SetLineColor(kGray);
    multi->Add(random, "c");

    TLegend *legend = new TLegend(0.2, 0.25, 0.45, 0.35);
    legend->SetTextFont(42);
    legend->SetFillStyle(0);
    legend->SetFillColor(0);
    legend->SetBorderSize(0);

    for(unsigned int hist_iter = 0; hist_iter < _signal_hists.size(); hist_iter++) {
        TH1F *signal_hist = _signal_hists[hist_iter];
        TH1F *background_hist = _background_hists[hist_iter];
        signal_hist->Scale(1./signal_hist->Integral(-1, signal_hist->GetNbinsX()+1));
        background_hist->Scale(1./background_hist->Integral(-1, background_hist->GetNbinsX()+1));

        float integrated_signal = 0;
        float integrated_background = 0;
        std::vector<float> vec_signal;
        std::vector<float> vec_background_rejection;
        std::vector<CustomSortPair> to_sort;
        for (unsigned int bin_iter = 0; bin_iter < _n_points; bin_iter++) {
            to_sort.push_back(CustomSortPair(signal_hist->GetBinContent(bin_iter),
                                             background_hist->GetBinContent(bin_iter)));

        }
        std::sort(to_sort.begin(), to_sort.end());
        for (unsigned int iter = 0; iter < _n_points; iter++) {
            integrated_signal += to_sort.at(iter).signal;
            integrated_background += to_sort.at(iter).background;
            vec_signal.push_back(integrated_signal);
            vec_background_rejection.push_back(1 - integrated_background);
        }
        TGraph *efficiency_graph = new TGraph(_n_points, &vec_signal.at(0), &vec_background_rejection.at(0));

        efficiency_graph->SetMarkerColor(_colors[hist_iter]);
        efficiency_graph->SetLineColor(_colors[hist_iter]);
        efficiency_graph->SetMarkerStyle(20); // circle

        legend->AddEntry(efficiency_graph, _labels[hist_iter].c_str(), "l");
        multi->Add(efficiency_graph);
    }

    TCanvas canvas("temporary", "", 0, 0, _canvas_x, _canvas_y);
    multi->Draw("ac");
    multi->GetXaxis()->SetTitle(_x_label.c_str());
    multi->GetYaxis()->SetTitle(_y_label.c_str());

    multi->GetXaxis()->SetLimits(0, 1);
    multi->GetYaxis()->SetLimits(0, 1);

    multi->GetXaxis()->SetNdivisions(_ticks);

    canvas.Update();

    legend->Draw();

    DrawAtlasLabel(_title, 0.2, 0.48);

    std::stringstream ss;
    ss << plot_prefix << _outfile_name;
    std::string out = ss.str();
    canvas.Print(out.c_str());

    // multi owns the efficiency graphs
    delete multi;
}

void StackedHistogramBase::Finish(__attribute__((unused)) EventManager const* event_manager) {
    SetupATLASStyle();

    TCanvas canvas("temporary", "", 0, 0, _canvas_x, _canvas_y);
    THStack hist_stack("temporary_stack", "");

    for (unsigned int hist_iter = 0; hist_iter < _hists.size(); hist_iter++) {
        TH1F *current_hist = _hists[hist_iter];
        hist_stack.Add(current_hist);
        
        // do other things like set the color
        current_hist->SetLineColor(_colors[hist_iter]);
        current_hist->SetLineStyle(_styles[hist_iter]);
        current_hist->SetLineWidth(2);
        
        assert(current_hist->Integral(-1, current_hist->GetNbinsX() + 1) > 0);
        if (_normalized) {
            current_hist->Scale(1./current_hist->Integral(-1, current_hist->GetNbinsX()+1));
        }
    }


    hist_stack.Draw("nostack");

    hist_stack.GetHistogram()->GetXaxis()->SetTitle(_x_label.c_str());
    hist_stack.GetHistogram()->GetYaxis()->SetTitle(_y_label.c_str());
    hist_stack.GetHistogram()->GetXaxis()->SetNdivisions(_ticks);

    canvas.Update();
    TLegend *legend = new TLegend(0.65, 0.75, 0.9, 0.85);
    for (unsigned int hist_iter = 0; hist_iter < _hists.size(); hist_iter++) {
        legend->AddEntry(_hists[hist_iter], _labels[hist_iter].c_str(), "l");
    }

    legend->SetTextFont(42);
    legend->SetFillStyle(0);
    legend->SetFillColor(0);
    legend->SetBorderSize(0);

    legend->Draw();

    DrawAtlasLabel(_title);
    //myText(0.18, 0.955, kBlack, _title.c_str());

    std::stringstream ss;
    ss << plot_prefix << _outfile_name;
    std::string out = ss.str();
    canvas.Print(out.c_str());

    delete legend;
}

void CorrelationBase::Finish(__attribute__((unused)) EventManager const* event_manager) {
    SetupATLASStyle();

    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    _hist->SetFillColor(46);

    TCanvas canvas("temporary", "", 0, 0, _canvas_x, _canvas_y);

    _hist->SetLabelSize(0.04, "X");
    _hist->SetLabelSize(0.04, "Y");
    _hist->SetXTitle(_x_label.c_str());
    _hist->SetYTitle(_y_label.c_str());

    _hist->Draw("colz");


    std::stringstream ss;
    ss.str(std::string());
    ss << _title;
    if (_correlation_in_title) {
        ss << "      Correlation = " << std::setprecision(2) << std::fixed << _hist->GetCorrelationFactor();
    }
    std::string canvas_title = ss.str();
    myText(0.18, 0.955, kBlack, canvas_title.c_str());

    
    ss.str(std::string());
    ss << plot_prefix << _outfile_name;
    std::string out = ss.str();
    canvas.Print(out.c_str());
}

void SigmaJetSizeCorrelationPoster::Update(EventManager const* event_manager) {
    float reweight = event_manager->Reweight("zprime");
    _hist->Fill(event_manager->_zprime_event.antikt_m/event_manager->_zprime_event.antikt_pt,
                event_manager->_zprime_event.mGMMc_r,
                reweight);
}

void WeightDistanceCorrelation::Update(EventManager const* event_manager) {
    std::stringstream ss;

    ss.str(std::string());
    ss << _alg_label << "_weight_vec";
    std::vector<float> *alg_weight_vec = (*event_manager)[_event_label].Get<std::vector<float> *>(ss.str());

    ss.str(std::string());
    ss << _alg_label << "_distance_vec";
    std::vector<float> *alg_distance_vec = (*event_manager)[_event_label].Get<std::vector<float> *>(ss.str());

    //ss.str(std::string());
    //ss << "is_lead_antikt_constituent_vec";
    //std::vector<bool> *is_constituent_vec = (*event_manager)[_event_label].Get<std::vector<bool> *>(ss.str());

    float reweight = event_manager->Reweight(_event_label);
    for (unsigned int particle_iter = 0; particle_iter < alg_weight_vec->size(); particle_iter++) {
        _hist->Fill(alg_weight_vec->at(particle_iter), alg_distance_vec->at(particle_iter), reweight);
    }
}

void PtCorrelation::Update(EventManager const* event_manager) {
    std::stringstream ss;

    ss.str(std::string());
    ss << _alg_label << "_pt";
    Float_t alg_val = (*event_manager)[_event_label].Get<Float_t>(ss.str());
    
    ss.str(std::string());
    ss << _other_alg_label << "_pt";
    Float_t other_alg_val = (*event_manager)[_event_label].Get<Float_t>(ss.str());

    float reweight = event_manager->Reweight(_event_label);
    _hist->Fill(alg_val, other_alg_val, reweight);
}

void SigmaEfficiencyPosterPlot::Update(EventManager const* event_manager) {
    float qcd_reweight = event_manager->Reweight("qcd");
    float wprime_reweight = event_manager->Reweight("wprime");

    _signal_hists.at(0)->Fill(event_manager->_wprime_event.mGMMc_r, wprime_reweight);
    _background_hists.at(0)->Fill(event_manager->_qcd_event.mGMMc_r, qcd_reweight);

    _signal_hists.at(1)->Fill(event_manager->_wprime_event.mGMMc_m, wprime_reweight);
    _background_hists.at(1)->Fill(event_manager->_qcd_event.mGMMc_m, qcd_reweight);

    _signal_hists.at(2)->Fill(event_manager->_wprime_event.antikt_m, wprime_reweight);
    _background_hists.at(2)->Fill(event_manager->_qcd_event.antikt_m, qcd_reweight);
}

void SigmaEfficiencyPlot::Update(EventManager const* event_manager) {
    float signal_reweight = event_manager->Reweight(_event_signal);
    float background_reweight = event_manager->Reweight(_event_background);

    float pT_temp = (*event_manager)[_event_signal].antikt_pt;
    if (_cut_low <= pT_temp && pT_temp < _cut_high)
        _signal_hists.at(0)->Fill((*event_manager)[_event_signal].mGMMc_r, signal_reweight);

    pT_temp = (*event_manager)[_event_background].antikt_pt;
    if (_cut_low <= pT_temp && pT_temp < _cut_high)
        _background_hists.at(0)->Fill((*event_manager)[_event_background].mGMMc_r, background_reweight);

    pT_temp = (*event_manager)[_event_signal].antikt_pt;
    if (_cut_low <= pT_temp && pT_temp < _cut_high)
        _signal_hists.at(1)->Fill((*event_manager)[_event_signal].antikt_m, signal_reweight);

    pT_temp = (*event_manager)[_event_background].antikt_pt;
    if (_cut_low <= pT_temp && pT_temp < _cut_high)
        _background_hists.at(1)->Fill((*event_manager)[_event_background].antikt_m, background_reweight);
}

void SkewEfficiencyPlot::Update(EventManager const* event_manager) {
    float qcd_reweight = event_manager->Reweight("qcd");
    float zprime_reweight = event_manager->Reweight("zprime");
    
    _signal_hists.at(0)->Fill(event_manager->_zprime_event.mGMMc_m_skew, zprime_reweight);
    _background_hists.at(0)->Fill(event_manager->_qcd_event.mGMMc_m_skew, qcd_reweight);
}

void FuzzyJetMassEfficiencyPlot::Update(EventManager const* event_manager) {
    float signal_reweight = event_manager->Reweight(_event_signal);
    float background_reweight = event_manager->Reweight(_event_background);

    std::stringstream ss;
    ss << _alg << "_m";

    float pT_temp = (*event_manager)[_event_signal].antikt_pt;
    if (_cut_low <= pT_temp && pT_temp < _cut_high) {
        Float_t mass = (*event_manager)[_event_signal].Get<Float_t>(ss.str());
        _signal_hists.at(0)->Fill(mass, signal_reweight);
    }

    pT_temp = (*event_manager)[_event_background].antikt_pt;
    if (_cut_low <= pT_temp && pT_temp < _cut_high) {
        Float_t mass = (*event_manager)[_event_background].Get<Float_t>(ss.str());
        _signal_hists.at(0)->Fill(mass, background_reweight);
    }

    pT_temp = (*event_manager)[_event_signal].antikt_pt;
    if (_cut_low <= pT_temp && pT_temp < _cut_high)
        _signal_hists.at(1)->Fill((*event_manager)[_event_signal].antikt_m, signal_reweight);

    pT_temp = (*event_manager)[_event_background].antikt_pt;
    if (_cut_low <= pT_temp && pT_temp < _cut_high)
        _background_hists.at(1)->Fill((*event_manager)[_event_background].antikt_m, background_reweight);
}

void MassCorrelation::Update(EventManager const* event_manager) {
    std::stringstream ss;

    ss.str(std::string());
    ss << _alg_label << "_m";
    Float_t alg_val = (*event_manager)[_event_label].Get<Float_t>(ss.str());
    
    ss.str(std::string());
    ss << _other_alg_label << "_m";
    Float_t other_alg_val = (*event_manager)[_event_label].Get<Float_t>(ss.str());

    float reweight = event_manager->Reweight(_event_label);
    _hist->Fill(alg_val, other_alg_val, reweight);
}

void SigmaJetSizeCorrelation::Update(EventManager const* event_manager) {
    std::stringstream ss;
    
    ss.str(std::string());
    ss << _alg_label << "_m";
    Float_t mass = (*event_manager)[_event_label].Get<Float_t>(ss.str());

    ss.str(std::string());
    ss << _alg_label << "_pt";
    Float_t pT = (*event_manager)[_event_label].Get<Float_t>(ss.str());

    Float_t sigma = (*event_manager)[_event_label].Get<Float_t>("mGMMc_r");

    float reweight = event_manager->Reweight(_event_label);
    _hist->Fill(mass/pT, sigma, reweight);
}

void SkewHistogram::Update(EventManager const* event_manager) {
    float zprime_reweight = event_manager->Reweight("zprime");
    float qcd_reweight = event_manager->Reweight("qcd");
    _hists[0]->Fill(event_manager->_zprime_event.mGMMc_m_skew, zprime_reweight);
    _hists[1]->Fill(event_manager->_qcd_event.mGMMc_m_skew, qcd_reweight);
}

void DeltaRHistogram::Update(EventManager const* event_manager) {
    std::stringstream ss;
    
    ss.str(std::string());
    ss << _alg_label << "_dr";
    Float_t zprime_dr = (*event_manager)["zprime"].Get<Float_t>(ss.str());
    Float_t wprime_dr = (*event_manager)["wprime"].Get<Float_t>(ss.str());
    Float_t qcd_dr = (*event_manager)["qcd"].Get<Float_t>(ss.str());

    float zprime_reweight = event_manager->Reweight("zprime");
    float wprime_reweight = event_manager->Reweight("wprime");
    float qcd_reweight = event_manager->Reweight("qcd");

    _hists[0]->Fill(zprime_dr, zprime_reweight);
    _hists[1]->Fill(wprime_dr, wprime_reweight);
    _hists[2]->Fill(qcd_dr, qcd_reweight);
}

void FuzzyAntiktPtCorrelation::Update(EventManager const* event_manager) {
    float reweight = event_manager->Reweight("zprime");
    _hist->Fill(event_manager->_zprime_event.antikt_pt, event_manager->_zprime_event.mGMMc_pt, reweight);
}

void RadiusComparisonHistogram::Update(EventManager const* event_manager) {
    float zprime_reweight = event_manager->Reweight("zprime");
    float wprime_reweight = event_manager->Reweight("wprime");
    float qcd_reweight = event_manager->Reweight("qcd");    

    _hists[0]->Fill(event_manager->_zprime_event.mGMMc_r, zprime_reweight);
    _hists[1]->Fill(event_manager->_wprime_event.mGMMc_r, wprime_reweight);
    _hists[2]->Fill(event_manager->_qcd_event.mGMMc_r, qcd_reweight);
}

void AverageRadiusComparisonHistogram::Update(EventManager const* event_manager) {
    float zprime_reweight = event_manager->Reweight("zprime");
    float wprime_reweight = event_manager->Reweight("wprime");
    float qcd_reweight = event_manager->Reweight("qcd");    

    _hists[0]->Fill(event_manager->_zprime_event.mGMMc_r_avg, zprime_reweight);
    _hists[1]->Fill(event_manager->_wprime_event.mGMMc_r_avg, wprime_reweight);
    _hists[2]->Fill(event_manager->_qcd_event.mGMMc_r_avg, qcd_reweight);
}
