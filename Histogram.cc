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

void StackedEfficiencyHistogramGen::Finish(__attribute__((unused)) EventManager
                                           const* event_manager) {
    SetupATLASStyle();

    TMultiGraph *multi = new TMultiGraph();
    TMultiGraph *inv_multi = new TMultiGraph();

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

    TLegend *inv_legend = new TLegend(0.2, 0.25, 0.45, 0.35);
    inv_legend->SetTextFont(42);
    inv_legend->SetFillStyle(0);
    inv_legend->SetFillColor(0);
    inv_legend->SetBorderSize(0);

    for (unsigned int attr_iter = 0; attr_iter < _colors.size(); attr_iter++) {
        float integral_sig = std::accumulate(_signals.at(attr_iter).begin(),
                                             _signals.at(attr_iter).end(),
                                             0.);
        float integral_background =
            std::accumulate(_backgrounds.at(attr_iter).begin(),
                            _backgrounds.at(attr_iter).end(),
                            0.);

        std::vector<CustomSortPair> to_sort;
        CustomSortPair pair(0,0);
        for (unsigned int iter = 0; iter < _signals.at(attr_iter).size(); iter++) {
            float p_sig = _signals.at(attr_iter).at(iter) / integral_sig;
            float p_back = _backgrounds.at(attr_iter).at(iter) / integral_background;
            pair.signal += p_sig;
            pair.background += p_back;
            if (pair.background > 0.005) {
                to_sort.push_back(CustomSortPair(pair.signal, pair.background));
                pair.signal = 0;
                pair.background = 0;
            }
        }
        if (pair.signal + pair.background > 0) {
            to_sort.push_back(CustomSortPair(pair.signal, pair.background));
        }
        std::sort(to_sort.begin(), to_sort.end());

        float integrated_signal = 0;
        float integrated_background = 0;
        std::vector<float> vec_signal;
        std::vector<float> vec_background_rejection;
        std::vector<float> vec_inv_background_efficiency;

        vec_signal.push_back(0);
        vec_background_rejection.push_back(1);
        for (unsigned int iter = 0; iter < to_sort.size(); iter++) {
            integrated_signal += to_sort.at(iter).signal;
            integrated_background += to_sort.at(iter).background;
            vec_signal.push_back(integrated_signal);
            vec_background_rejection.push_back(1 - integrated_background);
            if (integrated_background > 0)
                vec_inv_background_efficiency.push_back(1/integrated_background);
        }
        vec_signal.push_back(1);
        vec_background_rejection.push_back(0);
        vec_inv_background_efficiency.push_back(1);

        TGraph *efficiency_graph = new TGraph(vec_signal.size(),
                                              &vec_signal.at(0),
                                              &vec_background_rejection.at(0));
        TGraph *inv_efficiency_graph =
            new TGraph(vec_inv_background_efficiency.size(),
                       &vec_signal.at(vec_signal.size() -
                                      vec_inv_background_efficiency.size()),
                       &vec_inv_background_efficiency.at(0));

        efficiency_graph->SetMarkerColor(_colors[attr_iter]);
        efficiency_graph->SetLineColor(_colors[attr_iter]);
        efficiency_graph->SetMarkerStyle(20); // circle

        inv_efficiency_graph->SetMarkerColor(_colors[attr_iter]);
        inv_efficiency_graph->SetLineColor(_colors[attr_iter]);
        inv_efficiency_graph->SetMarkerStyle(20); // circle

        legend->AddEntry(efficiency_graph, _labels[attr_iter].c_str(), "l");
        inv_legend->AddEntry(inv_efficiency_graph,
                             _labels[attr_iter].c_str(), "l");

        multi->Add(efficiency_graph);
        inv_multi->Add(inv_efficiency_graph);
    }
    TCanvas canvas("temporary", "", 0, 0, _canvas_x, _canvas_y);
    canvas.cd();

    multi->Draw("ac");
    multi->GetXaxis()->SetTitle(_x_label.c_str());
    multi->GetYaxis()->SetTitle(_y_label.c_str());

    multi->GetXaxis()->SetLimits(0, 1);
    multi->GetYaxis()->SetLimits(0, 1);

    multi->GetXaxis()->SetNdivisions(_ticks);

    canvas.Update();

    legend->Draw();

    DrawAtlasLabel(_title, 0.2, 0.48);

    std::string out = plot_prefix + _outfile_name;
    canvas.Print(out.c_str());

    // multi owns the efficiency graphs
    delete multi;

    TCanvas inv_canvas("temporary2", "", 0,0, _canvas_x, _canvas_y);
    inv_canvas.cd();
    inv_canvas.SetLogy();

    inv_multi->Draw("ac");
    inv_multi->GetXaxis()->SetTitle(_x_label.c_str());
    inv_multi->GetYaxis()->SetTitle(_y_inv_label.c_str());

    inv_multi->GetXaxis()->SetLimits(0, 1);
    inv_multi->GetYaxis()->SetLimits(0, 1);

    inv_multi->GetXaxis()->SetNdivisions(_ticks);

    inv_canvas.Update();

    inv_legend->Draw();

    DrawAtlasLabel(_title, 0.2, 0.48);

    out = plot_prefix + "inv_" + _outfile_name;
    inv_canvas.Print(out.c_str());

    delete inv_multi;
}

void StackedEfficiencyHistogramBase::Finish(__attribute__((unused)) EventManager
                                            const* event_manager) {
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

        float norm = signal_hist->Integral(-1, signal_hist->GetNbinsX() + 1);
        signal_hist->Scale(1./norm);
        norm = background_hist->Integral(-1, background_hist->GetNbinsX() + 1);
        background_hist->Scale(1./norm);

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
        TGraph *efficiency_graph = new TGraph(_n_points, &vec_signal.at(0),
                                              &vec_background_rejection.at(0));

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

void StackedHistogramBase::Finish(__attribute__((unused)) EventManager
                                  const* event_manager) {
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
            float norm = current_hist->Integral(-1, current_hist->GetNbinsX() + 1);
            current_hist->Scale(1./norm);
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

void ScatterBase::Finish(__attribute__((unused)) EventManager const* event_manager) {
    SetupATLASStyle();

    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);

    TCanvas canvas("temporary", "", 0, 0, _canvas_x, _canvas_y);

    TMultiGraph *multi = new TMultiGraph();
    TGraph *g = new TGraph(_xs.size(), &_xs.at(0), &_ys.at(0));
    g->SetMarkerColor(_color);
    g->SetMarkerStyle(_style);
    multi->Add(g);

    TLegend *legend = new TLegend(0.65, 0.75, 0.9, 0.85);
    legend->SetTextFont(42);
    legend->SetFillStyle(0);
    legend->SetFillColor(0);
    legend->SetBorderSize(0);

    multi->Draw("ap");
    multi->GetXaxis()->SetTitle(_x_label.c_str());
    multi->GetYaxis()->SetTitle(_y_label.c_str());

    canvas.Update();
    //legend->Draw();

    //DrawAtlasLabel(_title, 0.2, 0.48);

    std::stringstream ss;
    ss << plot_prefix << _outfile_name;
    std::string out = ss.str();
    canvas.Print(out.c_str());

    delete legend;
    delete multi;
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
        float correlation = _hist->GetCorrelationFactor();
        ss << "      Correlation = " << std::setprecision(2)
           << std::fixed << correlation;
    }
    std::string canvas_title = ss.str();
    myText(0.18, 0.955, kBlack, canvas_title.c_str());


    ss.str(std::string());
    ss << plot_prefix << _outfile_name;
    std::string out = ss.str();
    canvas.Print(out.c_str());
}

void SigmaJetSizeCorrelationPoster::Update(EventManager const* event_manager) {
    float reweight = event_manager->Reweight("zprime_5");
    _hist->Fill((*event_manager)["zprime_5"].antikt_m/
                (*event_manager)["zprime_5"].antikt_pt,

                (*event_manager)["zprime_5"].mGMMc_r,
                reweight);
}

void WeightDistanceCorrelation::Update(EventManager const* event_manager) {
    std::stringstream ss;

    ss.str(std::string());
    ss << _alg_label << "_weight_vec";
    std::vector<float> *alg_weight_vec = (*event_manager)[_event_label].
        Get<std::vector<float> *>(ss.str());

    ss.str(std::string());
    ss << _alg_label << "_distance_vec";
    std::vector<float> *alg_distance_vec = (*event_manager)[_event_label].
        Get<std::vector<float> *>(ss.str());

    float reweight = event_manager->Reweight(_event_label);
    for (unsigned int particle_iter = 0; particle_iter < alg_weight_vec->size();
         particle_iter++) {

        _hist->Fill(alg_weight_vec->at(particle_iter),
                    alg_distance_vec->at(particle_iter),
                    reweight);
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

void SigmaImprovementEfficiencyMultiTau::Update(EventManager const* event_manager) {
    float background_reweight = event_manager->Reweight(_event_background);
    float signal_reweight = event_manager->Reweight(_event_signal);

    float pT_sig = (*event_manager)[_event_signal].antikt_pt;
    float pT_back = (*event_manager)[_event_background].antikt_pt;

    std::vector<float> *n_sub_sig = (*event_manager)[_event_signal].antikt_nsubjettiness;
    std::vector<float> *n_sub_back = (*event_manager)[_event_background].antikt_nsubjettiness;
    float tau1_sig;
    float tau2_sig;
    float tau3_sig;

    float tau1_back;
    float tau2_back;
    float tau3_back;

    if(n_sub_sig->size() > 0)
        tau1_sig = n_sub_sig->at(0);
    if(n_sub_back->size() > 0)
        tau1_back = n_sub_back->at(0);
    if(n_sub_sig->size() > 1)
        tau2_sig = n_sub_sig->at(1);
    if(n_sub_back->size() > 1)
        tau2_back = n_sub_back->at(1);
    if(n_sub_sig->size() > 2)
        tau3_sig = n_sub_sig->at(2);
    if(n_sub_back->size() > 2)
        tau3_back = n_sub_back->at(2);

    int max_sig = n_sub_sig->size();
    int max_back = n_sub_back->size();

    float sigma_sig = (*event_manager)[_event_signal].mGMMc_r;
    float sigma_back = (*event_manager)[_event_background].mGMMc_r;

    if (_cut_low <= pT_sig && pT_sig < _cut_high) {
        if (max_sig > 2) {
            fillSignal(0, {tau3_sig, tau2_sig}, signal_reweight);
            fillSignal(3, {tau3_sig, sigma_sig}, signal_reweight);
        }
        if (max_sig > 1) {
            fillSignal(1, {tau2_sig, tau1_sig}, signal_reweight);
            fillSignal(4, {tau2_sig, sigma_sig}, signal_reweight);
        }
        if (max_sig > 0) {
            fillSignal(5, {tau1_sig, sigma_sig}, signal_reweight);
        }
        fillSignal(2, {sigma_sig}, signal_reweight);
    }
    if (_cut_low <= pT_back && pT_back < _cut_high) {
        if (max_back > 2) {
            fillBackground(0, {tau3_back, tau2_back}, background_reweight);
            fillBackground(3, {tau3_back, sigma_back}, background_reweight);
        }
        if (max_back > 1) {
            fillBackground(1, {tau2_back, tau1_back}, background_reweight);
            fillBackground(4, {tau2_back, sigma_back}, background_reweight);
        }
        if (max_back > 0) {
            fillBackground(5, {tau1_back, sigma_back}, background_reweight);
        }

        fillBackground(2, {sigma_back}, background_reweight);
    }
}

void SigmaImprovementEfficiencyTau21::Update(EventManager const* event_manager) {
    float background_reweight = event_manager->Reweight(_event_background);
    float signal_reweight = event_manager->Reweight(_event_signal);

    float pT_sig = (*event_manager)[_event_signal].antikt_pt;
    float pT_back = (*event_manager)[_event_background].antikt_pt;

    std::vector<float> *n_sub_sig = (*event_manager)[_event_signal].antikt_nsubjettiness;
    std::vector<float> *n_sub_back = (*event_manager)[_event_background].antikt_nsubjettiness;
    float tau1_sig;
    float tau2_sig;
    float tau3_sig;
    float tau1_back;
    float tau2_back;
    float tau3_back;
    if(n_sub_sig->size() > 0)
        tau1_sig = n_sub_sig->at(0);
    if(n_sub_back->size() > 0)
        tau1_back = n_sub_back->at(0);
    if(n_sub_sig->size() > 1)
        tau2_sig = n_sub_sig->at(1);
    if(n_sub_back->size() > 1)
        tau2_back = n_sub_back->at(1);
    if(n_sub_sig->size() > 2)
        tau3_sig = n_sub_sig->at(2);
    if(n_sub_back->size() > 2)
        tau3_back = n_sub_back->at(2);

    int max_sig = n_sub_sig->size();
    int max_back = n_sub_back->size();

    float sigma_sig = (*event_manager)[_event_signal].mGMMc_r;
    float sigma_back = (*event_manager)[_event_background].mGMMc_r;

    if (_cut_low <= pT_sig && pT_sig < _cut_high) {
        if (max_sig > 2) {
            fillSignal(0, {tau3_sig/tau2_sig}, signal_reweight);
        }
        if (max_sig > 1) {
            fillSignal(1, {tau2_sig/tau1_sig}, signal_reweight);
            fillSignal(3, {tau2_sig/tau1_sig, sigma_sig}, signal_reweight);
        }
        fillSignal(2, {sigma_sig}, signal_reweight);
    }
    if (_cut_low <= pT_back && pT_back < _cut_high) {
        if (max_back > 2) {
            fillBackground(0, {tau3_back/tau2_back}, background_reweight);
        }
        if (max_back > 1) {
            fillBackground(1, {tau2_back/tau1_back}, background_reweight);
            fillBackground(3, {tau2_back/tau1_back, sigma_back}, background_reweight);
        }
        fillBackground(2, {sigma_back}, background_reweight);
    }
}

void AreaEfficiency::Update(EventManager const* event_manager) {
    float background_reweight = event_manager->Reweight(_event_background);
    float signal_reweight = event_manager->Reweight(_event_signal);

    float pT_sig = (*event_manager)[_event_signal].antikt_pt;
    float pT_back = (*event_manager)[_event_background].antikt_pt;

    float trim_pT_sig = (*event_manager)[_event_signal].antikt_pt_trimmed_three;
    float trim_pT_back = (*event_manager)[_event_background].antikt_pt_trimmed_three;

    float trim_mass_sig = (*event_manager)[_event_signal].antikt_m_trimmed_three;
    float trim_mass_back = (*event_manager)[_event_background].antikt_m_trimmed_three;

    float area2_sig = (*event_manager)[_event_signal].antikt_area_trimmed_two;
    float area2_back = (*event_manager)[_event_background].antikt_area_trimmed_two;

    float area3_sig = (*event_manager)[_event_signal].antikt_area_trimmed_three;
    float area3_back = (*event_manager)[_event_background].antikt_area_trimmed_three;

    float sigma_sig = (*event_manager)[_event_signal].mGMMc_r;
    float sigma_back = (*event_manager)[_event_background].mGMMc_r;

    if (_cut_low <= pT_sig && pT_sig < _cut_high) {
        fillSignal(0, {area2_sig}, signal_reweight);
        fillSignal(1, {area3_sig}, signal_reweight);
        fillSignal(2, {trim_mass_sig/trim_pT_sig}, signal_reweight);
        fillSignal(3, {area3_sig,trim_mass_sig/trim_pT_sig}, signal_reweight);
        fillSignal(4, {area3_sig,sigma_sig}, signal_reweight);
        fillSignal(5, {trim_mass_sig/trim_pT_sig, sigma_sig}, signal_reweight);
    }

    if (_cut_low <= pT_back && pT_back < _cut_high) {
        fillBackground(0, {area2_back}, background_reweight);
        fillBackground(1, {area3_back}, background_reweight);
        fillBackground(2, {trim_mass_back/trim_pT_back}, background_reweight);
        fillBackground(3, {area3_back,trim_mass_back/trim_pT_back}, background_reweight);
        fillBackground(4, {area3_back,sigma_back}, background_reweight);
        fillBackground(5, {trim_mass_back/trim_pT_back, sigma_back}, background_reweight);
    }
}

void SigmaImprovementEfficiencyTau32::Update(EventManager const* event_manager) {
    float background_reweight = event_manager->Reweight(_event_background);
    float signal_reweight = event_manager->Reweight(_event_signal);

    float pT_sig = (*event_manager)[_event_signal].antikt_pt;
    float pT_back = (*event_manager)[_event_background].antikt_pt;

    std::vector<float> *n_sub_sig = (*event_manager)[_event_signal].antikt_nsubjettiness;
    std::vector<float> *n_sub_back = (*event_manager)[_event_background].antikt_nsubjettiness;
    float tau1_sig;
    float tau2_sig;
    float tau3_sig;
    float tau1_back;
    float tau2_back;
    float tau3_back;

    if(n_sub_sig->size() > 0)
        tau1_sig = n_sub_sig->at(0);
    if(n_sub_back->size() > 0)
        tau1_back = n_sub_back->at(0);
    if(n_sub_sig->size() > 1)
        tau2_sig = n_sub_sig->at(1);
    if(n_sub_back->size() > 1)
        tau2_back = n_sub_back->at(1);
    if(n_sub_sig->size() > 2)
        tau3_sig = n_sub_sig->at(2);
    if(n_sub_back->size() > 2)
        tau3_back = n_sub_back->at(2);

    int max_sig = n_sub_sig->size();
    int max_back = n_sub_back->size();

    float sigma_sig = (*event_manager)[_event_signal].mGMMc_r;
    float sigma_back = (*event_manager)[_event_background].mGMMc_r;

    if (_cut_low <= pT_sig && pT_sig < _cut_high) {
        if (max_sig > 2) {
            fillSignal(0, {tau3_sig/tau2_sig}, signal_reweight);
            fillSignal(3, {tau3_sig/tau2_sig, sigma_sig}, signal_reweight);
        }
        if (max_sig > 1)
            fillSignal(1, {tau2_sig/tau1_sig}, signal_reweight);
        fillSignal(2, {sigma_sig}, signal_reweight);

    }
    if (_cut_low <= pT_back && pT_back < _cut_high) {
        if (max_back > 2) {
            fillBackground(0, {tau3_back/tau2_back}, background_reweight);
            fillBackground(3, {tau3_back/tau2_back, sigma_back}, background_reweight);
        }
        if (max_back > 1) {
            fillBackground(1, {tau2_back/tau1_back}, background_reweight);
        }
        fillBackground(2, {sigma_back}, background_reweight);
    }
}

void EfficiencyGenTest::Update(EventManager const* event_manager) {
    float qcd_reweight = event_manager->Reweight("qcd_5");
    float wprime_reweight = event_manager->Reweight("wprime_5");

    fillSignal(0, {(*event_manager)["wprime_5"].mGMMc_r}, wprime_reweight);
    fillBackground(0, {(*event_manager)["qcd_5"].mGMMc_r}, qcd_reweight);

    fillSignal(1, {(*event_manager)["wprime_5"].mGMMc_m}, wprime_reweight);
    fillBackground(1, {(*event_manager)["qcd_5"].mGMMc_m}, qcd_reweight);

    fillSignal(2, {(*event_manager)["wprime_5"].antikt_m}, wprime_reweight);
    fillBackground(2, {(*event_manager)["qcd_5"].antikt_m}, qcd_reweight);
}

void SigmaEfficiencyPosterPlot::Update(EventManager const* event_manager) {
    float qcd_reweight = event_manager->Reweight("qcd_5");
    float wprime_reweight = event_manager->Reweight("wprime_5");

    _signal_hists.at(0)->Fill((*event_manager)["wprime_5"].mGMMc_r, wprime_reweight);
    _background_hists.at(0)->Fill((*event_manager)["qcd_5"].mGMMc_r, qcd_reweight);

    _signal_hists.at(1)->Fill((*event_manager)["wprime_5"].mGMMc_m, wprime_reweight);
    _background_hists.at(1)->Fill((*event_manager)["qcd_5"].mGMMc_m, qcd_reweight);

    _signal_hists.at(2)->Fill((*event_manager)["wprime_5"].antikt_m, wprime_reweight);
    _background_hists.at(2)->Fill((*event_manager)["qcd_5"].antikt_m, qcd_reweight);
}

void SigmaEfficiencyPlot::Update(EventManager const* event_manager) {
    float signal_reweight = event_manager->Reweight(_event_signal);
    float background_reweight = event_manager->Reweight(_event_background);

    float pT_temp = (*event_manager)[_event_signal].antikt_pt;
    if (_cut_low <= pT_temp && pT_temp < _cut_high)
        _signal_hists.at(0)->Fill((*event_manager)[_event_signal].mGMMc_r,
                                  signal_reweight);

    pT_temp = (*event_manager)[_event_background].antikt_pt;
    if (_cut_low <= pT_temp && pT_temp < _cut_high)
        _background_hists.at(0)->Fill((*event_manager)[_event_background].mGMMc_r,
                                      background_reweight);

    pT_temp = (*event_manager)[_event_signal].antikt_pt;
    if (_cut_low <= pT_temp && pT_temp < _cut_high)
        _signal_hists.at(1)->Fill((*event_manager)[_event_signal].antikt_m,
                                  signal_reweight);

    pT_temp = (*event_manager)[_event_background].antikt_pt;
    if (_cut_low <= pT_temp && pT_temp < _cut_high)
        _background_hists.at(1)->Fill((*event_manager)[_event_background].antikt_m,
                                      background_reweight);
}

void SkewEfficiencyPlot::Update(EventManager const* event_manager) {
    float qcd_reweight = event_manager->Reweight("qcd_5");
    float zprime_reweight = event_manager->Reweight("zprime_5");

    _signal_hists.at(0)->Fill((*event_manager)["zprime_5"].mGMMc_m_skew, zprime_reweight);
    _background_hists.at(0)->Fill((*event_manager)["qcd_5"].mGMMc_m_skew, qcd_reweight);
}

void SigmaNSubjettinessEfficiencyPlot::Update(EventManager const* event_manager) {
    float signal_reweight = event_manager->Reweight(_event_signal);
    float background_reweight = event_manager->Reweight(_event_background);

    float pT_signal = (*event_manager)[_event_signal].antikt_pt;
    float pT_background = (*event_manager)[_event_background].antikt_pt;

    std::vector<float> *n_sub_sig = (*event_manager)[_event_signal].antikt_nsubjettiness;
    std::vector<float> *n_sub_back = (*event_manager)[_event_background].antikt_nsubjettiness;

    if (_cut_low <= pT_signal && pT_signal < _cut_high) {
        if (n_sub_sig->size() > 1) {
            _signal_hists.at(1)->Fill(n_sub_sig->at(1), signal_reweight);
            _signal_hists.at(2)->Fill(n_sub_sig->at(1)/n_sub_sig->at(0), signal_reweight);
        }

        if (n_sub_sig->size() > 0) {
            _signal_hists.at(0)->Fill(n_sub_sig->at(0), signal_reweight);
        }

        _signal_hists.at(3)->Fill((*event_manager)[_event_signal].mGMMc_r, signal_reweight);
    }
    if (_cut_low <= pT_background && pT_background < _cut_high) {
        if (n_sub_back->size() > 1) {
            _background_hists.at(1)->Fill(n_sub_back->at(1), background_reweight);
            _background_hists.at(2)->Fill(n_sub_back->at(1)/n_sub_back->at(0), background_reweight);
        }
        if (n_sub_back->size() > 0) {
            _background_hists.at(0)->Fill(n_sub_back->at(0), background_reweight);
        }

        _background_hists.at(3)->Fill((*event_manager)[_event_background].mGMMc_r, background_reweight);
    }
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
    float zprime_reweight = event_manager->Reweight("zprime_5");
    float qcd_reweight = event_manager->Reweight("qcd_5");
    _hists[0]->Fill((*event_manager)["zprime_5"].mGMMc_m_skew, zprime_reweight);
    _hists[1]->Fill((*event_manager)["qcd_5"].mGMMc_m_skew, qcd_reweight);
}

void DeltaRHistogram::Update(EventManager const* event_manager) {
    std::stringstream ss;

    ss.str(std::string());
    ss << _alg_label << "_dr";
    Float_t zprime_dr = (*event_manager)["zprime_5"].Get<Float_t>(ss.str());
    Float_t wprime_dr = (*event_manager)["wprime_5"].Get<Float_t>(ss.str());
    Float_t qcd_dr = (*event_manager)["qcd_5"].Get<Float_t>(ss.str());

    float zprime_reweight = event_manager->Reweight("zprime_5");
    float wprime_reweight = event_manager->Reweight("wprime_5");
    float qcd_reweight = event_manager->Reweight("qcd_5");

    _hists[0]->Fill(zprime_dr, zprime_reweight);
    _hists[1]->Fill(wprime_dr, wprime_reweight);
    _hists[2]->Fill(qcd_dr, qcd_reweight);
}

void FuzzyAntiktPtCorrelation::Update(EventManager const* event_manager) {
    float reweight = event_manager->Reweight("zprime_5");
    _hist->Fill((*event_manager)["zprime_5"].antikt_pt, (*event_manager)["zprime_5"].mGMMc_pt, reweight);
}

void RadiusComparisonHistogram::Update(EventManager const* event_manager) {
    float zprime_reweight = event_manager->Reweight("zprime_5");
    float wprime_reweight = event_manager->Reweight("wprime_5");
    float qcd_reweight = event_manager->Reweight("qcd_5");

    _hists[0]->Fill((*event_manager)["zprime_5"].mGMMc_r, zprime_reweight);
    _hists[1]->Fill((*event_manager)["wprime_5"].mGMMc_r, wprime_reweight);
    _hists[2]->Fill((*event_manager)["qcd_5"].mGMMc_r, qcd_reweight);
}

void AverageRadiusComparisonHistogram::Update(EventManager const* event_manager) {
    float zprime_reweight = event_manager->Reweight("zprime_5");
    float wprime_reweight = event_manager->Reweight("wprime_5");
    float qcd_reweight = event_manager->Reweight("qcd_5");

    _hists[0]->Fill((*event_manager)["zprime_5"].mGMMc_r_avg, zprime_reweight);
    _hists[1]->Fill((*event_manager)["wprime_5"].mGMMc_r_avg, wprime_reweight);
    _hists[2]->Fill((*event_manager)["qcd_5"].mGMMc_r_avg, qcd_reweight);
}

void SigmaAreaCorrelation::Update(EventManager const* event_manager) {
    float e_reweight = event_manager->Reweight(_event_label);
    float sigma = (*event_manager)[_event_label].mGMMc_r;
    float area = (*event_manager)[_event_label].Get<Float_t>(_alg_branch);
    _hist->Fill(sigma * sigma, area, e_reweight);
}

void JetMultiplicityPtCut::Update(EventManager const* event_manager) {
    float multiplicity = (*event_manager)[_event_label_base + "_5"].mGMM_etas->size();
    float reweight = event_manager->Reweight(_event_label_base + "_5");
    _ys[0] = (reweight * multiplicity + _ns[0]*_ys[0]) / (reweight + _ns[0]);
    _ns[0] += reweight;
    if(_ns[0] == 0) _ys[0] = 0;

    multiplicity = (*event_manager)[_event_label_base + "_15"].mGMM_etas->size();
    reweight = event_manager->Reweight(_event_label_base + "_15");
    _ys[1] = (reweight * multiplicity + _ns[1]*_ys[1]) / (reweight + _ns[1]);
    _ns[1] += reweight;
    if(_ns[1] == 0) _ys[1] = 0;

    multiplicity = (*event_manager)[_event_label_base + "_25"].mGMM_etas->size();
    reweight = event_manager->Reweight(_event_label_base + "_25");
    _ys[2] = (reweight * multiplicity + _ns[2]*_ys[2]) / (reweight + _ns[2]);
    _ns[2] += reweight;
    if(_ns[2] == 0) _ys[2] = 0;
}

void RadiusPtSeedHistogram::Update(EventManager const* event_manager) {
    float reweight = event_manager->Reweight(_event_label_base + "_5");
    float sigma = (*event_manager)[_event_label_base + "_5"].mGMMc_r;
    _hists[0]->Fill(sigma, reweight);

    reweight = event_manager->Reweight(_event_label_base + "_15");
    sigma = (*event_manager)[_event_label_base + "_15"].mGMMc_r;
    _hists[1]->Fill(sigma, reweight);

    reweight = event_manager->Reweight(_event_label_base + "_25");
    sigma = (*event_manager)[_event_label_base + "_25"].mGMMc_r;
    _hists[2]->Fill(sigma, reweight);
}

void AverageRadiusPtSeedHistogram::Update(EventManager const* event_manager) {
    float reweight = event_manager->Reweight(_event_label_base + "_5");
    float sigma = (*event_manager)[_event_label_base + "_5"].mGMMc_r_avg;
    _hists[0]->Fill(sigma, reweight);

    reweight = event_manager->Reweight(_event_label_base + "_15");
    sigma = (*event_manager)[_event_label_base + "_15"].mGMMc_r_avg;
    _hists[1]->Fill(sigma, reweight);

    reweight = event_manager->Reweight(_event_label_base + "_25");
    sigma = (*event_manager)[_event_label_base + "_25"].mGMMc_r_avg;
    _hists[2]->Fill(sigma, reweight);
}
