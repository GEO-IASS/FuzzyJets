#include <iostream>
#include <sstream>
#include <assert.h>
#include <math.h>

#include "TFile.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TH1F.h"
#include "THStack.h"
#include "TLeaf.h"
#include "TTree.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TList.h"
#include "TKey.h"
#include "TObject.h"
#include "TLegend.h"

#include "ComputationManager.h"
#include "boost/foreach.hpp"
#include "AtlasStyle.h"
//#include "Util.h"

std::vector<Color_t> color_prefs = {kBlue, kRed, kBlack};
std::string out_directory = "/u/at/chstan/nfs/summer_2014/ForConrad/results/tmp/";

class ReadsEventTree : public ComputationBase {
protected:
    TFile *_file;
    TTree *_tree;
    std::map<label, void*> _vals;
    size_t _current_event;
    size_t _n_events;

public:
    ReadsEventTree(std::string filename, std::string var_prefix) {
        _provides = {};
        _requires = {};
        _name = "ReadsFile:" + filename;
        _file = TFile::Open(filename.c_str(), "READ");
        _tree = (TTree *) _file->Get("EventTree");
        _current_event = 0;
        _n_events = _tree->GetEntries();

        {
            TObjArray *arr = _tree->GetListOfLeaves();
            TLeaf *leaf;

            TIter next(arr);

            std::stringstream ss;
            while ((leaf = (TLeaf*) next())) {
                ss.str(std::string());
                ss << var_prefix << leaf->GetName();
                _provides.push_back(ss.str());
                _vals[ss.str()] = new char[leaf->GetLenType()];
                _tree->SetBranchAddress(leaf->GetName(), _vals[ss.str()]);
            }
        }
    }

    void *how_to_provide(__attribute__((unused)) label l) {
        assert(_vals.find(l) != _vals.end());
        return _vals[l];
    }

    void compute(__attribute__((unused)) ComputationManager *manager) {
        assert(_current_event < _n_events);
        _tree->GetEvent(_current_event);
        _current_event++;
    }
};

class HistogramComputation : public ComputationBase {
protected:
    TH1F *_hist;
    label _to_histogram;

public:
    HistogramComputation(label to_histogram, size_t n_bins, float l, float h) {
        _hist = new TH1F((to_histogram + "_hh").c_str(), "", n_bins, l, h);
        _to_histogram = to_histogram;
        _provides = {"H(" + to_histogram + ")"};
        _requires = {to_histogram};
        _name = "histograms:" + to_histogram;
    }

    void *how_to_provide(label l) {
        assert(l == _provides.at(0));
        return _hist;
    }

    void compute(ComputationManager *manager) {
        float v = *(manager->Get<float *>(_to_histogram));
        _hist->Fill(v);
    }

    void publish(__attribute__((unused)) ComputationManager *manager) {
    }
};

class DiffSquareComputation : public ComputationBase {
protected:
    float _diff_sq;
    label _histAl;
    label _histBl;

public:
    DiffSquareComputation(label histAl, label histBl, label l) {
        _histAl = histAl;
        _histBl = histBl;
        _provides = {l};
        _requires = {histAl, histBl};
        _name = "dsq:" + histAl + "," + histBl;
    }

    void *how_to_provide(label l) {
        assert(l == _provides.at(0));
        return &_diff_sq;
    }

    void compute(__attribute__((unused)) ComputationManager *manager) {
        // nothing
    }

    void publish(__attribute__((unused)) ComputationManager *manager) {
        TH1F *_histA = manager->Get<TH1F *>(_histAl);
        TH1F *_histB = manager->Get<TH1F *>(_histBl);

        float normA = _histA->Integral(-1, _histA->GetNbinsX() + 1);
        float normB = _histB->Integral(-1, _histB->GetNbinsX() + 1);

        _diff_sq = 0;
        for (int idx = 0; idx <= _histA->GetXaxis()->GetNbins() + 1; idx++) {
            float v = (_histA->GetBinContent(idx) / normA)
                    - (_histB->GetBinContent(idx) / normB);
            _diff_sq += v * v / (1.5/40);
        }
    }
};

class GraphDiffSquares : public ComputationBase {
protected:
    std::vector<size_t> _NPVs;

public:
    GraphDiffSquares(std::vector<size_t> NPVs, label l) {
        _provides = {l};
        _NPVs = NPVs;

        _requires = {};
        std::stringstream ss;
        BOOST_FOREACH(size_t NPV, NPVs) {
            ss.str(std::string());
            ss << "C:dsq_sigma1_" << NPV;
            _requires.push_back(ss.str());
            ss.str(std::string());
            ss << "UC:dsq_sigma1_" << NPV;
            _requires.push_back(ss.str());
        }

        _name = "gdsq:" + l;
    }

    void *how_to_provide(label l) {
        assert(l == _provides.at(0));
        return NULL;
    }

    void compute(__attribute__((unused)) ComputationManager *manager) {
        // nothing
    }

    void publish(__attribute__((unused)) ComputationManager *manager) {
        std::vector<float> xs;
        std::vector<float> yCs;
        std::vector<float> yUCs;

        std::stringstream ss;
        BOOST_FOREACH(size_t NPV, _NPVs) {
            xs.push_back(static_cast<float>(NPV));

            ss.str(std::string());
            ss << "C:dsq_sigma1_" << NPV;
            yCs.push_back(*(manager->Get<float *>(ss.str())));

            ss.str(std::string());
            ss << "UC:dsq_sigma1_" << NPV;
            yUCs.push_back(*(manager->Get<float *>(ss.str())));
        }

        TCanvas canvas("canv", "", 600, 600);
        canvas.cd();

        TMultiGraph mg;
        TGraph *g_corrected = new TGraph(xs.size(), &xs.at(0), &yCs.at(0));
        g_corrected->SetMarkerColor(color_prefs.at(0));
        g_corrected->SetMarkerStyle(22);
        TGraph *g_uncorrected = new TGraph(xs.size(), &xs.at(0), &yUCs.at(0));
        g_uncorrected->SetMarkerColor(color_prefs.at(1));
        g_uncorrected->SetMarkerStyle(23);

        mg.Add(g_corrected, "p");
        mg.Add(g_uncorrected, "p");

        mg.Draw("ap");
        mg.GetXaxis()->SetTitle("");
        mg.GetYaxis()->SetTitle("");

        TLegend legend(0.2, 0.76, 0.45, 0.91);
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
        ss << out_directory << _name << ".pdf";
        canvas.Print(ss.str().c_str(), "pdf");
        canvas.Clear();
    }
};

class HistStackComputation : public ComputationBase {
protected:
    std::vector<std::string> _legend_labels;

public:
    HistStackComputation(std::vector<label> vars, std::vector<std::string> legend_labels) {
        _name = "stacks:";
        BOOST_FOREACH(label label, vars) {
            _name = _name + label + ",";
        }
        _legend_labels = legend_labels;
        std::stringstream ss;
        ss.str(std::string());
        ss << "S(";
        for (size_t i = 0; i < vars.size(); i++) {
            ss << vars.at(i);
            if (i < vars.size() - 1)
                ss << ",";
        }
        ss << ")";
        _provides = {ss.str()};
        _requires = vars;
    }

    void *how_to_provide(__attribute__((unused)) label l) {
        return NULL;
    }

    void compute(__attribute__((unused)) ComputationManager *manager) {
    }

    void publish(__attribute__((unused)) ComputationManager *manager) {
        TCanvas canvas("canv", "", 600, 600);
        canvas.cd();
        THStack hist_stack("stack", "");


        for (size_t idx = 0; idx < _requires.size(); idx++) {
            label l = _requires.at(idx);
            TH1F *h = manager->Get<TH1F *>(l);
            h->SetLineColor(color_prefs.at(idx));
            TColor *col = gROOT->GetColor(idx + 3);
            TColor *col_ref = gROOT->GetColor(color_prefs.at(idx));
            col->SetRGB(col_ref->GetRed(), col_ref->GetGreen(), col_ref->GetBlue());
            col->SetAlpha(0.28);

            h->SetLineStyle(1);
            h->SetFillColor(idx + 3);
            h->SetFillStyle(1001);
            h->SetLineWidth(2);

            float norm = h->Integral(-1, h->GetNbinsX() + 1);
            h->Scale(1./norm);
            hist_stack.Add(h);
        }

        hist_stack.Draw("nostack");
        hist_stack.GetHistogram()->GetYaxis()->SetTitle("");
        hist_stack.GetHistogram()->GetXaxis()->SetTitle("");
        hist_stack.GetHistogram()->GetXaxis()->SetNdivisions(505);

        hist_stack.SetMaximum(0.35);
        hist_stack.Draw("nostack");

        //TLegend legend(0.65, 0.75, 0.95, 0.92);
        //ss.str(std::string());
        //ss << "Corrected at NPV = " << NPV;
        //legend.AddEntry(ns1m[NPV], ss.str().c_str(), "f");
        //ss.str(std::string());
        //ss << "Uncorrected at NPV = " << NPV;
        //legend.AddEntry(ns1um[NPV], ss.str().c_str(), "f");
        //
        //legend.SetTextFont(42);
        //legend.SetFillStyle(0);
        //legend.SetFillColor(0);
        //legend.SetBorderSize(0);
        //legend.SetTextSize(1.3 * legend.GetTextSize());
        //
        //legend.Draw();

        canvas.Update();
        std::stringstream ss;
        ss.str(std::string());
        std::string out_directory = "~/nfs/summer_2014/ForConrad/results/tmp/";
        ss << out_directory << "test.pdf";
        canvas.Print(ss.str().c_str(), "pdf");
        canvas.Clear();
    }
};

void SetupATLASStyle() {
    gStyle->SetOptStat(0);
    gROOT->Reset();
    AtlasStyle();
    gROOT->SetStyle("ATLAS");
    gROOT->ForceStyle();
    gStyle->SetPadLeftMargin(0.16);
}

int main(__attribute__((unused)) int argc, __attribute__((unused)) char **argv) {
    SetupATLASStyle();

    std::string file_c = "~/nfs/summer_2014/ForConrad/files/20k_zprime_pileup_effect_corrected/2015_03_15_17h04m17s/fin.root";
    std::string file_uc = "~/nfs/summer_2014/ForConrad/files/20k_zprime_pileup_effect_uncorrected/2015_03_15_15h54m31s/fin.root";

    ComputationManager manager;

    manager.install_computation(new ReadsEventTree(file_c, "C:"));
    manager.install_computation(new ReadsEventTree(file_uc, "UC:"));

    std::vector<size_t> NPVs = {0, 2, 5, 10, 15, 20, 25, 30, 35, 40};
    BOOST_FOREACH(size_t NPV, NPVs) {
        std::string base_UC = "H(UC:sigma1_0)";
        std::string base_C = "H(C:sigma1_0)";

        std::stringstream ss;
        ss.str(std::string());
        ss << "C:sigma1_" << NPV;
        ComputationBase *c1 = new HistogramComputation(ss.str(), 40, 0, 1.5);
        manager.install_computation(c1);

        ss.str(std::string());
        ss << "UC:sigma1_" << NPV;
        ComputationBase *c2 = new HistogramComputation(ss.str(), 40, 0, 1.5);
        manager.install_computation(c2);

        ss.str(std::string());
        ss << "C:dsq_sigma1_" << NPV;
        manager.install_computation(new DiffSquareComputation(c1->provided().at(0), base_C, ss.str()));

        ss.str(std::string());
        ss << "UC:dsq_sigma1_" << NPV;
        manager.install_computation(new DiffSquareComputation(c2->provided().at(0), base_UC, ss.str()));
    }

    manager.push_to_final_state();

    std::vector<std::string> vars = {"H(C:sigma1_40)", "H(C:sigma1_0)"};
    std::vector<std::string> lls = {};

    manager.install_computation(new HistStackComputation(vars, lls));
    manager.install_computation(new GraphDiffSquares(NPVs, "test2"));

    manager.compile();
    manager.compute(20000);
    manager.publish();
    return 0;
}
