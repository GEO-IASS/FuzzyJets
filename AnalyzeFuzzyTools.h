#ifndef ANALYZEFUZZYTOOLS_H // header guard
#define ANALYZEFUZZYTOOLS_H

#include <vector>
#include <sstream>
#include <math.h>
#include <string>
#include <stdint.h>
#include <assert.h>

#include <TROOT.h>
#include <TColor.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TAxis.h>
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLegend.h"


#include "AtlasUtils.h"

namespace StyleTypes {
    enum HistOptions {
        NONE = 0,
        LINE = 1,
        FILL = 2,
        DASHED = 4,
        STRIPED = 8
    };
}

class CanvasHelper {
public:
    std::string x_label;
    std::string y_label;
    std::string title;
    std::string base;

    size_t width;
    size_t height;

    bool diff_scale;

    CanvasHelper(std::string x_label, std::string y_label, std::string title,
                 std::string base, size_t width, size_t height) {
        this->x_label = x_label;
        this->y_label = y_label;
        this->title = title;
        this->base = base;

        this->width = width;
        this->height = height;
        this->diff_scale = false;
    }
};

class HistHelper {
public:
    std::string file_name;
    std::string branch_name;
    std::string legend_label;

    size_t      ticks;

    double      min_edge, max_edge;
    size_t      n_bins;
    
    Style_t style;
    Color_t color;
    std::string draw_opts;

    StyleTypes::HistOptions options;

    HistHelper(std::string file_name, std::string branch_name,
               std::string legend_label, size_t ticks, double min_edge,
               double max_edge, size_t n_bins, StyleTypes::HistOptions options,
               Color_t color, Style_t style, std::string draw_opts) {
        this->file_name = file_name;
        this->branch_name = branch_name;
        this->legend_label = legend_label;
        this->ticks = ticks;
        this->options = options;
        this->min_edge = min_edge;
        this->max_edge = max_edge;
        this->n_bins = n_bins;
        this->color = color;
        this->style = style;
        this->draw_opts = draw_opts;
    }
};

// load a single branch from file as a vector, pretty efficient
template<typename T>
std::vector<T>
loadSingleBranch(std::string const& file, std::string const& branch) {
    TFile *f = new TFile(TString(file));
    assert(f && "Wasn't provided a valid file.");

    TTree *t = (TTree *) (f->Get("EventTree"));
    assert(t && "File contains no TTree by the name 'EventTree'");

    T branch_target;
    TBranch *b = t->GetBranch(TString(branch));
    b->SetAddress(&branch_target);

    std::vector<T> out;

    size_t n_entries = t->GetEntries();
    for (size_t entry_iter = 0; entry_iter < n_entries; entry_iter++) {
        b->GetEvent(entry_iter);
        out.push_back(branch_target);
    }

    //delete b;
    //delete f;
    delete f;

    return out;
}

// load a bunch of branches from file as a vector, not terribly efficient
// because it uses loadSingleBranch
template<typename T>
std::vector<std::vector<T> >
loadBranches(std::vector<std::string> const& files,
             std::vector<std::string> const& branches) {
    assert(files.size() == branches.size());

    std::vector<std::vector<T> > out;
    const size_t n_branches = branches.size();

    for (size_t branch_iter = 0; branch_iter < n_branches; branch_iter++) {
        out.push_back(loadSingleBranch<T>(files[branch_iter], branches[branch_iter]));
    }

    // want to make sure that we don't run into any issues with different
    // numbers of events, this might invalidate any results!
    for (size_t branch_iter = 0; branch_iter < n_branches; branch_iter++) {
        assert(out[branch_iter].size() == out[0].size());
    }
    return out;
}

// helpers needed for prettyHist
double HistMaximum(std::vector<TH1D *> const& hists);

size_t RoundDoubleUp(double d, size_t q);

std::string fileify(std::string text);

void SetupATLASStyle();

void DrawAtlasLabel(std::string title);


// pretty mass histograms
template<typename T>
void prettyHist(std::vector<HistHelper> const& hist_decs,
                CanvasHelper const& c_dec) {
std::vector<TH1D *> v_hists;
    std::vector<std::vector<T> > v_vals;
    const size_t n_hists = hist_decs.size();
    for (size_t hist_iter = 0; hist_iter < n_hists; hist_iter++) {
        HistHelper hist_dec = hist_decs[hist_iter];
        TH1D *h = new TH1D(TString(hist_dec.branch_name), TString(hist_dec.branch_name),
                           hist_dec.n_bins, hist_dec.min_edge, hist_dec.max_edge);
        v_hists.push_back(h);
        std::cout << hist_dec.file_name << std::endl;
        std::cout << hist_dec.branch_name << std::endl;
        v_vals.push_back(loadSingleBranch<T>(hist_dec.file_name, hist_dec.branch_name));
    }

    // ... if no branches provided
    if (!v_vals.size()) return;

    const size_t n_events = v_vals[0].size();
    for (size_t event_iter = 0; event_iter < n_events; event_iter++) {
        // fill each histogram with its next value
        for (size_t hist_iter = 0; hist_iter < n_hists; hist_iter++) {
            const double v = v_vals[hist_iter][event_iter];
            //if(hist_decs[hist_iter].min_edge <= v && v <= hist_decs[hist_iter].max_edge) {
            if (true) {
                v_hists[hist_iter]->Fill(v);
            }
        }
    }

    // Now that we have filled histograms, time to draw them
    SetupATLASStyle();
    TCanvas *canv = new TCanvas("c", c_dec.title.c_str(), 0, 0, c_dec.width, c_dec.height);
    double all_max = HistMaximum(v_hists);

    for (size_t hist_iter = 0; hist_iter < n_hists; hist_iter++) {
        TH1D *current_hist = v_hists[hist_iter];
        HistHelper hist_dec = hist_decs[hist_iter];

        // axes and labels
        current_hist->GetXaxis()->SetNdivisions(hist_dec.ticks);
        current_hist->GetXaxis()->SetTitleOffset(1.4);

        // draw options
        if (hist_dec.options & StyleTypes::DASHED) {
            current_hist->SetLineStyle(7); // dashed style
        }
        current_hist->SetLineColor(hist_dec.color);
        current_hist->SetMarkerColor(hist_dec.color);
        current_hist->SetMarkerStyle(hist_dec.style);

        if (hist_dec.options == StyleTypes::STRIPED) {
            current_hist->SetFillStyle(3004);
            current_hist->SetFillColor(hist_dec.color);
        }
        std::string postfix;
        if (c_dec.diff_scale) {
            current_hist->Scale(1./current_hist->Integral(-1, current_hist->GetNbinsX()+1));
        } else {
            current_hist->GetYaxis()->SetRangeUser(0, RoundDoubleUp(all_max*1.2, 10));
        }
        if (hist_iter == 0) {
            current_hist->GetXaxis()->SetTitle(c_dec.x_label.c_str());
            postfix = "";
        } else {
            postfix = " same";
        }
        //std::cout << current_hist->Integral(-1, current_hist->GetNbinsX() + 1) << std::endl;
        if (hist_dec.draw_opts == "") {
            if (postfix != "") {
                current_hist->Draw("same");
            } else {
                current_hist->Draw(TString::Format("%s",postfix.c_str()));
            }
        } else {
            current_hist->Draw(TString::Format("%s%s", hist_dec.draw_opts.c_str(), postfix.c_str()));
        }
    }

    TLegend *leggaa = new TLegend(0.6, 0.7, 0.9, 0.8);
    leggaa->SetTextFont(42);
    for (size_t hist_iter = 0; hist_iter < n_hists; hist_iter++) {
        if (hist_decs[hist_iter].draw_opts == "p") {
            leggaa->AddEntry(v_hists[hist_iter], hist_decs[hist_iter].legend_label.c_str(), "p");
        } else {
            leggaa->AddEntry(v_hists[hist_iter], hist_decs[hist_iter].legend_label.c_str(), "f");
        }
    }
    leggaa->SetFillStyle(0);
    leggaa->SetFillColor(0);
    leggaa->SetBorderSize(0);
    leggaa->Draw();
    DrawAtlasLabel(c_dec.title);

    std::string out_file_name_base = fileify(c_dec.title);
    std::stringstream ss;
    ss << c_dec.base << out_file_name_base << ".pdf";
    std::string full_file_name = ss.str();

    canv->Print(full_file_name.c_str());

    delete leggaa;
    delete canv;
    for (size_t hist_iter = 0; hist_iter < n_hists; hist_iter++) {
        delete v_hists[hist_iter];
    }
}


#endif
