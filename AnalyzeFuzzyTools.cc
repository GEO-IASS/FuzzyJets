#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>

#include <TROOT.h>
#include <TColor.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TAxis.h>
#include "TH1F.h"
#include "TTree.h"
#include "TFile.h"
#include "TLegend.h"

#include "AnalyzeFuzzyTools.h"
#include "AtlasUtils.h"

// load a single branch from file as a vector, pretty efficient
template<typename T>
std::vector<T>
loadSingleBranch(std::string const& file, std::string const& branch) {
    TFile *f = new TFile(TString(file));
    assert(("Wasn't provided a valid file.", f));

    TTree *t = (TTree *) (f->Get("EventTree"));
    assert(("File contains no TTree by the name 'EventTree'", t));

    T branchTarget;
    TBranch *b = t->GetBranch(TString(branch));
    b->SetAddress(&branchTarget);

    std::vector<T> out;

    size_t nEntries = t->GetEntries();
    for (size_t iEntry = 0; iEntry < nEntries; iEntry++) {
        b->GetEvent(iEntry);
        out.push_back(branchTarget);
    }
    delete b;
    delete f;
    delete t;

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
    const size_t nBranches = branches.size();

    for (size_t iBranch = 0; iBranch < nBranches; iBranch++) {
        out.push_back(loadSingleBranch<T>(files[iBranch], branches[iBranch]));
    }

    // want to make sure that we don't run into any issues with different
    // numbers of events, this might invalidate any results!
    for (size_t iBranch = 0; iBranch < nBranches; iBranch++) {
        assert(out[iBranch].size() == out[0].size());
    }
    return out;
}

double HistMaximum(std::vector<TH1D *> const& hists) {
    double max = -1;

    const size_t nHists = hists.size();
    for (size_t iHist = 0; iHist < nHists; iHist++) {
        double nmax = hists[iHist]->GetMaximum();
        if (nmax > max) {
            max = nmax;
        }
    }
    return max;
}

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
    gROOT->SetStyle("ATLAS");
    gROOT->ForceStyle();
    gStyle->SetPadLeftMargin(0.16);
}

void DrawAtlasLabel(std::string title) {
    ATLAS_LABEL(0.2, 0.88, kBlack);
    myText(0.38, 0.88, kBlack, "Internal");
    myText(0.2, 0.8, kBlack, "#sqrt{s} = 8000 GeV");
    std::stringstream ss;
    std::string temp = ss.str();

    ss << "#scale[0.5]{" << title << "}";
    myText(0.2, 0.74, kBlack, temp.c_str());
}

// pretty mass histograms
template<typename T>
void prettyHist(std::vector<HistHelper> const& histdecs,
                CanvasHelper const& cdec) {
std::vector<TH1D *> vHists;
    std::vector<std::vector<double> > vVals;
    const size_t nHists = histdecs.size();
    for (size_t iHist = 0; iHist < nHists; iHist++) {
        HistHelper histdec = histdecs[iHist];
        TH1D *h = new TH1D(TString(histdec.branchName), TString(histdec.branchName),
                           histdec.nBins, histdec.minEdge, histdec.maxEdge);
        vHists.push_back(h);
        vVals.push_back(loadSingleBranch<T>(histdec.fileName, histdec.branchName));
    }
    const size_t nEvents = vVals[0].size();
    for (size_t iEvent = 0; iEvent < nEvents; iEvent++) {
        // fill each histogram with its next value
        for (size_t iHist = 0; iHist < nHists; iHist++) {
            const double v = vVals[iHist][iEvent];
            if(histdecs[iHist].minEdge <= v && v <= histdecs[iHist].maxEdge) {
                vHists[iHist]->Fill(v);
            }
        }
    }

    // Now that we have filled histograms, time to draw them
    SetupATLASStyle();
    TCanvas *canv = new TCanvas("c", cdec.title.c_str(), 0, 0, cdec.width, cdec.height);
    double allMax = HistMaximum(vHists);

    for (size_t iHist = 0; iHist < nHists; iHist++) {
        TH1D *currentHist = vHists[iHist];
        HistHelper histdec = histdecs[iHist];

        // axes and labels
        currentHist->GetXaxis()->SetNdivisions(histdec.ticks);
        currentHist->GetXaxis()->SetTitleOffset(1.4);

        // draw options
        currentHist->SetLineColor(histdec.color);
        if (iHist == 0) {
            currentHist->GetYaxis()->SetRangeUser(0, RoundDoubleUp(allMax, 10));
            currentHist->GetXaxis()->SetTitle(cdec.xLabel.c_str());
            currentHist->Draw();
        } else {
            currentHist->Draw("same");
        }
    }

    TLegend *leggaa = new TLegend(0.6, 0.7, 0.9, 0.8);
    leggaa->SetTextFont(42);
    for (size_t iHist = 0; iHist < nHists; iHist++) {
        leggaa->AddEntry(vHists[iHist], histdecs[iHist].legendLabel.c_str(), "f");
    }
    leggaa->SetFillStyle(0);
    leggaa->SetFillColor(0);
    leggaa->SetBorderSize(0);
    leggaa->Draw();
    DrawAtlasLabel(cdec.title);

    std::string outfilenamebase = fileify(cdec.title);
    std::stringstream ss;
    ss << cdec.base << outfilenamebase << ".pdf";
    std::string fullfilename = ss.str();

    canv->Print(fullfilename.c_str());

    delete leggaa;
    delete canv;
    for (size_t iHist = 0; iHist < nHists; iHist++) {
        delete vHists[iHist];
    }
}
