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
    std::string xLabel;
    std::string yLabel;
    std::string title;
    std::string base;

    size_t width;
    size_t height;

    CanvasHelper(std::string xLabel, std::string yLabel, std::string title,
                 std::string base, size_t width, size_t height) {
        this->xLabel = xLabel;
        this->yLabel = yLabel;
        this->title = title;
        this->base = base;

        this->width = width;
        this->height = height;
    }
};

class HistHelper {
public:
    std::string fileName;
    std::string branchName;
    std::string legendLabel;

    size_t      ticks;

    double      minEdge, maxEdge;
    size_t      nBins;

    Color_t color;

    StyleTypes::HistOptions options;

    HistHelper(std::string fileName, std::string branchName,
               std::string legendLabel, size_t ticks, double minEdge,
               double maxEdge, size_t nBins, StyleTypes::HistOptions options,
               Color_t color) {
        this->fileName = fileName;
        this->branchName = branchName;
        this->legendLabel = legendLabel;
        this->ticks = ticks;
        this->options = options;
        this->minEdge = minEdge;
        this->maxEdge = maxEdge;
        this->nBins = nBins;
        this->color = color;
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

    T branchTarget;
    TBranch *b = t->GetBranch(TString(branch));
    b->SetAddress(&branchTarget);

    std::vector<T> out;

    size_t nEntries = t->GetEntries();
    for (size_t iEntry = 0; iEntry < nEntries; iEntry++) {
        b->GetEvent(iEntry);
        out.push_back(branchTarget);
    }
    //delete b;
    //delete f;
    //delete t;

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

// helpers needed for prettyHist
double HistMaximum(std::vector<TH1D *> const& hists);

size_t RoundDoubleUp(double d, size_t q);

std::string fileify(std::string text);

void SetupATLASStyle();

void DrawAtlasLabel(std::string title);


// pretty mass histograms
template<typename T>
void prettyHist(std::vector<HistHelper> const& histdecs,
                CanvasHelper const& cdec) {
std::vector<TH1D *> vHists;
    std::vector<std::vector<T> > vVals;
    const size_t nHists = histdecs.size();
    for (size_t iHist = 0; iHist < nHists; iHist++) {
        HistHelper histdec = histdecs[iHist];
        TH1D *h = new TH1D(TString(histdec.branchName), TString(histdec.branchName),
                           histdec.nBins, histdec.minEdge, histdec.maxEdge);
        vHists.push_back(h);
        std::cout << histdec.fileName << std::endl;
        std::cout << histdec.branchName << std::endl;
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
        if (histdec.options & StyleTypes::DASHED) {
            currentHist->SetLineStyle(7); // dashed style
        }
        currentHist->SetLineColor(histdec.color);

        if (iHist == 0) {
            currentHist->GetYaxis()->SetRangeUser(0, RoundDoubleUp(allMax*1.2, 10));
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


#endif
