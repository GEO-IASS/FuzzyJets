#include <iostream>
#include <string>
#include <vector>
#include <stdint.h>
#include <assert.h>
#include <ctype.h>

#include <TROOT.h>
#include "TH1F.h"

#include "AnalyzeFuzzyTools.h"

double HistMaximum(std::vector<TH1D *> const& hists) {
    double max = -1;

    const size_t n_hists = hists.size();
    for (size_t hist_iter = 0; hist_iter < n_hists; hist_iter++) {
        double n_max = hists[hist_iter]->GetMaximum();
        if (n_max > max) {
            max = n_max;
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
