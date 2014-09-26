#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <TH2F.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>


#include <sstream>
#include <iostream>
#include <stdlib.h>

class EventManager;

const std::string plot_prefix = "/u/at/chstan/nfs/summer_2014/ForConrad/results/plots/testing/";

class UpdatesOnEvent {
public:
    UpdatesOnEvent() {}
    virtual ~UpdatesOnEvent() {}

    virtual void Start(__attribute__((unused)) EventManager const* event_manager) {}
    virtual void Update(EventManager const* event_manager) = 0;
    virtual void Finish(__attribute__((unused)) EventManager const* event_manager) {}
};

void gen_string(char* s, const unsigned int len);

class StackedEfficiencyHistogramBase : public UpdatesOnEvent {
protected:    
    const static unsigned int _root_name_length = 20;
    std::vector<char *> _root_names_signal;
    std::vector<char *> _root_names_background;
    std::vector<TH1F *> _signal_hists;
    std::vector<TH1F *> _background_hists;
    std::vector<Color_t> _colors;
    std::vector<std::string> _labels;
    std::vector<int> _signs;

    std::vector<float> _lows;
    std::vector<float> _highs;

    std::string _title;
    std::string _outfile_name;

    unsigned int _canvas_x, _canvas_y;
    unsigned int _n_points;
    

    
    size_t _ticks;

    std::string _x_label;
    std::string _y_label;

    StackedEfficiencyHistogramBase() {
        _outfile_name = "UNNAMED.PDF";

        _title = "UNTITLED";

        _x_label = "X_LABEL";
        _y_label = "Y_LABEL";
        
        _canvas_x = _canvas_y = 800;
        _n_points = 200;
        
        _ticks = 505;
    }

    ~StackedEfficiencyHistogramBase() {
        for (unsigned int iter = 0; iter < _signal_hists.size(); iter++) {
            delete _signal_hists[iter];
            delete _background_hists[iter];
            delete[] _root_names_signal[iter];
            delete[] _root_names_background[iter];
        }
    }

public:
    void Start(__attribute__((unused)) EventManager const* event_manager) {
        for (unsigned int iter = 0; iter < _colors.size(); iter++) {
            _root_names_signal.push_back(new char[_root_name_length + 1]);
            _root_names_background.push_back(new char[_root_name_length + 1]);
            gen_string(_root_names_signal.at(iter), _root_name_length);
            gen_string(_root_names_background.at(iter), _root_name_length);
            _signal_hists.push_back(new TH1F(_root_names_signal[iter], _title.c_str(), _n_points-1, _lows[iter], _highs[iter]));
            _background_hists.push_back(new TH1F(_root_names_background[iter], _title.c_str(), _n_points-1, _lows[iter], _highs[iter]));
        }
    }

    virtual void Finish(__attribute__((unused)) EventManager const* event_manager);
};

class SigmaEfficiencyPosterPlot : public StackedEfficiencyHistogramBase {
public:
    SigmaEfficiencyPosterPlot() {
        _title = "";
        _x_label = "W' Efficiency";
        _y_label = "1 - QCD Efficiency";
        _outfile_name = "SigmaEfficiencyPoster.pdf";
        _n_points = 500;
        
        _colors = {kRed, kBlue, kBlack};
        _labels = {"Fuzzy #sigma", "Fuzzy Mass", "Anti-k_{T} Mass"};
        _signs = {-1, -1, -1};

        _ticks = 405;
        
        _lows = {0, 0, 0};
        _highs = {1, 400, 400};
    }
    void Update(EventManager const* event_manager);
};

class SigmaEfficiencyPlot : public StackedEfficiencyHistogramBase {
protected:
    std::string _event_signal;
    std::string _event_background;
    
    float _cut_low, _cut_high;

public:
    SigmaEfficiencyPlot(std::string event_signal, std::string event_background, 
                        float cut_low, float cut_high) {
        _event_signal = event_signal;
        _event_background = event_background;

        _cut_low = cut_low;
        _cut_high = cut_high;

        _title = "";

        std::stringstream ss;
        ss.str(std::string());
        ss << "SigmaEfficiency_" << event_signal << "_" << event_background << "_" 
           << (int) cut_low << "to" << (int) cut_high << "pT.pdf";
        _outfile_name = ss.str();

        ss.str(std::string());
        ss << _event_signal << " Efficiency";
        _x_label = ss.str();

        ss.str(std::string());
        ss << _event_background << "1 - Efficiency";
        _y_label = ss.str();

        _n_points = 1000;
        
        _colors = {kRed, kBlue};
        _labels = {"Fuzzy Sigma", "Anti-kt Mass"};
        _signs = {1, 1};
        
        _ticks = 405;

        _lows = {0, 0};
        _highs = {1, 400};
    }

    void Update(EventManager const* event_manager);
};

class SkewEfficiencyPlot : public StackedEfficiencyHistogramBase {
public:
    SkewEfficiencyPlot() {
        _title = "";
        _outfile_name = "SkewEfficiency.pdf";

        _x_label = "Z'#rightarrow tt' Efficiency";
        _y_label = "1 - QCD Efficiency";

        _n_points = 300;
        _colors = {kRed};
        _labels = {"Fuzzy Jet Skew"};
        _signs = {1};
        
        _ticks = 405;
        
        _lows = {-1};
        _highs = {1};
    }

    void Update(EventManager const* event_manager);
};

class FuzzyJetMassEfficiencyPlot  : public StackedEfficiencyHistogramBase {
protected:
    std::string _event_signal;
    std::string _event_background;
    std::string _alg;
    
    float _cut_low, _cut_high;

public:
    FuzzyJetMassEfficiencyPlot(std::string event_signal, std::string event_background, std::string alg,
                               float cut_low, float cut_high) {
        _event_signal = event_signal;
        _event_background = event_background;
        _alg = alg;

        _cut_low = cut_low;
        _cut_high = cut_high;

        _title = "";

        std::stringstream ss;
        ss.str(std::string());
        ss << "FuzzyJetMassEfficiency_" << event_signal << "_" << event_background << "_" 
           << alg << "_" << (int) cut_low << "to" << (int) cut_high << "pT.pdf";
        _outfile_name = ss.str();

        ss.str(std::string());
        ss << _event_signal << " Efficiency";
        _x_label = ss.str();

        ss.str(std::string());
        ss << _event_background << "1 - Efficiency";
        _y_label = ss.str();

        _n_points = 1000;
        
        _colors = {kRed, kBlue};

        ss.str(std::string());
        ss << _alg << " Mass";
        _labels = {ss.str(), "Anti-kt Mass"};
        _signs = {1, 1};
        
        _ticks = 405;

        _lows = {0, 0};
        _highs = {400, 400};
    }

    void Update(EventManager const* event_manager);
};

class ScatterBase : public UpdatesOnEvent {
protected:
    Color_t _color;
    Style_t _style;
    
    std::string _title;
    std::string _outfile_name;
    
    unsigned int _canvas_x, _canvas_y;

    std::vector<float> _xs, _ys;

    std::string _x_label, _y_label;

    ScatterBase() {
        _outfile_name = "UNNAMED.PDF";
        _title = "UNTITLED";
        
        _x_label = "X_LABEL";
        _y_label = "Y_LABEL";
        
        _canvas_x = _canvas_y = 800;        
    }
public:
    virtual void Finish(EventManager const* event_manager);
};

class JetMultiplicityPtCut : public ScatterBase {
protected:
    std::string _event_label_base;
    std::vector<float> _ns;
public:
    JetMultiplicityPtCut(std::string event_label_base) {
        _event_label_base = event_label_base;
        _xs = {5, 15, 25, 50};
        _ys = {0, 0, 0, 0};
        _ns = {0, 0, 0, 0}; // need to keep track of moving averages
        std::stringstream ss;
        ss.str(std::string());
        ss << _event_label_base << "_JetMultiplicityPtCut.pdf";
        _outfile_name = ss.str();
        _x_label = "p_{T} Cut";
        _y_label = "Average CA Jet Multiplicity";
        _color = kBlack;
        _style = 20;
    }
    void Update(EventManager const* event_manager);
};

class StackedHistogramBase : public UpdatesOnEvent {
protected:
    const static unsigned int _root_name_length = 20;
    std::vector<char *> _root_names;
    std::vector<TH1F *> _hists;
    std::vector<Color_t> _colors;
    std::vector<Style_t> _styles;
    std::vector<std::string> _labels;

    std::string _title;
    std::string _outfile_name;
    
    unsigned int _canvas_x, _canvas_y;
    unsigned int _n_bins;
    float _low, _high;

    size_t _ticks;
    
    std::string _x_label;
    std::string _y_label;

    bool _normalized;

    StackedHistogramBase() {
        _outfile_name = "UNNAMED.PDF";

        _title = "UNTITLED";

        _x_label = "X_LABEL";
        _y_label = "Y_LABEL";

        _canvas_x = _canvas_y = 800;
        _n_bins = 50;
        _low = 0;

        _normalized = true;

        _ticks = 505;
    }

    ~StackedHistogramBase() {
        for (unsigned int iter = 0; iter < _hists.size(); iter++) {
            delete _hists[iter];
            delete[] _root_names[iter];
        }
    }
    
public:
    void Start(__attribute__((unused)) EventManager const* event_manager) {
        for (unsigned int iter = 0; iter < _colors.size(); iter++) {
            _root_names.push_back(new char[_root_name_length + 1]);
            gen_string(_root_names.at(iter), _root_name_length);
            
            _hists.push_back(new TH1F(_root_names.at(iter), _title.c_str(), _n_bins, _low, _high));
        }
    }

    virtual void Finish(EventManager const* event_manager);
};

class CorrelationBase : public UpdatesOnEvent {
protected:
    const static unsigned int _root_name_length = 20;
    std::string _title;
    std::string _outfile_name;
    char _root_name[_root_name_length + 1]; // 20 + null terminator

    unsigned int _canvas_x, _canvas_y;
    unsigned int _n_bins_x, _n_bins_y;
    float _x_low, _x_high, _y_low, _y_high;

    bool _correlation_in_title;

    std::string _x_label;
    std::string _y_label;

    TH2F *_hist;

    CorrelationBase() {
        // prevent root from forgetting our damn histograms, 
        // would be better to use an UUID though
        gen_string(_root_name, _root_name_length);

        // sensible defaults
        _outfile_name = "UNNAMED.PDF";
        _title = "UNTITLED";
        
        _x_label = "X_LABEL";
        _y_label = "Y_LABEL";

        _canvas_x = _canvas_y = 800;

        _correlation_in_title = true;

        _n_bins_x = _n_bins_y = 50;
        _x_low = _y_low = 0;
    }

    
    ~CorrelationBase() {
        delete _hist;
    }

public:
    void Start(__attribute__((unused)) EventManager const* event_mangaer) {
        _hist = new TH2F(_root_name, _title.c_str(), _n_bins_x, _x_low, _x_high, _n_bins_y,
                         _y_low, _y_high);
    }
    virtual void Finish(__attribute__((unused)) EventManager const* event_manager);
};

class WeightDistanceCorrelation : public CorrelationBase {
protected:
    std::string _event_label;
    std::string _alg_label;

public:    
    WeightDistanceCorrelation(std::string event_label, std::string alg_label) {
        _x_low = 0.05;
        _x_high = 1.0;
        _y_high = 2.0;

        _n_bins_x = 19;
        _n_bins_y = 20;
        
        _title = "Weight & Distance";
        _event_label = event_label;
        _alg_label = alg_label;

        _x_label = "Weight";
        _y_label = "Distance to Leading Jet";
        
        std::stringstream ss;
        ss.str(std::string());
        ss << "Weight_Distance_Correlation_" << event_label << "_" << alg_label << ".pdf"; 
        _outfile_name = ss.str();
    }
    void Update(EventManager const* event_manager);
};

class PtCorrelation : public CorrelationBase {
protected:
    std::string _event_label;
    std::string _alg_label;
    std::string _other_alg_label;
    
public:
    PtCorrelation(std::string event_label, std::string alg_label, std::string other_alg_label) {
        _event_label = event_label;
        _alg_label = alg_label;
        _other_alg_label = other_alg_label;

        std::stringstream ss;

        ss.str(std::string());
        ss << alg_label << " Jet pT [GeV]";
        _x_label = ss.str();

        ss.str(std::string());
        ss << other_alg_label << " Jet pT [GeV]";
        _y_label = ss.str();

        ss.str(std::string());
        ss << "Pythia8 " << event_label;
        _title = ss.str();

        ss.str(std::string());
        ss << "pT_correlation_" << event_label << "_" 
           << alg_label << "_" << other_alg_label << ".pdf";
        _outfile_name = ss.str();

        _x_high = _y_high = 500;
        _n_bins_x = _n_bins_y = 20;
    }

    void Update(EventManager const* event_manager);
};

class MassCorrelation : public CorrelationBase {
protected:
    std::string _event_label;
    std::string _alg_label;
    std::string _other_alg_label;
    
public:
    MassCorrelation(std::string event_label, std::string alg_label, std::string other_alg_label) {
        _event_label = event_label;
        _alg_label = alg_label;
        _other_alg_label = other_alg_label;

        std::stringstream ss;

        ss.str(std::string());
        ss << alg_label << " Jet Mass [GeV]";
        _x_label = ss.str();

        ss.str(std::string());
        ss << other_alg_label << " Jet Mass [GeV]";
        _y_label = ss.str();

        ss.str(std::string());
        ss << "Pythia8 " << event_label;
        _title = ss.str();

        ss.str(std::string());
        ss << "Mass_correlation_" << event_label << "_" 
           << alg_label << "_" << other_alg_label << ".pdf";
        _outfile_name = ss.str();

        _x_high = _y_high = 250;
        _n_bins_x = _n_bins_y = 20;
    }

    void Update(EventManager const* event_manager);
};

class SigmaJetSizeCorrelation : public CorrelationBase {
protected:    
    std::string _event_label;
    std::string _alg_label;
public:
    SigmaJetSizeCorrelation(std::string event_label, std::string alg_label) {
        _x_high = 1;
        _y_high = 1;
        _n_bins_x = _n_bins_y = 20;

        _event_label = event_label;
        _alg_label = alg_label;

        std::stringstream ss;
        
        ss.str(std::string());
        ss << "Leading " << _alg_label << " Jet Mass/Jet pT";
        _x_label = ss.str();

        ss.str(std::string());
        ss << "Learned sigma";
        _y_label = ss.str();
        
        ss.str(std::string());
        ss << "Pythia8 " << _event_label;
        _title = ss.str();
        
        ss.str(std::string());
        ss << "SigmaJetSizeCorrelation_" << _event_label << "_" << _alg_label << ".pdf";
        _outfile_name = ss.str();
    }

    void Update(EventManager const* event_manager);
};

class SigmaJetSizeCorrelationPoster :  public CorrelationBase {
public:
    SigmaJetSizeCorrelationPoster() {
        _x_high = _y_high = 1;
        _canvas_x = 1200;
        _canvas_y = 1200;
        _n_bins_x = _n_bins_y = 20;
        _x_label = "Leading anti-k_{t} Jet Mass/Jet p_{T}";
        _y_label = "Learned #sigma";

        _title = "Z'#rightarrow t#bar{t}";
        _outfile_name = "SigmaSizeCorrelationPoster.pdf";
    }
    void Update(EventManager const* event_manager);
};

class SigmaAreaCorrelation : public CorrelationBase {
protected:
    std::string _event_label;
    std::string _alg_branch;
public:
    SigmaAreaCorrelation(std::string event_label, std::string alg) {
        _event_label = event_label;
        
        if (alg == "antikt") {
            _x_low = 0;
            _x_high = 1.2;
            _y_low = 2.8;
            _y_high = 3.6;
        } else {
            _x_low = 0;
            _x_high = 1;
            _y_low = 0.5;
            _y_high = 4.5;
        }

        _canvas_x = 1200;
        _canvas_y = 1200;
        
        _n_bins_x = _n_bins_y = 20;
        _x_label = "Leading jet #sigma^{2}";
        _y_label = alg + " jet area";

        std::stringstream ss;
        ss.str(std::string());
        ss << _event_label << "_" << alg << "_SigmaAreaCorrelation.pdf";
        _outfile_name = ss.str();
        _title = _event_label;

        ss.str(std::string());
        ss << alg << "_area";
        _alg_branch = ss.str();
    }

    void Update(EventManager const* event_manager);
};

class FuzzyAntiktPtCorrelation : public CorrelationBase {
public:
    FuzzyAntiktPtCorrelation() {
        _x_high = _y_high = 500;
        _n_bins_x = _n_bins_y = 20;
        _x_label = "Leading anti-kt Jet pT [GeV]";
        _y_label = "Leading Fuzzy Jet pT [GeV]";
        _title = "Z' -> ttbar";
        _outfile_name = "mGMMc_antikt_pt_correlation_zprime.pdf";
    }

    void Update(EventManager const* event_manager);
};

class SkewHistogram : public StackedHistogramBase {
public:
    SkewHistogram() {
        _high = 0.5;
        _low = -0.8;
        _n_bins = 30;

        _x_label = "Fuzzy Jet Mass Skew [GeV]";
        _y_label = "";
        
        _title = "";
        _outfile_name = "SkewHistogram.pdf";

        _colors = {kBlue, kRed};
        _styles = {1, 1};
        _labels = {"Z'#rightarrow tt'", "QCD"};
        
        _ticks = 210;
    }
    void Update(EventManager const* event_manager);
};

class DeltaRHistogram : public StackedHistogramBase {
protected:
    std::string _alg_label;
public:
    DeltaRHistogram(std::string alg_label) {
        _alg_label = alg_label;
         
        _high = 6;
        _n_bins = 20;
        
        std::stringstream ss;
        ss.str(std::string());
        ss << "#Delta r between " << _alg_label << " jets";
        _x_label = ss.str();
        _y_label = "";

        _title = "";

        ss.str(std::string());
        ss << "DeltaRHistogram_" << _alg_label << ".pdf";
        _outfile_name = ss.str();

        _colors = {kBlack, kBlue, kRed};
        _styles = {1, 1, 1};
        _labels = {"Z'#rightarrow tt'", "W'#rightarrow qq'", "QCD"};

        _ticks = 505;
    }
    
    void Update(EventManager const* event_manager);
};

class RadiusComparisonHistogram : public StackedHistogramBase {
public:
    RadiusComparisonHistogram() {
        _high = 0.8;
        _n_bins = 20;

        _x_label = "Learned #sigma";
        _y_label = "";
        _title = "";

        _outfile_name = "RadiusComparison.pdf";

        _colors = {kBlack, kBlue, kRed};
        _styles = {1, 1, 1};
        _labels = {"Z'#rightarrow tt'","W'#rightarrow qq'", "QCD"};

        _ticks = 502;
    }
    
    void Update(EventManager const* event_manager);
};

class AverageRadiusComparisonHistogram : public StackedHistogramBase {
public:
    AverageRadiusComparisonHistogram() {
        _high = 1;
        _n_bins = 20;

        _x_label = "Learned sigma";
        _y_label = "";

        _outfile_name = "AverageRadiusComparison.pdf";

        _colors = {kBlack, kBlue, kRed};
        _styles = {1, 1, 1};
        _labels = {"Z' -> ttbar","W'->qq'", "QCD"};

        _ticks = 502;
    }
    
    void Update(EventManager const* event_manager);
};

struct CustomSortPair {
public:
    float signal;
    float background;
    CustomSortPair(float s, float b) 
        : signal(s), background(b) {}
    
    bool operator < (const CustomSortPair& rhs) const {
        //        if (background == 0) return true;
        //        if (rhs.background == 0) return false;
        return ((signal / background) > (rhs.signal / rhs.background));
    }
};



#endif
