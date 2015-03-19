#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <TH2F.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>

#include <math.h>
#include <numeric>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <assert.h>

#include "Util.h"
#include "boost/foreach.hpp"

class EventManager;

const std::string plot_prefix = "/u/at/chstan/nfs/summer_2014/ForConrad/results/plots/event_jet_test/";

class UpdatesOnEvent {
public:
    UpdatesOnEvent() {}
    virtual ~UpdatesOnEvent() {}

    virtual void Start(__attribute__((unused)) EventManager const* event_manager) {}
    virtual void Update(EventManager const* event_manager) = 0;
    virtual void Finish(__attribute__((unused)) EventManager const* event_manager) {}
};

void gen_string(char* s, const unsigned int len);

class StackedEfficiencyHistogramGen : public UpdatesOnEvent {
protected:
    std::vector<std::vector<float>> _signals;
    std::vector<std::vector<float>> _backgrounds;

    std::vector<Color_t> _colors;
    std::vector<std::string> _labels;

    std::vector<std::vector<float>> _lows;
    std::vector<std::vector<float>> _highs;

    std::string _title;
    std::string _outfile_name;

    unsigned int _canvas_x, _canvas_y;
    std::vector<std::vector<unsigned int> > _n_points;

    size_t _ticks;

    std::string _x_label;
    std::string _y_label;
    std::string _y_inv_label;

    std::vector<unsigned int> _dimensions; // inferred

    StackedEfficiencyHistogramGen() {
        _outfile_name = "UNNAMED.PDF";

        _title = "UNTITLED";

        _x_label = "X_LABEL";
        _y_label = "Y_LABEL";
        _y_inv_label = "Y_LABEL";

        _canvas_x = _canvas_y = 800;
        _ticks = 505;
    }

    unsigned int determineBin(float val, float low, float high, unsigned int bin_count) {
        float n_val = (val - low) / (high - low);
        int bin = (int) floor(n_val * bin_count);
        if (bin < 0) return 0;
        if (bin >= (int) bin_count) return bin_count - 1;
        return (unsigned int) bin;
    }

    void fill(unsigned int attr_idx, std::vector<float> point, float weight, bool sig) {
        unsigned int raw_index = 0;
        unsigned int scalar = 1;
        unsigned int current_dimensions = _dimensions.at(attr_idx);

        std::vector<unsigned int> current_n_points = _n_points.at(attr_idx);

        assert(point.size() == current_dimensions);
        assert(point.size() == current_n_points.size());

        for (unsigned int dimension_iter = 0; dimension_iter < current_dimensions; dimension_iter++) {
            unsigned int index = determineBin(point.at(dimension_iter),
                                              _lows.at(attr_idx).at(dimension_iter),
                                              _highs.at(attr_idx).at(dimension_iter),
                                              current_n_points.at(dimension_iter));
            raw_index += index * scalar;
            scalar *= current_n_points.at(dimension_iter);
        }
        if (sig) {
            _signals.at(attr_idx).at(raw_index) += weight;
        } else {
            _backgrounds.at(attr_idx).at(raw_index) += weight;
        }
    }

    void fillSignal(unsigned int attr_idx, std::vector<float> point, float weight) {
        fill(attr_idx, point, weight, true);
    }
    void fillBackground(unsigned int attr_idx, std::vector<float> point, float weight) {
        fill(attr_idx, point, weight, false);
    }

public:
    void Start(__attribute__((unused)) EventManager const* event_manager) {
        unsigned int attribute_count = _colors.size();
        for(unsigned int attribute_iter = 0; attribute_iter < attribute_count; attribute_iter++) {
            _dimensions.push_back(_lows.at(attribute_iter).size());
            unsigned int current_dimensions = _dimensions.at(attribute_iter);

            std::vector<unsigned int> current_n_points = _n_points.at(attribute_iter);

            std::vector<float> v;
            std::vector<float> w;
            unsigned int count_max = 1;
            for(unsigned int iter = 0; iter < current_dimensions; iter++) {
                count_max *= current_n_points.at(iter);
            }
            for(unsigned int iter = 0; iter < count_max; iter++) {
                v.push_back(0);
                w.push_back(0);
            }
            _signals.push_back(v);
            _backgrounds.push_back(w);
        }
    }

    virtual void Finish(__attribute__((unused)) EventManager const* event_manager);
};

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

class AreaEfficiency : public StackedEfficiencyHistogramGen {
protected:
    std::string _event_signal;
    std::string _event_background;

    float _cut_low;
    float _cut_high;

public:
    AreaEfficiency (std::string event_signal,
                    std::string event_background,
                    float cut_low,
                    float cut_high) {
        _event_signal = event_signal;
        _event_background = event_background;

        _cut_low = cut_low;
        _cut_high = cut_high;

        _title = "";
        _x_label = _event_signal + " Efficiency";
        _y_label = "1 - " + _event_background + " Efficiency";
        _y_inv_label = "Inverse " + _event_background + " Efficiency";

        _n_points = {{600}, {600}, {600}, {40, 15}, {15, 40}, {15, 40}};

        _colors = {kRed, kBlue, kGreen, kBlack, kOrange, kViolet};
        _labels = {"Antikt trimmed 0.2 area",
                   "Antikt trimmed 0.3 area",
                   "Leading antikt trimmed 0.3 mass/pT",
                   "0.3 area, mass/pT",
                   "0.3 area, #sigma",
                   "mass/pT, #sigma"};

        _ticks = 405;
        _lows = {{0}, {0}, {0}, {0, 0}, {0, 0}, {0,0}};
        _highs = {{1}, {1.4}, {1}, {1.4, 1}, {1.4, 1.5}, {1, 1.5}};

        _outfile_name = "AreaEfficiency_" + event_signal +
            "_" + event_background + "_" + std::to_string((long long int) _cut_low) +
            "_to_" + std::to_string((long long int) _cut_high) + "pT.pdf";
    }
    void Update(EventManager const* event_manager);
};

class SigmaImprovementEfficiencyTau32 : public StackedEfficiencyHistogramGen {
protected:
    std::string _event_signal;
    std::string _event_background;

    float _cut_low;
    float _cut_high;

public:
    SigmaImprovementEfficiencyTau32 (std::string event_signal,
                                     std::string event_background,
                                     float cut_low,
                                     float cut_high) {
        _event_signal = event_signal;
        _event_background = event_background;

        _cut_low = cut_low;
        _cut_high = cut_high;

        _title = "";
        _x_label = _event_signal + " Efficiency";
        _y_label = "1 - " + _event_background + " Efficiency";
        _y_inv_label = "Inverse " + _event_background + " Efficiency";

        _n_points = {{400}, {400}, {400}, {10, 40}};

        _colors = {kRed, kBlue, kGreen, kBlack};
        _labels = {"#tau_{3} / #tau_{2}",
                   "#tau_{2} / #tau_{1}",
                   "#sigma",
                   "#tau_{3} / #tau_{2}, #sigma"};

        _ticks = 405;
        _lows = {{0}, {0}, {0}, {0, 0}};
        _highs = {{1}, {1}, {1.5}, {1, 1.5}};

        _outfile_name = "SigmaImprovementEfficiencyTau32_" + event_signal +
            "_" + event_background + "_" + std::to_string((long long int) _cut_low) +
            "_to_" + std::to_string((long long int) _cut_high) + "pT.pdf";
    }
    void Update(EventManager const* event_manager);
};

class SigmaImprovementEfficiencyTau21 : public StackedEfficiencyHistogramGen {
protected:
    std::string _event_signal;
    std::string _event_background;

    float _cut_low;
    float _cut_high;

public:
    SigmaImprovementEfficiencyTau21 (std::string event_signal,
                                     std::string event_background,
                                     float cut_low,
                                     float cut_high) {
        _event_signal = event_signal;
        _event_background = event_background;

        _cut_low = cut_low;
        _cut_high = cut_high;

        _title = "";
        _x_label = _event_signal + " Efficiency";
        _y_label = "1 - " + _event_background + " Efficiency";
        _y_inv_label = "Inverse " + _event_background + " Efficiency";

        _n_points = {{400}, {400}, {400}, {10, 40}};

        _colors = {kRed, kBlue, kGreen, kBlack};
        _labels = {"#tau_{3} / #tau_{2}",
                   "#tau_{2} / #tau_{1}",
                   "#sigma",
                   "#tau_{2} / #tau_{1}, #sigma"};

        _ticks = 405;
        _lows = {{0}, {0}, {0}, {0, 0}};
        _highs = {{1}, {1}, {1.5}, {1, 1.5}};

        _outfile_name = "SigmaImprovementEfficiencyTau21_" + event_signal +
            "_" + event_background + "_" + std::to_string((long long int) _cut_low) +
            "_to_" + std::to_string((long long int) _cut_high) + "pT.pdf";
    }
    void Update(EventManager const* event_manager);
};

class SigmaImprovementEfficiencyMultiTau : public StackedEfficiencyHistogramGen {
protected:
    std::string _event_signal;
    std::string _event_background;

    float _cut_low;
    float _cut_high;

public:
    SigmaImprovementEfficiencyMultiTau (std::string event_signal,
                                        std::string event_background,
                                        float cut_low,
                                        float cut_high) {
        _event_signal = event_signal;
        _event_background = event_background;

        _cut_low = cut_low;
        _cut_high = cut_high;

        _title = "";
        _x_label = _event_signal + " Efficiency";
        _y_label = "1 - " + _event_background + " Efficiency";
        _y_inv_label = "Inverse " + _event_background + " Efficiency";

        _n_points = {{20, 20}, {20, 20}, {400}, {10, 40},
                     {10, 40}, {10, 40}};

        _colors = {kRed, kBlue, kGreen, kBlack, kOrange, kViolet};
        _labels = {"#tau_{3}, #tau_{2}",
                   "#tau_{2}, #tau_{1}",
                   "#sigma",
                   "#tau_{3}, #sigma",
                   "#tau_{2}, #sigma",
                   "#tau_{1}, #sigma"};

        _ticks = 405;
        _lows = {{0, 0}, {0, 0}, {0}, {0, 0}, {0, 0}, {0, 0}};
        _highs = {{0.5, 0.8}, {0.8, 0.8}, {1.5}, {0.5, 1.5},
                  {0.8, 1.5}, {0.8, 1.5}};

        _outfile_name = "SigmaImprovementEfficiencyMultiTau_" + event_signal +
            "_" + event_background + "_" + std::to_string((long long int) _cut_low) +
            "_to_" + std::to_string((long long int) _cut_high) + "pT.pdf";
    }
    void Update(EventManager const* event_manager);
};

class EfficiencyGenTest : public StackedEfficiencyHistogramGen {
public:
    EfficiencyGenTest() {
        _title = "";
        _x_label = "W' Efficiency";
        _y_label = "1 - QCD Efficiency";
        _y_inv_label = "Inverse QCD Efficiency";

        _outfile_name = "EfficiencyGenTest.pdf";
        _n_points = {{500}, {500}, {500}};

        _colors = {kRed, kBlue, kBlack};
        _labels = {"Fuzzy #sigma", "Fuzzy Mass", "Anti-k_{T} Mass"};

        _ticks = 405;

        _lows = {{0}, {0}, {0}};
        _highs = {{1}, {400}, {400}};
    }
    void Update(EventManager const* event_manager);
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

class SigmaNSubjettinessEfficiencyPlot : public StackedEfficiencyHistogramBase {
protected:
    std::string _event_signal;
    std::string _event_background;

    float _cut_low, _cut_high;
public:
    SigmaNSubjettinessEfficiencyPlot(std::string event_signal, std::string event_background,
                                     float cut_low, float cut_high) {
        _event_signal = event_signal;
        _event_background = event_background;

        _cut_low = cut_low;
        _cut_high = cut_high;

        _title = "";

        _outfile_name = "SigmaNSubjettinessEfficiency_" + event_signal + "_" + event_background
            + "_" + std::to_string((long long int) cut_low) + "_to_" + std::to_string((long long int) cut_high)
            + "pT.pdf";
        _x_label = _event_signal + " Efficiency";
        _y_label = "1 - " + _event_background + " Efficiency";

        _n_points = 1000;
        _colors = {kRed, kBlue, kBlack, kGreen};
        _labels = {"#tau_1", "#tau_2", "#tau_2 / #tau_1", "#sigma"};
        _signs = {1, 1, 1, 1};

        _ticks = 405;
        _lows = {0, 0, 0, 0};
        _highs = {1, 1, 1, 1};
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
           << alg << "_" << (int) cut_low << "_to_" << (int) cut_high << "pT.pdf";
        _outfile_name = ss.str();

        ss.str(std::string());
        ss << _event_signal << " Efficiency";
        _x_label = ss.str();

        ss.str(std::string());
        ss << "1 - " << _event_background << " Efficiency";
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
        _xs = {5, 15, 25};
        _ys = {0, 0, 0};
        _ns = {0, 0, 0}; // need to keep track of moving averages
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

class JetMultiplicityProcess : public ScatterBase {
protected:
    util::ParameterSet _parameters;
    std::vector<float> _ns;

public:
    JetMultiplicityProcess(util::ParameterSet parameters) {
        _parameters = parameters; // all are used except process

        for (unsigned int i = 1; i <= util::processes.size(); i++) {
            _xs.push_back(i);
        }
        _ys = std::vector<float>(_xs.size(), 0);
        _ns = std::vector<float>(_xs.size(), 0);
        std::stringstream ss;
        ss.str(std::string());
        ss << "JetMultiplicityProcess_"
            //<< _parameters._process << "_"
           << _parameters._NPV << "mu_"
           << _parameters._EJW << "ejw_"
           << _parameters._EJO << "off_"
           << _parameters._PP << "pp_"
           << _parameters._TBS << "tbs_"
           << _parameters._seed_pT_cut << "cut"
           << ".pdf";
        _outfile_name = ss.str();
        _color = kBlack;
        _style = 20;
        _x_label = "Processes ";
        BOOST_FOREACH(auto process, util::processes) {
            _x_label += process;
            if (process != util::processes.back())
                _x_label += ", ";
        }
        _y_label = "Average Jet Multiplicity";
    }
    void Update(EventManager const* event_manager);
};

class JetMultiplicityNPV : public ScatterBase {
protected:
    util::ParameterSet _parameters;
    std::vector<float> _ns;

public:
    JetMultiplicityNPV(util::ParameterSet parameters) {
        _parameters = parameters; // all are used except process

        BOOST_FOREACH(auto NPV, util::NPVs) {
            if (_parameters._process == "background" && NPV == 0)
                continue;
            _xs.push_back(NPV);
        }
        _ys = std::vector<float>(_xs.size(), 0);
        _ns = std::vector<float>(_xs.size(), 0);
        std::stringstream ss;
        ss.str(std::string());
        ss << "JetMultiplicityNPV_"
           << _parameters._process << "_"
            //<< _parameters._NPV << "mu_"
           << _parameters._EJW << "ejw_"
           << _parameters._EJO << "off_"
           << _parameters._PP << "pp_"
           << _parameters._TBS << "tbs_"
           << _parameters._seed_pT_cut << "cut"
           << ".pdf";
        _outfile_name = ss.str();
        _color = kBlack;
        _style = 20;
        _x_label = "NPVs";
        _y_label = "Average Jet Multiplicity";
    }
    void Update(EventManager const* event_manager);
};

class JetMultiplicityEJW : public ScatterBase {
protected:
    util::ParameterSet _parameters;
    std::vector<float> _ns;

public:
    JetMultiplicityEJW(util::ParameterSet parameters) {
        _parameters = parameters; // all are used except process

        BOOST_FOREACH(auto EJW, util::EJWs) {
            _xs.push_back(EJW);
        }
        _ys = std::vector<float>(_xs.size(), 0);
        _ns = std::vector<float>(_xs.size(), 0);
        std::stringstream ss;
        ss.str(std::string());
        ss << "JetMultiplicityEJW_"
           << _parameters._process << "_"
           << _parameters._NPV << "mu_"
            //<< _parameters._EJW << "ejw_"
           << _parameters._EJO << "off_"
           << _parameters._PP << "pp_"
           << _parameters._TBS << "tbs_"
           << _parameters._seed_pT_cut << "cut"
           << ".pdf";
        _outfile_name = ss.str();
        _color = kBlack;
        _style = 20;
        _x_label = "EJWs";
        _y_label = "Average Jet Multiplicity";
    }
    void Update(EventManager const* event_manager);
};

class JetMultiplicityEJO : public ScatterBase {
protected:
    util::ParameterSet _parameters;
    std::vector<float> _ns;

public:
    JetMultiplicityEJO(util::ParameterSet parameters) {
        _parameters = parameters; // all are used except process

        BOOST_FOREACH(auto EJO, util::EJOs) {
            _xs.push_back(EJO);
        }
        _ys = std::vector<float>(_xs.size(), 0);
        _ns = std::vector<float>(_xs.size(), 0);
        std::stringstream ss;
        ss.str(std::string());
        ss << "JetMultiplicityEJO_"
           << _parameters._process << "_"
           << _parameters._NPV << "mu_"
           << _parameters._EJW << "ejw_"
            //<< _parameters._EJO << "off_"
           << _parameters._PP << "pp_"
           << _parameters._TBS << "tbs_"
           << _parameters._seed_pT_cut << "cut"
           << ".pdf";
        _outfile_name = ss.str();
        _color = kBlack;
        _style = 20;
        _x_label = "EJOs";
        _y_label = "Average Jet Multiplicity";
    }
    void Update(EventManager const* event_manager);
};

class JetMultiplicityPP : public ScatterBase {
protected:
    util::ParameterSet _parameters;
    std::vector<float> _ns;

public:
    JetMultiplicityPP(util::ParameterSet parameters) {
        _parameters = parameters; // all are used except process

        for (unsigned int i = 1; i <= util::PPs.size(); i++) {
            _xs.push_back(util::PPs.at(i));
        }
        _ys = std::vector<float>(_xs.size(), 0);
        _ns = std::vector<float>(_xs.size(), 0);
        std::stringstream ss;
        ss.str(std::string());
        ss << "JetMultiplicityPP_"
           << _parameters._process << "_"
           << _parameters._NPV << "mu_"
           << _parameters._EJW << "ejw_"
           << _parameters._EJO << "off_"
            //<< _parameters._PP << "pp_"
           << _parameters._TBS << "tbs_"
           << _parameters._seed_pT_cut << "cut"
           << ".pdf";
        _outfile_name = ss.str();
        _color = kBlack;
        _style = 20;
        _x_label = "PPs: 0 = NONE, 1 = SINGLE MERGER";
        _y_label = "Average Jet Multiplicity";
    }
    void Update(EventManager const* event_manager);
};

class JetMultiplicitySeedCut : public ScatterBase {
protected:
    util::ParameterSet _parameters;
    std::vector<float> _ns;
public:
    JetMultiplicitySeedCut(util::ParameterSet parameters) {
        _parameters = parameters; // all are used except process

        BOOST_FOREACH(auto seed_pT_cut, util::seed_pT_cuts) {
            _xs.push_back(seed_pT_cut);
        }
        _ys = std::vector<float>(_xs.size(), 0);
        _ns = std::vector<float>(_xs.size(), 0);
        std::stringstream ss;
        ss.str(std::string());
        ss << "JetMultiplicitySeedCut_"
           << _parameters._process << "_"
           << _parameters._NPV << "mu_"
           << _parameters._EJW << "ejw_"
           << _parameters._EJO << "off_"
           << _parameters._PP << "pp_"
           << _parameters._TBS << "tbs_"
            //<< _parameters._seed_pT_cut << "cut"
           << ".pdf";
        _outfile_name = ss.str();
        _color = kBlack;
        _style = 20;
        _x_label = "Seed pT Cuts";
        _y_label = "Average Jet Multiplicity";
    }
    void Update(EventManager const* event_manager);
};

class MeanFloatVarProcess : public ScatterBase {
protected:
    util::ParameterSet _parameters;
    std::vector<float> _ns;
    std::string _branch_label;

public:
    MeanFloatVarProcess(util::ParameterSet parameters, std::string branch_label) {
        _parameters = parameters; // all are used except process
        _branch_label = branch_label;

        for (unsigned int i = 1; i <= util::processes.size(); i++) {
            _xs.push_back(i);
        }
        _ys = std::vector<float>(_xs.size(), 0);
        _ns = std::vector<float>(_xs.size(), 0);
        std::stringstream ss;
        ss.str(std::string());
        ss << "MeanFloatVarProcess_"
           << _branch_label << "_"
            //<< _parameters._process << "_"
           << _parameters._NPV << "mu_"
           << _parameters._EJW << "ejw_"
           << _parameters._EJO << "off_"
           << _parameters._PP << "pp_"
           << _parameters._TBS << "tbs_"
           << _parameters._seed_pT_cut << "cut"
           << ".pdf";
        _outfile_name = ss.str();
        _color = kBlack;
        _style = 20;
        _x_label = "Processes ";
        BOOST_FOREACH(auto process, util::processes) {
            _x_label += process;
            if (process != util::processes.back())
                _x_label += ", ";
        }
        _y_label = "Average " + _branch_label;
    }
    void Update(EventManager const* event_manager);
};

class MeanFloatVarNPV : public ScatterBase {
protected:
    util::ParameterSet _parameters;
    std::vector<float> _ns;
    std::string _branch_label;

public:
    MeanFloatVarNPV(util::ParameterSet parameters, std::string branch_label) {
        _parameters = parameters; // all are used except process
        _branch_label = branch_label;

        BOOST_FOREACH(auto NPV, util::NPVs) {
            if (_parameters._process == "background" && NPV == 0)
                continue;
            _xs.push_back(NPV);
        }
        _ys = std::vector<float>(_xs.size(), 0);
        _ns = std::vector<float>(_xs.size(), 0);
        std::stringstream ss;
        ss.str(std::string());
        ss << "MeanFloatVarNPV_"
           << _branch_label << "_"
           << _parameters._process << "_"
            //<< _parameters._NPV << "mu_"
           << _parameters._EJW << "ejw_"
           << _parameters._EJO << "off_"
           << _parameters._PP << "pp_"
           << _parameters._TBS << "tbs_"
           << _parameters._seed_pT_cut << "cut"
           << ".pdf";
        _outfile_name = ss.str();
        _color = kBlack;
        _style = 20;
        _x_label = "NPVs";
        _y_label = "Average " + _branch_label;
    }
    void Update(EventManager const* event_manager);
};

class MeanFloatVarEJW : public ScatterBase {
protected:
    util::ParameterSet _parameters;
    std::vector<float> _ns;
    std::string _branch_label;

public:
    MeanFloatVarEJW(util::ParameterSet parameters, std::string branch_label) {
        _parameters = parameters; // all are used except process
        _branch_label = branch_label;

        BOOST_FOREACH(auto EJW, util::EJWs) {
            _xs.push_back(EJW);
        }
        _ys = std::vector<float>(_xs.size(), 0);
        _ns = std::vector<float>(_xs.size(), 0);
        std::stringstream ss;
        ss.str(std::string());
        ss << "MeanFloatVarEJW_"
           << _branch_label << "_"
           << _parameters._process << "_"
           << _parameters._NPV << "mu_"
            //<< _parameters._EJW << "ejw_"
           << _parameters._EJO << "off_"
           << _parameters._PP << "pp_"
           << _parameters._TBS << "tbs_"
           << _parameters._seed_pT_cut << "cut"
           << ".pdf";
        _outfile_name = ss.str();
        _color = kBlack;
        _style = 20;
        _x_label = "EJWs";
        _y_label = "Average " + _branch_label;
}
    void Update(EventManager const* event_manager);
};

class MeanFloatVarEJO : public ScatterBase {
protected:
    util::ParameterSet _parameters;
    std::vector<float> _ns;
    std::string _branch_label;

public:
    MeanFloatVarEJO(util::ParameterSet parameters, std::string branch_label) {
        _parameters = parameters; // all are used except process
        _branch_label = branch_label;

        BOOST_FOREACH(auto EJO, util::EJOs) {
            _xs.push_back(EJO);
        }
        _ys = std::vector<float>(_xs.size(), 0);
        _ns = std::vector<float>(_xs.size(), 0);
        std::stringstream ss;
        ss.str(std::string());
        ss << "MeanFloatVarEJO_"
           << _branch_label << "_"
           << _parameters._process << "_"
           << _parameters._NPV << "mu_"
           << _parameters._EJW << "ejw_"
            //<< _parameters._EJO << "off_"
           << _parameters._PP << "pp_"
           << _parameters._TBS << "tbs_"
           << _parameters._seed_pT_cut << "cut"
           << ".pdf";
        _outfile_name = ss.str();
        _color = kBlack;
        _style = 20;
        _x_label = "EJOs";
        _y_label = "Average " + _branch_label;
    }
    void Update(EventManager const* event_manager);
};

class MeanFloatVarPP : public ScatterBase {
protected:
    util::ParameterSet _parameters;
    std::vector<float> _ns;
    std::string _branch_label;

public:
    MeanFloatVarPP(util::ParameterSet parameters, std::string branch_label) {
        _parameters = parameters; // all are used except process
        _branch_label = branch_label;

        for (unsigned int i = 1; i <= util::PPs.size(); i++) {
            _xs.push_back(util::PPs.at(i));
        }
        _ys = std::vector<float>(_xs.size(), 0);
        _ns = std::vector<float>(_xs.size(), 0);
        std::stringstream ss;
        ss.str(std::string());
        ss << "MeanFloatVarPP_"
           << _branch_label << "_"
           << _parameters._process << "_"
           << _parameters._NPV << "mu_"
           << _parameters._EJW << "ejw_"
           << _parameters._EJO << "off_"
            //<< _parameters._PP << "pp_"
           << _parameters._TBS << "tbs_"
           << _parameters._seed_pT_cut << "cut"
           << ".pdf";
        _outfile_name = ss.str();
        _color = kBlack;
        _style = 20;
        _x_label = "PPs: 0 = NONE, 1 = SINGLE MERGER";
        _y_label = "Average " + _branch_label;
    }
    void Update(EventManager const* event_manager);
};

class MeanFloatVarSeedCut : public ScatterBase {
protected:
    util::ParameterSet _parameters;
    std::vector<float> _ns;
    std::string _branch_label;

public:
    MeanFloatVarSeedCut(util::ParameterSet parameters, std::string branch_label) {
        _parameters = parameters; // all are used except process
        _branch_label = branch_label;

        BOOST_FOREACH(auto seed_pT_cut, util::seed_pT_cuts) {
            _xs.push_back(seed_pT_cut);
        }
        _ys = std::vector<float>(_xs.size(), 0);
        _ns = std::vector<float>(_xs.size(), 0);
        std::stringstream ss;
        ss.str(std::string());
        ss << "MeanFloatVarSeedCut_"
           << _branch_label << "_"
           << _parameters._process << "_"
           << _parameters._NPV << "mu_"
           << _parameters._EJW << "ejw_"
           << _parameters._EJO << "off_"
           << _parameters._PP << "pp_"
           << _parameters._TBS << "tbs_"
            //<< _parameters._seed_pT_cut << "cut"
           << ".pdf";
        _outfile_name = ss.str();
        _color = kBlack;
        _style = 20;
        _x_label = "Seed pT Cuts";
        _y_label = "Average " + _branch_label;
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
    bool _log_x;

    StackedHistogramBase() {
        _log_x = false;
        _normalized = false;

        _outfile_name = "UNNAMED.PDF";

        _title = "UNTITLED";

        _x_label = "X_LABEL";
        _y_label = "Y_LABEL";

        _canvas_x = _canvas_y = 800;
        _n_bins = 50;
        _low = 0;

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

            TH1F *new_hist = new TH1F(_root_names.at(iter), _title.c_str(), _n_bins, _low, _high);
            if (_log_x) {
                // setup log x binning
                TAxis *axis = new_hist->GetXaxis();
                int bins = axis->GetNbins();

                Axis_t from = axis->GetXmin();
                Axis_t to = axis->GetXmax();
                Axis_t width = (to - from) / bins;
                Axis_t *new_bins = new Axis_t[bins + 1];

                for (int i = 0; i <= bins; i++) {
                    new_bins[i] = TMath::Power(10, from + i * width);
                }
                axis->Set(bins, new_bins);
                delete new_bins;
            }
            _hists.push_back(new_hist);
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

class GeneralCorrelation : public CorrelationBase {
protected:
    util::ParameterSet _param;
    std::string _var_a;
    std::string _var_b;
    float _pT_low;
    float _pT_high;

public:
    GeneralCorrelation(util::ParameterSet param,
                       std::string var_a, std::string var_b,
                       float a_low, float a_high, float b_low, float b_high,
                       float pT_low, float pT_high) {
        _pT_low = pT_low;
        _pT_high = pT_high;

        _canvas_x = 1200;
        _canvas_y = 1000;

        _param = param;
        _var_a = var_a;
        _var_b = var_b;
        _x_low = a_low;
        _y_low = b_low;
        _x_high = a_high;
        _y_high = b_high;

        _n_bins_x = 50;
        _n_bins_y = 50;

        std::stringstream ss;
        ss.str(std::string());
        _x_label = util::fancy_var_name(var_a);
        _y_label = util::fancy_var_name(var_b);
        ss << "Pythia 8  " << util::fancy_var_name(_param._process);
        _title = ss.str();

        ss.str(std::string());
        ss << "General_Correlation_"
           << static_cast<unsigned int>(pT_low) << "_"
           << static_cast<unsigned int>(pT_high) << "_"
           << util::sanitize(_var_a) << "_"
           << util::sanitize(_var_b) << "_"
           << _param.Print() << ".pdf";
        _outfile_name = ss.str();
    }
    void Update(EventManager const* event_manager);
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
            _alg_branch = "antikt_area";
        } else if (alg == "CA") {
            _x_low = 0;
            _x_high = 1;
            _y_low = 0.5;
            _y_high = 4.5;
            _alg_branch = "CA_area";
        } else if (alg == "antikt_trimmed_two") {
            _x_low = 0;
            _x_high = 1.2;
            _y_low = 0;
            _y_high = 2;
            _alg_branch = "antikt_area_trimmed_two";
        } else if (alg == "antikt_trimmed_three") {
            _x_low = 0;
            _x_high = 1.2;
            _y_low = 0;
            _y_high = 2;
            _alg_branch = "antikt_area_trimmed_three";
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

class GeneralHistogram : public StackedHistogramBase {
protected:
    std::vector<util::ParameterSet> _params;

    // distinguished_param can be any of
    // process, EJW, EJO, NPV, PP, TBS, cut
    std::string _distinguished_param;
    std::string _var;
    float _pT_low;
    float _pT_high;

public:
    GeneralHistogram(std::vector<util::ParameterSet> params, std::string distinguished_param,
                     std::string var, float low, float high, float pT_low, float pT_high,
                     const bool log_x = false) {
        _log_x = log_x;
        _distinguished_param = distinguished_param;
        _params = params;
        _var = var;
        _low = low;
        _high = high;
        _pT_low = pT_low;
        _normalized = true;
        _pT_high = pT_high;

        std::stringstream ss;
        ss.str(std::string());
        ss << "GeneralHistogram_" << _distinguished_param << "_"
           << _pT_low << "_" << _pT_high << "_"
           << util::sanitize(_var) << "_"
           << _params.at(0).PrintExceptDistinguished(_distinguished_param) << ".pdf";
        _outfile_name = ss.str();
        _y_label = "Normalized to Unity";
        _x_label = util::fancy_var_name(_var);
        _ticks = 505;

        _title = util::latexed_pT_range(_pT_low, _pT_high);

        // generate colors, styles, and labels
        _labels = {};
        for (unsigned int i = 0; i < _params.size(); i++) {
            _labels.push_back(util::fancy_distinguished_param_name(_params.at(i), _distinguished_param));
        }
        _styles = std::vector<Style_t>(_labels.size(), 1);
        _colors = util::get_colors(_labels.size());
    }

    void Update(EventManager const* event_manager);
};

class GeneralDiffgram : public StackedHistogramBase {
protected:
    std::vector<util::ParameterSet> _params;
    std::vector<util::ParameterSet> _mparams;

    // distinguished_param can be any of
    // process, EJW, EJO, NPV, PP, TBS, cut
    std::string _distinguished_param;

    std::string _var;
    float _pT_low;
    float _pT_high;

public:
    GeneralDiffgram(std::vector<util::ParameterSet> params,
                     std::vector<util::ParameterSet> mparams,
                     std::string distinguished_param,
                     std::string var,
                     float low, float high, float pT_low, float pT_high,
                     const bool log_x = false) {
        _log_x = log_x;
        _normalized = false;
        _distinguished_param = distinguished_param;
        _params = params;
        _mparams = mparams;
        _var = var;
        _low = low;
        _high = high;
        _pT_low = pT_low;
        _pT_high = pT_high;

        std::stringstream ss;
        ss.str(std::string());
        ss << "GeneralDiffgram_" << _distinguished_param << "_"
           << _pT_low << "_" << _pT_high << "_"
           << util::sanitize(_var) << "_"
           << _params.at(0).PrintExceptDistinguished(_distinguished_param)
           << "_minus_"
           << _mparams.at(0).PrintExceptDistinguished(_distinguished_param) << ".pdf";
        _outfile_name = ss.str();
        _y_label = "Arbitrary Units";
        _x_label = util::fancy_var_name(_var);
        _ticks = 505;

        _title = util::latexed_pT_range(_pT_low, _pT_high);

        // generate colors, styles, and labels
        _labels = {};
        for (unsigned int i = 0; i < _params.size(); i++) {
            std::string a = util::fancy_distinguished_param_name(_params.at(i), _distinguished_param);
            std::string b = util::fancy_distinguished_param_name(_mparams.at(i), _distinguished_param);
            _labels.push_back(a + " - " + b);
        }
        _styles = std::vector<Style_t>(_labels.size(), 1);
        _colors = util::get_colors(_labels.size());
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

class SigmaEventJetStrength : public StackedHistogramBase {
protected:
    std::string _process;
    int _NPV;
    // int _EJW;
    int _EJO;
    int _PP;
    int _TBS;
    int _seed_pT_cut;
    int _seed_noise;

public:
    SigmaEventJetStrength(std::string process, int NPV, int EJO,
                          int PP, int TBS, int seed_pT_cut, int seed_noise) {
        _process = process;
        _NPV = NPV;
        _EJO = EJO;
        _PP = PP;
        _TBS = TBS;
        _seed_pT_cut = seed_pT_cut;
        _seed_noise = seed_noise;

        _high = 1;
        _n_bins = 20;

        _x_label = "Learned #sigma";
        _y_label = "Arb. units";
        _title = "";

        std::stringstream ss;
        ss.str(std::string());
        ss << "SigmaEventJetStrength_" << _process << "_"
           << _NPV << "mu_"
           << _EJO << "off_"
           << _PP  << "PP_"
           << _TBS << "TBS_"
           << _seed_pT_cut << "cut_"
           << _seed_noise << "noise"
           << ".pdf";
        _outfile_name = ss.str();

        // use EJWs from util
        _labels = {};
        BOOST_FOREACH(auto EJW, util::EJWs) {
            float EJW_f = static_cast<float>(EJW)/1000;
            _labels.push_back(util::format_float(EJW_f, 3));
        }
        _colors = util::get_colors(_labels.size());
        _styles = std::vector<Style_t>(_labels.size(), 1);

        _ticks = 502;
    }

    void Update(EventManager const* event_manager);
};

class SigmaEventJetOffset : public StackedHistogramBase {
protected:
    std::string _process;
    int _NPV;
    int _EJW;
    //int _EJO;
    int _PP;
    int _TBS;
    int _seed_pT_cut;
    int _seed_noise;

public:
    SigmaEventJetOffset(std::string process, int NPV, int EJW,
                        int PP, int TBS, int seed_pT_cut, int seed_noise) {
        _process = process;
        _NPV = NPV;
        _EJW = EJW;
        //_EJO = EJO;
        _PP = PP;
        _TBS = TBS;
        _seed_pT_cut = seed_pT_cut;
        _seed_noise = seed_noise;

        _high = 1;
        _n_bins = 20;

        _x_label = "Learned #sigma";
        _y_label = "Arb. units";
        _title = "";

        std::stringstream ss;
        ss.str(std::string());
        ss << "SigmaEventJetOffset_" << _process << "_"
           << _NPV << "mu_"
           << _EJW << "ejw_"
            //<< _EJO << "off_"
           << _PP  << "PP_"
           << _TBS  << "TBS_"
           << _seed_pT_cut << "cut_"
           << _seed_noise << "noise"
           << ".pdf";
        _outfile_name = ss.str();

        // use EJOs from util
        _labels = {};
        BOOST_FOREACH(auto EJO, util::EJOs) {
            ss.str(std::string());
            ss << EJO;
            _labels.push_back(ss.str());
        }
        _colors = util::get_colors(_labels.size());
        _styles = std::vector<Style_t>(_labels.size(), 1);

        _ticks = 502;
    }

    void Update(EventManager const* event_manager);
};

class SigmaEventJetPP : public StackedHistogramBase {
protected:
    std::string _process;
    int _NPV;
    int _EJW;
    int _EJO;
    //int _PP;
    int _TBS;
    int _seed_pT_cut;
    int _seed_noise;

public:
    SigmaEventJetPP(std::string process, int NPV, int EJW,
                    int EJO, int TBS, int seed_pT_cut, int seed_noise) {
        _process = process;
        _NPV = NPV;
        _EJW = EJW;
        _EJO = EJO;
        //_PP = PP;
        _TBS = TBS;
        _seed_pT_cut = seed_pT_cut;
        _seed_noise = seed_noise;

        _high = 1;
        _n_bins = 20;

        _x_label = "Learned #sigma";
        _y_label = "Arb. units";
        _title = "";

        std::stringstream ss;
        ss.str(std::string());
        ss << "SigmaEventJetPP_" << _process << "_"
           << _NPV << "mu_"
           << _EJW << "ejw_"
           << _EJO << "off_"
            //<< _PP  << "PP_"
           << _TBS  << "TBS_"
           << _seed_pT_cut << "cut_"
           << _seed_noise << "noise"
           << ".pdf";
        _outfile_name = ss.str();

        // use PPs from util
        _labels = {};
        BOOST_FOREACH(auto PP, util::PPs) {
            ss.str(std::string());
            ss << PP;
            _labels.push_back(ss.str());
        }
        _colors = util::get_colors(_labels.size());
        _styles = std::vector<Style_t>(_labels.size(), 1);

        _ticks = 502;
    }

    void Update(EventManager const* event_manager);
};

class SigmaEventJetSeedCut : public StackedHistogramBase {
protected:
    std::string _process;
    int _NPV;
    int _EJW;
    int _EJO;
    int _PP;
    int _TBS;
    //int _seed_pT_cut;
    int _seed_noise;

public:
    SigmaEventJetSeedCut(std::string process, int NPV, int EJW,
                         int EJO, int PP, int TBS, int seed_noise) {
        _process = process;
        _NPV = NPV;
        _EJW = EJW;
        _EJO = EJO;
        _PP = PP;
        _TBS = TBS;
        //_seed_pT_cut = seed_pT_cut;
        _seed_noise = seed_noise;

        _high = 1;
        _n_bins = 20;

        _x_label = "Learned #sigma";
        _y_label = "Arb. units";
        _title = "";

        std::stringstream ss;
        ss.str(std::string());
        ss << "SigmaEventJetSeedCut_" << _process << "_"
           << _NPV << "mu_"
           << _EJW << "ejw_"
           << _EJO << "off_"
           << _PP  << "PP_"
           << _TBS  << "TBS_"
            //<< _seed_pT_cut << "cut"
           << _seed_noise << "noise"
           << ".pdf";
        _outfile_name = ss.str();

        // use cuts from util
        _labels = {};
        BOOST_FOREACH(auto seed_pT_cut, util::seed_pT_cuts) {
            ss.str(std::string());
            ss << seed_pT_cut;
            _labels.push_back(ss.str());
        }
        _colors = util::get_colors(_labels.size());
        _styles = std::vector<Style_t>(_labels.size(), 1);

        _ticks = 502;
    }

    void Update(EventManager const* event_manager);
};

class SigmaEventJetNPV : public StackedHistogramBase {
protected:
    std::string _process;
    //int _NPV;
    int _EJW;
    int _EJO;
    int _PP;
    int _TBS;
    int _seed_pT_cut;
    int _seed_noise;

public:
    SigmaEventJetNPV(std::string process, int EJW, int EJO,
                     int PP, int TBS, int seed_pT_cut, int seed_noise) {
        _process = process;
        _EJW = EJW;
        _EJO = EJO;
        _PP = PP;
        _TBS = TBS;
        _seed_pT_cut = seed_pT_cut;
        _seed_noise = seed_noise;

        _high = 1;
        _n_bins = 20;

        _x_label = "Learned #sigma";
        _y_label = "Arb. units";
        _title = "";

        std::stringstream ss;
        ss.str(std::string());
        ss << "SigmaEventJetStrength_" << _process << "_"
           << _EJW << "ejw_"
           << _EJO << "off_"
           << _PP  << "PP_"
           << _TBS  << "TBS_"
           << _seed_pT_cut << "cut_"
           << _seed_noise << "noise"
           << ".pdf";
        _outfile_name = ss.str();

        _labels = {};
        BOOST_FOREACH(auto NPV, util::NPVs) {
            if (_process == "background" && NPV == 0)
                continue;
            ss.str(std::string());
            ss << "#mu = " << NPV;
            _labels.push_back(ss.str());
        }
        _colors = util::get_colors(_labels.size());
        _styles = std::vector<Style_t>(_labels.size(), 1);

        _ticks = 502;
    }

    void Update(EventManager const* event_manager);
};

class RadiusPtSeedHistogram : public StackedHistogramBase {
protected:
    std::string _event_label_base;

public:
    RadiusPtSeedHistogram(std::string event_label_base) {
        _event_label_base = event_label_base;

        _high = 1;
        _n_bins = 50;

        _x_label = "Learned #sigma";
        _y_label = "Arb. units";
        _title = "";

        _outfile_name = "RadiusPtSeed_" + _event_label_base + ".pdf";
        _colors = {kBlack, kBlue, kRed};
        _styles = {1, 1, 1};
        _labels = {"5 GeV Cut", "15 GeV Cut", "25 GeV Cut"};

        _ticks = 502;
    }

    void Update(EventManager const* event_manager);
};

class AverageRadiusPtSeedHistogram : public StackedHistogramBase {
protected:
    std::string _event_label_base;

public:
    AverageRadiusPtSeedHistogram(std::string event_label_base) {
        _event_label_base = event_label_base;

        _high = 1;
        _n_bins = 50;

        _x_label = "Learned Event Average #sigma";
        _y_label = "Arb. units";
        _title = "";

        _outfile_name = "AverageRadiusPtSeed_" + _event_label_base + ".pdf";
        _colors = {kBlack, kBlue, kRed};
        _styles = {1, 1, 1};
        _labels = {"5 GeV Cut", "15 GeV Cut", "25 GeV Cut"};

        _ticks = 502;
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
        if (background == 0 && rhs.background == 0) {
            // both nan candidates, better do something consistent
            return signal > rhs.signal;
        }
        if (background == 0) return true;
        if (rhs.background == 0) return false;
        return ((signal / background) > (rhs.signal / rhs.background));
    }
};



#endif
