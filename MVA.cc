#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <map>
#include <iostream>
#include <assert.h>

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TROOT.h"
#include "TColor.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLegend.h"

#include "MVA.h"
#include "Histogram.h"
#include "AtlasUtils.h"
#include "AnalyzeFuzzyTools.h"


void MVATest() {
    // load the TMVA library
    TMVA::Tools::Instance();

    // instantiate TTrees and set up friendships
    // that the the MVA will have event weights
    std::string qcd_base_5 = "/u/at/chstan/nfs/summer_2014/ForConrad/files/150kevts_qcd_area/2014_10_14_11h23m44s/10s_0mu_0lw_0rec_5cut";
    std::string wprime_base_5 = "/u/at/chstan/nfs/summer_2014/ForConrad/files/150kevts_wprime_area/2014_10_15_00h40m54s/10s_0mu_0lw_0rec_5cut";
    std::string zprime_base_5 = "/u/at/chstan/nfs/summer_2014/ForConrad/files/150kevts_zprime_area/2014_10_13_22h17m12s/10s_0mu_0lw_0rec_5cut";

    std::string qcd_location_5 = qcd_base_5 + ".root";
    std::string wprime_location_5 = wprime_base_5 + ".root";
    std::string zprime_location_5 = zprime_base_5 + ".root";

    std::string qcd_addendum_location_5 = qcd_base_5 + "_w.root";
    std::string wprime_addendum_location_5 = wprime_base_5 + "_w.root";
    std::string zprime_addendum_location_5 = zprime_base_5 + "_w.root";

    std::map<std::string, TFile *> base_file_map;
    std::map<std::string, TFile *> addendum_file_map;
    std::map<std::string, TTree *> base_ttree_map;
    std::map<std::string, TTree *> addendum_ttree_map;

    base_file_map["qcd_5"] = new TFile(qcd_location_5.c_str());
    base_file_map["wprime_5"] = new TFile(wprime_location_5.c_str());
    base_file_map["zprime_5"] = new TFile(zprime_location_5.c_str());

    addendum_file_map["qcd_5"] = new TFile(qcd_addendum_location_5.c_str());
    addendum_file_map["wprime_5"] = new TFile(wprime_addendum_location_5.c_str());
    addendum_file_map["zprime_5"] = new TFile(zprime_addendum_location_5.c_str());

    for (auto it = base_file_map.begin(); it != base_file_map.end(); it++) {
        base_ttree_map[it->first] = (TTree *) it->second->Get("EventTree");
    }
    for (auto it = addendum_file_map.begin(); it != addendum_file_map.end(); it++) {
        addendum_ttree_map[it->first] = (TTree *) it->second->Get("EventTree");
    }

    // marry ttree data temporarily with a tree friendship
    for (auto it = base_ttree_map.begin(); it != base_ttree_map.end(); it++) {
        it->second->AddFriend(addendum_ttree_map[it->first]);
    }

    std::vector<std::string> processes = {"qcd", "wprime", "zprime"};
    std::vector<float> pT_cuts = {0, 300, 350, 400, 1000};
    for (unsigned int it = 0; it < processes.size(); it++) {
        for (unsigned int o_it = 0; o_it < it; o_it++) {
            std::string signal_process = processes.at(it) + "_5";
            std::string background_process = processes.at(o_it) + "_5";
            for (unsigned int pT_it = 0; pT_it < pT_cuts.size() - 1; pT_it++) {
                float pT_low = pT_cuts.at(pT_it);
                float pT_high = pT_cuts.at(pT_it + 1);
                std::cout << "Generating plots with " << background_process
                          << " background and "  << signal_process << " signal \n"
                          << "between pT GeV \t" << pT_low << " and " << pT_high
                          << std::endl;

                MVAPlots(base_ttree_map[signal_process],
                         base_ttree_map[background_process],
                         pT_low, pT_high,
                         signal_process,
                         background_process);

                std::cout << std::endl << std::endl;
            }
        }
    }

    // cleanup
    for (auto it = base_file_map.begin(); it != base_file_map.end(); it++) {
        delete it->second;
    }
    for (auto it = addendum_file_map.begin(); it != addendum_file_map.end(); it++) {
        delete it->second;
    }
}

void MVAPlots(TTree *signal_tree, TTree *background_tree,
              float cut_low, float cut_high,
              std::string signal_label,
              std::string background_label) {
    // map of variable set labels to efficiency vectors
    std::map<std::string, std::vector<std::vector<float> > > efficiencies_map;

    // first thing is to add all the needed data to this efficiencies map
    std::vector<std::vector<std::string>> var_set_labels = {
        std::vector<std::string>{"AkT 0.2 Area"},
        std::vector<std::string>{"AkT 0.3 Area"},
        std::vector<std::string>{"AkT 0.3 m/pT"},
        std::vector<std::string>{"AkT 0.3 m/pT", "AkT 0.3 Area"},
        std::vector<std::string>{"AkT 0.3 Area", "#sigma"},
        std::vector<std::string>{"AkT 0.3 m/pT", "#sigma"},
        std::vector<std::string>{"#tau_{3}/#tau_{2}"},
        std::vector<std::string>{"#tau_{2}/#tau_{1}"},
        std::vector<std::string>{"#sigma"},
        std::vector<std::string>{"#tau_{3}/#tau_{2}", "#sigma"},
        std::vector<std::string>{"#tau_{2}/#tau_{1}", "#sigma"},
        std::vector<std::string>{"#tau_{3}", "#tau_{2}"},
        std::vector<std::string>{"#tau_{2}", "#tau_{1}"},
        std::vector<std::string>{"#tau_{3}", "#sigma"},
        std::vector<std::string>{"#tau_{2}", "#sigma"},
        std::vector<std::string>{"#tau_{1}", "#sigma"}
    };
    std::vector<std::vector<std::string>> var_set_vars = {
        std::vector<std::string>{"antikt_area_trimmed_two"},
        std::vector<std::string>{"antikt_area_trimmed_three"},
        std::vector<std::string>{"antikt_m_trimmed_three/antikt_pt_trimmed_three"},
        std::vector<std::string>{"antikt_m_trimmed_three/antikt_pt_trimmed_three",
                                 "antikt_area_trimmed_three"},
        std::vector<std::string>{"antikt_area_trimmed_three", "mGMMc_r"},
        std::vector<std::string>{"antikt_m_trimmed_three/antikt_pt_trimmed_three",
                                 "mGMMc_r"},
        std::vector<std::string>{"antikt_nsubjettiness[2]/antikt_nsubjettiness[1]"},
        std::vector<std::string>{"antikt_nsubjettiness[1]/antikt_nsubjettiness[0]"},
        std::vector<std::string>{"mGMMc_r"},
        std::vector<std::string>{"antikt_nsubjettiness[2]/antikt_nsubjettiness[1]",
                                 "mGMMc_r"},
        std::vector<std::string>{"antikt_nsubjettiness[1]/antikt_nsubjettiness[0]",
                                 "mGMMc_r"},
        std::vector<std::string>{"antikt_nsubjettiness[2]", "antikt_nsubjettiness[1]"},
        std::vector<std::string>{"antikt_nsubjettiness[1]", "antikt_nsubjettiness[0]"},
        std::vector<std::string>{"antikt_nsubjettiness[2]", "mGMMc_r"},
        std::vector<std::string>{"antikt_nsubjettiness[1]", "mGMMc_r"},
        std::vector<std::string>{"antikt_nsubjettiness[0]", "mGMMc_r"}
    };
    std::vector<std::vector<std::string>> var_set_branches = {
        std::vector<std::string>{"antikt_area_trimmed_two"},
        std::vector<std::string>{"antikt_area_trimmed_three"},
        std::vector<std::string>{"antikt_m_trimmed_three", "antikt_pt_trimmed_three"},
        std::vector<std::string>{"antikt_area_trimmed_three",
                                 "antikt_m_trimmed_three",
                                 "antikt_pt_trimmed_three"},
        std::vector<std::string>{"antikt_area_trimmed_three", "mGMMc_r"},
        std::vector<std::string>{"antikt_m_trimmed_three", "antikt_pt_trimmed_three",
                                 "mGMMc_r"},
        std::vector<std::string>{"antikt_nsubjettiness"},
        std::vector<std::string>{"antikt_nsubjettiness"},
        std::vector<std::string>{"mGMMc_r"},
        std::vector<std::string>{"antikt_nsubjettiness", "mGMMc_r"},
        std::vector<std::string>{"antikt_nsubjettiness", "mGMMc_r"},
        std::vector<std::string>{"antikt_nsubjettiness"},
        std::vector<std::string>{"antikt_nsubjettiness"},
        std::vector<std::string>{"antikt_nsubjettiness", "mGMMc_r"},
        std::vector<std::string>{"antikt_nsubjettiness", "mGMMc_r"},
        std::vector<std::string>{"antikt_nsubjettiness", "mGMMc_r"}
    };

    assert(var_set_labels.size() == var_set_branches.size());
    assert(var_set_labels.size() == var_set_vars.size());

    for (unsigned int var_set_it = 0; var_set_it < var_set_labels.size(); var_set_it++) {
        std::vector<float> sig_eff;
        std::vector<float> back_eff;
        MVAEfficiency(signal_tree, background_tree,
                      cut_low, cut_high,
                      var_set_vars.at(var_set_it),
                      var_set_labels.at(var_set_it),
                      var_set_branches.at(var_set_it),
                      sig_eff, back_eff);

        // build the full label to index the results
        std::ostringstream ss;
        std::copy(var_set_labels.at(var_set_it).begin(),
                  var_set_labels.at(var_set_it).end() - 1,
                  std::ostream_iterator<std::string>(ss, ", "));
        ss << var_set_labels.at(var_set_it).back();

        // add it to the index
        std::vector<std::vector<float> > t;
        t.push_back(sig_eff);
        t.push_back(back_eff);
        efficiencies_map[ss.str()] = t;
    }

    for (auto it = efficiencies_map.begin(); it != efficiencies_map.end(); it++) {
        std::cout << it->first << std::endl;
    }

    // make plots
    std::vector<std::vector<std::string> > plot_var_set = {
        std::vector<std::string>{"AkT 0.2 Area", "AkT 0.3 Area",
                                 "AkT 0.3 m/pT", "AkT 0.3 m/pT, AkT 0.3 Area",
                                 "AkT 0.3 Area, #sigma",
                                 "AkT 0.3 m/pT, #sigma"},
        std::vector<std::string>{"#tau_{3}/#tau_{2}",
                                 "#tau_{2}/#tau_{1}",
                                 "#sigma",
                                 "#tau_{3}/#tau_{2}, #sigma"},
        std::vector<std::string>{"#tau_{3}/#tau_{2}",
                                 "#tau_{2}/#tau_{1}",
                                 "#sigma",
                                 "#tau_{2}/#tau_{1}, #sigma"},
        std::vector<std::string>{"#tau_{3}/#tau_{2}",
                                 "#tau_{2}/#tau_{1}",
                                 "#sigma",
                                 "#tau_{3}, #sigma",
                                 "#tau_{2}, #sigma",
                                 "#tau_{1}, #sigma"}
    };
    std::vector<std::string> plot_titles = {"area", "tau32", "tau21", "multitau"};
    for (unsigned int plot_it = 0; plot_it < plot_var_set.size(); plot_it++) {
        std::ostringstream ss;
        ss << "/u/at/chstan/nfs/summer_2014/ForConrad/results/plots/MVA/"
           << plot_titles.at(plot_it) << "_" << signal_label << "_"
           << background_label << "_" << (int) cut_low << "_"
           << (int) cut_high;

        SinglePlot(ss.str(), plot_var_set.at(plot_it), efficiencies_map);
    }
}

void SinglePlot(std::string partial_path, std::vector<std::string> vars,
                std::map<std::string, std::vector<std::vector<float> > > eff_map) {
    SetupATLASStyle();

    TMultiGraph *multi = new TMultiGraph();
    TMultiGraph *inv_multi = new TMultiGraph();

    std::vector<float> random_x = {0, 1};
    std::vector<float> random_y = {1, 0};

    TGraph *random = new TGraph(2, &random_x.at(0), &random_y.at(0));
    random->SetLineColor(kGray);
    multi->Add(random, "c");

    std::vector<float> inv_random_x;
    std::vector<float> inv_random_y;
    for (unsigned int i = 1; i <= 100; i++) {
        inv_random_x.push_back(1.*i/100);
        inv_random_y.push_back(1./inv_random_x.back());
    }
    TGraph *inv_random = new TGraph(inv_random_x.size(), &inv_random_x.at(0),
                                    &inv_random_y.at(0));
    random->SetLineColor(kGray);
    inv_multi->Add(inv_random, "c");

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

    std::vector<Color_t> available_colors = {kBlack, kBlue,
                                             kRed, kGreen, kOrange,
                                             kViolet, kAzure};
    // add graphs according to variables
    for (unsigned int var_iter = 0; var_iter < vars.size(); var_iter++) {
        std::vector<float> vec_signal;
        std::vector<float> vec_background_rej;
        std::vector<float> vec_inv_background_eff;

        std::string label = vars.at(var_iter);
        std::vector<std::vector<float> > data = eff_map[label];
        vec_signal.push_back(0);
        vec_background_rej.push_back(1);
        for (unsigned int data_iter = 0; data_iter < data.at(0).size(); data_iter++) {
            float s_e = data.at(0).at(data_iter);
            float b_e = data.at(1).at(data_iter);
            vec_signal.push_back(s_e);
            vec_background_rej.push_back(1.-b_e);
            vec_inv_background_eff.push_back(1./b_e);
        }
        vec_signal.push_back(1);
        vec_background_rej.push_back(0);
        vec_inv_background_eff.push_back(1);

        TGraph *efficiency_graph = new TGraph(vec_signal.size(),
                                              &vec_signal.at(0),
                                              &vec_background_rej.at(0));
        TGraph *inv_efficiency_graph = new TGraph(vec_inv_background_eff.size(),
                                                  &vec_signal.at(vec_signal.size() -
                                                                 vec_inv_background_eff.size()),
                                                  &vec_inv_background_eff.at(0));

        efficiency_graph->SetMarkerColor(available_colors[var_iter]);
        efficiency_graph->SetLineColor(available_colors[var_iter]);
        efficiency_graph->SetMarkerStyle(20); // circle

        inv_efficiency_graph->SetMarkerColor(available_colors[var_iter]);
        inv_efficiency_graph->SetLineColor(available_colors[var_iter]);
        inv_efficiency_graph->SetMarkerStyle(20); // circle

        legend->AddEntry(efficiency_graph, label.c_str(), "l");
        inv_legend->AddEntry(inv_efficiency_graph,
                             label.c_str(), "l");

        multi->Add(efficiency_graph);
        inv_multi->Add(inv_efficiency_graph);
    }

    TCanvas canvas("temporary", "", 0, 0, 1000, 1000);
    canvas.cd();

    multi->Draw("ac");
    multi->GetXaxis()->SetTitle("Signal Efficiency");
    multi->GetYaxis()->SetTitle("Background Rejection");

    multi->GetXaxis()->SetLimits(0, 1);
    multi->GetYaxis()->SetLimits(0, 1);

    multi->GetXaxis()->SetNdivisions(405);

    canvas.Update();

    legend->Draw();

    DrawAtlasLabel("", 0.2, 0.48);

    std::string out = partial_path + ".pdf";
    canvas.Print(out.c_str());

    // multi owns the efficiency graphs
    delete multi;

    TCanvas inv_canvas("temporary2", "", 0,0, 1000, 1000);
    inv_canvas.cd();
    inv_canvas.SetLogy();

    inv_multi->Draw("ac");
    inv_multi->GetXaxis()->SetTitle("Signal Efficiency");
    inv_multi->GetYaxis()->SetTitle("Inverse Background Efficiency");

    inv_multi->GetXaxis()->SetLimits(0, 1);
    inv_multi->GetYaxis()->SetLimits(0, 1);

    inv_multi->GetXaxis()->SetNdivisions(405);

    inv_canvas.Update();

    inv_legend->Draw();

    DrawAtlasLabel("", 0.2, 0.48);

    out = partial_path + "_inv.pdf";

    inv_canvas.Print(out.c_str());

    delete inv_multi;
}

void MVAEfficiency(TTree *signal_tree, TTree *background_tree,
                   float cut_low, float cut_high,
                   std::vector<std::string> const& variables,
                   std::vector<std::string> const& variable_labels,
                   std::vector<std::string> const& required_branches,
                   std::vector<float> &sig_eff,
                   std::vector<float> &back_eff) {
    std::string output_file_location = "/u/at/chstan/nfs/summer_2014/ForConrad/results/TMVA/TMVA.root";

    // disable all unnecessary branches
    signal_tree->SetBranchStatus("*", 0);
    background_tree->SetBranchStatus("*", 0);

    signal_tree->SetBranchStatus("antikt_pt", 1);
    background_tree->SetBranchStatus("antikt_pt", 1);

    signal_tree->SetBranchStatus("pT_reweight", 1);
    background_tree->SetBranchStatus("pT_reweight", 1);
    for (unsigned int it = 0; it < required_branches.size(); it++) {
        signal_tree->SetBranchStatus(required_branches.at(it).c_str(), 1);
        background_tree->SetBranchStatus(required_branches.at(it).c_str(), 1);
    }

    std::map<std::string, int> Methods;
    // --- Cut optimisation
    Methods["Cuts"]            = 0;
    Methods["CutsD"]           = 0;
    Methods["CutsPCA"]         = 0;
    Methods["CutsGA"]          = 0;
    Methods["CutsSA"]          = 0;

    // --- 1-dimensional likelihood ("naive Bayes estimator")
    Methods["Likelihood"]      = 0;
    Methods["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
    Methods["LikelihoodPCA"]   = 0; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
    Methods["LikelihoodKDE"]   = 0;
    Methods["LikelihoodMIX"]   = 0;

    // --- Mutidimensional likelihood and Nearest-Neighbour methods
    Methods["PDERS"]           = 0;
    Methods["PDERSD"]          = 0;
    Methods["PDERSPCA"]        = 0;
    Methods["PDEFoam"]         = 0;
    Methods["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
    Methods["KNN"]             = 1; // k-nearest neighbour method

    // --- Linear Discriminant Analysis
    Methods["LD"]              = 0; // Linear Discriminant identical to Fisher
    Methods["Fisher"]          = 0;
    Methods["FisherG"]         = 0;
    Methods["BoostedFisher"]   = 0; // uses generalised MVA method boosting
    Methods["HMatrix"]         = 0;

    // --- Function Discriminant analysis
    Methods["FDA_GA"]          = 0; // minimisation of user-defined function using Genetics Algorithm
    Methods["FDA_SA"]          = 0;
    Methods["FDA_MC"]          = 0;
    Methods["FDA_MT"]          = 0;
    Methods["FDA_GAMT"]        = 0;
    Methods["FDA_MCMT"]        = 0;

    // --- Neural Networks (all are feed-forward Multilayer Perceptrons)
    Methods["MLP"]             = 0; // Recommended ANN
    Methods["MLPBFGS"]         = 0; // Recommended ANN with optional training method
    Methods["MLPBNN"]          = 0; // Recommended ANN with BFGS training method and bayesian regulator
    Methods["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
    Methods["TMlpANN"]         = 0; // ROOT's own ANN

    // --- Support Vector Machine
    Methods["SVM"]             = 0;

    // --- Boosted Decision Trees
    Methods["BDT"]             = 0; // uses Adaptive Boost
    Methods["BDTG"]            = 0; // uses Gradient Boost
    Methods["BDTB"]            = 0; // uses Bagging
    Methods["BDTD"]            = 0; // decorrelation + Adaptive Boost

    // --- Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
    Methods["RuleFit"]         = 0;

    TFile *outputFile = TFile::Open(output_file_location.c_str(), "RECREATE");
    const char *factory_options = "!V:"
                                  "Silent:"
                                  "Color:"
                                  "DrawProgressBar:"
                                  "Transformations=I:"
                                  "AnalysisType=Classification";

    TMVA::Factory *factory = new TMVA::Factory("TMVAClassification", outputFile, factory_options);

    assert(variables.size() == variable_labels.size());
    for (unsigned int it = 0; it < variables.size(); it++) {
        factory->AddVariable(variables.at(it).c_str(),
                             variable_labels.at(it).c_str(), "", 'F');
    }

    factory->AddSignalTree(signal_tree, 1.);
    factory->AddBackgroundTree(background_tree, 1.);

    factory->SetWeightExpression("pT_reweight");

    std::stringstream ss;
    ss.str(std::string());
    ss << "antikt_pt >= " << cut_low << " && antikt_pt < " << cut_high;
    std::string cut_string = ss.str();
    TCut pT_cut = TCut(cut_string.c_str());
    const char *preparation_options = "nTrain_Signal=0:"
                                      "nTrain_Background=0:"
                                      "nTest_Signal=0:"
                                      "nTest_Background=0:"
                                      "SplitMode=Random:"
                                      "NormMode=NumEvents:"
                                      "!V";
    factory->PrepareTrainingAndTestTree(pT_cut, pT_cut,
                                        preparation_options);

    // Book MVA Methods

    const char *PDEFoam_opts = "VolFrac=0.0333:"
                               "nActiveCells=500:"
                               "nSampl=1000:";
    if (Methods["PDEFoam"])
        factory->BookMethod(TMVA::Types::kPDEFoam, "PDEFoam", PDEFoam_opts);

    const char *KNN_opts = "nkNN=20";
    if (Methods["KNN"])
        factory->BookMethod( TMVA::Types::kKNN, "KNN", KNN_opts);

    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();

    outputFile->Close();

    getMVAData(output_file_location, sig_eff, back_eff);

    delete factory;
    delete outputFile; // ??
}

void getMVAData(std::string location,
                std::vector<float> &s_eff,
                std::vector<float> &b_eff) {
    s_eff.clear();
    b_eff.clear();

    TFile *fMVA = TFile::Open(location.c_str());
    fMVA->cd("Method_KNN");
    gDirectory->cd("KNN");

    TH1F *hist = (TH1F *) gDirectory->Get("MVA_KNN_effBvsS");

    Int_t n_bins = hist->GetNbinsX();
    for (Int_t bin_iter = 2; bin_iter < n_bins + 1; bin_iter++) {
        s_eff.push_back(hist->GetBinCenter(bin_iter));
        b_eff.push_back(hist->GetBinContent(bin_iter));
    }
}
