#ifndef MVA_H
#define MVA_H

void MVATest();
void MVAEventJetTest();

void MVAEfficiency(TTree *signal_tree, TTree *background_tree,
                   float cut_low, float cut_high,
                   std::vector<std::string> const& variables,
                   std::vector<std::string> const& variable_labels,
                   std::vector<std::string> const& required_branches,
                   std::vector<float> &sig_eff,
                   std::vector<float> &back_eff);

void MVAPlots(TTree *signal_tree, TTree *background_tree,
              float cut_low, float cut_high,
              std::string signal_label,
              std::string background_label);

void getMVAData(std::string location, std::vector<float> &s_eff,
                std::vector<float> &b_eff);

void SinglePlot(std::string partial_path, std::vector<std::string> vars,
                std::map<std::string, std::vector<std::vector<float> > > eff_map);


#endif
