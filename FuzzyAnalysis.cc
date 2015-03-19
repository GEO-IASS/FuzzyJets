#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <ctime>
#include <numeric>
#include <float.h>
#include <set>
#include <assert.h>
#include <exception>
#include <time.h>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

#include "myFastJetBase.h"
#include "fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/ClusterSequenceActiveArea.hh"
#include "fastjet/RangeDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/Selector.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/Njettiness.hh"

#include "Pythia8/Pythia.h"

#include "FuzzyAnalysis.h"
#include "FuzzyTools.h"
#include "FuzzyUtil.h"
#include "ROOTConf.h"


#include "boost/foreach.hpp"

#ifdef WITHROOT
#include "TFile.h"
#include "THStack.h"
#include "TColor.h"
#include "TROOT.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TMath.h"
#include "TVector3.h"
#include "TMatrix.h"
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include "TStyle.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"
#endif

// Privately separate the logic of different analysis modes from using them
namespace {
    float randb() {
        float random = ((float) rand()) / (float) RAND_MAX;
        return (random >= 0.5) ? 1 : -1;
    }
    float randf(float low, float high) {
        float random = ((float) rand()) / (float) RAND_MAX;
        float d = high - low;
        return low + d*random;
    }

    void multi_with_errors(std::vector<float> const &x, std::vector<float> const &y1,
                           std::vector<float> const &ey1,
                           std::vector<float> const &y2,
                           std::vector<float> const &ey2,
                           std::string var, std::string label,
                           std::string prefix) {
        TCanvas *c = new TCanvas("scatter_canv_mwe", "", 1100, 1000);
        c->cd();

        TMultiGraph *mgr = new TMultiGraph();
        TGraphErrors *gr1 = new TGraphErrors(x.size(), &x.at(0), &y1.at(0), 0, &ey1.at(0));
        TGraphErrors *gr2 = new TGraphErrors(x.size(), &x.at(0), &y2.at(0), 0, &ey2.at(0));

        // use colors 5 and 6 because they are crap
        float trans = 0.2;
        TColor *colA = gROOT->GetColor(5);
        colA->SetRGB(0, 0, 1);
        colA->SetAlpha(trans);

        gr1->SetMarkerColor(kBlue);
        gr1->SetLineColor(kBlue);
        gr1->SetFillColor(5);
        gr1->SetMarkerStyle(22);

        mgr->Add(gr1);

        TColor *colB = gROOT->GetColor(6);
        colB->SetRGB(1, 0, 0);
        colB->SetAlpha(trans);

        gr2->SetMarkerColor(kRed);
        gr2->SetLineColor(kRed);
        gr2->SetFillColor(6);
        gr2->SetMarkerStyle(23);


        mgr->Add(gr2);

        if (var == "Sigma") {
            mgr->SetMinimum(0.4);
            mgr->SetMaximum(1.1);
        }

        std::stringstream ss;
        ss.str(std::string());
        ss << var << " for a Fixed Event;NPV;" << var << " + Standard Deviation";
        mgr->SetTitle(ss.str().c_str());
        mgr->Draw("a3");
        mgr->Draw("xp");

        gPad->Modified();
        //mgr->GetXaxis()->SetLimits(-2., 80.1);
        mgr->GetYaxis()->SetTitleOffset(1.4);

        c->Update();

        ss.str(std::string());
        ss << prefix << label << "DepOnMu.pdf";
        c->Print(ss.str().c_str(), "pdf");
        delete mgr;
        delete c;
    }

    void calc_stats(std::vector<float> const& v, float &mean, float &stdev) {
        mean = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
        float sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
        stdev = std::sqrt(sq_sum / v.size() - mean * mean);
    }

    vecPseudoJet randomSubvector(vecPseudoJet const& origin, unsigned int how_many) {
        assert(origin.size() >= how_many);
        vecPseudoJet out;

        std::vector<unsigned int> indices;
        for(unsigned int iter = 0; iter < origin.size(); iter++) {
            indices.push_back(iter);
        }
        std::random_shuffle(indices.begin(), indices.end());

        for(unsigned int iter = 0; iter < how_many; iter++) {
            unsigned int index = indices.at(iter);
            out.push_back(origin[index]);
        }

        for(unsigned int iter = 0; iter < out.size(); iter++) {
            float rap = out.at(iter).rapidity() + randb() * randf(0.8, 1);
            float phi = out.at(iter).phi() + randb() * randf(0.8, 1);
            out.at(iter).reset_PtYPhiM(1.0, rap, phi, 0);

        }

        return out;
    }

    void WeightDistribution(__attribute__((unused)) vector<vector<double> > const& weights,
                            __attribute__((unused)) int which,
                            __attribute__((unused)) std::string out,
                            __attribute__((unused)) int iter) {
#ifdef WITHROOT
        TTree aux(TString::Format("WT_%s_%d", out.c_str(), iter), "");
        double w;
        int p_count = weights.size();
        aux.Branch("w", &w, "w/D");
        for (int p_idx = 0; p_idx < p_count; p_idx++) {
            w = weights[p_idx][which];
            aux.Fill();
        }
        aux.Write();
#endif
    }

    // is the reference particle given by index p_idx clustered at all?
    bool isClustered(vector<vector<double> > const& weights,
                     int p_idx) {
        bool clustered = false;
        int n_clusters = weights[0].size();
        for (int j = 0; j < n_clusters; j++) {
            double w = weights[p_idx][j];
            if (!std::isnan(w) && w > 0) {
                clustered = true;
            }
        }
        return clustered;
    }

    // what is the index of the fuzzy jet that the particle belongs to?
    int belongsTo(vector<vector<double> > const& weights,
                  int p_idx) {
        double w_max = -1;
        int best_idx = -1;
        const int n_clusters = weights[0].size();
        for (int j = 0; j < n_clusters; j++) {
            double w = weights[p_idx][j];
            if(!std::isnan(w) && w > 0 && w > w_max) {
                w_max = w;
                best_idx = j;
            }
        }
        return best_idx;
    }

    // Compute the mass of a Fuzzy Jet which comes from pileup
    // Please note that this is somewhat ill defined, due to how particles add
    double JetPuMassHard(vecPseudoJet const& particles,
                         vector<vector<double> > const& weights,
                         int jet_idx) {
        fastjet::PseudoJet pu_jet;
        const int n_particles = particles.size();
        for (int i = 0; i < n_particles; i++) {
            if (belongsTo(weights, i) == jet_idx &&
                particles[i].user_info<MyUserInfo>().isPU()) {
                pu_jet += particles[i];
            }
        }
        return pu_jet.m();
    }

    // Compute the soft mass due to pileup
    // Please note that this is somewhat ill defined, due to how particles add
    double JetPuMassSoft(vecPseudoJet const& particles,
                         vector<vector<double> > const& weights,
                         int jet_idx) {
        fastjet::PseudoJet pu_jet;
        const int n_particles = particles.size();
        for (int i = 0; i < n_particles; i++) {
            if (particles[i].user_info<MyUserInfo>().isPU()) {
                pu_jet += particles[i] * weights[i][jet_idx];
            }
        }
        return pu_jet.m();
    }

    // Compute the mass of a fastjet jet which is due to pileup
    // Please note that this is somewhat ill defined, due to how particles add
    double JetPuMassFastjet(fastjet::PseudoJet const& jet) {
        const vecPseudoJet c = jet.constituents();
        const int n_particles = c.size();
        fastjet::PseudoJet my_jet;
        for (int i = 0; i < n_particles; i++) {
            if(c[i].user_info<MyUserInfo>().isPU()) {
                my_jet += c[i];
            }
        }
        return my_jet.m();
    }

    // Compute the fraction of particles in a jet which are pileup
    double JetPuFracHard(vecPseudoJet const& particles,
                         vector<vector<double> > const& weights,
                         int jet_idx) {
        const int n_particles = particles.size();
        int clustered_particles = 0;
        int clustered_pu = 0;
        for (int i = 0; i < n_particles; i++) {
            if (belongsTo(weights, i) == jet_idx) {
                clustered_particles++;
                if(particles[i].user_info<MyUserInfo>().isPU()) {
                    clustered_pu++;
                }
            }
        }
        return (double)clustered_pu / clustered_particles;
    }

    // Compute the fraction of particles in a jet which are pileup by soft assignment
    double JetPuFracSoft(vecPseudoJet const& particles,
                         vector<vector<double> > const& weights,
                         int jet_idx) {
        const int n_particles = particles.size();
        const int n_clusters = weights[0].size();

        assert(n_particles == (int) weights.size());
        assert(jet_idx >= 0);
        assert(jet_idx < n_clusters);

        double clustered_particles = 0;
        double clustered_pu = 0;
        for (int i = 0; i < n_particles; i++) {
            if (std::isnan(weights[i][jet_idx])) continue;
            clustered_particles += weights[i][jet_idx];
            if (particles[i].user_info<MyUserInfo>().isPU()) {
                clustered_pu += weights[i][jet_idx];
            }
        }
        return clustered_pu / clustered_particles;
    }

    // Compute the fraction of particles in a jet which are pileup for standard jets
    double JetPuFracFastjet(fastjet::PseudoJet const& jet) {
        const vecPseudoJet c = jet.constituents();
        const int n_particles = c.size();
        int pu = 0;
        for (int i = 0; i < n_particles; i++) {
            if (c[i].user_info<MyUserInfo>().isPU()) {
                pu++;
            }
        }
        return (double)pu / n_particles;
    }

    // Compute (#clustered pileup particles) / (#pileup particles)
    double ClusteredPileupFrac(vecPseudoJet const& particles,
                               vector<vector<double> > const& weights) {
        int pileup = 0;
        int clustered_pileup = 0;
        for(unsigned int i = 0; i < particles.size(); i++) {
            fastjet::PseudoJet p = particles[i];
            if (p.user_info<MyUserInfo>().isPU()) {
                pileup++;
                if(isClustered(weights, i))
                    clustered_pileup++;
            }
        }
        return (double)clustered_pileup / pileup;
    }

    // Compute (#unclustered pileup particles) / (#unclustered particles)
    double UnclusteredPileupComposition(vecPseudoJet particles,
                                        vector<vector<double> > const& weights) {
        int unclustered_particles = 0;
        int unclustered_pileup = 0;
        for(unsigned int i = 0; i < particles.size(); i++) {
            if(!isClustered(weights, i)) {
                unclustered_particles++;
                fastjet::PseudoJet p = particles[i];
                if(p.user_info<MyUserInfo>().isPU()) {
                    unclustered_pileup++;
                }
            }
        }
        return (double)unclustered_pileup / unclustered_particles;
    }

    void DoFastJetFinding(vecPseudoJet particles,
                          fastjet::JetDefinition *p_jet_def,
                          double pT_min,
                          vecPseudoJet& jets) {
        fastjet::ClusterSequence cs_large_r_ca(particles, *p_jet_def);
        jets = fastjet::sorted_by_pt(cs_large_r_ca.inclusive_jets(pT_min));
    }

    struct SortablePt {
    public:
        float pT;
        unsigned int index;

        SortablePt(int input_index, float input_pT) : pT(input_pT), index(input_index) {}

        bool operator < (const SortablePt& rhs) const {
            return (pT < rhs.pT);
        }
    };

    void FindLeadingJet(vecPseudoJet& particles,
                        vecPseudoJet& jets,
                        vector<vector<double> >& particle_weights,
                        FuzzyTools *tool,
                        vector<int>& indices,
                        double& pT) {

        pT = -1;
        indices.clear();
        vector<SortablePt> pTs;
        int cluster_count = particle_weights[0].size();
        for (unsigned int i=0; i < jets.size(); i++) {
            double holdpT = tool->MLpT(particles, particle_weights, i,
                                       cluster_count, 0, true);
            pTs.push_back(SortablePt(i, holdpT));
        }

        sort(pTs.begin(), pTs.end());
        // reverse the list of indices and push it back onto the output
        for (unsigned int iter = jets.size(); iter-- > 0;) {
            indices.push_back((signed int) pTs.at(iter).index);
            //std::cout << pTs.at(iter).index << " : " << pTs.at(iter).pT << std::endl;
        }
        if (jets.size()) {
            pT = pTs.at(pTs.size() - 1).pT;
        }
    }

    double estimateRho(vecPseudoJet arg_particles) {
        fastjet::JetDefinition jet_def_estimator =
            fastjet::JetDefinition(fastjet::kt_algorithm, 0.4);
        fastjet::AreaDefinition area_def_estimator = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts);
        fastjet::Selector rap_sel = fastjet::SelectorAbsRapMax(1.5);
        fastjet::JetMedianBackgroundEstimator bge(rap_sel, jet_def_estimator,
                                                  area_def_estimator);
        bge.set_particles(arg_particles);
        double rho = bge.rho();
        return rho;
    }

    vecPseudoJet groomedInitialLocations(vecPseudoJet const& particles,
                                         float akt_jet_size, float filter_size,
                                         float filter_pT_frac, float pT_min,
                                         double rho, float seed_location_noise) {
        // build a list of input jet locations on the basis of groomed anti-kt
        // jet locations
        fastjet::JetDefinition m_jet_def(fastjet::antikt_algorithm, akt_jet_size);
        fastjet::AreaDefinition a_def = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts);
        fastjet::ClusterSequenceArea cs_seed_akt_ungroomed(particles, m_jet_def, a_def);
        fastjet::Filter seed_filter(filter_size,
                                    fastjet::SelectorPtFractionMin(filter_pT_frac));

        vecPseudoJet original_jets = cs_seed_akt_ungroomed.inclusive_jets(0);
        vecPseudoJet final_locations;
        final_locations.clear();

        for (unsigned int j_i = 0; j_i < original_jets.size(); j_i++) {
            fastjet::PseudoJet trimmed = seed_filter(original_jets.at(j_i));
            if (trimmed.pt() > pT_min + rho * trimmed.area()) {
                fastjet::PseudoJet accepted_location;

                double noise_amp = seed_location_noise * sqrt(randf(0, 1));
                if (seed_location_noise <= 0) noise_amp = 0;
                double noise_angle = randf(0, 2 * M_PI);
                double x = trimmed.rapidity() + noise_amp * cos(noise_angle);
                double y = trimmed.phi() + noise_amp * sin(noise_angle);
                while (y < 0) y += 2 * M_PI;
                while (y > 2 * M_PI) y -= 2 * M_PI;
                accepted_location.
                    reset_PtYPhiM(1.0, x, y, 0);
                final_locations.push_back(accepted_location);
            }
        }

        return final_locations;
    }

    void DoMUMMJetFinding(vecPseudoJet& particles,
                          vecPseudoJet& seeds,
                          bool learn_weights,
                          bool do_recombination,
                          double size,
                          vector<int>& lead_indices,
                          double& lead_pT,
                          FuzzyTools *tool,
                          vecPseudoJet& jets,
                          vector<vector<double> >& particle_weights,
                          vector<double>& jet_weights,
                          unsigned int &iters,
                          FuzzyTools::EventJet event_jet_type,
                          FuzzyTools::PostProcess post_processing_method,
                          double event_jet_weight) {
        vecPseudoJet parts = particles;
        tool->SetEventJetType(event_jet_type);
        tool->SetPostProcess(post_processing_method);
        if (event_jet_type == FuzzyTools::FLAT) {
            tool->SetEventJetWeight(event_jet_weight);
        }
        tool->SetKernelType(FuzzyTools::UNIFORM);
        tool->SetLearnWeights(learn_weights);
        if (do_recombination) {
            tool->SetClusteringMode(FuzzyTools::RECOMBINATION);
            tool->SetSeeds(parts);
        } else {
            tool->SetClusteringMode(FuzzyTools::FIXED);
            tool->SetSeeds(seeds);
        }
        tool->SetR(size);
        jets = tool->ClusterFuzzyUniform(particles,
                                         &particle_weights,
                                         &jet_weights,
                                         iters);
        FindLeadingJet(particles, jets, particle_weights, tool, lead_indices, lead_pT);
    }

    void DoMTGMMJetFinding(vecPseudoJet& particles,
                           vecPseudoJet& seeds,
                           bool learn_weights,
                           bool learn_shape,
                           bool do_recombination,
                           double size,
                           vector<int>& lead_indices,
                           double& lead_pT,
                           FuzzyTools *tool,
                           vecPseudoJet& jets,
                           vector<vector<double> >& particle_weights,
                           vector<MatTwo>& parameters,
                           vector<double>& jet_weights,
                           unsigned int &iters,
                           FuzzyTools::EventJet event_jet_type,
                           FuzzyTools::PostProcess post_processing_method,
                           double event_jet_weight) {
        tool->SetKernelType(FuzzyTools::TRUNCGAUSSIAN);
        tool->SetLearnWeights(learn_weights);
        tool->SetEventJetType(event_jet_type);
        tool->SetPostProcess(post_processing_method);
        if (event_jet_type == FuzzyTools::FLAT) {
            tool->SetEventJetWeight(event_jet_weight);
        }
        if(learn_shape) {
            tool->SetClusteringMode(FuzzyTools::FIXED);
            tool->SetSeeds(seeds);
        } else {
            if (do_recombination) {
                tool->SetClusteringMode(FuzzyTools::RECOMBINATION);
                tool->SetSeeds(particles);
            } else {
                tool->SetClusteringMode(FuzzyTools::FIXED);
                tool->SetSeeds(seeds);
            }
        }
        tool->SetLearnShape(learn_shape);
        tool->SetR(size);
        jets = tool->ClusterFuzzyTruncGaus(particles,
                                           &particle_weights,
                                           &parameters,
                                           &jet_weights,
                                           iters);
        FindLeadingJet(particles, jets, particle_weights, tool, lead_indices, lead_pT);
    }

    void DoMGMMCJetFinding(vecPseudoJet& particles,
                           vecPseudoJet& seeds,
                           bool learn_weights,
                           bool do_recombination,
                           vector<int>& lead_indices,
                           double& lead_pT,
                           FuzzyTools *tool,
                           vecPseudoJet& jets,
                           vector<vector<double> >& particle_weights,
                           vector<MatTwo>& parameters,
                           vector<double>& jet_weights,
                           unsigned int &iters,
                           FuzzyTools::EventJet event_jet_type,
                           FuzzyTools::PostProcess post_processing_method,
                           double event_jet_weight) {
        tool->SetKernelType(FuzzyTools::GAUSSIAN);
        tool->SetLearnWeights(learn_weights);
        tool->SetEventJetType(event_jet_type);
        tool->SetPostProcess(post_processing_method);
        if (event_jet_type == FuzzyTools::FLAT) {
            tool->SetEventJetWeight(event_jet_weight);
        }
        if (do_recombination) {
            tool->SetClusteringMode(FuzzyTools::RECOMBINATION);
            tool->SetSeeds(particles);
        } else {
            tool->SetClusteringMode(FuzzyTools::FIXED);
            tool->SetSeeds(seeds);
        }
        jets = tool->ClusterFuzzyGaussianC(particles,
                                           &particle_weights,
                                           &parameters,
                                           &jet_weights,
                                           iters);
        FindLeadingJet(particles, jets, particle_weights, tool, lead_indices, lead_pT);
    }


    void DoMGMMJetFinding(vecPseudoJet& particles,
                          vecPseudoJet& seeds,
                          bool learn_weights,
                          bool learn_shape,
                          bool do_recombination,
                          vector<int>& lead_indices,
                          double& lead_pT,
                          FuzzyTools *tool,
                          vecPseudoJet& jets,
                          vector<vector<double> >& particle_weights,
                          vector<MatTwo>& parameters,
                          vector<double>& jet_weights,
                          unsigned int &iters,
                          FuzzyTools::EventJet event_jet_type,
                          FuzzyTools::PostProcess post_processing_method,
                          double event_jet_weight) {
        tool->SetKernelType(FuzzyTools::GAUSSIAN);
        tool->SetLearnWeights(learn_weights);
        tool->SetEventJetType(event_jet_type);
        tool->SetPostProcess(post_processing_method);
        if (event_jet_type == FuzzyTools::FLAT) {
            tool->SetEventJetWeight(event_jet_weight);
        }
        if(learn_shape) {
            tool->SetClusteringMode(FuzzyTools::FIXED);
            tool->SetSeeds(seeds);
        } else {
            if (do_recombination) {
                tool->SetClusteringMode(FuzzyTools::RECOMBINATION);
                tool->SetSeeds(particles);
            } else {
                tool->SetClusteringMode(FuzzyTools::FIXED);
                tool->SetSeeds(seeds);
            }
        }
        tool->SetLearnShape(learn_shape);
        jets = tool->ClusterFuzzyGaussian(particles,
                                          &particle_weights,
                                          &parameters,
                                          &jet_weights,
                                          iters);
        FindLeadingJet(particles, jets, particle_weights, tool, lead_indices, lead_pT);
    }
}

// Constructor
FuzzyAnalysis::FuzzyAnalysis(){
    f_debug = false;
    if(f_debug) cout << "FuzzyAnalysis::FuzzyAnalysis Start " << endl;
    f_test = 0;
    f_out_name = "test.root";
    directory_prefix = "results/";

    stop_before_clustering = false;
    _non_standard_study = false;
    batched = false;
    should_print = false;

    pT_min = 5;

    tool = new FuzzyTools();

    tool->SetClusteringMode(FuzzyTools::RECOMBINATION);


    // jet def
    m_jet_def                 = new fastjet::JetDefinition(fastjet::antikt_algorithm, 0.4);
    m_jet_def_large_r_antikt  = new fastjet::JetDefinition(fastjet::antikt_algorithm, 1.0);
    m_jet_def_large_r_ca      = new fastjet::JetDefinition(fastjet::cambridge_algorithm, 1.0);
    m_jet_def_vlarge_r_antikt  = new fastjet::JetDefinition(fastjet::antikt_algorithm, 1.5);
    m_jet_def_vlarge_r_ca      = new fastjet::JetDefinition(fastjet::cambridge_algorithm, 1.5);

    if(f_debug) cout << "FuzzyAnalysis::FuzzyAnalysis End " << endl;
}

// Destructor
FuzzyAnalysis::~FuzzyAnalysis(){
    delete tool;
    delete m_jet_def;
    delete m_jet_def_large_r_antikt;
    delete m_jet_def_large_r_ca;
    delete m_jet_def_vlarge_r_antikt;
    delete m_jet_def_vlarge_r_ca;
}

// Begin method
void FuzzyAnalysis::Begin(){
    // Declare TTree
    string fullName = directory_prefix + f_out_name;
    t_f = new TFile(fullName.c_str(), "RECREATE");
    t_t = new TTree("EventTree", "Event Tree for Fuzzy");
    if (!_non_standard_study) {
#ifdef WITHROOT
        if(!batched) {
            SetupHistosMap();
        }
#endif
        DeclareBranches();
        ResetBranches();
    }

    tool->SetPrefix(directory_prefix);

    return;
}

vecPseudoJet FuzzyAnalysis::DiscretizeParticles(vecPseudoJet const& arg_particles,
                                                int n_phi, int n_eta,
                                                float eta_cutoff, bool keep_zero_cells) {
    TH2F *grid_E = new TH2F("grid_E", "", n_eta, -eta_cutoff, eta_cutoff, n_phi,
                            0, 2*M_PI);
    TH2F *grid_charge = new TH2F("grid_charge", "", n_eta, -eta_cutoff, eta_cutoff, n_phi,
                                 0, 2*M_PI);
    for (unsigned int p_i = 0; p_i < arg_particles.size(); p_i++) {
        fastjet::PseudoJet p = arg_particles[p_i];
        double eta = p.eta();
        double phi = p.phi();
        if (phi < 0) phi += (2*M_PI);
        double charge = p.user_info<MyUserInfo>().charge();
        double E = p.E();
        grid_E->Fill(eta, phi, E);
        grid_charge->Fill(eta, phi, charge);
    }

    // build the discretized particles back into jets
    vecPseudoJet out;
    for (int x_i = 1; x_i < grid_E->GetNbinsX() + 1; x_i++) {
        for (int y_i = 1; y_i < grid_E->GetNbinsY() + 1; y_i++) {
            float eta = grid_E->GetXaxis()->GetBinCenter(x_i);
            float phi = grid_E->GetXaxis()->GetBinCenter(y_i);
            float E = grid_E->GetBinContent(x_i, y_i);
            float charge = grid_charge->GetBinContent(x_i, y_i);
            fastjet::PseudoJet p;
            // eta = phi because m = 0
            // don't have a particle ID anymore
            if (E) {
                p.reset_PtYPhiM(E/cosh(eta), eta, phi, 0); // equivalent to energy with m=0
                p.set_user_info(new MyUserInfo(0, 0, charge, false, false));
                out.push_back(p);
            } else if (keep_zero_cells) {
                // "ghost" the cell to prevent fastjet screwing with its four vector
                p.reset_PtYPhiM(FLT_EPSILON, eta, phi, 0);
                p.set_user_info(new MyUserInfo(0, 0, charge, false, false));
                out.push_back(p);
            }
        }
    }
    delete grid_E;
    delete grid_charge;

    return out;
}

vecPseudoJet FuzzyAnalysis::ChargedParticles(vecPseudoJet const& arg_particles,
                                             double eta_max) {
    vecPseudoJet charged_particles;
    for (unsigned int particle_iter = 0; particle_iter < arg_particles.size(); particle_iter++) {
        double charge = arg_particles[particle_iter].user_info<MyUserInfo>().charge();
        double abs_eta = fabs(arg_particles[particle_iter].eta());
        if ((fabs(charge) > 1000 * DBL_EPSILON) && abs_eta < eta_max) {
            // is charged
            charged_particles.push_back(arg_particles[particle_iter]);
        }
    }
    return charged_particles;
}

vecPseudoJet FuzzyAnalysis::UnchargedParticles(vecPseudoJet const& arg_particles,
                                               double eta_max) {
    vecPseudoJet uncharged_particles;
    for (unsigned int particle_iter = 0; particle_iter < arg_particles.size(); particle_iter++) {
        double charge = arg_particles[particle_iter].user_info<MyUserInfo>().charge();
        double abs_eta = fabs(arg_particles[particle_iter].eta());
        if ((fabs(charge) <= 1000 * DBL_EPSILON) || abs_eta >= eta_max) {
            // is uncharged
            uncharged_particles.push_back(arg_particles[particle_iter]);
        }
    }
    return uncharged_particles;
}

void FuzzyAnalysis::SigmaMassStudy(Pythia8::Pythia *pythia_8, Pythia8::Pythia *pythia_MB) {
    batched = true; // don't do EDs
    gStyle->SetOptStat(0);
    SetEventJetType(FuzzyTools::NONE);

    SetToolAlpha(1);
    SetDoTowerSubtraction(false);

    TCanvas canv("massinsigma", "", 1400, 800);
    TH1F * mass_hist = new TH1F("H1", "", 20, 0, 250);
    TH1F * sigma_mass_hist = new TH1F("H2", "", 20, 0, 250);
    TH1F * antikt_mass_hist = new TH1F("H3", "", 20, 0, 250);

    mass_hist->SetLineColor(kRed);
    sigma_mass_hist->SetLineColor(kBlue);
    antikt_mass_hist->SetLineColor(kBlack);

    mass_hist->SetFillStyle(0);
    sigma_mass_hist->SetFillStyle(0);
    antikt_mass_hist->SetFillStyle(0);

    mass_hist->SetLineStyle(1);
    sigma_mass_hist->SetLineStyle(1);
    antikt_mass_hist->SetLineStyle(1);

    mass_hist->SetLineWidth(2);
    sigma_mass_hist->SetLineWidth(2);
    antikt_mass_hist->SetLineWidth(2);

    THStack stack("tstack", "");

    for (int event_iter = 0; event_iter < 1000; event_iter++) {
        if (event_iter % 50 == 0)
            std::cout << event_iter << std::endl;
        AnalyzeEvent(event_iter, pythia_8, pythia_MB, 0);
        mass_hist->Fill(fTmGMMc_m);
        antikt_mass_hist->Fill(fTantikt_m);
        sigma_mass_hist->Fill(fTmGMMc_m_sigma);
    }

    double n = mass_hist->Integral(-1, mass_hist->GetNbinsX() + 1);
    mass_hist->Scale(1./n);

    n = antikt_mass_hist->Integral(-1, antikt_mass_hist->GetNbinsX() + 1);
    antikt_mass_hist->Scale(1./n);

    n = sigma_mass_hist->Integral(-1, sigma_mass_hist->GetNbinsX() + 1);
    sigma_mass_hist->Scale(1./n);

    stack.Add(mass_hist);
    stack.Add(sigma_mass_hist);
    stack.Add(antikt_mass_hist);

    stack.Draw("nostack");

    stack.GetHistogram()->GetXaxis()->SetTitle("Mass [GeV]");
    stack.GetHistogram()->GetYaxis()->SetTitle("Arb. Units");
    stack.GetHistogram()->GetXaxis()->SetNdivisions(505);

    canv.Update();

    TLegend legend(0.65, 0.75, 0.9, 0.85);
    legend.AddEntry(mass_hist, "Fuzzy mass", "l");
    legend.AddEntry(antikt_mass_hist, "Anti-k_{t} mass", "l");
    legend.AddEntry(sigma_mass_hist, "Fuzzy mass in 1-#sigma", "l");

    legend.SetTextFont(42);
    legend.SetFillStyle(0);
    legend.SetFillColor(0);
    legend.SetBorderSize(0);

    legend.Draw();

    canv.Print("test2.pdf");
}

void FuzzyAnalysis::SigmaDependenceOnMuCompleteStudy(Pythia8::Pythia *pythia_8,
                                                     Pythia8::Pythia *pythia_MB) {
    // basic initialization
    gStyle->SetOptStat(0);
    bool CORRECTED = true;
    SetToolAlpha(1);
    if (CORRECTED) {
        SetEventJetType(FuzzyTools::FLAT);
        SetEventJetStrength(0.01);
        SetEventJetRhoOffset(0);
        SetDoTowerSubtraction(true);
        pT_min = 5;
    } else {
        SetEventJetType(FuzzyTools::NONE);
        SetEventJetStrength(0);
        SetDoTowerSubtraction(false);
        pT_min = 30;
    }

    // don't use HS (pythia_8)!!!
    // this Pythia instantiation will intentionally cause a crash
    // if using it is attempted because init has not been called
    Pythia8::Pythia *pythia_fail = new Pythia8::Pythia();


    // for lead radius and subleading radius
    // average deviation from \sigma_NPV=0 as a function of NPV
    // average variance of \sigma_NPV=0 as a function of NPV
    // histogram of deviation from \sigma_NPV=0
    // histogram of variance as a function of NPV

    // setup containers for data
    size_t n_events = 200;//200
    static const size_t NPVs[10] = {0, 2, 5, 10, 15, 20, 25, 30, 35, 40};

    size_t npv_iter_max = sizeof(NPVs) / sizeof(NPVs[0]);

    std::vector<float> sigma1;
    std::vector<float> sigma2;
    std::vector<float> eta1;
    std::vector<float> phi1;
    std::vector<float> mass1;
    std::vector<float> mass2;
    std::vector<float> mass_trimmed1;
    std::vector<float> mass_trimmed2;

    for (size_t r_i = 0; r_i < npv_iter_max; r_i++) {
        sigma1.push_back(0);
        sigma2.push_back(0);
        mass1.push_back(0);
        mass2.push_back(0);
        eta1.push_back(0);
        phi1.push_back(0);
        mass_trimmed1.push_back(0);
        mass_trimmed2.push_back(0);
    }

    for (size_t npv_iter = 0; npv_iter < npv_iter_max; npv_iter++) {
        size_t npv = NPVs[npv_iter];
        std::stringstream ss;

        ss.str(std::string());
        ss << "sigma1_" << npv;
        t_t->Branch(ss.str().c_str(), &sigma1.at(npv_iter), (ss.str() + "/F").c_str());

        ss.str(std::string());
        ss << "sigma2_" << npv;
        t_t->Branch(ss.str().c_str(), &sigma2.at(npv_iter), (ss.str() + "/F").c_str());

        ss.str(std::string());
        ss << "mass1_" << npv;
        t_t->Branch(ss.str().c_str(), &mass1.at(npv_iter), (ss.str() + "/F").c_str());

        ss.str(std::string());
        ss << "mass2_" << npv;
        t_t->Branch(ss.str().c_str(), &mass2.at(npv_iter), (ss.str() + "/F").c_str());

        ss.str(std::string());
        ss << "eta1_" << npv;
        t_t->Branch(ss.str().c_str(), &eta1.at(npv_iter), (ss.str() + "/F").c_str());

        ss.str(std::string());
        ss << "phi1_" << npv;
        t_t->Branch(ss.str().c_str(), &phi1.at(npv_iter), (ss.str() + "/F").c_str());

        ss.str(std::string());
        ss << "mass_trimmed1_" << npv;
        t_t->Branch(ss.str().c_str(), &mass_trimmed1.at(npv_iter), (ss.str() + "/F").c_str());

        ss.str(std::string());
        ss << "mass_trimmed2_" << npv;
        t_t->Branch(ss.str().c_str(), &mass_trimmed2.at(npv_iter), (ss.str() + "/F").c_str());
    }

    for (size_t event_idx = 0; event_idx < n_events; event_idx++) {
        vecPseudoJet hs_particles;
        vecPseudoJet tops;
        float max_rap = -1;

        generate_particles(hs_particles, tops, max_rap, false, pythia_8);

        batched = true;
        stop_before_clustering = false;

        for (size_t npv_iter = 0; npv_iter < npv_iter_max; npv_iter++) {
            unsigned int NPV = NPVs[npv_iter];
            AnalyzeEvent(event_idx, pythia_fail, pythia_MB, static_cast<int>(NPV), hs_particles, tops);
            // sigma1
            sigma1.at(npv_iter) = fTmGMMc_r;

            // sigma2
            sigma2.at(npv_iter) = fTmGMMc_r_second;

            // m1
            mass1.at(npv_iter) = fTantikt_m;

            // m2
            mass2.at(npv_iter) = fTantikt_m_second;

            // eta1
            eta1.at(npv_iter) = fTmGMMc_eta;

            // phi1
            phi1.at(npv_iter) = fTmGMMc_phi;

            // m_tr1
            mass_trimmed1.at(npv_iter) = fTantikt_m_trimmed_three;

            // m_tr2
            mass_trimmed2.at(npv_iter) = fTantikt_m_second_trimmed_three;
        }
        t_t->Fill();
    }

    // save and quit
    delete pythia_fail;
}

void FuzzyAnalysis::NoiseStudy(Pythia8::Pythia *pythia_8,
                               Pythia8::Pythia *pythia_MB, int NPV, size_t n_events) {
    // study to look at what happens when we add noise to the seeds

    // basic initialization
    gStyle->SetOptStat(0);
    SetToolAlpha(1);

    // don't use HS (pythia_8)!!!
    // this Pythia instantiation will intentionally cause a crash
    // if using it is attempted because init has not been called
    Pythia8::Pythia *pythia_fail = new Pythia8::Pythia();


    // for lead radius and subleading radius
    // average deviation from \sigma_NPV=0 as a function of NPV
    // average variance of \sigma_NPV=0 as a function of NPV
    // histogram of deviation from \sigma_NPV=0
    // histogram of variance as a function of NPV

    // setup containers for data
    static const float noise_amts[5] = {0, 0.5, 1.0, 1.5, 2.0};
    size_t noise_iter_max = sizeof(noise_amts)/sizeof(noise_amts[0]);

    std::vector<float> sigma1;
    std::vector<float> eta1;
    std::vector<float> phi1;

    for (size_t noise_iter = 0; noise_iter < noise_iter_max; noise_iter++) {
        sigma1.push_back(0);
        eta1.push_back(0);
        phi1.push_back(0);
    }
    for (size_t noise_iter = 0; noise_iter < noise_iter_max; noise_iter++) {
        int noise_print = static_cast<int>(noise_amts[noise_iter] * 100);

        std::stringstream ss;

        ss.str(std::string());
        ss << "sigma1_" << noise_print;
        t_t->Branch(ss.str().c_str(), &sigma1.at(noise_iter), (ss.str() + "/F").c_str());

        ss.str(std::string());
        ss << "eta1_" << noise_print;
        t_t->Branch(ss.str().c_str(), &eta1.at(noise_iter), (ss.str() + "/F").c_str());

        ss.str(std::string());
        ss << "phi1_" << noise_print;
        t_t->Branch(ss.str().c_str(), &phi1.at(noise_iter), (ss.str() + "/F").c_str());
    }

    for (size_t event_idx = 0; event_idx < n_events; event_idx++) {
        vecPseudoJet hs_particles;
        vecPseudoJet tops;
        float max_rap = -1;

        generate_particles(hs_particles, tops, max_rap, false, pythia_8);

        batched = true;
        stop_before_clustering = false;
        for (size_t noise_iter = 0; noise_iter < noise_iter_max; noise_iter++) {
            SetSeedLocationNoise(noise_amts[noise_iter]);
            AnalyzeEvent(event_idx, pythia_fail, pythia_MB, NPV, hs_particles, tops);

            sigma1.at(noise_iter) = fTmGMMc_r;
            phi1.at(noise_iter) = fTmGMMc_phi;
            eta1.at(noise_iter) = fTmGMMc_eta;
        }
        t_t->Fill();
    }

    // save and quit
    delete pythia_fail;
}

void FuzzyAnalysis::PtCutWanderingStudy(__attribute__((unused)) Pythia8::Pythia *pythia_8,
                                        __attribute__((unused)) Pythia8::Pythia *pythia_MB,
                                        int NPV, size_t n_events) {
    SetToolAlpha(1);

    float pT_baseline = 5;
    pT_min = pT_baseline;

    Pythia8::Pythia *pythia_fail = new Pythia8::Pythia();

    size_t cut_iter_max = 6; // 2 <= cut_iter <= cut_iter_max

    std::vector<float> sigma_deviation;
    std::vector<float> phi_deviation;
    std::vector<float> eta_deviation;
    for (size_t r_i = 0; r_i < cut_iter_max - 1; r_i++) {
        sigma_deviation.push_back(0);
        phi_deviation.push_back(0);
        eta_deviation.push_back(0);
    }

    for (size_t cut_iter = 0; cut_iter < cut_iter_max - 1; cut_iter++) {
        pT_min = (cut_iter + 2) * pT_baseline;
        std::stringstream ss;

        ss.str(std::string());
        ss << "sigma_norm_deviation_" << pT_min;
        t_t->Branch(ss.str().c_str(), &sigma_deviation.at(cut_iter), (ss.str()+"/F").c_str());

        ss.str(std::string());
        ss << "phi_deviation_" << pT_min;
        t_t->Branch(ss.str().c_str(), &phi_deviation.at(cut_iter), (ss.str()+"/F").c_str());

        ss.str(std::string());
        ss << "eta_deviation_" << pT_min;
        t_t->Branch(ss.str().c_str(), &eta_deviation.at(cut_iter), (ss.str()+"/F").c_str());
    }

    for (size_t event_idx = 0; event_idx < n_events; event_idx++) {
        vecPseudoJet hs_particles;
        vecPseudoJet tops;
        float max_rap = -1;

        generate_particles(hs_particles, tops, max_rap, false, pythia_8);

        batched = true;
        stop_before_clustering = false;

        float pT_baseline_sigma = -999999;
        float pT_baseline_eta = -999999;
        float pT_baseline_phi = -999999;
        for (size_t cut_iter = 0; cut_iter < cut_iter_max; cut_iter++) {
            pT_min = (cut_iter + 1) * pT_baseline;
            AnalyzeEvent(0, pythia_fail, pythia_MB, NPV, hs_particles, tops);
            if (cut_iter == 0) {
                pT_baseline_sigma = fTmGMMc_r;
                pT_baseline_phi = fTmGMMc_phi;
                pT_baseline_eta = fTmGMMc_eta;
                continue;
            }
            sigma_deviation.at(cut_iter - 1) = (fTmGMMc_r - pT_baseline_sigma)/pT_baseline_sigma;
            phi_deviation.at(cut_iter - 1) = fTmGMMc_phi - pT_baseline_phi;
            eta_deviation.at(cut_iter - 1) = fTmGMMc_eta - pT_baseline_eta;
        }
        t_t->Fill();
    }

    delete pythia_fail;
}

void FuzzyAnalysis::SigmaDependenceOnMuStudy(Pythia8::Pythia *pythia_8, Pythia8::Pythia *pythia_MB) {
    // basic initialization
    gStyle->SetOptStat(0);
    SetEventJetType(FuzzyTools::NONE);
    //SetEventJetStrength(0.01);
    //SetEventJetRhoOffset(0);
    SetToolAlpha(1);
    SetDoTowerSubtraction(false);
    pT_min = 30;

    // don't use HS (pythia_8)!!!
    // this Pythia instantiation will intentionally cause a crash
    // if using it is attempted because init has not been called
    Pythia8::Pythia *pythia_fail = new Pythia8::Pythia();
    vecPseudoJet hs_particles;
    vecPseudoJet tops;
    float max_rap = -1;

    generate_particles(hs_particles, tops, max_rap, false, pythia_8);

    batched = true;
    stop_before_clustering = false;

    unsigned int n_samples = 150;
    // sigma
    std::vector<float> avg_sigma1;
    std::vector<float> avg_sigma2;
    std::vector<float> avg_sigma_stdevs1;
    std::vector<float> avg_sigma_stdevs2;

    // antikt_m
    std::vector<float> avg_m1;
    std::vector<float> avg_m_stdevs1;
    std::vector<float> avg_m2;
    std::vector<float> avg_m_stdevs2;

    // antikt_m_trimmed
    std::vector<float> avg_m_tr1;
    std::vector<float> avg_m_tr_stdevs1;
    std::vector<float> avg_m_tr2;
    std::vector<float> avg_m_tr_stdevs2;

    std::vector<float> NPVs_float;
    for (unsigned int NPV = 0; NPV <= 80; NPV += 5) {
        NPVs_float.push_back(static_cast<float>(NPV));
        std::cout << NPV << std::flush;
        std::vector<float> samples_sigma1;
        std::vector<float> samples_sigma2;
        std::vector<float> samples_m1;
        std::vector<float> samples_m2;
        std::vector<float> samples_m_tr1;
        std::vector<float> samples_m_tr2;
        for (unsigned int s_i = 0; s_i < n_samples; s_i++) {
            std::cout << "." << std::flush;
            AnalyzeEvent(s_i, pythia_fail, pythia_MB, static_cast<int>(NPV), hs_particles, tops);
            samples_sigma1.push_back(fTmGMMc_r);
            samples_sigma2.push_back(fTmGMMc_r_second);

            samples_m1.push_back(fTantikt_m);
            samples_m2.push_back(fTantikt_m_second);

            samples_m_tr1.push_back(fTantikt_m_trimmed_three);
            samples_m_tr2.push_back(fTantikt_m_second_trimmed_three);
        }
        std::cout << std::endl;
        float mean, stdev;
        // sigma1
        calc_stats(samples_sigma1, mean, stdev);
        avg_sigma1.push_back(mean);
        avg_sigma_stdevs1.push_back(stdev);

        // sigma2
        calc_stats(samples_sigma2, mean, stdev);
        avg_sigma2.push_back(mean);
        avg_sigma_stdevs2.push_back(stdev);

        // m1
        calc_stats(samples_m1, mean, stdev);
        avg_m1.push_back(mean);
        avg_m_stdevs1.push_back(stdev);

        // m2
        calc_stats(samples_m2, mean, stdev);
        avg_m2.push_back(mean);
        avg_m_stdevs2.push_back(stdev);

        // m_tr1
        calc_stats(samples_m_tr1, mean, stdev);
        avg_m_tr1.push_back(mean);
        avg_m_tr_stdevs1.push_back(stdev);

        // m_tr2
        calc_stats(samples_m_tr2, mean, stdev);
        avg_m_tr2.push_back(mean);
        avg_m_tr_stdevs2.push_back(stdev);
    }

    // draw sigma canvas
    multi_with_errors(NPVs_float, avg_sigma1, avg_sigma_stdevs1,
                      avg_sigma2, avg_sigma_stdevs2, "Sigma", "mGMMc_r_",
                      directory_prefix);

    // draw m canvas
    multi_with_errors(NPVs_float, avg_m1, avg_m_stdevs1,
                      avg_m2, avg_m_stdevs2, "Mass", "antikt_m_",
                      directory_prefix);

    // draw m_trimmed canvas
    multi_with_errors(NPVs_float, avg_m_tr1, avg_m_tr_stdevs1,
                      avg_m_tr2, avg_m_tr_stdevs2,
                      "Trimmed Mass", "antikt_m_trimmed_three_",
                      directory_prefix);

    //TLegend *legend = new TLegend(0.2, 0.25, 0.45, 0.35);
    //legend->SetTextFont(42);
    //legend->SetFillStyle(0);
    //legend->SetFillColor(0);
    //legend->SetBorderSize(0);

    // save and quit
    delete pythia_fail;
}

void FuzzyAnalysis::GroomingStudy(vecPseudoJet leading_particles,
                                  int event_iter) {
    // recluster the particles belonging to the leading jet
    // and investigate a number of techniques for grooming

    // 1) Recluster with 0.2 radius mGMMc jets with seeds from the location
    //    of kt jets with similar size. Remove jets which are in the bottom Xth
    //    percentile by way of pT.
    gStyle->SetOptStat(0);
    fastjet::JetDefinition *m_jet_def_seed = new fastjet::JetDefinition(fastjet::antikt_algorithm, 0.2);
    fastjet::ClusterSequence cs_seed_antikt(leading_particles, *m_jet_def_seed);
    double minimum_seed_pT = 1.0;
    vecPseudoJet seeds = fastjet::sorted_by_pt(cs_seed_antikt.inclusive_jets(minimum_seed_pT));
    if (seeds.size() == 0)
        return;

    bool mUMM_on = true;
    bool mGMM_on = true;
    bool mGMMc_on = true;
    bool mGMMs_on = true;
    bool mTGMM_on = true;
    bool mTGMMs_on = true;

    double subs_sigma = 0.15;
    tool->SetDefaultSigma(MatTwo(subs_sigma*subs_sigma, 0, 0, subs_sigma*subs_sigma));
    tool->SetR(0.2);
    tool->SetMergeDistance(0.05);

    unsigned int dummy;
    FuzzyTools::EventJet event_jet_type = FuzzyTools::NONE;
    FuzzyTools::PostProcess post_processing_method = FuzzyTools::NO_POST_PROCESS;
    double event_jet_weight = 0;

    // Fuzzy Jets: mGMMs --------------------
    vector<vector<double> > mGMMs_particle_weights;
    vector<MatTwo> mGMMs_jets_params;
    vector<double> mGMMs_weights;
    vecPseudoJet mGMMs_jets;
    vector<int> mGMMs_indices;
    double max_pT_mGMMs;
    if(mGMMs_on) {
        DoMGMMJetFinding(leading_particles, seeds,
                         f_learn_weights, true, false,
                         mGMMs_indices, max_pT_mGMMs,
                         tool, mGMMs_jets, mGMMs_particle_weights,
                         mGMMs_jets_params, mGMMs_weights, dummy,
                         event_jet_type, post_processing_method,
                         event_jet_weight);
    }
    __attribute__((unused)) int lead_mGMMs_index = mGMMs_indices.size() ? mGMMs_indices.at(0) : -1;


    // Fuzzy Jets: mGMMc --------------------
    vector<vector<double> > mGMMc_particle_weights;
    vector<MatTwo> mGMMc_jets_params;
    vector<double> mGMMc_weights;
    vecPseudoJet mGMMc_jets;
    vector<int> mGMMc_indices;
    double max_pT_mGMMc;
    if(mGMMc_on) {
        DoMGMMCJetFinding(leading_particles, seeds,
                          f_learn_weights, false,
                          mGMMc_indices, max_pT_mGMMc,
                          tool, mGMMc_jets, mGMMc_particle_weights,
                          mGMMc_jets_params, mGMMc_weights, dummy,
                          event_jet_type, post_processing_method,
                          event_jet_weight);
    }
    __attribute__((unused)) int lead_mGMMc_index = mGMMc_indices.size() ? mGMMc_indices.at(0) : -1;

    // Fuzzy Jets: mTGMMs -------------------
    vector<vector<double > > mTGMMs_particle_weights;
    vector<MatTwo> mTGMMs_jets_params;
    vector<double> mTGMMs_weights;
    vecPseudoJet mTGMMs_jets;
    vector<int> mTGMMs_indices;
    double max_pT_mTGMMs;
    if(mTGMMs_on) {
        DoMTGMMJetFinding(leading_particles, seeds,
                          f_learn_weights, true, false,
                          f_size, mTGMMs_indices, max_pT_mTGMMs,
                          tool, mTGMMs_jets, mTGMMs_particle_weights,
                          mTGMMs_jets_params, mTGMMs_weights, dummy,
                          event_jet_type, post_processing_method,
                          event_jet_weight);
    }
    __attribute__((unused)) int lead_mTGMMs_index = mTGMMs_indices.size() ? mTGMMs_indices.at(0) : -1;

    // Fuzzy Jets: mGMM ---------------------
    vector<vector<double> >mGMM_particle_weights;
    vector<MatTwo> mGMM_jets_params;
    vector<double> mGMM_weights;
    vecPseudoJet mGMM_jets;
    vector<int> mGMM_indices;
    double max_pT_mGMM;
    if(mGMM_on) {
        DoMGMMJetFinding(leading_particles, seeds,
                         f_learn_weights, false, false,
                         mGMM_indices, max_pT_mGMM,
                         tool, mGMM_jets, mGMM_particle_weights,
                         mGMM_jets_params, mGMM_weights, dummy, event_jet_type,
                         post_processing_method, event_jet_weight);
    }
    __attribute__((unused)) int lead_mGMM_index = mGMM_indices.size() ? mGMM_indices.at(0) : -1;

    // Fuzzy Jets: mUMM ---------------------
    vector<vector<double> > mUMM_particle_weights;
    vector<double> mUMM_weights;
    vecPseudoJet mUMM_jets;
    vector<int> mUMM_indices;
    double max_pT_mUMM;
    if(mUMM_on) {
        DoMUMMJetFinding(leading_particles, seeds,
                         f_learn_weights, f_size, false,
                         mUMM_indices, max_pT_mUMM, tool, mUMM_jets,
                         mUMM_particle_weights, mUMM_weights, dummy, event_jet_type,
                         post_processing_method, event_jet_weight);
    }
    __attribute__((unused)) int lead_mUMM_index = mUMM_indices.size() ? mUMM_indices.at(0) : -1;

    // Fuzzy Jets: mTGMM --------------------
    vector<vector<double> > mTGMM_particle_weights;
    vector<double> mTGMM_weights;
    vecPseudoJet mTGMM_jets;
    vector<MatTwo> mTGMM_jets_params;
    vector<int> mTGMM_indices;
    double max_pT_mTGMM;
    if(mTGMM_on) {
        DoMTGMMJetFinding(leading_particles, seeds,
                          f_learn_weights, false, false, f_size,
                          mTGMM_indices, max_pT_mTGMM, tool, mTGMM_jets,
                          mTGMM_particle_weights, mTGMM_jets_params,
                          mTGMM_weights, dummy, event_jet_type,
                          post_processing_method, event_jet_weight);
    }
    __attribute__((unused)) int lead_mTGMM_index = mTGMM_indices.size() ? mTGMM_indices.at(0) : -1;

    if(!mGMM_weights.size()) {
        mGMM_on = false;
    }
    if(!mGMMs_weights.size()) {
        mGMMs_on = false;
    }
    if(!mGMMc_weights.size()) {
        mGMMc_on = false;
    }
    if(!mTGMM_weights.size()) {
        mTGMM_on = false;
    }
    if(!mTGMMs_weights.size()) {
        mTGMMs_on = false;
    }
    if(!mUMM_weights.size()) {
        mUMM_on = false;
    }

    double total_pT = 0;
    for (unsigned int p_i = 0; p_i < leading_particles.size(); p_i++) {
        total_pT += leading_particles.at(p_i).pt();
    }
    // remove jets up to Xth percentile in pT
    float X = 0.25;
    if (mGMM_on) {
        // find groomed components for mGMM subjets
        vector<int> indices;
        double pT;
        vecPseudoJet groomed_particles;
        FindLeadingJet(leading_particles,
                       mGMM_jets, mGMM_particle_weights,
                       tool, indices, pT);
        for (unsigned int p_i = 0; p_i < leading_particles.size(); p_i++) {
            double max_weight = -1;
            int belongs_idx = -1;
            for (unsigned int j_i = 0; j_i < mGMM_jets.size(); j_i++) {
                double current_weight = mGMM_particle_weights.at(p_i).at(j_i);
                if (current_weight > max_weight) {
                    max_weight = current_weight;
                    belongs_idx = j_i;
                }
            }
            if (max_weight >= 0 && ((1.0*indices.at(belongs_idx)/mGMM_jets.size()) < (1-X))) {
                groomed_particles.push_back(leading_particles.at(p_i));
            }
        }

        // UNFINISHED DISPLAY ALL JETS AND PARTICLES AND GROOMED ONES
        tool->NewEventDisplay(leading_particles, leading_particles,
                              leading_particles, mGMM_jets, mGMM_particle_weights,
                              0, mGMM_jets_params, mGMM_weights, "mGMM_ungroomed", event_iter);
        vecPseudoJet groomed_jets;
        //vector<vector<double>> groomed_weights; not necessary it turns out,
        // since these just get logged
        vector<MatTwo> groomed_jets_params;
        vector<double> groomed_jet_weights;
        for (unsigned int j_i = 0; j_i < mGMM_jets.size(); j_i++) {
            if ((1.0*indices.at(j_i)/mGMM_jets.size()) < (1-X)) {
                groomed_jets.push_back(mGMM_jets.at(j_i));
                groomed_jets_params.push_back(mGMM_jets_params.at(j_i));
                groomed_jet_weights.push_back(mGMM_weights.at(j_i));
            }
        }
        tool->NewEventDisplay(groomed_particles, groomed_particles,
                              groomed_particles, groomed_jets, mGMM_particle_weights,
                              0, groomed_jets_params, groomed_jet_weights, "mGMM_groomed", event_iter);
    }
    if (mGMMs_on) {
        // find groomed components for mGMMs subjets
        vector<int> indices;
        double pT;
        vecPseudoJet groomed_particles;
        FindLeadingJet(leading_particles,
                       mGMMs_jets, mGMMs_particle_weights,
                       tool, indices, pT);
        for (unsigned int p_i = 0; p_i < leading_particles.size(); p_i++) {
            double max_weight = -1;
            int belongs_idx = -1;
            for (unsigned int j_i = 0; j_i < mGMMs_jets.size(); j_i++) {
                double current_weight = mGMMs_particle_weights.at(p_i).at(j_i);
                if (current_weight > max_weight) {
                    max_weight = current_weight;
                    belongs_idx = j_i;
                }
            }
            if (max_weight >= 0 && ((1.0*indices.at(belongs_idx)/mGMMs_jets.size()) < (1-X))) {
                groomed_particles.push_back(leading_particles.at(p_i));
            }
        }

        // UNFINISHED DISPLAY ALL JETS AND PARTICLES AND GROOMED ONES
        tool->NewEventDisplay(leading_particles, leading_particles,
                              leading_particles, mGMMs_jets, mGMMs_particle_weights,
                              0, mGMMs_jets_params, mGMMs_weights, "mGMMs_ungroomed", event_iter);
        vecPseudoJet groomed_jets;
        //vector<vector<double>> groomed_weights; not necessary it turns out,
        // since these just get logged
        vector<MatTwo> groomed_jets_params;
        vector<double> groomed_jet_weights;
        for (unsigned int j_i = 0; j_i < mGMMs_jets.size(); j_i++) {
            if ((1.0*indices.at(j_i)/mGMMs_jets.size()) < (1-X)) {
                groomed_jets.push_back(mGMMs_jets.at(j_i));
                groomed_jets_params.push_back(mGMMs_jets_params.at(j_i));
                groomed_jet_weights.push_back(mGMMs_weights.at(j_i));
            }
        }
        tool->NewEventDisplay(groomed_particles, groomed_particles,
                              groomed_particles, groomed_jets, mGMMs_particle_weights,
                              0, groomed_jets_params, groomed_jet_weights, "mGMMs_groomed", event_iter);
    }
    if (mGMMc_on) {
        // find groomed components for mGMMc subjets
        vector<int> indices;
        double pT;
        vecPseudoJet groomed_particles;
        FindLeadingJet(leading_particles,
                       mGMMc_jets, mGMMc_particle_weights,
                       tool, indices, pT);
        for (unsigned int p_i = 0; p_i < leading_particles.size(); p_i++) {
            double max_weight = -1;
            int belongs_idx = -1;
            for (unsigned int j_i = 0; j_i < mGMMc_jets.size(); j_i++) {
                double current_weight = mGMMc_particle_weights.at(p_i).at(j_i);
                if (current_weight > max_weight) {
                    max_weight = current_weight;
                    belongs_idx = j_i;
                }
            }
            if (max_weight >= 0 && ((1.0*indices.at(belongs_idx)/mGMMc_jets.size()) < (1-X))) {
                groomed_particles.push_back(leading_particles.at(p_i));
            }
        }

        // UNFINISHED DISPLAY ALL JETS AND PARTICLES AND GROOMED ONES
        tool->NewEventDisplay(leading_particles, leading_particles,
                              leading_particles, mGMMc_jets, mGMMc_particle_weights,
                              0, mGMMc_jets_params, mGMMc_weights, "mGMMc_ungroomed", event_iter);
        vecPseudoJet groomed_jets;
        //vector<vector<double>> groomed_weights; not necessary it turns out,
        // since these just get logged
        vector<MatTwo> groomed_jets_params;
        vector<double> groomed_jet_weights;
        for (unsigned int j_i = 0; j_i < mGMMc_jets.size(); j_i++) {
            if ((1.0*indices.at(j_i)/mGMMc_jets.size()) < (1-X)) {
                groomed_jets.push_back(mGMMc_jets.at(j_i));
                groomed_jets_params.push_back(mGMMc_jets_params.at(j_i));
                groomed_jet_weights.push_back(mGMMc_weights.at(j_i));
            }
        }
        tool->NewEventDisplay(groomed_particles, groomed_particles,
                              groomed_particles, groomed_jets, mGMMc_particle_weights,
                              0, groomed_jets_params, groomed_jet_weights, "mGMMc_groomed", event_iter);
    }
    if (mTGMM_on) {
        // find groomed components for mTGMM subjets
        vector<int> indices;
        double pT;
        vecPseudoJet groomed_particles;
        FindLeadingJet(leading_particles,
                       mTGMM_jets, mTGMM_particle_weights,
                       tool, indices, pT);
        for (unsigned int p_i = 0; p_i < leading_particles.size(); p_i++) {
            double max_weight = -1;
            int belongs_idx = -1;
            for (unsigned int j_i = 0; j_i < mTGMM_jets.size(); j_i++) {
                double current_weight = mTGMM_particle_weights.at(p_i).at(j_i);
                if (current_weight > max_weight) {
                    max_weight = current_weight;
                    belongs_idx = j_i;
                }
            }
            if (max_weight >= 0 && ((1.0*indices.at(belongs_idx)/mTGMM_jets.size()) < (1-X))) {
                groomed_particles.push_back(leading_particles.at(p_i));
            }
        }

        // UNFINISHED DISPLAY ALL JETS AND PARTICLES AND GROOMED ONES
        tool->NewEventDisplay(leading_particles, leading_particles,
                              leading_particles, mTGMM_jets, mTGMM_particle_weights,
                              0, mTGMM_jets_params, mTGMM_weights, "mTGMM_ungroomed", event_iter);
        vecPseudoJet groomed_jets;
        //vector<vector<double>> groomed_weights; not necessary it turns out,
        // since these just get logged
        vector<MatTwo> groomed_jets_params;
        vector<double> groomed_jet_weights;
        for (unsigned int j_i = 0; j_i < mTGMM_jets.size(); j_i++) {
            if ((1.0*indices.at(j_i)/mTGMM_jets.size()) < (1-X)) {
                groomed_jets.push_back(mTGMM_jets.at(j_i));
                groomed_jets_params.push_back(mTGMM_jets_params.at(j_i));
                groomed_jet_weights.push_back(mTGMM_weights.at(j_i));
            }
        }
        tool->NewEventDisplay(groomed_particles, groomed_particles,
                              groomed_particles, groomed_jets, mTGMM_particle_weights,
                              0, groomed_jets_params, groomed_jet_weights, "mTGMM_groomed", event_iter);
    }
    if (mTGMMs_on) {
        // find groomed components for mTGMMs subjets
        vector<int> indices;
        double pT;
        vecPseudoJet groomed_particles;
        FindLeadingJet(leading_particles,
                       mTGMMs_jets, mTGMMs_particle_weights,
                       tool, indices, pT);
        for (unsigned int p_i = 0; p_i < leading_particles.size(); p_i++) {
            double max_weight = -1;
            int belongs_idx = -1;
            for (unsigned int j_i = 0; j_i < mTGMMs_jets.size(); j_i++) {
                double current_weight = mTGMMs_particle_weights.at(p_i).at(j_i);
                if (current_weight > max_weight) {
                    max_weight = current_weight;
                    belongs_idx = j_i;
                }
            }
            if (max_weight >= 0 && ((1.0*indices.at(belongs_idx)/mTGMMs_jets.size()) < (1-X))) {
                groomed_particles.push_back(leading_particles.at(p_i));
            }
        }

        // UNFINISHED DISPLAY ALL JETS AND PARTICLES AND GROOMED ONES
        tool->NewEventDisplay(leading_particles, leading_particles,
                              leading_particles, mTGMMs_jets, mTGMMs_particle_weights,
                              0, mTGMMs_jets_params, mTGMMs_weights, "mTGMMs_ungroomed", event_iter);
        vecPseudoJet groomed_jets;
        //vector<vector<double>> groomed_weights; not necessary it turns out,
        // since these just get logged
        vector<MatTwo> groomed_jets_params;
        vector<double> groomed_jet_weights;
        for (unsigned int j_i = 0; j_i < mTGMMs_jets.size(); j_i++) {
            if ((1.0*indices.at(j_i)/mTGMMs_jets.size()) < (1-X)) {
                groomed_jets.push_back(mTGMMs_jets.at(j_i));
                groomed_jets_params.push_back(mTGMMs_jets_params.at(j_i));
                groomed_jet_weights.push_back(mTGMMs_weights.at(j_i));
            }
        }
        tool->NewEventDisplay(groomed_particles, groomed_particles,
                              groomed_particles, groomed_jets, mTGMMs_particle_weights,
                              0, groomed_jets_params, groomed_jet_weights, "mTGMMs_groomed", event_iter);
    }
    if (mUMM_on) {
        // find groomed components for mUMM subjets
        vector<int> indices;
        double pT;
        vecPseudoJet groomed_particles;
        FindLeadingJet(leading_particles,
                       mUMM_jets, mUMM_particle_weights,
                       tool, indices, pT);
        for (unsigned int p_i = 0; p_i < leading_particles.size(); p_i++) {
            double max_weight = -1;
            int belongs_idx = -1;
            for (unsigned int j_i = 0; j_i < mUMM_jets.size(); j_i++) {
                double current_weight = mUMM_particle_weights.at(p_i).at(j_i);
                if (current_weight > max_weight) {
                    max_weight = current_weight;
                    belongs_idx = j_i;
                }
            }
            if (max_weight >= 0 && ((1.0*indices.at(belongs_idx)/mUMM_jets.size()) < (1-X))) {
                groomed_particles.push_back(leading_particles.at(p_i));
            }
        }

        // UNFINISHED DISPLAY ALL JETS AND PARTICLES AND GROOMED ONES
        tool->NewEventDisplayUniform(leading_particles, leading_particles,
                                     leading_particles, mUMM_jets, mUMM_particle_weights,
                                     0, mUMM_weights, "mUMM_ungroomed", event_iter);
        vecPseudoJet groomed_jets;
        //vector<vector<double>> groomed_weights; not necessary it turns out,
        // since these just get logged
        vector<double> groomed_jet_weights;
        for (unsigned int j_i = 0; j_i < mUMM_jets.size(); j_i++) {
            if ((1.0*indices.at(j_i)/mUMM_jets.size()) < (1-X)) {
                groomed_jets.push_back(mUMM_jets.at(j_i));
                groomed_jet_weights.push_back(mUMM_weights.at(j_i));
            }
        }
        tool->NewEventDisplayUniform(groomed_particles, groomed_particles,
                                     groomed_particles, groomed_jets, mUMM_particle_weights,
                                     0, groomed_jet_weights, "mUMM_groomed", event_iter);
    }

    delete m_jet_def_seed;
}

void FuzzyAnalysis::SubstructureStudy(vecPseudoJet ca_jets,
                                      __attribute__((unused)) vecPseudoJet antikt_jets,
                                      __attribute__((unused)) int event_iter) {
    vecPseudoJet particles_for_jets = ca_jets[0].constituents();
    cout << "Constituents: " << particles_for_jets.size();

    unsigned int max_num_clusters = 5;
    for (unsigned int num_clusters = 2; num_clusters <= max_num_clusters; num_clusters++) {
        vecPseudoJet random_seeds = randomSubvector(particles_for_jets, num_clusters);

        bool mUMM_on = true;
        bool mGMM_on = true;
        bool mGMMc_on = true;
        bool mGMMs_on = true;
        bool mTGMM_on = true;
        bool mTGMMs_on = true;

        double subs_sigma = 0.3;
        tool->SetDefaultSigma(MatTwo(subs_sigma*subs_sigma, 0, 0, subs_sigma*subs_sigma));
        tool->SetMergeDistance(0.05);

        unsigned int dummy;
        FuzzyTools::EventJet event_jet_type = FuzzyTools::NONE;
        FuzzyTools::PostProcess post_processing_method = FuzzyTools::NO_POST_PROCESS;
        double event_jet_weight = 0;

        // Fuzzy Jets: mGMMs --------------------
        vector<vector<double> > mGMMs_particle_weights;
        vector<MatTwo> mGMMs_jets_params;
        vector<double> mGMMs_weights;
        vecPseudoJet mGMMs_jets;
        vector<int> mGMMs_indices;
        double max_pT_mGMMs;
        if(mGMMs_on) {
            DoMGMMJetFinding(particles_for_jets, random_seeds,
                             f_learn_weights, true, false,
                             mGMMs_indices, max_pT_mGMMs,
                             tool, mGMMs_jets, mGMMs_particle_weights,
                             mGMMs_jets_params, mGMMs_weights, dummy,
                             event_jet_type, post_processing_method, event_jet_weight);
        }
        int lead_mGMMs_index = mGMMs_indices.size() ? mGMMs_indices.at(0) : -1;


        // Fuzzy Jets: mGMMc --------------------
        vector<vector<double> > mGMMc_particle_weights;
        vector<MatTwo> mGMMc_jets_params;
        vector<double> mGMMc_weights;
        vecPseudoJet mGMMc_jets;
        vector<int> mGMMc_indices;
        double max_pT_mGMMc;
        if(mGMMc_on) {
            DoMGMMCJetFinding(particles_for_jets, random_seeds,
                              f_learn_weights, false,
                              mGMMc_indices, max_pT_mGMMc,
                              tool, mGMMc_jets, mGMMc_particle_weights,
                              mGMMc_jets_params, mGMMc_weights, dummy,
                              event_jet_type, post_processing_method, event_jet_weight);
        }
        int lead_mGMMc_index = mGMMc_indices.size() ? mGMMc_indices.at(0) : -1;

        // Fuzzy Jets: mTGMMs -------------------
        vector<vector<double > > mTGMMs_particle_weights;
        vector<MatTwo> mTGMMs_jets_params;
        vector<double> mTGMMs_weights;
        vecPseudoJet mTGMMs_jets;
        vector<int> mTGMMs_indices;
        double max_pT_mTGMMs;
        if(mTGMMs_on) {
            DoMTGMMJetFinding(particles_for_jets, random_seeds,
                              f_learn_weights, true, false,
                              f_size, mTGMMs_indices, max_pT_mTGMMs,
                              tool, mTGMMs_jets, mTGMMs_particle_weights,
                              mTGMMs_jets_params, mTGMMs_weights, dummy,
                              event_jet_type, post_processing_method, event_jet_weight);
        }
        int lead_mTGMMs_index = mTGMMs_indices.size() ? mTGMMs_indices.at(0) : -1;

        // Fuzzy Jets: mGMM ---------------------
        vector<vector<double> >mGMM_particle_weights;
        vector<MatTwo> mGMM_jets_params;
        vector<double> mGMM_weights;
        vecPseudoJet mGMM_jets;
        vector<int> mGMM_indices;
        double max_pT_mGMM;
        if(mGMM_on) {
            DoMGMMJetFinding(particles_for_jets, random_seeds,
                             f_learn_weights, false, false,
                             mGMM_indices, max_pT_mGMM,
                             tool, mGMM_jets, mGMM_particle_weights,
                             mGMM_jets_params, mGMM_weights, dummy,
                             event_jet_type, post_processing_method, event_jet_weight);
        }
        int lead_mGMM_index = mGMM_indices.size() ? mGMM_indices.at(0) : -1;

        // Fuzzy Jets: mUMM ---------------------
        vector<vector<double> > mUMM_particle_weights;
        vector<double> mUMM_weights;
        vecPseudoJet mUMM_jets;
        vector<int> mUMM_indices;
        double max_pT_mUMM;
        if(mUMM_on) {
            DoMUMMJetFinding(particles_for_jets, random_seeds,
                             f_learn_weights, f_size, false,
                             mUMM_indices, max_pT_mUMM, tool, mUMM_jets,
                             mUMM_particle_weights, mUMM_weights, dummy,
                             event_jet_type, post_processing_method, event_jet_weight);
        }
        __attribute__((unused)) int lead_mUMM_index = mUMM_indices.size() ? mUMM_indices.at(0) : -1;

        // Fuzzy Jets: mTGMM --------------------
        vector<vector<double> > mTGMM_particle_weights;
        vector<double> mTGMM_weights;
        vecPseudoJet mTGMM_jets;
        vector<MatTwo> mTGMM_jets_params;
        vector<int> mTGMM_indices;
        double max_pT_mTGMM;
        if(mTGMM_on) {
            DoMTGMMJetFinding(particles_for_jets, random_seeds,
                              f_learn_weights, false, false, f_size,
                              mTGMM_indices, max_pT_mTGMM, tool, mTGMM_jets,
                              mTGMM_particle_weights, mTGMM_jets_params,
                              mTGMM_weights, dummy, event_jet_type,
                              post_processing_method, event_jet_weight);
        }
        int lead_mTGMM_index = mTGMM_indices.size() ? mTGMM_indices.at(0) : -1;

        if(!mGMM_weights.size()) {
            mGMM_on = false;
        }
        if(!mGMMs_weights.size()) {
            mGMMs_on = false;
        }
        if(!mGMMc_weights.size()) {
            mGMMc_on = false;
        }
        if(!mTGMM_weights.size()) {
            mTGMM_on = false;
        }
        if(!mTGMMs_weights.size()) {
            mTGMMs_on = false;
        }
        if(!mUMM_weights.size()) {
            mUMM_on = false;
        }

        std::stringstream ss;

        // do event displays
        if (mGMM_on) {
            ss.str(std::string());
            ss << directory_prefix << "SubstructureEvent" << "_mGMM_" << event_iter;
            if (num_clusters == 2) {
                ss << ".pdf(";
            } else {
                if (num_clusters == max_num_clusters) {
                    ss << ".pdf)";
                } else {
                    ss << ".pdf";
                }
            }
            std::string file = ss.str();
            ss.str(std::string());
            ss << "Substructure mGMM k=" << num_clusters;
            std::string title = ss.str();

            tool->SubsEventDisplay(particles_for_jets, mGMM_jets, mGMM_particle_weights,
                                   lead_mGMM_index, mGMM_jets_params, file, title);
        }
        if (mGMMs_on) {
            ss.str(std::string());
            ss << directory_prefix << "SubstructureEvent" << "_mGMMs_" << event_iter;
            if (num_clusters == 2) {
                ss << ".pdf(";
            } else {
                if (num_clusters == max_num_clusters) {
                    ss << ".pdf)";
                } else {
                    ss << ".pdf";
                }
            }
            std::string file = ss.str();
            ss.str(std::string());
            ss << "Substructure mGMMs k=" << num_clusters;
            std::string title = ss.str();
            tool->SubsEventDisplay(particles_for_jets, mGMMs_jets, mGMMs_particle_weights,
                                   lead_mGMMs_index, mGMMs_jets_params, file, title);
        }
        if (mGMMc_on) {
            ss.str(std::string());
            ss << directory_prefix << "SubstructureEvent" << "_mGMMc_" << event_iter;
            if (num_clusters == 2) {
                ss << ".pdf(";
            } else {
                if (num_clusters == max_num_clusters) {
                    ss << ".pdf)";
                } else {
                    ss << ".pdf";
                }
            }
            std::string file = ss.str();
            ss.str(std::string());
            ss << "Substructure mGMMc k=" << num_clusters;
            std::string title = ss.str();
            tool->SubsEventDisplay(particles_for_jets, mGMMc_jets, mGMMc_particle_weights,
                                   lead_mGMMc_index, mGMMc_jets_params, file, title);
        }
        if (mTGMM_on) {
            ss.str(std::string());
            ss << directory_prefix << "SubstructureEvent" << "_mTGMM_" << event_iter;
            if (num_clusters == 2) {
                ss << ".pdf(";
            } else {
                if (num_clusters == max_num_clusters) {
                    ss << ".pdf)";
                } else {
                    ss << ".pdf";
                }
            }
            std::string file = ss.str();
            ss.str(std::string());
            ss << "Substructure mTGMM k=" << num_clusters;
            std::string title = ss.str();
            tool->SubsEventDisplay(particles_for_jets, mTGMM_jets, mTGMM_particle_weights,
                                   lead_mTGMM_index, mTGMM_jets_params, file, title);
        }
        if (mTGMMs_on) {
            ss.str(std::string());
            ss << directory_prefix << "SubstructureEvent" << "_mTGMMs_" << event_iter;
            if (num_clusters == 2) {
                ss << ".pdf(";
            } else {
                if (num_clusters == max_num_clusters) {
                    ss << ".pdf)";
                } else {
                    ss << ".pdf";
                }
            }
            std::string file = ss.str();
            ss.str(std::string());
            ss << "Substructure mTGMMs k=" << num_clusters;
            std::string title = ss.str();
            tool->SubsEventDisplay(particles_for_jets, mTGMMs_jets, mTGMMs_particle_weights,
                                   lead_mTGMMs_index, mTGMMs_jets_params, file, title);
        }
    }
}

void FuzzyAnalysis::SetupHistosMap() {
    static const std::string algs_arr[] =
        {"mGMM", "mGMMs", "mGMMc", "mTGMM", "mTGMMs", "mUMM"};
    std::vector<std::string> algs(algs_arr, algs_arr+sizeof(algs_arr) / sizeof(algs_arr[0]));
    for (unsigned int alg_iter = 0; alg_iter < algs.size(); alg_iter++) {
        std::stringstream ss;

        ss.str(std::string());
        ss << algs[alg_iter] << "_hs";
        map_weight_vecs[ss.str()] = std::vector<float>();

        ss.str(std::string());
        ss << algs[alg_iter] << "_pu";
        map_weight_vecs[ss.str()] = std::vector<float>();

    }
}

void FuzzyAnalysis::WriteHistosMap() {
#ifdef WITHROOT
    static const std::string algs_arr[] =
        {"mGMM", "mGMMs", "mGMMc", "mTGMM", "mTGMMs", "mUMM"};
    std::vector<std::string> algs(algs_arr, algs_arr+sizeof(algs_arr) / sizeof(algs_arr[0]));
    for (unsigned int alg_iter = 0; alg_iter < algs.size(); alg_iter++) {
        std::stringstream ss;

        ss.str(std::string());
        ss << algs[alg_iter] << "_hs";
        std::string hs_name = ss.str();
        std::vector<float> const& hs_weights = map_weight_vecs[hs_name];
        TH1F hs_hist(TString::Format("Hard Scatter Weights %s", algs[alg_iter].c_str()),
                     TString::Format("Hard Scatter Weights for %s", algs[alg_iter].c_str()),
                     48, 0.04, 1);
        for (unsigned int weight_iter = 0; weight_iter < hs_weights.size(); weight_iter++) {
            hs_hist.Fill(hs_weights[weight_iter]);
        }
        hs_hist.Write();

        ss.str(std::string());
        ss << algs[alg_iter] << "_pu";
        std::string pu_name = ss.str();
        std::vector<float> const& pu_weights = map_weight_vecs[pu_name];
        TH1F pu_hist(TString::Format("Pileup Weights %s", algs[alg_iter].c_str()),
                     TString::Format("Pileup Weights for %s", algs[alg_iter].c_str()),
                     48, 0.04, 1);
        for (unsigned int weight_iter = 0; weight_iter < pu_weights.size(); weight_iter++) {
            pu_hist.Fill(pu_weights[weight_iter]);
        }
        pu_hist.Write();

        // don't write canvases if we are on the batch,
        // we can do that later
        if (batched) continue;

        TCanvas *c = new TCanvas(TString::Format("HSPU Weight Comparison %s", algs[alg_iter].c_str()),
                                 "Hard Scatter and Pileup Weight Comparison", 800, 800);
        pu_hist.SetLineColor(kBlue);
        pu_hist.SetFillStyle(3004);
        pu_hist.SetFillColor(kBlue);

        hs_hist.SetLineColor(kRed);
        hs_hist.SetFillStyle(3004);
        hs_hist.SetFillColor(kRed);

        if(hs_hist.Integral(-1, hs_hist.GetNbinsX()+1) > 0) {
            hs_hist.Scale(1./hs_hist.Integral(-1, hs_hist.GetNbinsX()+1));
        }
        if(pu_hist.Integral(-1, pu_hist.GetNbinsX()+1) > 0) {
            pu_hist.Scale(1./pu_hist.Integral(-1, pu_hist.GetNbinsX()+1));
        }

        hs_hist.Draw("");
        pu_hist.Draw("same");

        TLegend *leggaa = new TLegend(0.6, 0.7, 0.9, 0.8);
        leggaa->SetTextFont(42);
        leggaa->AddEntry(&pu_hist, "Pileup weights", "l");
        leggaa->AddEntry(&hs_hist, "Hard scatter weights", "l");
        leggaa->SetFillColor(0);
        leggaa->SetFillStyle(0);
        leggaa->SetBorderSize(0);
        leggaa->Draw();

        c->Write();
        delete c;
        delete leggaa;

    }
#endif
}

// End
void FuzzyAnalysis::End(){
#ifdef WITHROOT
    t_t->Write();
    delete t_f;

    if (!_non_standard_study) {
        WriteHistosMap();
    }
#endif
    return;
}


#ifdef NEW
// Test new fuzzy tools
void FuzzyAnalysis::AnalyzeEventNew(int event_iter, Pythia8::Pythia* pythia8, Pythia8::Pythia* pythia_MB, int NPV){
    // generate a new event
    if (!pythia8->next()) return;

    // lazy determine if should generate new pythia pileup event
    if (NPV && !pythia_MB->next()) return;
    if(f_debug) cout << "FuzzyAnalysis::AnalyzeEvent Event Number " << event_iter << endl;

    // reset branches
    ResetBranches();

    // new event-----------------------
    fTEventNumber = event_iter;
    std::vector <fastjet::PseudoJet>           particles_for_jets;
    vector<fastjet::PseudoJet> tops;
    fastjet::PseudoJet dl;
    tops.push_back(dl);
    fastjet::PseudoJet dl2;
    tops.push_back(dl2);

    fTNPV = NPV;
    // Pileup loop -------------------------------------------------------------
    double px, py, pz, e;
    for (int pileup_idx = 0; pileup_idx < NPV; ++pileup_idx) {
        for (unsigned int particle_idx = 0; particle_idx < (unsigned)pythia_MB->event.size(); ++particle_idx) {
            if (!pythia_MB->event[particle_idx].isFinal()) continue;
            if (fabs(pythia_MB->event[particle_idx].id())==11) continue;
            if (fabs(pythia_MB->event[particle_idx].id())==12) continue;
            if (fabs(pythia_MB->event[particle_idx].id())==14) continue;
            if (fabs(pythia_MB->event[particle_idx].id())==13) continue;
            if (fabs(pythia_MB->event[particle_idx].id())==16) continue;
            //if (pythia_MB->event[particle_idx].pT() < 0.5)     continue;
            px = pythia_MB->event[particle_idx].px();
            py = pythia_MB->event[particle_idx].py();
            pz = pythia_MB->event[particle_idx].pz();
            e  = pythia_MB->event[particle_idx].e();

            fastjet::PseudoJet p(px, py, pz, e);
            p.reset_PtYPhiM(p.pt(), p.rapidity(), p.phi(), 0.);
            // note that we don't really keep the particle number! only store particle_idx, should really store particle_idx and pileup_idx
            p.set_user_info(new MyUserInfo(pythia_MB->event[particle_idx].id(),particle_idx,pythia_MB->event[particle_idx].charge(),true,false));
            particles_for_jets.push_back(p);
        }
        if (!pythia_MB->next()) continue;
    }

    // Particle loop -----------------------------------------------------------
    // The Pythia event listing contains a lot more than we want to process,
    // we prune out certain particles (muons / neutrinos) and only add final
    // state particles
    for (unsigned int particle_idx = 0; particle_idx < (unsigned) pythia8->event.size(); ++particle_idx){
        px = pythia8->event[particle_idx].px();
        py = pythia8->event[particle_idx].py();
        pz = pythia8->event[particle_idx].pz();
        e  = pythia8->event[particle_idx].e();
        fastjet::PseudoJet p(px, py, pz, e);
        p.set_user_info(new MyUserInfo(pythia8->event[particle_idx].id(),particle_idx,pythia8->event[particle_idx].charge(),false,false));

        // In reality we should be more careful about finding tops,
        // but this will do for now. In the future consider refactoring
        // and tracing to find a top quark with no daughters
        if (pythia8->event[particle_idx].id()  ==6) tops[0]=p;
        if (pythia8->event[particle_idx].id()  ==-6) tops[1]=p;

        // prune uninteresting particles
        if (!pythia8->event[particle_idx].isFinal() )      continue; // only final state
        if (fabs(pythia8->event[particle_idx].id())  ==11) continue; // ...   electron
        if (fabs(pythia8->event[particle_idx].id())  ==12) continue; // prune nu-e
        if (fabs(pythia8->event[particle_idx].id())  ==13) continue; // ...   mu
        if (fabs(pythia8->event[particle_idx].id())  ==14) continue; // ...   nu-mu
        if (fabs(pythia8->event[particle_idx].id())  ==16) continue; // ...   nu-tau
        if (pythia8->event[particle_idx].pT()       < 0.5) continue; // ...   low pT

        particles_for_jets.push_back(p);

    } // end particle loop -----------------------------------------------
    // have particles now, time to try clustering

    tool->SetMergeDistance(0.05);
    double sigma_squared = 0.5;
    tool->SetDefaultSigma(MatTwo(sigma_squared, 0, 0, sigma_squared));
    tool->SetClusteringMode(FuzzyTools::RECOMBINATION);
    tool->SetSeeds(particles_for_jets);
    tool->SetR(1.0);
    tool->SetLogLogLikelihoodLimit(-10);
    vector<GaussianKernel *> gaussian_jets = MakeGaussianKernels(*tool);
}
#endif

void FuzzyAnalysis::JetMultiplicityStudy(Pythia8::Pythia *pythia8, Pythia8::Pythia *pythia_MB) {
    batched = true;
    stop_before_clustering = true;
    static const int NPVs_arr[] = {0, 10, 20, 30, 40, 50, 60, 70, 80};
    static const float NPVs_float_arr[] = {0, 10, 20, 30, 40, 50, 60, 70, 80};
    static const int cuts_arr[] = {15, 20, 25, 30, 35, 40, 45, 50};

    std::vector<float> NPVs_float(NPVs_float_arr, NPVs_float_arr+sizeof(NPVs_float_arr)/sizeof(NPVs_float_arr[0]));
    std::vector<int> NPVs(NPVs_arr, NPVs_arr+sizeof(NPVs_arr)/sizeof(NPVs_arr[0]));
    std::vector<int> cuts(cuts_arr, cuts_arr+sizeof(cuts_arr)/sizeof(cuts_arr[0]));

    unsigned int n_samples = 30;
    BOOST_FOREACH(int cut, cuts) {
        pT_min = cut;
        TCanvas *c = new TCanvas("scatter_canv", "", 600, 400);
        c->cd();
        std::vector<float> avg_mults;
        std::vector<float> avg_mults_stdevs;

        BOOST_FOREACH(int NPV, NPVs) {
            std::cout << NPV << std::endl;
            std::vector<float> samples;
            for (unsigned int s_i = 0; s_i < n_samples; s_i++) {
                AnalyzeEvent(0, pythia8, pythia_MB, NPV);
                samples.push_back(fTn_jet_seeds);

            }
            float sum = std::accumulate(samples.begin(), samples.end(), 0.0);
            float mean = sum / samples.size();

            float sq_sum = std::inner_product(samples.begin(), samples.end(), samples.begin(), 0.0);
            float stdev = std::sqrt(sq_sum / samples.size() - mean * mean);
            avg_mults.push_back(mean);
            avg_mults_stdevs.push_back(stdev);
        }

        TGraphErrors *gr = new TGraphErrors(avg_mults.size(), &NPVs_float.at(0), &avg_mults.at(0), 0, &avg_mults_stdevs.at(0));
        gr->SetTitle(TString::Format("Average number of seed jets with a cut of %d GeV",cut));
        gr->SetMarkerStyle(22);
        gr->GetXaxis()->SetTitle("NPV");
        gr->GetYaxis()->SetTitle("Average # jets");
        gr->Draw("AP");
        std::stringstream ss;
        ss.str(std::string());
        ss << directory_prefix << "JetMultiplicity_" << cut << "cut.pdf";
        c->Print(ss.str().c_str(), "pdf");
        delete c;
    }
}

void FuzzyAnalysis::PiFixStudy() {
    gStyle->SetOptStat(0);
    // set up a few test particles and initial jets
    vecPseudoJet test_particles;
    vecPseudoJet test_jets;

    fastjet::PseudoJet p(0,0,0,0);

    p.reset_PtYPhiM(40, 0, 0.1, 0);
    p.set_user_info(new MyUserInfo(0, 0, 0, false, false));
    test_particles.push_back(p);
    test_jets.push_back(p);


    p.reset_PtYPhiM(20, 0, 6.2, 0);
    p.set_user_info(new MyUserInfo(0, 0, 0, false, false));
    test_particles.push_back(p);

    // Separated particle
    p.reset_PtYPhiM(10, -1, 5.7, 0);
    p.set_user_info(new MyUserInfo(0, 0, 0, false, false));
    test_particles.push_back(p);
    test_jets.push_back(p);

    bool mGMM_on = true;
    bool mGMMc_on = true;

    // Fuzzy Jets: mGMMc --------------------
    vector<vector<double> > mGMMc_particle_weights;
    vector<MatTwo> mGMMc_jets_params;
    vector<double> mGMMc_weights;
    vecPseudoJet mGMMc_jets;
    vector<int> mGMMc_indices;
    double max_pT_mGMMc;
    if(mGMMc_on) {
        DoMGMMCJetFinding(test_particles, test_jets,
                          f_learn_weights, do_recombination,
                          mGMMc_indices, max_pT_mGMMc,
                          tool, mGMMc_jets, mGMMc_particle_weights,
                          mGMMc_jets_params, mGMMc_weights, fTmGMMc_iter,
                          FuzzyTools::NONE, FuzzyTools::NO_POST_PROCESS, 0);
    }
    int lead_mGMMc_index = mGMMc_indices.size() ? mGMMc_indices.at(0) : -1;

    // Fuzzy Jets: mGMM ---------------------
    vector<vector<double> >mGMM_particle_weights;
    vector<MatTwo> mGMM_jets_params;
    vector<double> mGMM_weights;
    vector<fastjet::PseudoJet> mGMM_jets;
    vector<int> mGMM_indices;
    double max_pT_mGMM;
    if(mGMM_on) {
        DoMGMMJetFinding(test_particles, test_jets,
                         f_learn_weights, false, do_recombination,
                         mGMM_indices, max_pT_mGMM,
                         tool, mGMM_jets, mGMM_particle_weights,
                         mGMM_jets_params, mGMM_weights, fTmGMM_iter,
                         FuzzyTools::NONE, FuzzyTools::NO_POST_PROCESS, 0);
    }
    int lead_mGMM_index = mGMM_indices.size() ? mGMM_indices.at(0) : -1;

    tool->NewEventDisplay(test_particles,
                          test_jets, test_particles,
                          mGMM_jets,
                          mGMM_particle_weights,
                          lead_mGMM_index,
                          mGMM_jets_params,
                          mGMM_weights,
                          "mGMM_mod_pi",
                          0);

    tool->NewEventDisplay(test_particles,
                          test_jets, test_particles,
                          mGMMc_jets,
                          mGMMc_particle_weights,
                          lead_mGMMc_index,
                          mGMMc_jets_params,
                          mGMMc_weights,
                          "mGMMc_mod_pi",
                          0);
}

// Analyze
void FuzzyAnalysis::AnalyzeEvent(int event_iter, Pythia8::Pythia* pythia8, Pythia8::Pythia* pythia_MB, int NPV, vecPseudoJet const & hs_particles_in, vecPseudoJet const & tops_in){
    bool use_arg_particles = false;
    if (hs_particles_in.size() != 0) {
        // don't touch pythia8
        use_arg_particles = true;
    }

    if(f_debug) cout << "FuzzyAnalysis::AnalyzeEvent Begin " << endl;

    bool hard_scatter_on = true;
    bool pileup_on = true;
    bool remove_negative_towers = true;

    // generate a new event
    if (!use_arg_particles) {
        if (!pythia8->next()) return;
    }

    // lazy determine if should generate new pythia pileup event
    if(f_debug) cout << "FuzzyAnalysis::AnalyzeEvent Event Number " << event_iter << endl;

    // reset branches
    ResetBranches();

    // new event-----------------------
    fTEventNumber = event_iter;
    vecPseudoJet particles_for_jets;
    vecPseudoJet pileup_for_jets;
    vecPseudoJet tops;
    vecPseudoJet ftops;
    float max_rap = -1;

    fTNPV = NPV;

    if (pileup_on) {
        generate_particles_multi(pileup_for_jets, ftops, max_rap, true,
                                 pythia_MB, static_cast<size_t>(NPV));
    }
    if (!use_arg_particles && hard_scatter_on) {
        generate_particles(particles_for_jets, tops, max_rap, false, pythia8);
    } else {
        particles_for_jets = hs_particles_in;
        assert(tops_in.size());
        tops = tops_in;
        // have to calculate max_rap from the input particles
        for (unsigned int p_i = 0; p_i < particles_for_jets.size(); p_i++) {
            fastjet::PseudoJet p = particles_for_jets.at(p_i);
            max_rap = max_rap < fabs(p.rap()) ? fabs(p.rap()) : max_rap;
        }
    }

    vecPseudoJet charged_hs = ChargedParticles(particles_for_jets, 2.5);
    vecPseudoJet uncharged_hs = UnchargedParticles(particles_for_jets, 2.5);
    vecPseudoJet charged_pileup = ChargedParticles(pileup_for_jets, 2.5);
    vecPseudoJet uncharged_pileup = UnchargedParticles(pileup_for_jets, 2.5);

    vecPseudoJet all_uncharged;
    for (unsigned int p_i = 0; p_i < uncharged_hs.size(); p_i++) {
        all_uncharged.push_back(uncharged_hs.at(p_i));
    }
    for (unsigned int p_i = 0; p_i < uncharged_pileup.size(); p_i++) {
        all_uncharged.push_back(uncharged_pileup.at(p_i));
    }

    vecPseudoJet pre_disc_for_jets = all_uncharged;
    for (unsigned int p_i = 0; p_i < charged_hs.size(); p_i++) {
        pre_disc_for_jets.push_back(charged_hs.at(p_i));
    }

    double disc_max_eta = 5;
    int bins_eta = 100;
    int bins_phi = 60;
    vecPseudoJet discretized_all_for_jets =
        DiscretizeParticles(pre_disc_for_jets, bins_phi, bins_eta, disc_max_eta,
                            _do_tower_subtraction);

    vecPseudoJet all_particles_for_jets = discretized_all_for_jets;
    //for (unsigned int p_i = 0; p_i < charged_hs.size(); p_i++) {
    //    all_particles_for_jets.push_back(charged_hs.at(p_i));
    //}

    particles_for_jets = all_particles_for_jets;
    vecPseudoJet clean_particles_for_jets;
    for (unsigned int p_i = 0; p_i < particles_for_jets.size(); p_i++) {
        if (particles_for_jets.at(p_i).pt() > 3*FLT_EPSILON) {
            clean_particles_for_jets.push_back(particles_for_jets.at(p_i));
        }
    }

    // ======================================
    // Calculate rho ------------------------
    // ======================================
    // use charged pileup and all uncharged for background estimation
    vecPseudoJet rho_estimator_particles;
    for (unsigned int p_i = 0; p_i < all_uncharged.size(); p_i++) {
        rho_estimator_particles.push_back(all_uncharged.at(p_i));
    }
    //for (unsigned int p_i = 0; p_i < charged_pileup.size(); p_i++) {
    //    rho_estimator_particles.push_back(charged_pileup.at(p_i));
    //}

    vecPseudoJet discretized_rho_estimator_particles =
        DiscretizeParticles(rho_estimator_particles, bins_phi, bins_eta, disc_max_eta,
                            false); // don't keep the zero pT particles!
                                    // fastjet will die on area computation

    fTrho = estimateRho(discretized_rho_estimator_particles);

    // Area Definition ---------------------
    fastjet::GhostedAreaSpec area_spec(max_rap);
    fastjet::AreaDefinition area_def(fastjet::active_area, area_spec);
    fastjet::AreaDefinition area_def_exp = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts);

    // N Subjettiness
    fastjet::contrib::Nsubjettiness n_subjettiness_1(1, fastjet::contrib::Njettiness::kt_axes, 1, 1, 1);
    fastjet::contrib::Nsubjettiness n_subjettiness_2(2, fastjet::contrib::Njettiness::kt_axes, 1, 1, 1);
    fastjet::contrib::Nsubjettiness n_subjettiness_3(3, fastjet::contrib::Njettiness::kt_axes, 1, 1, 1);

    // large-R jets: C/A --------------------
    {
        fastjet::ClusterSequenceArea cs_large_r_ca(clean_particles_for_jets, *m_jet_def_large_r_ca, area_def);
        vecPseudoJet my_jets_large_r_ca = fastjet::sorted_by_pt(cs_large_r_ca.inclusive_jets(pT_min));
        if (my_jets_large_r_ca.size() != 0) {
            fTCA_pufrac = JetPuFracFastjet(my_jets_large_r_ca[0]);
            fTCA_m_pu = JetPuMassFastjet(my_jets_large_r_ca[0]);
        }
    }
    fastjet::ClusterSequenceArea cs_large_r_ca(clean_particles_for_jets, *m_jet_def_large_r_ca, area_def_exp);
    vecPseudoJet my_jets_large_r_ca = fastjet::sorted_by_pt(cs_large_r_ca.inclusive_jets(pT_min));

    // this is a very temporary fix, it appears that there are no jets sometimes with pT at least
    if (my_jets_large_r_ca.size() != 0) {
        fTCA_m = my_jets_large_r_ca[0].m();
        fTCA_pt = my_jets_large_r_ca[0].pt();
        fTCA_area = my_jets_large_r_ca[0].area();
        if (my_jets_large_r_ca.size() >= 2) {
            fTCA_dr = my_jets_large_r_ca[0].delta_R(my_jets_large_r_ca[1]);
        }
        CA_nsubjettiness.push_back(n_subjettiness_1(my_jets_large_r_ca[0]));
        CA_nsubjettiness.push_back(n_subjettiness_2(my_jets_large_r_ca[0]));
        CA_nsubjettiness.push_back(n_subjettiness_3(my_jets_large_r_ca[0]));
    }

    // anti-kt R:1.0 trimmed ----------------
    fastjet::Filter filter_two(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2),
                               fastjet::SelectorPtFractionMin(0.05));
    fastjet::Filter filter_three(fastjet::JetDefinition(fastjet::kt_algorithm, 0.3),
                                 fastjet::SelectorPtFractionMin(0.05));

    {
        fastjet::ClusterSequenceArea cs_large_r_antikt(clean_particles_for_jets, *m_jet_def_large_r_antikt, area_def);
        vecPseudoJet my_jets_large_r_antikt = fastjet::sorted_by_pt(cs_large_r_antikt.inclusive_jets(pT_min));
        fastjet::PseudoJet lead_akt = my_jets_large_r_antikt[0];
        fastjet::PseudoJet lead_akt_filter_two = filter_two(lead_akt);
        fastjet::PseudoJet lead_akt_filter_three = filter_three(lead_akt);

        fTantikt_m_pu = JetPuMassFastjet(lead_akt);
        fTantikt_pufrac = JetPuFracFastjet(lead_akt);
        fTantikt_m_pu_trimmed_two = JetPuMassFastjet(lead_akt_filter_two);
        fTantikt_pufrac_trimmed_two = JetPuFracFastjet(lead_akt_filter_two);
        fTantikt_m_pu_trimmed_three = JetPuMassFastjet(lead_akt_filter_three);
        fTantikt_pufrac_trimmed_three = JetPuFracFastjet(lead_akt_filter_three);
    }

    fastjet::ClusterSequenceArea cs_large_r_antikt(clean_particles_for_jets, *m_jet_def_large_r_antikt, area_def_exp);
    vecPseudoJet my_jets_large_r_antikt = fastjet::sorted_by_pt(cs_large_r_antikt.inclusive_jets(pT_min));

    if (my_jets_large_r_antikt.size() != 0) {
        fastjet::PseudoJet lead_akt = my_jets_large_r_antikt[0];
        fastjet::PseudoJet lead_akt_filter_two = filter_two(lead_akt);
        fastjet::PseudoJet lead_akt_filter_three = filter_three(lead_akt);

        fTantikt_area = lead_akt.area();
        fTantikt_m = lead_akt.m();
        if (my_jets_large_r_antikt.size() >= 2) {
            fTantikt_m_second = my_jets_large_r_antikt.at(1).m();
        }

        fTantikt_pt = lead_akt.pt();

        fTantikt_area_trimmed_two = lead_akt_filter_two.area();
        fTantikt_m_trimmed_two = lead_akt_filter_two.m();
        fTantikt_pt_trimmed_two = lead_akt_filter_two.pt();

        fTantikt_area_trimmed_three = lead_akt_filter_three.area();
        fTantikt_m_trimmed_three = lead_akt_filter_three.m();
        if (my_jets_large_r_antikt.size() >= 2) {
            fTantikt_m_second_trimmed_three = filter_three(my_jets_large_r_antikt.at(1)).m();
        }
        fTantikt_pt_trimmed_three = lead_akt_filter_three.pt();

        antikt_nsubjettiness.push_back(n_subjettiness_1(lead_akt));
        antikt_nsubjettiness.push_back(n_subjettiness_2(lead_akt));
        antikt_nsubjettiness.push_back(n_subjettiness_3(lead_akt));
        if (my_jets_large_r_antikt.size() >= 2) {
            fTantikt_dr = my_jets_large_r_antikt[0].delta_R(my_jets_large_r_antikt[1]);
        }

        for (unsigned int particle_iter = 0; particle_iter < clean_particles_for_jets.size(); particle_iter++) {
            fastjet::PseudoJet const& particle = clean_particles_for_jets[particle_iter];
            is_lead_antikt_constituent_vec.push_back(lead_akt.contains(particle));
        }
        is_lead_antikt_constituent_vec.clear(); // don't need this info, doesn't make sense
    }
    for (unsigned int particle_iter = 0; particle_iter < clean_particles_for_jets.size(); particle_iter++) {
        fastjet::PseudoJet const& particle = clean_particles_for_jets[particle_iter];
        is_pileup_vec.push_back(particle.user_info<MyUserInfo>().isPU());
    }
    is_pileup_vec.clear(); // this doesn't make sense anymore after discretizing

    bool do_substructure_study = false;
    if (do_substructure_study && !batched && event_iter < 10) {
        // generate particles for use in the substructure stuff
        fastjet::ClusterSequence cs_vlarge_r_antikt(clean_particles_for_jets, *m_jet_def_vlarge_r_antikt);
        fastjet::ClusterSequence cs_vlarge_r_ca(clean_particles_for_jets, *m_jet_def_vlarge_r_ca);
        vecPseudoJet my_jets_vlarge_r_antikt = fastjet::sorted_by_pt(cs_vlarge_r_antikt.inclusive_jets(pT_min));
        vecPseudoJet my_jets_vlarge_r_ca = fastjet::sorted_by_pt(cs_vlarge_r_ca.inclusive_jets(pT_min));
        if (!my_jets_vlarge_r_antikt.size() || !my_jets_vlarge_r_ca.size()) return;

        SubstructureStudy(my_jets_vlarge_r_ca, my_jets_vlarge_r_antikt, event_iter);
        return;
    }

    bool do_grooming_study = false;
    if (do_grooming_study && !batched) {
        if (!my_jets_large_r_antikt.size()) return;

        GroomingStudy(my_jets_large_r_antikt.at(0).constituents(), event_iter);
        return;
    }

    vecPseudoJet seeds = groomedInitialLocations(clean_particles_for_jets,
                                                 1.0, 0.2, 0.05, pT_min, fTrho,
                                                 fTseed_location_noise);
    fTn_jet_seeds = seeds.size(); // used for jet multiplicity studies
    if (stop_before_clustering)
        return;

    // ======================================
    // Tower based subtraction --------------
    // ======================================
    //std::cout << fTrho << std::endl;
    double area_per_cell = (2*M_PI) * (2*disc_max_eta) / (bins_phi * bins_eta);
    //std::cout << area_per_cell << " : " << bins_phi * bins_eta << " : "
    //          << bins_phi * bins_eta * area_per_cell << std::endl;
    double pt_before_subs = 0;
    double total_pt = 0;
    double positive_pt = 0;
    double pt_before_subs_inner = 0;
    double pt_before_subs_outer = 0;
    double total_pt_inner = 0;
    double total_pt_outer = 0;
    vecPseudoJet tower_subtracted_particles_for_jets;
    // particles_for_jets includes empty towers
    for (unsigned int p_i = 0; p_i < particles_for_jets.size(); p_i++) {
        fastjet::PseudoJet p_orig = particles_for_jets.at(p_i);
        if (_do_tower_subtraction) {
            double npt = p_orig.pt() - (area_per_cell * fTrho);
            pt_before_subs += p_orig.pt();
            total_pt += npt;
            if (fabs(p_orig.eta()) > 2.5) {
                total_pt_outer += npt;
                pt_before_subs_outer += p_orig.pt();
            } else {
                total_pt_inner += npt;
                pt_before_subs_inner += p_orig.pt();
            }
            if (npt > 0) {
                positive_pt += npt;
            }
            double eta = p_orig.eta();
            double phi = p_orig.phi();
            double charge = p_orig.user_info<MyUserInfo>().charge();
            fastjet::PseudoJet p_new;
            p_new.reset_PtYPhiM(fabs(npt), eta, phi, 0);
            p_new.set_user_info(new MyUserInfo(0, 0, charge, false,
                                               npt > 0 ? false : true));
            // negative tower removal
            if (!remove_negative_towers || npt > 0) {
                if (fabs(npt) > 3*FLT_EPSILON)
                    tower_subtracted_particles_for_jets.push_back(p_new);
            }
        } else {
            if (p_orig.pt() > 3*FLT_EPSILON)
                tower_subtracted_particles_for_jets.push_back(p_orig);
        }
    }

    //std::cout << total_pt << " : " << positive_pt << " : " << pt_before_subs << std::endl;
    //std::cout << fTrho * (2*disc_max_eta) * (2*M_PI) << std::endl;
    //std::cout << "Inner : " << total_pt_inner << " : " << pt_before_subs_inner << std::endl;
    //std::cout << "Outer : " << total_pt_outer << " : " << pt_before_subs_outer << std::endl;
    // ======================================
    // Various mixture models ---------------
    // ======================================
    tool->SetMergeDistance(0.05);
    double sigma_squared = 0.5;
    tool->SetDefaultSigma(MatTwo(sigma_squared, 0, 0, sigma_squared));

    // which jets to run
    bool mUMM_on = false;
    bool mGMM_on = false;
    bool mGMMc_on = true;
    bool mGMMs_on = false;
    bool mTGMM_on = false;
    bool mTGMMs_on = false;

    bool make_trace_diagrams = true && !batched;
    tool->SetMakeTraceDiagrams(make_trace_diagrams);

    fTevent_jet_strength = _event_jet_strength;
    fTevent_jet_rho_offset = _event_jet_rho_offset;

    // convention is that offsets should be positive, so that a reasonable
    // offset is <hs_rho>.
    double event_jet_weight = (fTrho - fTevent_jet_rho_offset) * fTevent_jet_strength;
    // don't allow negative likelihood backgrounds
    if (event_jet_weight < 0)
        event_jet_weight = 0;

    // Fuzzy Jets: mGMMs --------------------
    vector<vector<double> > mGMMs_particle_weights;
    vector<MatTwo> mGMMs_jets_params;
    vector<double> mGMMs_weights;
    vecPseudoJet mGMMs_jets;
    vector<int> mGMMs_indices;
    double max_pT_mGMMs;
    if(mGMMs_on) {
        DoMGMMJetFinding(tower_subtracted_particles_for_jets, seeds,
                         f_learn_weights, true, do_recombination,
                         mGMMs_indices, max_pT_mGMMs,
                         tool, mGMMs_jets, mGMMs_particle_weights,
                         mGMMs_jets_params, mGMMs_weights, fTmGMMs_iter,
                         _event_jet_type, _post_processing_method, event_jet_weight);
    }
    int lead_mGMMs_index = mGMMs_indices.size() ? mGMMs_indices.at(0) : -1;
    if (make_trace_diagrams && mGMMs_on) {
        tool->MakeTraceDiagram(tower_subtracted_particles_for_jets,
                               lead_mGMMs_index, "mGMMs", event_iter);
    }



    // Fuzzy Jets: mGMMc --------------------
    vector<vector<double> > mGMMc_particle_weights;
    vector<MatTwo> mGMMc_jets_params;
    vector<double> mGMMc_weights;
    vecPseudoJet mGMMc_jets;
    vector<int> mGMMc_indices;
    double max_pT_mGMMc;
    if(mGMMc_on) {
        DoMGMMCJetFinding(tower_subtracted_particles_for_jets, seeds,
                          f_learn_weights, do_recombination,
                          mGMMc_indices, max_pT_mGMMc,
                          tool, mGMMc_jets, mGMMc_particle_weights,
                          mGMMc_jets_params, mGMMc_weights, fTmGMMc_iter,
                          _event_jet_type, _post_processing_method, event_jet_weight);
    }
    int lead_mGMMc_index = mGMMc_indices.size() ? mGMMc_indices.at(0) : -1;
    if (make_trace_diagrams && mGMMc_on) {
        tool->MakeTraceDiagram(tower_subtracted_particles_for_jets,
                               lead_mGMMc_index, "mGMMc", event_iter);
    }


    // Fuzzy Jets: mTGMMs -------------------
    vector<vector<double > > mTGMMs_particle_weights;
    vector<MatTwo> mTGMMs_jets_params;
    vector<double> mTGMMs_weights;
    vecPseudoJet mTGMMs_jets;
    vector<int> mTGMMs_indices;
    double max_pT_mTGMMs;
    if(mTGMMs_on) {
        DoMTGMMJetFinding(tower_subtracted_particles_for_jets, seeds,
                          f_learn_weights, true, do_recombination,
                          f_size, mTGMMs_indices, max_pT_mTGMMs,
                          tool, mTGMMs_jets, mTGMMs_particle_weights,
                          mTGMMs_jets_params, mTGMMs_weights, fTmTGMMs_iter,
                          _event_jet_type, _post_processing_method, event_jet_weight);
    }
    int lead_mTGMMs_index = mTGMMs_indices.size() ? mTGMMs_indices.at(0) : -1;
    if (make_trace_diagrams && mTGMMs_on) {
        tool->MakeTraceDiagram(tower_subtracted_particles_for_jets,
                               lead_mTGMMs_index, "mTGMMs", event_iter);
    }

    // Fuzzy Jets: mGMM ---------------------
    vector<vector<double> >mGMM_particle_weights;
    vector<MatTwo> mGMM_jets_params;
    vector<double> mGMM_weights;
    vector<fastjet::PseudoJet> mGMM_jets;
    vector<int> mGMM_indices;
    double max_pT_mGMM;
    if(mGMM_on) {
        DoMGMMJetFinding(tower_subtracted_particles_for_jets, seeds,
                         f_learn_weights, false, do_recombination,
                         mGMM_indices, max_pT_mGMM,
                         tool, mGMM_jets, mGMM_particle_weights,
                         mGMM_jets_params, mGMM_weights, fTmGMM_iter,
                         _event_jet_type, _post_processing_method, event_jet_weight);
    }
    int lead_mGMM_index = mGMM_indices.size() ? mGMM_indices.at(0) : -1;
    if (make_trace_diagrams && mGMM_on) {
        tool->MakeTraceDiagram(tower_subtracted_particles_for_jets,
                               lead_mGMM_index, "mGMM", event_iter);
    }

    // Fuzzy Jets: mUMM ---------------------
    vector<vector<double> > mUMM_particle_weights;
    vector<double> mUMM_weights;
    vecPseudoJet mUMM_jets;
    vector<int> mUMM_indices;
    double max_pT_mUMM;
    if(mUMM_on) {
        DoMUMMJetFinding(tower_subtracted_particles_for_jets, seeds,
                         f_learn_weights, f_size, do_recombination,
                         mUMM_indices, max_pT_mUMM, tool, mUMM_jets,
                         mUMM_particle_weights, mUMM_weights, fTmUMM_iter,
                         _event_jet_type, _post_processing_method, event_jet_weight);
    }
    int lead_mUMM_index = mUMM_indices.size() ? mUMM_indices.at(0) : -1;
    if (make_trace_diagrams && mUMM_on) {
        tool->MakeTraceDiagram(tower_subtracted_particles_for_jets,
                               lead_mUMM_index, "mUMM", event_iter);
    }

    // Fuzzy Jets: mTGMM --------------------
    vector<vector<double> > mTGMM_particle_weights;
    vector<double> mTGMM_weights;
    vecPseudoJet mTGMM_jets;
    vector<MatTwo> mTGMM_jets_params;
    vector<int> mTGMM_indices;
    double max_pT_mTGMM;
    if(mTGMM_on) {
        DoMTGMMJetFinding(tower_subtracted_particles_for_jets, seeds,
                          f_learn_weights, false, do_recombination, f_size,
                          mTGMM_indices, max_pT_mTGMM, tool, mTGMM_jets,
                          mTGMM_particle_weights, mTGMM_jets_params,
                          mTGMM_weights, fTmTGMM_iter, _event_jet_type,
                          _post_processing_method, event_jet_weight);
    }
    int lead_mTGMM_index = mTGMM_indices.size() ? mTGMM_indices.at(0) : -1;
    if (make_trace_diagrams && mTGMM_on) {
        tool->MakeTraceDiagram(tower_subtracted_particles_for_jets,
                               lead_mTGMM_index, "mTGMM", event_iter);
    }

    // Having found jets, now do a bit of data logging and analysis
    // turn off any degenerate algs for this particular event
    if (!mGMM_weights.size()) {
        mGMM_on = false;
    }
    if (!mGMMc_weights.size()) {
        mGMMc_on = false;
    }
    if (!mGMMs_weights.size()) {
        mGMMs_on = false;
    }
    if (!mUMM_weights.size()) {
        mUMM_on = false;
    }
    if (!mTGMM_weights.size()) {
        mTGMM_on = false;
    }
    if (!mTGMMs_weights.size()) {
        mTGMMs_on = false;
    }

    bool do_weight_distributions = false;
    if (do_weight_distributions && event_iter < 10 && !batched) {

        if(mGMM_on) {
            WeightDistribution(mGMM_particle_weights, lead_mGMM_index,
                               "mGMM", event_iter);
        }
        if(mGMMc_on) {
            WeightDistribution(mGMMc_particle_weights, lead_mGMMc_index,
                               "mGMMc", event_iter);
        }
        if(mTGMM_on) {
            WeightDistribution(mTGMM_particle_weights, lead_mTGMM_index,
                               "mTGMM", event_iter);
        }
        if(mUMM_on) {
            WeightDistribution(mUMM_particle_weights, lead_mUMM_index,
                               "mUMM", event_iter);
        }
        if(mGMMs_on) {
            WeightDistribution(mGMMs_particle_weights, lead_mGMMs_index,
                               "mGMMs", event_iter);
        }
        if(mTGMMs_on) {
            WeightDistribution(mTGMMs_particle_weights, lead_mTGMMs_index,
                               "mTGMMs", event_iter);
        }
    }

    // Some sporadic values for specific algorithms
    if (mGMMc_on) {
        fTmGMMc_r = sqrt(mGMMc_jets_params[lead_mGMMc_index].xx);
        fTmGMMc_phi = mGMMc_jets[lead_mGMMc_index].phi();
        fTmGMMc_eta = mGMMc_jets[lead_mGMMc_index].eta();

        double sum = 0;
        double weighted_sum = 0;
        for (unsigned int cluster_iter = 0; cluster_iter < mGMMc_jets.size(); cluster_iter++) {
            sum += sqrt(mGMMc_jets_params[cluster_iter].xx);
            weighted_sum += sqrt(mGMMc_jets_params[cluster_iter].xx) * mGMMc_weights[cluster_iter];
        }
        fTmGMMc_r_avg = sum / mGMMc_jets.size();
        fTmGMMc_r_weighted_avg = weighted_sum; // already normalized because weights sum to one
        if (mGMMc_indices.size() >= 2) {
            fTmGMMc_r_second = sqrt(mGMMc_jets_params[mGMMc_indices[1]].xx);
        }
        if (mGMMc_indices.size() >= 3) {
            fTmGMMc_r_third = sqrt(mGMMc_jets_params[mGMMc_indices[2]].xx);
        }
    }

    // Event displays
    bool do_event_displays = true;
    if (do_event_displays && !batched) {
        if(event_iter < 10 && mGMM_on) {
            tool->EventDisplay(tower_subtracted_particles_for_jets,
                               my_jets_large_r_ca,tops,
                               mGMM_jets,
                               mGMM_particle_weights,
                               lead_mGMM_index,
                               mGMM_jets_params,
                               "",
                               event_iter);
        }
        if(mGMM_on && static_cast<unsigned int>(event_iter) < _n_events) {
            tool->NewEventDisplay(tower_subtracted_particles_for_jets,
                                  my_jets_large_r_ca,tops,
                                  mGMM_jets,
                                  mGMM_particle_weights,
                                  lead_mGMM_index,
                                  mGMM_jets_params,
                                  mGMM_weights,
                                  "mGMM",
                                  event_iter);
            tool->JetContributionDisplay(tower_subtracted_particles_for_jets,
                                         mGMM_particle_weights,
                                         lead_mGMM_index,
                                         1,
                                         "m_mGMM",
                                         event_iter);
            tool->JetContributionDisplay(tower_subtracted_particles_for_jets,
                                         mGMM_particle_weights,
                                         lead_mGMM_index,
                                         0,
                                         "pt_mGMM",
                                         event_iter);
            if (_event_jet_type != FuzzyTools::NONE) {
                tool->EventJetDisplay(tower_subtracted_particles_for_jets,
                                      mGMM_jets,
                                      mGMM_particle_weights,
                                      mGMM_jets_params,
                                      mGMM_weights,
                                      "mGMM", event_iter);
            }
        }
        if(mTGMM_on && static_cast<unsigned int>(event_iter) < _n_events) {
            tool->NewEventDisplay(tower_subtracted_particles_for_jets,
                                  my_jets_large_r_ca,tops,
                                  mTGMM_jets,
                                  mTGMM_particle_weights,
                                  lead_mTGMM_index,
                                  mTGMM_jets_params,
                                  mTGMM_weights,
                                  "mTGMM",
                                  event_iter);
            if (_event_jet_type != FuzzyTools::NONE) {
                tool->EventJetDisplay(tower_subtracted_particles_for_jets,
                                      mTGMM_jets,
                                      mTGMM_particle_weights,
                                      mTGMM_jets_params,
                                      mTGMM_weights,
                                      "mTGMM", event_iter);
            }
        }
        if(mUMM_on && static_cast<unsigned int>(event_iter) < _n_events) {
            tool->NewEventDisplayUniform(tower_subtracted_particles_for_jets,
                                         my_jets_large_r_ca,tops,
                                         mUMM_jets,
                                         mUMM_particle_weights,
                                         lead_mUMM_index,
                                         mUMM_weights,
                                         "mUMM",
                                         event_iter);
            if (_event_jet_type != FuzzyTools::NONE) {
                //tool->EventJetDisplay(tower_subtracted_particles_for_jets,
                //                      mUMM_jets,
                //                      mUMM_particle_weights,
                //                      mUMM_jets_params,
                //                      mUMM_weights,
                //                      "mUMM", event_iter);
            }
        }
        if(mGMMs_on && static_cast<unsigned int>(event_iter) < _n_events) {
            tool->NewEventDisplay(tower_subtracted_particles_for_jets,
                                  my_jets_large_r_ca, tops,
                                  mGMMs_jets,
                                  mGMMs_particle_weights,
                                  lead_mGMMs_index,
                                  mGMMs_jets_params,
                                  mGMMs_weights,
                                  "mGMMs",
                                  event_iter);
            if (_event_jet_type != FuzzyTools::NONE) {
                tool->EventJetDisplay(tower_subtracted_particles_for_jets,
                                      mGMMs_jets,
                                      mGMMs_particle_weights,
                                      mGMMs_jets_params,
                                      mGMMs_weights,
                                      "mGMMs", event_iter);
            }
        }
        if(mGMMc_on && static_cast<unsigned int>(event_iter) < _n_events) {
            tool->ComparisonED(tower_subtracted_particles_for_jets,
                               my_jets_large_r_antikt, tops,
                               mGMMc_jets,
                               mGMMc_particle_weights,
                               lead_mGMMc_index,
                               mGMMc_jets_params,
                               mGMMc_weights,
                               "mGMMc",
                               event_iter);
            tool->NewEventDisplay(tower_subtracted_particles_for_jets,
                                  my_jets_large_r_ca, tops,
                                  mGMMc_jets,
                                  mGMMc_particle_weights,
                                  lead_mGMMc_index,
                                  mGMMc_jets_params,
                                  mGMMc_weights,
                                  "mGMMc",
                                  event_iter);
            if (_event_jet_type != FuzzyTools::NONE) {
                tool->EventJetDisplay(tower_subtracted_particles_for_jets,
                                      mGMMc_jets,
                                      mGMMc_particle_weights,
                                      mGMMc_jets_params,
                                      mGMMc_weights,
                                      "mGMMc", event_iter);
            }
            tool->NewEventDisplayPoster(tower_subtracted_particles_for_jets,
                                        my_jets_large_r_ca, tops,
                                        mGMMc_jets,
                                        mGMMc_particle_weights,
                                        lead_mGMMc_index,
                                        mGMMc_jets_params,
                                        mGMMc_weights,
                                        "mGMMc",
                                        event_iter);
        }
        if(mTGMMs_on && static_cast<unsigned int>(event_iter) < _n_events) {
            tool->NewEventDisplay(tower_subtracted_particles_for_jets,
                                  my_jets_large_r_ca,tops,
                                  mTGMMs_jets,
                                  mTGMMs_particle_weights,
                                  lead_mTGMMs_index,
                                  mTGMMs_jets_params,
                                  mTGMMs_weights,
                                  "mTGMMs",
                                  event_iter);
            if (_event_jet_type != FuzzyTools::NONE) {
                tool->EventJetDisplay(tower_subtracted_particles_for_jets,
                                      mTGMMs_jets,
                                      mTGMMs_particle_weights,
                                      mTGMMs_jets_params,
                                      mTGMMs_weights,
                                      "mTGMMs", event_iter);
            }
        }
    }

    if (mGMM_indices.size() >= 2) {
        fTmGMM_dr = mGMM_jets[mGMM_indices[0]].delta_R(mGMM_jets[mGMM_indices[1]]);
    }
    if (mGMMs_indices.size() >= 2) {
        fTmGMMs_dr = mGMMs_jets[mGMMs_indices[0]].delta_R(mGMMs_jets[mGMMs_indices[1]]);
    }
    if (mGMMc_indices.size() >= 2) {
        fTmGMMc_dr = mGMMc_jets[mGMMc_indices[0]].delta_R(mGMMc_jets[mGMMc_indices[1]]);
    }
    if (mTGMM_indices.size() >= 2) {
        fTmTGMM_dr = mTGMM_jets[mTGMM_indices[0]].delta_R(mTGMM_jets[mTGMM_indices[1]]);
    }
    if (mTGMMs_indices.size() >= 2) {
        fTmTGMMs_dr = mTGMMs_jets[mTGMMs_indices[0]].delta_R(mTGMMs_jets[mTGMMs_indices[1]]);
    }
    if (mUMM_indices.size() >= 2) {
        fTmUMM_dr = mUMM_jets[mUMM_indices[0]].delta_R(mUMM_jets[mUMM_indices[1]]);
    }
    for(unsigned int cluster_iter = 0; cluster_iter < mGMM_jets.size(); cluster_iter++) {
        fTmGMM_etas.push_back(mGMM_jets.at(cluster_iter).eta());
        fTmGMM_phis.push_back(mGMM_jets.at(cluster_iter).eta());
    }
    for(unsigned int cluster_iter = 0; cluster_iter < mGMMs_jets.size(); cluster_iter++) {
        fTmGMMs_etas.push_back(mGMMs_jets.at(cluster_iter).eta());
        fTmGMMs_phis.push_back(mGMMs_jets.at(cluster_iter).eta());
    }
    for(unsigned int cluster_iter = 0; cluster_iter < mGMMc_jets.size(); cluster_iter++) {
        fTmGMMc_etas.push_back(mGMMc_jets.at(cluster_iter).eta());
        fTmGMMc_phis.push_back(mGMMc_jets.at(cluster_iter).eta());
    }
    for(unsigned int cluster_iter = 0; cluster_iter < mUMM_jets.size(); cluster_iter++) {
        fTmUMM_etas.push_back(mUMM_jets.at(cluster_iter).eta());
        fTmUMM_phis.push_back(mUMM_jets.at(cluster_iter).eta());
    }
    for(unsigned int cluster_iter = 0; cluster_iter < mTGMM_jets.size(); cluster_iter++) {
        fTmTGMM_etas.push_back(mTGMM_jets.at(cluster_iter).eta());
        fTmTGMM_phis.push_back(mTGMM_jets.at(cluster_iter).eta());
    }
    for(unsigned int cluster_iter = 0; cluster_iter < mTGMMs_jets.size(); cluster_iter++) {
        fTmTGMMs_etas.push_back(mTGMMs_jets.at(cluster_iter).eta());
        fTmTGMMs_phis.push_back(mTGMMs_jets.at(cluster_iter).eta());
    }
    for(unsigned int particle_iter = 0; particle_iter < tower_subtracted_particles_for_jets.size(); particle_iter++) {
        fastjet::PseudoJet const& p = tower_subtracted_particles_for_jets[particle_iter];
        if (mGMM_on) {
            mGMM_weight_vec.push_back(mGMM_particle_weights[particle_iter][lead_mGMM_index]);
            mGMM_distance_vec.push_back(p.delta_R(mGMM_jets[lead_mGMM_index]));
        } else {
            mGMM_weight_vec.push_back(-1);
            mGMM_distance_vec.push_back(-1);
        }
        if (mGMMs_on) {
            mGMMs_weight_vec.push_back(mGMMs_particle_weights[particle_iter][lead_mGMMs_index]);
            mGMMs_distance_vec.push_back(p.delta_R(mGMMs_jets[lead_mGMMs_index]));
        } else {
            mGMMs_weight_vec.push_back(-1);
            mGMMs_distance_vec.push_back(-1);
        }
        if (mGMMc_on) {
            mGMMc_weight_vec.push_back(mGMMc_particle_weights[particle_iter][lead_mGMMc_index]);
            mGMMc_distance_vec.push_back(p.delta_R(mGMMc_jets[lead_mGMMc_index]));
        } else {
            mGMMc_weight_vec.push_back(-1);
            mGMMc_distance_vec.push_back(-1);
        }
        if (mTGMM_on) {
            mTGMM_weight_vec.push_back(mTGMM_particle_weights[particle_iter][lead_mTGMM_index]);
            mTGMM_distance_vec.push_back(p.delta_R(mTGMM_jets[lead_mTGMM_index]));
        } else {
            mTGMM_weight_vec.push_back(-1);
            mTGMM_distance_vec.push_back(-1);
        }
        if (mTGMMs_on) {
            mTGMMs_weight_vec.push_back(mTGMMs_particle_weights[particle_iter][lead_mTGMMs_index]);
            mTGMMs_distance_vec.push_back(p.delta_R(mTGMMs_jets[lead_mTGMMs_index]));
        } else {
            mTGMMs_weight_vec.push_back(-1);
            mTGMMs_distance_vec.push_back(-1);
        }
        if (mUMM_on) {
            mUMM_weight_vec.push_back(mUMM_particle_weights[particle_iter][lead_mUMM_index]);
            mUMM_distance_vec.push_back(p.delta_R(mUMM_jets[lead_mUMM_index]));
        } else {
            mUMM_weight_vec.push_back(-1);
            mUMM_distance_vec.push_back(-1);
        }
        // some older weight histo stuff, can be rebuilt from above
        if (p.user_info<MyUserInfo>().isPU()) {
            if(mGMM_on) {
                map_weight_vecs["mGMM_pu"].push_back(mGMM_particle_weights[particle_iter][lead_mGMM_index]);
            }
            if(mGMMs_on) {
                map_weight_vecs["mGMMs_pu"].push_back(mGMMs_particle_weights[particle_iter][lead_mGMMs_index]);
            }
            if(mGMMc_on) {
                map_weight_vecs["mGMMc_pu"].push_back(mGMMc_particle_weights[particle_iter][lead_mGMMc_index]);
            }
            if(mTGMM_on) {
                map_weight_vecs["mTGMM_pu"].push_back(mTGMM_particle_weights[particle_iter][lead_mTGMM_index]);
            }
            if(mTGMMs_on) {
                map_weight_vecs["mTGMMs_pu"].push_back(mTGMMs_particle_weights[particle_iter][lead_mTGMMs_index]);
            }
            if(mUMM_on) {
                map_weight_vecs["mUMM_pu"].push_back(mUMM_particle_weights[particle_iter][lead_mUMM_index]);
            }
        } else {
            if(mGMM_on) {
                map_weight_vecs["mGMM_hs"].push_back(mGMM_particle_weights[particle_iter][lead_mGMM_index]);
            }
            if(mGMMs_on) {
                map_weight_vecs["mGMMs_hs"].push_back(mGMMs_particle_weights[particle_iter][lead_mGMMs_index]);
            }
            if(mGMMc_on) {
                map_weight_vecs["mGMMc_hs"].push_back(mGMMc_particle_weights[particle_iter][lead_mGMMc_index]);
            }
            if(mTGMM_on) {
                map_weight_vecs["mTGMM_hs"].push_back(mTGMM_particle_weights[particle_iter][lead_mTGMM_index]);
            }
            if(mTGMMs_on) {
                map_weight_vecs["mTGMMs_hs"].push_back(mTGMMs_particle_weights[particle_iter][lead_mTGMMs_index]);
            }
            if(mUMM_on) {
                map_weight_vecs["mUMM_hs"].push_back(mUMM_particle_weights[particle_iter][lead_mUMM_index]);
            }
        }
    }
    // weight information takes up a ton of room, until we need it, turn it off
    mUMM_weight_vec.clear();
    mGMM_weight_vec.clear();
    mGMMc_weight_vec.clear();
    mGMMs_weight_vec.clear();
    mTGMM_weight_vec.clear();
    mTGMMs_weight_vec.clear();

    // same goes for distance
    mUMM_distance_vec.clear();
    mGMM_distance_vec.clear();
    mGMMc_distance_vec.clear();
    mGMMs_distance_vec.clear();
    mTGMM_distance_vec.clear();
    mTGMMs_distance_vec.clear();

    if(mGMM_on) {
        fTmGMM_m = tool->MLpT(tower_subtracted_particles_for_jets,mGMM_particle_weights,
                              lead_mGMM_index,mGMM_particle_weights[0].size(),1, true);
        fTmGMM_pt = tool->MLpT(tower_subtracted_particles_for_jets,mGMM_particle_weights,
                               lead_mGMM_index,mGMM_particle_weights[0].size(),0, true);
        fTmGMM_ml = tool->MLlpTGaussian(tower_subtracted_particles_for_jets,mGMM_jets[lead_mGMM_index],
                                        mGMM_jets_params[lead_mGMM_index], mGMM_weights[lead_mGMM_index], 1);
        fTmGMM_ptl = tool->MLlpTGaussian(tower_subtracted_particles_for_jets,mGMM_jets[lead_mGMM_index],
                                         mGMM_jets_params[lead_mGMM_index], mGMM_weights[lead_mGMM_index], 0);
        fTmGMM_pufrac_soft = JetPuFracSoft(tower_subtracted_particles_for_jets, mGMM_particle_weights, lead_mGMM_index);
        fTmGMM_pufrac_hard = JetPuFracHard(tower_subtracted_particles_for_jets, mGMM_particle_weights, lead_mGMM_index);
        fTmGMM_m_pu_soft = JetPuMassSoft(tower_subtracted_particles_for_jets, mGMM_particle_weights, lead_mGMM_index);
        fTmGMM_m_pu_hard = JetPuMassHard(tower_subtracted_particles_for_jets, mGMM_particle_weights, lead_mGMM_index);
    }
    if(mUMM_on) {
        fTmUMM_m = tool->MLpT(tower_subtracted_particles_for_jets, mUMM_particle_weights,
                              lead_mUMM_index, mUMM_particle_weights[0].size(), 1, true);
        fTmUMM_pt = max_pT_mUMM;
        fTmUMM_ml = tool->MLlpTUniform(tower_subtracted_particles_for_jets,mUMM_jets[lead_mUMM_index],
                                       mUMM_weights[lead_mUMM_index], 1);
        fTmUMM_ptl = tool->MLlpTUniform(tower_subtracted_particles_for_jets,mUMM_jets[lead_mUMM_index],
                                        mUMM_weights[lead_mUMM_index], 0);
        fTmUMM_pufrac_soft = JetPuFracSoft(tower_subtracted_particles_for_jets, mUMM_particle_weights, lead_mUMM_index);
        fTmUMM_pufrac_hard = JetPuFracHard(tower_subtracted_particles_for_jets, mUMM_particle_weights, lead_mUMM_index);
        fTmUMM_m_pu_soft = JetPuMassSoft(tower_subtracted_particles_for_jets, mUMM_particle_weights, lead_mUMM_index);
        fTmUMM_m_pu_hard = JetPuMassHard(tower_subtracted_particles_for_jets, mUMM_particle_weights, lead_mUMM_index);
    }
    if(mTGMM_on) {
        fTmTGMM_m = tool->MLpT(tower_subtracted_particles_for_jets, mTGMM_particle_weights,
                               lead_mTGMM_index, mTGMM_particle_weights[0].size(), 1, true);
        fTmTGMM_pt = max_pT_mTGMM;
        fTmTGMM_ml = tool->MLlpTTruncGaus(tower_subtracted_particles_for_jets,mTGMM_jets[lead_mTGMM_index],
                                          mTGMM_jets_params[lead_mTGMM_index], mTGMM_weights[lead_mTGMM_index],1);
        fTmTGMM_ptl = tool->MLlpTTruncGaus(tower_subtracted_particles_for_jets,mTGMM_jets[lead_mTGMM_index],
                                           mTGMM_jets_params[lead_mTGMM_index], mTGMM_weights[lead_mTGMM_index],0);
        fTmTGMM_pufrac_soft = JetPuFracSoft(tower_subtracted_particles_for_jets, mTGMM_particle_weights, lead_mTGMM_index);
        fTmTGMM_pufrac_hard = JetPuFracHard(tower_subtracted_particles_for_jets, mTGMM_particle_weights, lead_mTGMM_index);
        fTmTGMM_m_pu_soft = JetPuMassSoft(tower_subtracted_particles_for_jets, mTGMM_particle_weights, lead_mTGMM_index);
        fTmTGMM_m_pu_hard = JetPuMassHard(tower_subtracted_particles_for_jets, mTGMM_particle_weights, lead_mTGMM_index);
    }
    if(mTGMMs_on) {
        fTmTGMMs_m = tool->MLpT(tower_subtracted_particles_for_jets, mTGMMs_particle_weights,
                                lead_mTGMMs_index, mTGMMs_particle_weights[0].size(), 1, true);
        fTmTGMMs_pt = max_pT_mTGMMs;
        fTmTGMMs_ml = tool->MLlpTTruncGaus(tower_subtracted_particles_for_jets,mTGMMs_jets[lead_mTGMMs_index],
                                           mTGMMs_jets_params[lead_mTGMMs_index], mTGMMs_weights[lead_mTGMMs_index],1);
        fTmTGMMs_ptl = tool->MLlpTTruncGaus(tower_subtracted_particles_for_jets,mTGMMs_jets[lead_mTGMMs_index],
                                            mTGMMs_jets_params[lead_mTGMMs_index], mTGMMs_weights[lead_mTGMMs_index],0);
        fTmTGMMs_pufrac_soft = JetPuFracSoft(tower_subtracted_particles_for_jets, mTGMMs_particle_weights, lead_mTGMMs_index);
        fTmTGMMs_pufrac_hard = JetPuFracHard(tower_subtracted_particles_for_jets, mTGMMs_particle_weights, lead_mTGMMs_index);
        fTmTGMMs_m_pu_soft = JetPuMassSoft(tower_subtracted_particles_for_jets, mTGMMs_particle_weights, lead_mTGMMs_index);
        fTmTGMMs_m_pu_hard = JetPuMassHard(tower_subtracted_particles_for_jets, mTGMMs_particle_weights, lead_mTGMMs_index);
    }
    if(mGMMs_on) {
        fTmGMMs_m = tool->MLpT(tower_subtracted_particles_for_jets, mGMMs_particle_weights,
                               lead_mGMMs_index, mGMMs_particle_weights[0].size(), 1, true);
        fTmGMMs_pt = max_pT_mGMMs;
        fTmGMMs_ml = tool->MLlpTGaussian(tower_subtracted_particles_for_jets,mGMMs_jets[lead_mGMMs_index],
                                         mGMMs_jets_params[lead_mGMMs_index], mGMMs_weights[lead_mGMMs_index],1);
        fTmGMMs_ptl = tool->MLlpTGaussian(tower_subtracted_particles_for_jets,mGMMs_jets[lead_mGMMs_index],
                                          mGMMs_jets_params[lead_mGMMs_index], mGMMs_weights[lead_mGMMs_index],0);
        fTmGMMs_pufrac_soft = JetPuFracSoft(tower_subtracted_particles_for_jets, mGMMs_particle_weights, lead_mGMMs_index);
        fTmGMMs_pufrac_hard = JetPuFracHard(tower_subtracted_particles_for_jets, mGMMs_particle_weights, lead_mGMMs_index);
        fTmGMMs_m_pu_soft = JetPuMassSoft(tower_subtracted_particles_for_jets, mGMMs_particle_weights, lead_mGMMs_index);
        fTmGMMs_m_pu_hard = JetPuMassHard(tower_subtracted_particles_for_jets, mGMMs_particle_weights, lead_mGMMs_index);
    }
    if(mGMMc_on) {
        // compute the mass of all jets inside 1 sigma
        // collect all jets inside 1 sigma
        std::vector<unsigned int> jets_in_sigma;
        for (unsigned int j_i = 0; j_i < mGMMc_jets.size(); j_i++) {
            if (mGMMc_jets.at(j_i).delta_R(mGMMc_jets[lead_mGMMc_index]) < 0.8*fTmGMMc_r)
                jets_in_sigma.push_back(j_i);
        }
        fastjet::PseudoJet p;
        for (unsigned int p_i = 0; p_i < tower_subtracted_particles_for_jets.size(); p_i++) {
            int idx = tool->belongs_idx(mGMMc_particle_weights,
                                        p_i,
                                        true);
            if (std::find(jets_in_sigma.begin(), jets_in_sigma.end(), idx) != jets_in_sigma.end()) {
                p += tower_subtracted_particles_for_jets.at(p_i);
            }
        }
        fTmGMMc_m_sigma = p.m();

        fTmGMMc_m = tool->MLpT(tower_subtracted_particles_for_jets, mGMMc_particle_weights,
                               lead_mGMMc_index, mGMMc_particle_weights[0].size(), 1, true);
        fTmGMMc_pt = max_pT_mGMMc;
        fTmGMMc_ml = tool->MLlpTGaussian(tower_subtracted_particles_for_jets,mGMMc_jets[lead_mGMMc_index],
                                         mGMMc_jets_params[lead_mGMMc_index], mGMMc_weights[lead_mGMMc_index],1);
        fTmGMMc_ptl = tool->MLlpTGaussian(tower_subtracted_particles_for_jets,mGMMc_jets[lead_mGMMc_index],
                                          mGMMc_jets_params[lead_mGMMc_index], mGMMc_weights[lead_mGMMc_index],0);
        fTmGMMc_pufrac_soft = JetPuFracSoft(tower_subtracted_particles_for_jets, mGMMc_particle_weights, lead_mGMMc_index);
        fTmGMMc_pufrac_hard = JetPuFracHard(tower_subtracted_particles_for_jets, mGMMc_particle_weights, lead_mGMMc_index);
        fTmGMMc_m_pu_soft = JetPuMassSoft(tower_subtracted_particles_for_jets, mGMMc_particle_weights, lead_mGMMc_index);
        fTmGMMc_m_pu_hard = JetPuMassHard(tower_subtracted_particles_for_jets, mGMMc_particle_weights, lead_mGMMc_index);
    }

    int my_top = 0;
    if (tops[1].pt()> tops[0].pt()){
        my_top = 1;
    }

    fTtoppt = tops[0].pt();

    if(mGMM_on) {
        fTdeltatop_mGMM = tops[0].delta_R(mGMM_jets[lead_mGMM_index]);

        if (tops[1].delta_R(mGMM_jets[lead_mGMM_index]) < fTdeltatop_mGMM){
            fTdeltatop_mGMM = tops[1].delta_R(mGMM_jets[lead_mGMM_index]);
            fTtoppt = tops[1].pt();
        }
        fTmGMM_m_soft = tool->SoftpT(tower_subtracted_particles_for_jets,
                                     mGMM_particle_weights,
                                     lead_mGMM_index,
                                     1);
        fTmGMM_pt_soft = tool->SoftpT(tower_subtracted_particles_for_jets,
                                      mGMM_particle_weights,
                                      lead_mGMM_index,
                                      0);
        fTmGMM_ucpu = UnclusteredPileupComposition(tower_subtracted_particles_for_jets,
                                                   mGMM_particle_weights);
        fTmGMM_clpu = ClusteredPileupFrac(tower_subtracted_particles_for_jets,
                                          mGMM_particle_weights);
    }
    if(mUMM_on) {
        fTdeltatop_mUMM = tops[0].delta_R(mUMM_jets[lead_mUMM_index]);
        if (tops[1].delta_R(mUMM_jets[lead_mUMM_index]) <  fTdeltatop_mUMM) {
            fTdeltatop_mUMM = tops[1].delta_R(mUMM_jets[lead_mUMM_index]);
        }
        fTmUMM_m_soft = tool->SoftpT(tower_subtracted_particles_for_jets,
                                     mUMM_particle_weights,
                                     lead_mUMM_index,
                                     1);
        fTmUMM_pt_soft = tool->SoftpT(tower_subtracted_particles_for_jets,
                                      mUMM_particle_weights,
                                      lead_mUMM_index,
                                      0);
        fTmUMM_ucpu = UnclusteredPileupComposition(tower_subtracted_particles_for_jets,
                                                   mUMM_particle_weights);
        fTmUMM_clpu = ClusteredPileupFrac(tower_subtracted_particles_for_jets,
                                          mUMM_particle_weights);
    }
    if(mTGMM_on) {
        fTdeltatop_mTGMM = tops[0].delta_R(mTGMM_jets[lead_mTGMM_index]);
        if (tops[1].delta_R(mTGMM_jets[lead_mTGMM_index]) < fTdeltatop_mTGMM) {
            fTdeltatop_mTGMM = tops[1].delta_R(mTGMM_jets[lead_mTGMM_index]);
        }
        fTmTGMM_m_soft = tool->SoftpT(tower_subtracted_particles_for_jets,
                                      mTGMM_particle_weights,
                                      lead_mTGMM_index,
                                      1);
        fTmTGMM_pt_soft = tool->SoftpT(tower_subtracted_particles_for_jets,
                                       mTGMM_particle_weights,
                                       lead_mTGMM_index,
                                       0);
        fTmTGMM_ucpu = UnclusteredPileupComposition(tower_subtracted_particles_for_jets,
                                                    mTGMM_particle_weights);
        fTmTGMM_clpu = ClusteredPileupFrac(tower_subtracted_particles_for_jets,
                                           mTGMM_particle_weights);
    }
    if(mTGMMs_on) {
        fTdeltatop_mTGMMs = tops[0].delta_R(mTGMMs_jets[lead_mTGMMs_index]);
        if (tops[1].delta_R(mTGMMs_jets[lead_mTGMMs_index]) < fTdeltatop_mTGMMs) {
            fTdeltatop_mTGMMs = tops[1].delta_R(mTGMMs_jets[lead_mTGMMs_index]);
        }
        fTmTGMMs_m_soft = tool->SoftpT(tower_subtracted_particles_for_jets,
                                       mTGMMs_particle_weights,
                                       lead_mTGMMs_index,
                                       1);
        fTmTGMMs_pt_soft = tool->SoftpT(tower_subtracted_particles_for_jets,
                                        mTGMMs_particle_weights,
                                        lead_mTGMMs_index,
                                        0);
        fTmTGMMs_ucpu = UnclusteredPileupComposition(tower_subtracted_particles_for_jets,
                                                     mTGMMs_particle_weights);
        fTmTGMMs_clpu = ClusteredPileupFrac(tower_subtracted_particles_for_jets,
                                            mTGMMs_particle_weights);
    }
    if(mGMMs_on) {
        fTdeltatop_mGMMs = tops[0].delta_R(mGMMs_jets[lead_mGMMs_index]);
        if (tops[1].delta_R(mGMMs_jets[lead_mGMMs_index]) < fTdeltatop_mGMMs) {
            fTdeltatop_mGMMs = tops[1].delta_R(mGMMs_jets[lead_mGMMs_index]);
        }
        fTmGMMs_m_soft = tool->SoftpT(tower_subtracted_particles_for_jets,
                                      mGMMs_particle_weights,
                                      lead_mGMMs_index,
                                      1);
        fTmGMMs_pt_soft = tool->SoftpT(tower_subtracted_particles_for_jets,
                                       mGMMs_particle_weights,
                                       lead_mGMMs_index,
                                       0);
        fTmGMMs_ucpu = UnclusteredPileupComposition(tower_subtracted_particles_for_jets,
                                                    mGMMs_particle_weights);
        fTmGMMs_clpu = ClusteredPileupFrac(tower_subtracted_particles_for_jets,
                                           mGMMs_particle_weights);
    }
    if(mGMMc_on) {
        fTdeltatop_mGMMc = tops[0].delta_R(mGMMc_jets[lead_mGMMc_index]);
        if (tops[1].delta_R(mGMMc_jets[lead_mGMMc_index]) < fTdeltatop_mGMMc) {
            fTdeltatop_mGMMc = tops[1].delta_R(mGMMc_jets[lead_mGMMc_index]);
        }
        fTmGMMc_m_soft = tool->SoftpT(tower_subtracted_particles_for_jets,
                                      mGMMc_particle_weights,
                                      lead_mGMMc_index,
                                      1);
        fTmGMMc_pt_soft = tool->SoftpT(tower_subtracted_particles_for_jets,
                                       mGMMc_particle_weights,
                                       lead_mGMMc_index,
                                       0);
        fTmGMMc_ucpu = UnclusteredPileupComposition(tower_subtracted_particles_for_jets,
                                                    mGMMc_particle_weights);
        fTmGMMc_clpu = ClusteredPileupFrac(tower_subtracted_particles_for_jets,
                                           mGMMc_particle_weights);
    }

    // Moments
    vector<double> moments_m;
    vector<double> moments_pt;
    double sig;
    if(mGMM_on) {
        moments_m = tool->CentralMoments(tower_subtracted_particles_for_jets, mGMM_particle_weights,
                                         lead_mGMM_index, 3, &totalMass);

        moments_pt = tool->CentralMoments(tower_subtracted_particles_for_jets, mGMM_particle_weights,
                                          lead_mGMM_index, 3, &totalpT);

        fTmGMM_m_mean = moments_m[0];
        fTmGMM_m_var  = moments_m[1];
        sig = sqrt(fTmGMM_m_var);
        fTmGMM_m_skew = moments_m[2];
        fTmGMM_pt_mean = moments_pt[0];
        fTmGMM_pt_var = moments_pt[1];
        sig = sqrt(fTmGMM_m_var);
        fTmGMM_pt_skew = moments_pt[2];
    }

    moments_m.clear();
    moments_pt.clear();

    if(mUMM_on) {
        moments_m = tool->CentralMoments(tower_subtracted_particles_for_jets,
                                         mUMM_particle_weights,
                                         lead_mUMM_index, 3, &totalMass);

        moments_pt = tool->CentralMoments(tower_subtracted_particles_for_jets, mUMM_particle_weights,
                                          lead_mUMM_index, 3, &totalpT);

        fTmUMM_m_mean = moments_m[0];
        fTmUMM_m_var  = moments_m[1];
        sig = sqrt(fTmUMM_m_var);
        fTmUMM_m_skew = moments_m[2];
        fTmUMM_pt_mean = moments_pt[0];
        fTmUMM_pt_var = moments_pt[1];
        sig = sqrt(fTmUMM_m_var);
        fTmUMM_pt_skew = moments_pt[2];
    }
    moments_m.clear();
    moments_pt.clear();

    if(mTGMM_on) {
        moments_m = tool->CentralMoments(tower_subtracted_particles_for_jets, mTGMM_particle_weights,
                                         lead_mTGMM_index, 3, &totalMass);

        moments_pt = tool->CentralMoments(tower_subtracted_particles_for_jets, mTGMM_particle_weights,
                                          lead_mTGMM_index, 3, &totalpT);

        fTmTGMM_m_mean = moments_m[0];
        fTmTGMM_m_var  = moments_m[1];
        sig = sqrt(fTmTGMM_m_var);
        fTmTGMM_m_skew = moments_m[2];
        fTmTGMM_pt_mean = moments_pt[0];
        fTmTGMM_pt_var = moments_pt[1];
        sig = sqrt(fTmTGMM_m_var);
        fTmTGMM_pt_skew = moments_pt[2];
    }

    moments_m.clear();
    moments_pt.clear();

    if(mTGMMs_on) {
        moments_m = tool->CentralMoments(tower_subtracted_particles_for_jets, mTGMMs_particle_weights,
                                         lead_mTGMMs_index, 3, &totalMass);

        moments_pt = tool->CentralMoments(tower_subtracted_particles_for_jets, mTGMMs_particle_weights,
                                          lead_mTGMMs_index, 3, &totalpT);

        fTmTGMMs_m_mean = moments_m[0];
        fTmTGMMs_m_var  = moments_m[1];
        sig = sqrt(fTmTGMMs_m_var);
        fTmTGMMs_m_skew = moments_m[2];
        fTmTGMMs_pt_mean = moments_pt[0];
        fTmTGMMs_pt_var = moments_pt[1];
        sig = sqrt(fTmTGMMs_m_var);
        fTmTGMMs_pt_skew = moments_pt[2];
    }

    moments_m.clear();
    moments_pt.clear();

    if(mGMMs_on) {
        moments_m = tool->CentralMoments(tower_subtracted_particles_for_jets, mGMMs_particle_weights,
                                         lead_mGMMs_index, 3, &totalMass);

        moments_pt = tool->CentralMoments(tower_subtracted_particles_for_jets, mGMMs_particle_weights,
                                          lead_mGMMs_index, 3, &totalpT);

        fTmGMMs_m_mean = moments_m[0];
        fTmGMMs_m_var  = moments_m[1];
        sig = sqrt(fTmGMMs_m_var);
        fTmGMMs_m_skew = moments_m[2];
        fTmGMMs_pt_mean = moments_pt[0];
        fTmGMMs_pt_var = moments_pt[1];
        sig = sqrt(fTmGMMs_m_var);
        fTmGMMs_pt_skew = moments_pt[2];
    }

    moments_m.clear();
    moments_pt.clear();

    if(mGMMc_on) {
        moments_m = tool->CentralMoments(tower_subtracted_particles_for_jets, mGMMc_particle_weights,
                                         lead_mGMMc_index, 3, &totalMass);

        moments_pt = tool->CentralMoments(tower_subtracted_particles_for_jets, mGMMc_particle_weights,
                                          lead_mGMMc_index, 3, &totalpT);

        fTmGMMc_m_mean = moments_m[0];
        fTmGMMc_m_var  = moments_m[1];
        sig = sqrt(fTmGMMc_m_var);
        fTmGMMc_m_skew = moments_m[2];
        fTmGMMc_pt_mean = moments_pt[0];
        fTmGMMc_pt_var = moments_pt[1];
        sig = sqrt(fTmGMMc_m_var);
        fTmGMMc_pt_skew = moments_pt[2];
    }

    moments_m.clear();
    moments_pt.clear();

#ifdef WITHROOT
    if (!_non_standard_study) {
        t_t->Fill();
    }

    if (!batched && should_print) {
#else
        if (true) {
#endif
            map<string, float*>::const_iterator iter;
            for (iter = tree_vars.begin(); iter != tree_vars.end(); iter++) {
                cout << iter->first << ": " << *(iter->second) << endl;
            }
        }

        if(f_debug) cout << "FuzzyAnalysis::AnalyzeEvent End " << endl;
        return;
    }



    // declate branches
    void FuzzyAnalysis::DeclareBranches(){
        tree_vars["rho"] = &fTrho;
        tree_vars["event_jet_strength"] = &fTevent_jet_strength;
        tree_vars["event_jet_rho_offset"] = &fTevent_jet_rho_offset;

        tree_vars["CA_m"] = &fTCA_m;
        tree_vars["CA_pt"] = &fTCA_pt;
        tree_vars["CA_pufrac"] = &fTCA_pufrac;
        tree_vars["CA_m_pu"] = &fTCA_m_pu;
        tree_vars["CA_dr"] = &fTCA_dr;
        tree_vars["CA_area"] = &fTCA_area;

        tree_vars["toppt"] = &fTtoppt;

        tree_vars["antikt_m"] = &fTantikt_m;
        tree_vars["antikt_m_pu"] = &fTantikt_m_pu;
        tree_vars["antikt_pt"] = &fTantikt_pt;
        tree_vars["antikt_pufrac"] = &fTantikt_pufrac;
        tree_vars["antikt_dr"] = &fTantikt_dr;
        tree_vars["antikt_area"] = &fTantikt_area;

        tree_vars["antikt_m_trimmed_two"] = &fTantikt_m_trimmed_two;
        tree_vars["antikt_m_pu_trimmed_two"] = &fTantikt_m_pu_trimmed_two;
        tree_vars["antikt_pt_trimmed_two"] = &fTantikt_pt_trimmed_two;
        tree_vars["antikt_pufrac_trimmed_two"] = &fTantikt_pufrac_trimmed_two;
        tree_vars["antikt_area_trimmed_two"] = &fTantikt_area_trimmed_two;

        tree_vars["antikt_m_trimmed_three"] = &fTantikt_m_trimmed_three;
        tree_vars["antikt_m_pu_trimmed_three"] = &fTantikt_m_pu_trimmed_three;
        tree_vars["antikt_pt_trimmed_three"] = &fTantikt_pt_trimmed_three;
        tree_vars["antikt_pufrac_trimmed_three"] = &fTantikt_pufrac_trimmed_three;
        tree_vars["antikt_area_trimmed_three"] = &fTantikt_area_trimmed_three;

        tree_vars["mGMMs_m"] = &fTmGMMs_m;
        tree_vars["mGMMs_pt"] = &fTmGMMs_pt;
        tree_vars["deltatop_mGMMs"] = &fTdeltatop_mGMMs;
        tree_vars["mGMMs_m_mean"] = &fTmGMMs_m_mean;
        tree_vars["mGMMs_m_var"] = &fTmGMMs_m_var;
        tree_vars["mGMMs_m_skew"] = &fTmGMMs_m_skew;
        tree_vars["mGMMs_pt_mean"] = &fTmGMMs_pt_mean;
        tree_vars["mGMMs_pt_var"] = &fTmGMMs_pt_var;
        tree_vars["mGMMs_pt_skew"] = &fTmGMMs_pt_skew;
        tree_vars["mGMMs_ptl"] = &fTmGMMs_ptl;
        tree_vars["mGMMs_ml"] = &fTmGMMs_ml;
        tree_vars["mGMMs_m_soft"] = &fTmGMMs_m_soft;
        tree_vars["mGMMs_pt_soft"] = &fTmGMMs_pt_soft;
        tree_vars["mGMMs_ucpu"] = &fTmGMMs_ucpu;
        tree_vars["mGMMs_clpu"] = &fTmGMMs_clpu;
        tree_vars["mGMMs_pufrac_soft"] = &fTmGMMs_pufrac_soft;
        tree_vars["mGMMs_pufrac_hard"] = &fTmGMMs_pufrac_hard;
        tree_vars["mGMMs_m_pu_soft"] = &fTmGMMs_m_pu_soft;
        tree_vars["mGMMs_m_pu_hard"] = &fTmGMMs_m_pu_hard;
        tree_vars["mGMMs_dr"] = &fTmGMMs_dr;

        tree_vars["mGMMc_m"] = &fTmGMMc_m;
        tree_vars["mGMMc_pt"] = &fTmGMMc_pt;
        tree_vars["deltatop_mGMMc"] = &fTdeltatop_mGMMc;
        tree_vars["mGMMc_m_mean"] = &fTmGMMc_m_mean;
        tree_vars["mGMMc_m_var"] = &fTmGMMc_m_var;
        tree_vars["mGMMc_m_skew"] = &fTmGMMc_m_skew;
        tree_vars["mGMMc_pt_mean"] = &fTmGMMc_pt_mean;
        tree_vars["mGMMc_pt_var"] = &fTmGMMc_pt_var;
        tree_vars["mGMMc_pt_skew"] = &fTmGMMc_pt_skew;
        tree_vars["mGMMc_ptl"] = &fTmGMMc_ptl;
        tree_vars["mGMMc_ml"] = &fTmGMMc_ml;
        tree_vars["mGMMc_m_soft"] = &fTmGMMc_m_soft;
        tree_vars["mGMMc_pt_soft"] = &fTmGMMc_pt_soft;
        tree_vars["mGMMc_ucpu"] = &fTmGMMc_ucpu;
        tree_vars["mGMMc_clpu"] = &fTmGMMc_clpu;
        tree_vars["mGMMc_pufrac_soft"] = &fTmGMMc_pufrac_soft;
        tree_vars["mGMMc_pufrac_hard"] = &fTmGMMc_pufrac_hard;
        tree_vars["mGMMc_m_pu_soft"] = &fTmGMMc_m_pu_soft;
        tree_vars["mGMMc_m_pu_hard"] = &fTmGMMc_m_pu_hard;
        tree_vars["mGMMc_dr"] = &fTmGMMc_dr;

        tree_vars["mGMMc_r"] = &fTmGMMc_r;
        tree_vars["mGMMc_r_second"] = &fTmGMMc_r_second;
        tree_vars["mGMMc_r_third"] = &fTmGMMc_r_third;
        tree_vars["mGMMc_r_avg"] = &fTmGMMc_r_avg;
        tree_vars["mGMMc_r_weighted_avg"] = &fTmGMMc_r_weighted_avg;

        tree_vars["mTGMMs_m"] = &fTmTGMMs_m;
        tree_vars["mTGMMs_pt"] = &fTmTGMMs_pt;
        tree_vars["deltatop_mTGMMs"] = &fTdeltatop_mTGMMs;
        tree_vars["mTGMMs_m_mean"] = &fTmTGMMs_m_mean;
        tree_vars["mTGMMs_m_var"] = &fTmTGMMs_m_var;
        tree_vars["mTGMMs_m_skew"] = &fTmTGMMs_m_skew;
        tree_vars["mTGMMs_pt_mean"] = &fTmTGMMs_pt_mean;
        tree_vars["mTGMMs_pt_var"] = &fTmTGMMs_pt_var;
        tree_vars["mTGMMs_pt_skew"] = &fTmTGMMs_pt_skew;
        tree_vars["mTGMMs_ptl"] = &fTmTGMMs_ptl;
        tree_vars["mTGMMs_ml"] = &fTmTGMMs_ml;
        tree_vars["mTGMMs_m_soft"] = &fTmTGMMs_m_soft;
        tree_vars["mTGMMs_pt_soft"] = &fTmTGMMs_pt_soft;
        tree_vars["mTGMMs_ucpu"] = &fTmTGMMs_ucpu;
        tree_vars["mTGMMs_clpu"] = &fTmTGMMs_clpu;
        tree_vars["mTGMMs_pufrac_soft"] = &fTmTGMMs_pufrac_soft;
        tree_vars["mTGMMs_pufrac_hard"] = &fTmTGMMs_pufrac_hard;
        tree_vars["mTGMMs_m_pu_soft"] = &fTmTGMMs_m_pu_soft;
        tree_vars["mTGMMs_m_pu_hard"] = &fTmTGMMs_m_pu_hard;
        tree_vars["mTGMMs_dr"] = &fTmTGMMs_dr;

        tree_vars["mUMM_m"] = &fTmUMM_m;
        tree_vars["mUMM_pt"] = &fTmUMM_pt;
        tree_vars["deltatop_mUMM"] = &fTdeltatop_mUMM;
        tree_vars["mUMM_m_mean"] = &fTmUMM_m_mean;
        tree_vars["mUMM_m_var"] = &fTmUMM_m_var;
        tree_vars["mUMM_m_skew"] = &fTmUMM_m_skew;
        tree_vars["mUMM_pt_mean"] = &fTmUMM_pt_mean;
        tree_vars["mUMM_pt_var"] = &fTmUMM_pt_var;
        tree_vars["mUMM_pt_skew"] = &fTmUMM_pt_skew;
        tree_vars["mUMM_ptl"] = &fTmUMM_ptl;
        tree_vars["mUMM_ml"] = &fTmUMM_ml;
        tree_vars["mUMM_m_soft"] = &fTmUMM_m_soft;
        tree_vars["mUMM_pt_soft"] = &fTmUMM_pt_soft;
        tree_vars["mUMM_ucpu"] = &fTmUMM_ucpu;
        tree_vars["mUMM_clpu"] = &fTmUMM_clpu;
        tree_vars["mUMM_pufrac_soft"] = &fTmUMM_pufrac_soft;
        tree_vars["mUMM_pufrac_hard"] = &fTmUMM_pufrac_hard;
        tree_vars["mUMM_m_pu_soft"] = &fTmUMM_m_pu_soft;
        tree_vars["mUMM_m_pu_hard"] = &fTmUMM_m_pu_hard;
        tree_vars["mUMM_dr"] = &fTmUMM_dr;

        tree_vars["mGMM_m"] = &fTmGMM_m;
        tree_vars["mGMM_pt"] = &fTmGMM_pt;
        tree_vars["deltatop_mGMM"] = &fTdeltatop_mGMM;
        tree_vars["mGMM_m_mean"] = &fTmGMM_m_mean;
        tree_vars["mGMM_m_var"] = &fTmGMM_m_var;
        tree_vars["mGMM_m_skew"] = &fTmGMM_m_skew;
        tree_vars["mGMM_pt_mean"] = &fTmGMM_pt_mean;
        tree_vars["mGMM_pt_var"] = &fTmGMM_pt_var;
        tree_vars["mGMM_pt_skew"] = &fTmGMM_pt_skew;
        tree_vars["mGMM_ptl"] = &fTmGMM_ptl;
        tree_vars["mGMM_ml"] = &fTmGMM_ml;
        tree_vars["mGMM_m_soft"] = &fTmUMM_m_soft;
        tree_vars["mGMM_pt_soft"] = &fTmUMM_pt_soft;
        tree_vars["mGMM_ucpu"] = &fTmGMM_ucpu;
        tree_vars["mGMM_clpu"] = &fTmGMM_clpu;
        tree_vars["mGMM_pufrac_soft"] = &fTmGMM_pufrac_soft;
        tree_vars["mGMM_pufrac_hard"] = &fTmGMM_pufrac_hard;
        tree_vars["mGMM_m_pu_soft"] = &fTmGMM_m_pu_soft;
        tree_vars["mGMM_m_pu_hard"] = &fTmGMM_m_pu_hard;
        tree_vars["mGMM_dr"] = &fTmGMM_dr;

        tree_vars["mTGMM_m"] = &fTmTGMM_m;
        tree_vars["mTGMM_pt"] = &fTmTGMM_pt;
        tree_vars["deltatop_mTGMM"] = &fTdeltatop_mTGMM;
        tree_vars["mTGMM_m_mean"] = &fTmTGMM_m_mean;
        tree_vars["mTGMM_m_var"] = &fTmTGMM_m_var;
        tree_vars["mTGMM_m_skew"] = &fTmTGMM_m_skew;
        tree_vars["mTGMM_pt_mean"] = &fTmTGMM_pt_mean;
        tree_vars["mTGMM_pt_var"] = &fTmTGMM_pt_var;
        tree_vars["mTGMM_pt_skew"] = &fTmTGMM_pt_skew;
        tree_vars["mTGMM_ptl"] = &fTmTGMM_ptl;
        tree_vars["mTGMM_ml"] = &fTmTGMM_ml;
        tree_vars["mTGMM_m_soft"] = &fTmTGMM_m_soft;
        tree_vars["mTGMM_pt_soft"] = &fTmTGMM_pt_soft;
        tree_vars["mTGMM_ucpu"] = &fTmTGMM_ucpu;
        tree_vars["mTGMM_clpu"] = &fTmTGMM_clpu;
        tree_vars["mTGMM_pufrac_soft"] = &fTmTGMM_pufrac_soft;
        tree_vars["mTGMM_pufrac_hard"] = &fTmTGMM_pufrac_hard;
        tree_vars["mTGMM_m_pu_soft"] = &fTmTGMM_m_pu_soft;
        tree_vars["mTGMM_m_pu_hard"] = &fTmTGMM_m_pu_hard;
        tree_vars["mTGMM_dr"] = &fTmTGMM_dr;

        t_t->Branch("seed_location_noise", &fTseed_location_noise, "seed_location_noise/F");
        // Event Properties
#ifdef WITHROOT
        t_t->Branch("EventNumber", &fTEventNumber, "EventNumber/I");
        t_t->Branch("NPV",         &fTNPV,         "NPV/I");
        t_t->Branch("n_jet_seeds", &fTn_jet_seeds, "n_jet_seeds/I");

        // declare non-float branches
        t_t->Branch("mGMM_iter", &fTmGMM_iter, "mGMM_iter/i");
        t_t->Branch("mGMMc_iter", &fTmGMMc_iter, "mGMMc_iter/i");
        t_t->Branch("mGMMs_iter", &fTmGMMs_iter, "mGMMs_iter/i");
        t_t->Branch("mTGMM_iter", &fTmTGMM_iter, "mTGMM_iter/i");
        t_t->Branch("mTGMMs_iter", &fTmTGMMs_iter, "mTGMMs_iter/i");
        t_t->Branch("mUMM_iter", &fTmUMM_iter, "mUMM_iter/i");

        // declare vector branches
        t_t->Branch("CA_nsubjettiness", "vector<float>", &CA_nsubjettiness);
        t_t->Branch("antikt_nsubjettiness", "vector<float>", &antikt_nsubjettiness);

        t_t->Branch("mGMM_etas", "vector<float>", &fTmGMM_etas);
        t_t->Branch("mGMM_phis", "vector<float>", &fTmGMM_phis);
        t_t->Branch("mGMMc_etas", "vector<float>", &fTmGMMc_etas);
        t_t->Branch("mGMMc_phis", "vector<float>", &fTmGMMc_phis);
        t_t->Branch("mGMMs_etas", "vector<float>", &fTmGMMs_etas);
        t_t->Branch("mGMMs_phis", "vector<float>", &fTmGMMs_phis);
        t_t->Branch("mUMM_etas", "vector<float>", &fTmUMM_etas);
        t_t->Branch("mUMM_phis", "vector<float>", &fTmUMM_phis);
        t_t->Branch("mTGMM_etas", "vector<float>", &fTmTGMM_etas);
        t_t->Branch("mTGMM_phis", "vector<float>", &fTmTGMM_phis);
        t_t->Branch("mTGMMs_etas", "vector<float>", &fTmTGMMs_etas);
        t_t->Branch("mTGMMs_phis", "vector<float>", &fTmTGMMs_phis);

        t_t->Branch("mGMM_weight_vec", "vector<float>", &mGMM_weight_vec);
        t_t->Branch("mGMMs_weight_vec", "vector<float>", &mGMMs_weight_vec);
        t_t->Branch("mGMMc_weight_vec", "vector<float>", &mGMMc_weight_vec);
        t_t->Branch("mTGMM_weight_vec", "vector<float>", &mTGMM_weight_vec);
        t_t->Branch("mTGMMs_weight_vec", "vector<float>", &mTGMMs_weight_vec);
        t_t->Branch("mUMM_weight_vec", "vector<float>", &mUMM_weight_vec);

        t_t->Branch("mGMM_distance_vec", "vector<float>", &mGMM_distance_vec);
        t_t->Branch("mGMMs_distance_vec", "vector<float>", &mGMMs_distance_vec);
        t_t->Branch("mGMMc_distance_vec", "vector<float>", &mGMMc_distance_vec);
        t_t->Branch("mTGMM_distance_vec", "vector<float>", &mTGMM_distance_vec);
        t_t->Branch("mTGMMs_distance_vec", "vector<float>", &mTGMMs_distance_vec);
        t_t->Branch("mUMM_distance_vec", "vector<float>", &mUMM_distance_vec);

        t_t->Branch("is_pileup_vec", &is_pileup_vec);
        t_t->Branch("is_lead_antikt_constituent_vec", &is_lead_antikt_constituent_vec);
        t_t->Branch("is_lead_antikt_trimmed_two_constituent_vec", &is_lead_antikt_trimmed_two_constituent_vec);

        stringstream ss;
        map<string, float*>::const_iterator iter;
        for (iter = tree_vars.begin(); iter != tree_vars.end(); iter++) {
            string branch_name = iter->first;
            ss.clear();
            ss.str("");
            ss << branch_name << "/F";
            string branch_typename = ss.str();
            t_t->Branch(branch_name.c_str(), iter->second, branch_typename.c_str());
        }
#endif

        return;
    }


    // resets vars
    void FuzzyAnalysis::ResetBranches(){
        // reset branches
        fTEventNumber = -999;
        fTNPV         = -999;
        fTn_jet_seeds = -999;

        map<string, float*>::const_iterator iter;
        for (iter = tree_vars.begin(); iter != tree_vars.end(); iter++) {
            *(iter->second) = -999999;
        }

        // reset non-float basic branches
        fTmGMM_iter = 0;
        fTmGMMc_iter = 0;
        fTmGMMs_iter = 0;
        fTmTGMM_iter = 0;
        fTmTGMMs_iter = 0;
        fTmUMM_iter = 0;

        // reset vector branches

        CA_nsubjettiness.clear();
        antikt_nsubjettiness.clear();

        fTmGMM_etas.clear();
        fTmGMM_phis.clear();
        fTmGMMc_etas.clear();
        fTmGMMc_phis.clear();
        fTmGMMs_etas.clear();
        fTmGMMs_phis.clear();
        fTmUMM_etas.clear();
        fTmUMM_phis.clear();
        fTmTGMM_etas.clear();
        fTmTGMM_phis.clear();
        fTmTGMMs_etas.clear();
        fTmTGMMs_phis.clear();

        mGMM_weight_vec.clear();
        mGMMs_weight_vec.clear();
        mGMMc_weight_vec.clear();
        mTGMM_weight_vec.clear();
        mTGMMs_weight_vec.clear();
        mUMM_weight_vec.clear();

        mGMM_distance_vec.clear();
        mGMMs_distance_vec.clear();
        mGMMc_distance_vec.clear();
        mTGMM_distance_vec.clear();
        mTGMMs_distance_vec.clear();
        mUMM_distance_vec.clear();

        is_pileup_vec.clear();
        is_lead_antikt_constituent_vec.clear();
        is_lead_antikt_trimmed_two_constituent_vec.clear();
    }
