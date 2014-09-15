#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <set>
#include <assert.h>

#include "myFastJetBase.h"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/Selector.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/Njettiness.hh"

#include "Pythia8/Pythia.h"

#include "FuzzyAnalysis.h"
#include "FuzzyTools.h"
#include "ROOTConf.h"

#ifdef WITHROOT
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TMath.h"
#include "TVector3.h"
#include "TMatrix.h"
#include "TH1F.h"
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
                                       cluster_count, 0);
            pTs.push_back(SortablePt(i, holdpT));
        }
        
        sort(pTs.begin(), pTs.end());
        // reverse the list of indices and push it back onto the output
        for (unsigned int iter = jets.size(); iter-- > 0;) {
            indices.push_back((signed int) pTs.at(iter).index);
        }
        if (jets.size()) {
            pT = pTs.at(pTs.size() - 1).pT;
        }
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
                          vector<double>& jet_weights) {
        vecPseudoJet parts = particles;
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
                                         &jet_weights);
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
                           vector<double>& jet_weights) {
        tool->SetKernelType(FuzzyTools::TRUNCGAUSSIAN);
        tool->SetLearnWeights(learn_weights);
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
                                           &jet_weights);
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
                           vector<double>& jet_weights) {
        tool->SetKernelType(FuzzyTools::GAUSSIAN);
        tool->SetLearnWeights(learn_weights);
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
                                          &jet_weights);
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
                          vector<double>& jet_weights) {
        tool->SetKernelType(FuzzyTools::GAUSSIAN);
        tool->SetLearnWeights(learn_weights);
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
                                          &jet_weights);
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

    batched = false;
    should_print = true;
    
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
    #ifdef WITHROOT
    t_f = new TFile(fullName.c_str(), "RECREATE");
    t_t = new TTree("EventTree", "Event Tree for Fuzzy");
    #endif
    
    if(!batched) {
        SetupHistosMap();
    }
    DeclareBranches();

    ResetBranches();

    tool->SetPrefix(directory_prefix);

    return;
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
                             mGMMs_jets_params, mGMMs_weights);
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
                              mGMMc_jets_params, mGMMc_weights);
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
                              mTGMMs_jets_params, mTGMMs_weights);
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
                             mGMM_jets_params, mGMM_weights);
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
                             mUMM_particle_weights, mUMM_weights);
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
                              mTGMM_particle_weights, mTGMM_jets_params, mTGMM_weights);
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
    WriteHistosMap();
    t_t->Write();
    delete t_f;
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
    tops.push_back(dl);

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
            p.set_user_info(new MyUserInfo(pythia_MB->event[particle_idx].id(),particle_idx,pythia_MB->event[particle_idx].charge(),true));
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
        p.set_user_info(new MyUserInfo(pythia8->event[particle_idx].id(),particle_idx,pythia8->event[particle_idx].charge(),false));

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

// Analyze
void FuzzyAnalysis::AnalyzeEvent(int event_iter, Pythia8::Pythia* pythia8, Pythia8::Pythia* pythia_MB, int NPV){
    if(f_debug) cout << "FuzzyAnalysis::AnalyzeEvent Begin " << endl;

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
    tops.push_back(dl);

    fTNPV = NPV;
    float max_rap = -1;
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
            p.set_user_info(new MyUserInfo(pythia_MB->event[particle_idx].id(),particle_idx,pythia_MB->event[particle_idx].charge(),true));

            max_rap = max_rap < fabs(p.rap()) ? fabs(p.rap()) : max_rap;
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
        p.set_user_info(new MyUserInfo(pythia8->event[particle_idx].id(),particle_idx,pythia8->event[particle_idx].charge(),false));

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

        max_rap = max_rap < fabs(p.rap()) ? fabs(p.rap()) : max_rap;
        particles_for_jets.push_back(p);

    } // end particle loop -----------------------------------------------

    // Area Definition ---------------------
    fastjet::GhostedAreaSpec area_spec(max_rap);
    fastjet::AreaDefinition area_def(fastjet::active_area, area_spec);

    // N Subjettiness
    fastjet::contrib::Nsubjettiness n_subjettiness_1(1, fastjet::contrib::Njettiness::kt_axes, 1, 1, 1);
    fastjet::contrib::Nsubjettiness n_subjettiness_2(2, fastjet::contrib::Njettiness::kt_axes, 1, 1, 1);
    fastjet::contrib::Nsubjettiness n_subjettiness_3(3, fastjet::contrib::Njettiness::kt_axes, 1, 1, 1);

    // large-R jets: C/A --------------------
    fastjet::ClusterSequenceArea cs_large_r_ca(particles_for_jets, *m_jet_def_large_r_ca, area_def);
    vecPseudoJet my_jets_large_r_ca = fastjet::sorted_by_pt(cs_large_r_ca.inclusive_jets(pT_min));

    // this is a very temporary fix, it appears that there are no jets sometimes with pT at least 
    if (my_jets_large_r_ca.size() != 0) {
        fTCA_m = my_jets_large_r_ca[0].m();
        fTCA_pt = my_jets_large_r_ca[0].pt();
        fTCA_pufrac = JetPuFracFastjet(my_jets_large_r_ca[0]);
        fTCA_m_pu = JetPuMassFastjet(my_jets_large_r_ca[0]);
        if (my_jets_large_r_ca.size() >= 2) {
            fTCA_dr = my_jets_large_r_ca[0].delta_R(my_jets_large_r_ca[1]);
        }
        CA_nsubjettiness.push_back(n_subjettiness_1(my_jets_large_r_ca[0]));
        CA_nsubjettiness.push_back(n_subjettiness_2(my_jets_large_r_ca[0]));
        CA_nsubjettiness.push_back(n_subjettiness_3(my_jets_large_r_ca[0]));
    }

    // anti-kt R:1.0 trimmed ----------------
    fastjet::Filter filter_two(0.2, fastjet::SelectorPtFractionMin(0.05));
    fastjet::Filter filter_three(0.3, fastjet::SelectorPtFractionMin(0.05));

    fastjet::ClusterSequenceArea cs_large_r_antikt(particles_for_jets, *m_jet_def_large_r_antikt, area_def);
    vecPseudoJet my_jets_large_r_antikt = fastjet::sorted_by_pt(cs_large_r_antikt.inclusive_jets(pT_min));

    if (my_jets_large_r_antikt.size() != 0) {
        fastjet::PseudoJet lead_akt = my_jets_large_r_antikt[0];
        fastjet::PseudoJet lead_akt_filter_two = filter_two(lead_akt);
        fastjet::PseudoJet lead_akt_filter_three = filter_three(lead_akt);
        fTantikt_m = lead_akt.m();
        fTantikt_m_pu = JetPuMassFastjet(lead_akt);
        fTantikt_pt = lead_akt.pt();
        fTantikt_pufrac = JetPuFracFastjet(lead_akt);

        fTantikt_m_trimmed_two = lead_akt_filter_two.m();
        fTantikt_m_pu_trimmed_two = JetPuMassFastjet(lead_akt_filter_two);
        fTantikt_pt_trimmed_two = lead_akt_filter_two.pt();
        fTantikt_pufrac_trimmed_two = JetPuFracFastjet(lead_akt_filter_two);

        fTantikt_m_trimmed_three = lead_akt_filter_three.m();
        fTantikt_m_pu_trimmed_three = JetPuMassFastjet(lead_akt_filter_three);
        fTantikt_pt_trimmed_three = lead_akt_filter_three.pt();
        fTantikt_pufrac_trimmed_three = JetPuFracFastjet(lead_akt_filter_three);

        antikt_nsubjettiness.push_back(n_subjettiness_1(lead_akt));
        antikt_nsubjettiness.push_back(n_subjettiness_2(lead_akt));
        antikt_nsubjettiness.push_back(n_subjettiness_3(lead_akt));

        if (my_jets_large_r_antikt.size() >= 2) {
            fTantikt_dr = my_jets_large_r_antikt[0].delta_R(my_jets_large_r_antikt[1]);
        }

        for (unsigned int particle_iter = 0; particle_iter < particles_for_jets.size(); particle_iter++) {
            fastjet::PseudoJet const& particle = particles_for_jets[particle_iter];
            is_lead_antikt_constituent_vec.push_back(lead_akt.contains(particle));
        }
    }
    for (unsigned int particle_iter = 0; particle_iter < particles_for_jets.size(); particle_iter++) {
        fastjet::PseudoJet const& particle = particles_for_jets[particle_iter];
        is_pileup_vec.push_back(particle.user_info<MyUserInfo>().isPU());
    }

    bool do_substructure_study = false;
    if (do_substructure_study && !batched && event_iter < 10) {
        // generate particles for use in the substructure stuff
        fastjet::ClusterSequence cs_vlarge_r_antikt(particles_for_jets, *m_jet_def_vlarge_r_antikt);
        fastjet::ClusterSequence cs_vlarge_r_ca(particles_for_jets, *m_jet_def_vlarge_r_ca);
        vecPseudoJet my_jets_vlarge_r_antikt = fastjet::sorted_by_pt(cs_vlarge_r_antikt.inclusive_jets(pT_min));
        vecPseudoJet my_jets_vlarge_r_ca = fastjet::sorted_by_pt(cs_vlarge_r_ca.inclusive_jets(pT_min));
        if (!my_jets_vlarge_r_antikt.size() || !my_jets_vlarge_r_ca.size()) return;

        SubstructureStudy(my_jets_vlarge_r_ca, my_jets_vlarge_r_antikt, event_iter);
        return;
    }
    
    // ======================================
    // Various mixture models ---------------
    // ======================================
    tool->SetMergeDistance(0.05);
    double sigma_squared = 0.5;
    tool->SetDefaultSigma(MatTwo(sigma_squared, 0, 0, sigma_squared));

    // which jets to run
    bool mUMM_on = true;
    bool mGMM_on = true;
    bool mGMMc_on = true;
    bool mGMMs_on = true;
    bool mTGMM_on = true;
    bool mTGMMs_on = true;


    // Fuzzy Jets: mGMMs --------------------
    vector<vector<double> > mGMMs_particle_weights;
    vector<MatTwo> mGMMs_jets_params;
    vector<double> mGMMs_weights;
    vecPseudoJet mGMMs_jets;
    vector<int> mGMMs_indices;
    double max_pT_mGMMs;
    if(mGMMs_on) {
        DoMGMMJetFinding(particles_for_jets, my_jets_large_r_ca,
                         f_learn_weights, true, do_recombination,
                         mGMMs_indices, max_pT_mGMMs,
                         tool, mGMMs_jets, mGMMs_particle_weights,
                         mGMMs_jets_params, mGMMs_weights);
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
        DoMGMMCJetFinding(particles_for_jets, my_jets_large_r_ca,
                          f_learn_weights, do_recombination,
                          mGMMc_indices, max_pT_mGMMc,
                          tool, mGMMc_jets, mGMMc_particle_weights,
                          mGMMc_jets_params, mGMMc_weights);
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
        DoMTGMMJetFinding(particles_for_jets, my_jets_large_r_ca,
                          f_learn_weights, true, do_recombination,
                          f_size, mTGMMs_indices, max_pT_mTGMMs,
                          tool, mTGMMs_jets, mTGMMs_particle_weights,
                          mTGMMs_jets_params, mTGMMs_weights);
    }
    int lead_mTGMMs_index = mTGMMs_indices.size() ? mTGMMs_indices.at(0) : -1;

    // Fuzzy Jets: mGMM ---------------------
    vector<vector<double> >mGMM_particle_weights;
    vector<MatTwo> mGMM_jets_params;
    vector<double> mGMM_weights;
    vector<fastjet::PseudoJet> mGMM_jets;
    vector<int> mGMM_indices;
    double max_pT_mGMM;
    if(mGMM_on) {
        DoMGMMJetFinding(particles_for_jets, my_jets_large_r_ca,
                         f_learn_weights, false, do_recombination,
                         mGMM_indices, max_pT_mGMM,
                         tool, mGMM_jets, mGMM_particle_weights,
                         mGMM_jets_params, mGMM_weights);
    }
    int lead_mGMM_index = mGMM_indices.size() ? mGMM_indices.at(0) : -1;

    // Fuzzy Jets: mUMM ---------------------
    vector<vector<double> > mUMM_particle_weights;
    vector<double> mUMM_weights;
    vecPseudoJet mUMM_jets;
    vector<int> mUMM_indices;
    double max_pT_mUMM;
    if(mUMM_on) {
        DoMUMMJetFinding(particles_for_jets, my_jets_large_r_ca,
                         f_learn_weights, f_size, do_recombination,
                         mUMM_indices, max_pT_mUMM, tool, mUMM_jets,
                         mUMM_particle_weights, mUMM_weights);
    }
    int lead_mUMM_index = mUMM_indices.size() ? mUMM_indices.at(0) : -1;

    // Fuzzy Jets: mTGMM --------------------
    vector<vector<double> > mTGMM_particle_weights;
    vector<double> mTGMM_weights;
    vecPseudoJet mTGMM_jets;
    vector<MatTwo> mTGMM_jets_params;
    vector<int> mTGMM_indices;
    double max_pT_mTGMM;
    if(mTGMM_on) {
        DoMTGMMJetFinding(particles_for_jets, my_jets_large_r_ca,
                          f_learn_weights, false, do_recombination, f_size,
                          mTGMM_indices, max_pT_mTGMM, tool, mTGMM_jets,
                          mTGMM_particle_weights, mTGMM_jets_params, mTGMM_weights);
    }
    int lead_mTGMM_index = mTGMM_indices.size() ? mTGMM_indices.at(0) : -1;

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
            tool->EventDisplay(particles_for_jets,
                               my_jets_large_r_ca,tops,
                               mGMM_jets,
                               mGMM_particle_weights,
                               lead_mGMM_index,
                               mGMM_jets_params,
                               "",
                               event_iter);
        }
        if(mGMM_on) {
            tool->NewEventDisplay(particles_for_jets,
                                  my_jets_large_r_ca,tops,
                                  mGMM_jets,
                                  mGMM_particle_weights,
                                  lead_mGMM_index,
                                  mGMM_jets_params,
                                  mGMM_weights,
                                  "mGMM",
                                  event_iter);
            tool->JetContributionDisplay(particles_for_jets,
                                         mGMM_particle_weights,
                                         lead_mGMM_index,
                                         1, 
                                         "m_mGMM",
                                         event_iter);
            tool->JetContributionDisplay(particles_for_jets,
                                         mGMM_particle_weights,
                                         lead_mGMM_index,
                                         0, 
                                         "pt_mGMM",
                                         event_iter);
        }
        if(mTGMM_on) {
            tool->NewEventDisplay(particles_for_jets,
                                  my_jets_large_r_ca,tops,
                                  mTGMM_jets,
                                  mTGMM_particle_weights,
                                  lead_mTGMM_index,
                                  mTGMM_jets_params,
                                  mTGMM_weights,
                                  "mTGMM",
                                  event_iter);
        }
        if(mUMM_on) {
            tool->NewEventDisplayUniform(particles_for_jets,
                                         my_jets_large_r_ca,tops,
                                         mUMM_jets,
                                         mUMM_particle_weights,
                                         lead_mUMM_index,
                                         mUMM_weights,
                                         "mUMM",
                                         event_iter);
        }
        if(mGMMs_on) {
            tool->NewEventDisplay(particles_for_jets,
                                  my_jets_large_r_ca, tops,
                                  mGMMs_jets,
                                  mGMMs_particle_weights,
                                  lead_mGMMs_index,
                                  mGMMs_jets_params,
                                  mGMMs_weights,
                                  "mGMMs",
                                  event_iter);
        }
        if(mGMMc_on) {
            tool->NewEventDisplay(particles_for_jets,
                                  my_jets_large_r_ca, tops,
                                  mGMMc_jets,
                                  mGMMc_particle_weights,
                                  lead_mGMMc_index,
                                  mGMMc_jets_params,
                                  mGMMc_weights,
                                  "mGMMc",
                                  event_iter);
            tool->NewEventDisplayPoster(particles_for_jets,
                                  my_jets_large_r_ca, tops,
                                  mGMMc_jets,
                                  mGMMc_particle_weights,
                                  lead_mGMMc_index,
                                  mGMMc_jets_params,
                                  mGMMc_weights,
                                  "mGMMc",
                                  event_iter);
        }
        if(mTGMMs_on) {
            tool->NewEventDisplay(particles_for_jets,
                                  my_jets_large_r_ca,tops,
                                  mTGMMs_jets,
                                  mTGMMs_particle_weights,
                                  lead_mTGMMs_index,
                                  mTGMMs_jets_params,
                                  mTGMMs_weights,
                                  "mTGMMs",
                                  event_iter);
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
    for(unsigned int particle_iter = 0; particle_iter < particles_for_jets.size(); particle_iter++) {
        fastjet::PseudoJet const& p = particles_for_jets[particle_iter];
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

    if(mGMM_on) {
        fTmGMM_m = tool->MLpT(particles_for_jets,mGMM_particle_weights,
                              lead_mGMM_index,mGMM_particle_weights[0].size(),1);
        fTmGMM_pt = tool->MLpT(particles_for_jets,mGMM_particle_weights,
                               lead_mGMM_index,mGMM_particle_weights[0].size(),0);
        fTmGMM_ml = tool->MLlpTGaussian(particles_for_jets,mGMM_jets[lead_mGMM_index],
                                        mGMM_jets_params[lead_mGMM_index], mGMM_weights[lead_mGMM_index], 1);
        fTmGMM_ptl = tool->MLlpTGaussian(particles_for_jets,mGMM_jets[lead_mGMM_index],
                                         mGMM_jets_params[lead_mGMM_index], mGMM_weights[lead_mGMM_index], 0);
        fTmGMM_pufrac_soft = JetPuFracSoft(particles_for_jets, mGMM_particle_weights, lead_mGMM_index);
        fTmGMM_pufrac_hard = JetPuFracHard(particles_for_jets, mGMM_particle_weights, lead_mGMM_index);
        fTmGMM_m_pu_soft = JetPuMassSoft(particles_for_jets, mGMM_particle_weights, lead_mGMM_index);
        fTmGMM_m_pu_hard = JetPuMassHard(particles_for_jets, mGMM_particle_weights, lead_mGMM_index);
    }
    if(mUMM_on) {
        fTmUMM_m = tool->MLpT(particles_for_jets, mUMM_particle_weights,
                              lead_mUMM_index, mUMM_particle_weights[0].size(), 1);
        fTmUMM_pt = max_pT_mUMM;
        fTmUMM_ml = tool->MLlpTUniform(particles_for_jets,mUMM_jets[lead_mUMM_index],
                                       mUMM_weights[lead_mUMM_index], 1);
        fTmUMM_ptl = tool->MLlpTUniform(particles_for_jets,mUMM_jets[lead_mUMM_index],
                                        mUMM_weights[lead_mUMM_index], 0);
        fTmUMM_pufrac_soft = JetPuFracSoft(particles_for_jets, mUMM_particle_weights, lead_mUMM_index);
        fTmUMM_pufrac_hard = JetPuFracHard(particles_for_jets, mUMM_particle_weights, lead_mUMM_index);
        fTmUMM_m_pu_soft = JetPuMassSoft(particles_for_jets, mUMM_particle_weights, lead_mUMM_index);
        fTmUMM_m_pu_hard = JetPuMassHard(particles_for_jets, mUMM_particle_weights, lead_mUMM_index);
    }
    if(mTGMM_on) {
        fTmTGMM_m = tool->MLpT(particles_for_jets, mTGMM_particle_weights,
                               lead_mTGMM_index, mTGMM_particle_weights[0].size(), 1);
        fTmTGMM_pt = max_pT_mTGMM;
        fTmTGMM_ml = tool->MLlpTTruncGaus(particles_for_jets,mTGMM_jets[lead_mTGMM_index],
                                          mTGMM_jets_params[lead_mTGMM_index], mTGMM_weights[lead_mTGMM_index],1);
        fTmTGMM_ptl = tool->MLlpTTruncGaus(particles_for_jets,mTGMM_jets[lead_mTGMM_index],
                                           mTGMM_jets_params[lead_mTGMM_index], mTGMM_weights[lead_mTGMM_index],0);
        fTmTGMM_pufrac_soft = JetPuFracSoft(particles_for_jets, mTGMM_particle_weights, lead_mTGMM_index);
        fTmTGMM_pufrac_hard = JetPuFracHard(particles_for_jets, mTGMM_particle_weights, lead_mTGMM_index);
        fTmTGMM_m_pu_soft = JetPuMassSoft(particles_for_jets, mTGMM_particle_weights, lead_mTGMM_index);
        fTmTGMM_m_pu_hard = JetPuMassHard(particles_for_jets, mTGMM_particle_weights, lead_mTGMM_index);
    }
    if(mTGMMs_on) {
        fTmTGMMs_m = tool->MLpT(particles_for_jets, mTGMMs_particle_weights,
                                lead_mTGMMs_index, mTGMMs_particle_weights[0].size(), 1);
        fTmTGMMs_pt = max_pT_mTGMMs;
        fTmTGMMs_ml = tool->MLlpTTruncGaus(particles_for_jets,mTGMMs_jets[lead_mTGMMs_index],
                                           mTGMMs_jets_params[lead_mTGMMs_index], mTGMMs_weights[lead_mTGMMs_index],1);
        fTmTGMMs_ptl = tool->MLlpTTruncGaus(particles_for_jets,mTGMMs_jets[lead_mTGMMs_index],
                                            mTGMMs_jets_params[lead_mTGMMs_index], mTGMMs_weights[lead_mTGMMs_index],0);
        fTmTGMMs_pufrac_soft = JetPuFracSoft(particles_for_jets, mTGMMs_particle_weights, lead_mTGMMs_index);
        fTmTGMMs_pufrac_hard = JetPuFracHard(particles_for_jets, mTGMMs_particle_weights, lead_mTGMMs_index);
        fTmTGMMs_m_pu_soft = JetPuMassSoft(particles_for_jets, mTGMMs_particle_weights, lead_mTGMMs_index);
        fTmTGMMs_m_pu_hard = JetPuMassHard(particles_for_jets, mTGMMs_particle_weights, lead_mTGMMs_index);
    }
    if(mGMMs_on) {
        fTmGMMs_m = tool->MLpT(particles_for_jets, mGMMs_particle_weights,
                               lead_mGMMs_index, mGMMs_particle_weights[0].size(), 1);
        fTmGMMs_pt = max_pT_mGMMs;
        fTmGMMs_ml = tool->MLlpTGaussian(particles_for_jets,mGMMs_jets[lead_mGMMs_index],
                                         mGMMs_jets_params[lead_mGMMs_index], mGMMs_weights[lead_mGMMs_index],1);
        fTmGMMs_ptl = tool->MLlpTGaussian(particles_for_jets,mGMMs_jets[lead_mGMMs_index],
                                          mGMMs_jets_params[lead_mGMMs_index], mGMMs_weights[lead_mGMMs_index],0);
        fTmGMMs_pufrac_soft = JetPuFracSoft(particles_for_jets, mGMMs_particle_weights, lead_mGMMs_index);
        fTmGMMs_pufrac_hard = JetPuFracHard(particles_for_jets, mGMMs_particle_weights, lead_mGMMs_index);
        fTmGMMs_m_pu_soft = JetPuMassSoft(particles_for_jets, mGMMs_particle_weights, lead_mGMMs_index);
        fTmGMMs_m_pu_hard = JetPuMassHard(particles_for_jets, mGMMs_particle_weights, lead_mGMMs_index);
    }
    if(mGMMc_on) {
        fTmGMMc_m = tool->MLpT(particles_for_jets, mGMMc_particle_weights,
                               lead_mGMMc_index, mGMMc_particle_weights[0].size(), 1);
        fTmGMMc_pt = max_pT_mGMMc;
        fTmGMMc_ml = tool->MLlpTGaussian(particles_for_jets,mGMMc_jets[lead_mGMMc_index],
                                         mGMMc_jets_params[lead_mGMMc_index], mGMMc_weights[lead_mGMMc_index],1);
        fTmGMMc_ptl = tool->MLlpTGaussian(particles_for_jets,mGMMc_jets[lead_mGMMc_index],
                                          mGMMc_jets_params[lead_mGMMc_index], mGMMc_weights[lead_mGMMc_index],0);
        fTmGMMc_pufrac_soft = JetPuFracSoft(particles_for_jets, mGMMc_particle_weights, lead_mGMMc_index);
        fTmGMMc_pufrac_hard = JetPuFracHard(particles_for_jets, mGMMc_particle_weights, lead_mGMMc_index);
        fTmGMMc_m_pu_soft = JetPuMassSoft(particles_for_jets, mGMMc_particle_weights, lead_mGMMc_index);
        fTmGMMc_m_pu_hard = JetPuMassHard(particles_for_jets, mGMMc_particle_weights, lead_mGMMc_index);
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
        fTmGMM_m_soft = tool->SoftpT(particles_for_jets,
                                     mGMM_particle_weights,
                                     lead_mGMM_index,
                                     1);
        fTmGMM_pt_soft = tool->SoftpT(particles_for_jets,
                                      mGMM_particle_weights,
                                      lead_mGMM_index,
                                      0);
        fTmGMM_ucpu = UnclusteredPileupComposition(particles_for_jets,
                                                   mGMM_particle_weights);
        fTmGMM_clpu = ClusteredPileupFrac(particles_for_jets,
                                          mGMM_particle_weights);
    }
    if(mUMM_on) {
        fTdeltatop_mUMM = tops[0].delta_R(mUMM_jets[lead_mUMM_index]);
        if (tops[1].delta_R(mUMM_jets[lead_mUMM_index]) <  fTdeltatop_mUMM) {
            fTdeltatop_mUMM = tops[1].delta_R(mUMM_jets[lead_mUMM_index]);
        }
        fTmUMM_m_soft = tool->SoftpT(particles_for_jets,
                                     mUMM_particle_weights,
                                     lead_mUMM_index,
                                     1);
        fTmUMM_pt_soft = tool->SoftpT(particles_for_jets,
                                      mUMM_particle_weights,
                                      lead_mUMM_index,
                                      0);
        fTmUMM_ucpu = UnclusteredPileupComposition(particles_for_jets,
                                                   mUMM_particle_weights);
        fTmUMM_clpu = ClusteredPileupFrac(particles_for_jets,
                                          mUMM_particle_weights);
    }
    if(mTGMM_on) {
        fTdeltatop_mTGMM = tops[0].delta_R(mTGMM_jets[lead_mTGMM_index]);
        if (tops[1].delta_R(mTGMM_jets[lead_mTGMM_index]) < fTdeltatop_mTGMM) {
            fTdeltatop_mTGMM = tops[1].delta_R(mTGMM_jets[lead_mTGMM_index]);
        }
        fTmTGMM_m_soft = tool->SoftpT(particles_for_jets,
                                      mTGMM_particle_weights,
                                      lead_mTGMM_index,
                                      1);
        fTmTGMM_pt_soft = tool->SoftpT(particles_for_jets,
                                       mTGMM_particle_weights,
                                       lead_mTGMM_index,
                                       0);
        fTmTGMM_ucpu = UnclusteredPileupComposition(particles_for_jets,
                                                    mTGMM_particle_weights);
        fTmTGMM_clpu = ClusteredPileupFrac(particles_for_jets,
                                           mTGMM_particle_weights);
    }
    if(mTGMMs_on) {
        fTdeltatop_mTGMMs = tops[0].delta_R(mTGMMs_jets[lead_mTGMMs_index]);
        if (tops[1].delta_R(mTGMMs_jets[lead_mTGMMs_index]) < fTdeltatop_mTGMMs) {
            fTdeltatop_mTGMMs = tops[1].delta_R(mTGMMs_jets[lead_mTGMMs_index]);
        }
        fTmTGMMs_m_soft = tool->SoftpT(particles_for_jets,
                                       mTGMMs_particle_weights,
                                       lead_mTGMMs_index,
                                       1);
        fTmTGMMs_pt_soft = tool->SoftpT(particles_for_jets,
                                        mTGMMs_particle_weights,
                                        lead_mTGMMs_index,
                                        0);
        fTmTGMMs_ucpu = UnclusteredPileupComposition(particles_for_jets,
                                                     mTGMMs_particle_weights);
        fTmTGMMs_clpu = ClusteredPileupFrac(particles_for_jets,
                                            mTGMMs_particle_weights);
    }
    if(mGMMs_on) {
        fTdeltatop_mGMMs = tops[0].delta_R(mGMMs_jets[lead_mGMMs_index]);
        if (tops[1].delta_R(mGMMs_jets[lead_mGMMs_index]) < fTdeltatop_mGMMs) {
            fTdeltatop_mGMMs = tops[1].delta_R(mGMMs_jets[lead_mGMMs_index]);
        }
        fTmGMMs_m_soft = tool->SoftpT(particles_for_jets,
                                      mGMMs_particle_weights,
                                      lead_mGMMs_index,
                                      1);
        fTmGMMs_pt_soft = tool->SoftpT(particles_for_jets,
                                       mGMMs_particle_weights,
                                       lead_mGMMs_index,
                                       0);
        fTmGMMs_ucpu = UnclusteredPileupComposition(particles_for_jets,
                                                    mGMMs_particle_weights);
        fTmGMMs_clpu = ClusteredPileupFrac(particles_for_jets,
                                           mGMMs_particle_weights);
    }
    if(mGMMc_on) {
        fTdeltatop_mGMMc = tops[0].delta_R(mGMMc_jets[lead_mGMMc_index]);
        if (tops[1].delta_R(mGMMc_jets[lead_mGMMc_index]) < fTdeltatop_mGMMc) {
            fTdeltatop_mGMMc = tops[1].delta_R(mGMMc_jets[lead_mGMMc_index]);
        }
        fTmGMMc_m_soft = tool->SoftpT(particles_for_jets,
                                      mGMMc_particle_weights,
                                      lead_mGMMc_index,
                                      1);
        fTmGMMc_pt_soft = tool->SoftpT(particles_for_jets,
                                       mGMMc_particle_weights,
                                       lead_mGMMc_index,
                                       0);
        fTmGMMc_ucpu = UnclusteredPileupComposition(particles_for_jets,
                                                    mGMMc_particle_weights);
        fTmGMMc_clpu = ClusteredPileupFrac(particles_for_jets,
                                           mGMMc_particle_weights);
    }

    // Moments
    vector<double> moments_m;
    vector<double> moments_pt;
    double sig;
    if(mGMM_on) {
        moments_m = tool->CentralMoments(particles_for_jets, mGMM_particle_weights,
                                         lead_mGMM_index, 3, &totalMass);

        moments_pt = tool->CentralMoments(particles_for_jets, mGMM_particle_weights,
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
        moments_m = tool->CentralMoments(particles_for_jets, mUMM_particle_weights,
                                         lead_mUMM_index, 3, &totalMass);

        moments_pt = tool->CentralMoments(particles_for_jets, mUMM_particle_weights,
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
        moments_m = tool->CentralMoments(particles_for_jets, mTGMM_particle_weights,
                                         lead_mTGMM_index, 3, &totalMass);

        moments_pt = tool->CentralMoments(particles_for_jets, mTGMM_particle_weights,
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
        moments_m = tool->CentralMoments(particles_for_jets, mTGMMs_particle_weights,
                                         lead_mTGMMs_index, 3, &totalMass);

        moments_pt = tool->CentralMoments(particles_for_jets, mTGMMs_particle_weights,
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
        moments_m = tool->CentralMoments(particles_for_jets, mGMMs_particle_weights,
                                         lead_mGMMs_index, 3, &totalMass);

        moments_pt = tool->CentralMoments(particles_for_jets, mGMMs_particle_weights,
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
        moments_m = tool->CentralMoments(particles_for_jets, mGMMc_particle_weights,
                                         lead_mGMMc_index, 3, &totalMass);

        moments_pt = tool->CentralMoments(particles_for_jets, mGMMc_particle_weights,
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
    t_t->Fill();

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

    // Event Properties
    #ifdef WITHROOT
    t_t->Branch("EventNumber",    &fTEventNumber,    "EventNumber/I");
    t_t->Branch("NPV",            &fTNPV,            "NPV/I");

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
    map<string, float*>::const_iterator iter;
    for (iter = tree_vars.begin(); iter != tree_vars.end(); iter++) {
        *(iter->second) = -999999;
    }
    
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
