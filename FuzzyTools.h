#ifndef FUZZYTOOLS_H
#define FUZZYTOOLS_H

#include <vector>
#include <set>
#include <math.h>
#include <string>

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

#include "Pythia8/Pythia.h"

#include "myFastJetBase.h"
#include "ROOTConf.h"

#ifdef WITHROOT
#include "TMatrix.h"
#endif

typedef std::vector<fastjet::PseudoJet> vecPseudoJet;

// declare helper functions to CentralMoments so they are accessible
// to any client code
double totalMass(vecPseudoJet particles, vector<unsigned int> indices);
double totalpT(vecPseudoJet particles, vector<unsigned int> indices);

// used to extract meaningful signed pTs
int pTSign(fastjet::PseudoJet p);
float signedPt(fastjet::PseudoJet p);

class MatTwo {
public:
    double xx;
    double xy;
    double yx;
    double yy;

    MatTwo(double xx_in, double xy_in, double yx_in, double yy_in)
        : xx(xx_in)
        , xy(xy_in)
        , yx(yx_in)
        , yy(yy_in)
    {}


    double determinant() const {
        return (xx * yy) - (xy * yx);
    }
};

#ifdef NEW
class FuzzyTools;

class AbstractKernel {
protected:
    std::vector<double> _weights;
    double _mu_x, _mu_y;
    double _weight;
public:
    AbstractKernel() {}

    virtual ~AbstractKernel() {}

    virtual double PDF(double x, double y) = 0;

    virtual void SetParticleWeight(int index, double value) {
        _weights.at(index) = value;
    }

    virtual void AdditionalUpdate(vecPseudoJet const& particles,
                                  FuzzyTools const& tool,
                                  double cluster_weighted_pt) = 0;

    double Weight() {
        return _weight;
    }

    void SetLocation(double mu_x, double mu_y) {
        _mu_x = mu_x;
        _mu_y = mu_y;
    }

    void SetWeight(double new_weight) {
        _weight = new_weight;
    }

    std::vector<double> const& ParticleWeights() {
        return _weights;
    }
};

class GaussianKernel : public AbstractKernel {
    double _mu_x, _mu_y;
    MatTwo _sigma;

public:
    GaussianKernel(double weight, double mu_x, double mu_y, MatTwo sigma)
        :_mu_x(mu_x), _mu_y(mu_y), _sigma(sigma) {
        _weight = weight;
    }

    ~GaussianKernel() {}

    virtual double PDF(double x, double y);

    virtual void AdditionalUpdate(vecPseudoJet const& particles,
                                  FuzzyTools const& tool,
                                  double cluster_weighted_pt);
};

void ClusterFuzzy(vecPseudoJet const& particles,
                  vector<AbstractKernel *> & jets,
                  FuzzyTools const& tool);

set<unsigned int> ClustersForRemoval(FuzzyTools const& tool,
                                     vector<AbstractKernel *>& jets);

vector<GaussianKernel *>
MakeGaussianKernels(FuzzyTools const& tool);
#endif

/* TODO:
 * There is no need to pass around constants giving cluster size etc,
 * as they can be extracted locally out of the C++ vector, which gives
 * O(1) size operation. At the moment though, everything is working,
 * so I won't bother to change this.
 *
 * Pass debug flags and print useful information.
 */
class FuzzyTools {
 public:
    enum Kernel {
        NOKERNEL,
        GAUSSIAN,
        TRUNCGAUSSIAN,
        UNIFORM
    };

    enum PostProcess {
        NO_POST_PROCESS,                  // Standard tool mode
        ONE_DISTANCE_MERGER    // Do one merger step at the end of clustering
    };

    enum EventJet {
        NONE,         // Standard tool mode
        FLAT,      // In this mode some likelihood is uniformly
        // smeared out as and attributed to pileup vertices.
        // Particle weights will *not* sum to 1 in this mode,
        // because every particle probabilistically belongs to this
        // pileup jet.
        ONE_CLUSTER   // Runs the tool on input particles with k=1 before running
                      // on k=k_0 + 1 where the last jet is the one found in the
                      // k=1 mode. Weights will once again not sum to 1 in this
                      // mode.
    };

    enum Mode {
        NOCLUSTERINGMODE,
        RECOMBINATION,
        FIXED
    };

 private:
    vector<vecPseudoJet> _historical_jets;
    vector<vector<MatTwo> > _historical_params;

    unsigned int _n_events;

    int m_test;
    double alpha;
    double R;
    int max_iters;

    MatTwo default_sigma;

    string directory_prefix;

    bool learn_weights;
    bool learn_shape;

    bool _make_trace_diagrams;

    double merge_dist;
    double min_weight;
    double min_sigma;

    double log_log_likelihood_limit;

    PostProcess post_process_method;

    EventJet event_jet_type;
    double event_jet_weight;

    int cluster_count;
    Mode clustering_mode;
    Kernel kernel_type;
    vecPseudoJet seeds;

 public:
    FuzzyTools();

    // methods
    vecPseudoJet Initialize(vecPseudoJet particles,
                            int k,
                            vecPseudoJet jets);

    void LogVectors(vecPseudoJet jets, vector<MatTwo> jet_parameters);

    void MakeTraceDiagram(vecPseudoJet const& particles,
                          int max_pT_idx,
                          std::string label, unsigned int iter);

    void SetEventJetWeight(double w) {
        event_jet_weight = w;
    }

    void SetDefaultSigma(MatTwo s) {
        default_sigma = s;
    }

    void SetNEvents (unsigned int n) {
        _n_events = n;
    }

    void SetLogLogLikelihoodLimit(double l) {
        log_log_likelihood_limit = l;
    }

    void SetMakeTraceDiagrams(bool b) {
        _make_trace_diagrams = b;
    }

    void SetLearnWeights(bool b) {
        learn_weights = b;
    }

    void SetR(double r) {
        R = r;
    }

    void SetPostProcess(PostProcess ppm) {
        post_process_method = ppm;
    }

    void SetEventJetType(EventJet type) {
        event_jet_type = type;
    }

    void SetKernelType(Kernel k) {
        kernel_type = k;
    }

    void SetMergeDistance(double d) {
        merge_dist = d;
    }

    void SetMinimumWeight(double m) {
        min_weight = m;
    }

    void SetLearnShape(bool s) {
        learn_shape = s;
    }

    void SetMinimumSigma(double s) {
        min_sigma = s;
    }

    void SetPrefix(string p) {
        directory_prefix = p;
    }

    vector<double> CentralMoments(vecPseudoJet const& particles,
                                  vector<vector<double> > const& weights,
                                  unsigned int cluster_id,
                                  unsigned int moment_count,
                                  double (*f) (vecPseudoJet, vector<unsigned int>));

    void SetClusteringMode(Mode m) {
        clustering_mode = m;
    }

    void SetAlpha(double a) {
        alpha = a;
    }

    PostProcess GetPostProcess() {
        return post_process_method;
    }

    double GetEventJetWeight() const {
        return event_jet_weight;
    }

    double GetAlpha() const {
        return alpha;
    }

    MatTwo GetDefaultSigma() const {
        return default_sigma;
    }

    vecPseudoJet const& GetSeeds() const {
        return seeds;
    }

    int GetMaxIters() const {
        return max_iters;
    }

    bool GetLearnWeights() const {
        return learn_weights;
    }

    Mode GetClusteringMode() const {
        return clustering_mode;
    }

    EventJet GetEventJetType() const {
        return event_jet_type;
    }

    set<unsigned int> ClustersForRemovalGaussian(vecPseudoJet const& mGMM_jets,
                                                 vector<MatTwo> const& mGMM_jets_params,
                                                 vector<double> const& mGMM_weights);

    set<unsigned int> ClustersForRemovalUniform(vecPseudoJet const& mUMM_jets,
                                                vector<double> const& mUMM_weights);

    void SetSeeds(vecPseudoJet s) {
        seeds = s;
    }

    vector<vector<double> > InitWeights(vecPseudoJet const& particles,int k);

    double MDist(double x1, double x2, double mu1, double mu2, MatTwo const& sigma);

    double doGaus(double x1, double x2, double mu1, double mu2, MatTwo const& sigma);

    double doTruncGaus(double x1, double x2, double mu1, double mu2, MatTwo const& sigma);

    vector<MatTwo> Initializeparams(vecPseudoJet const& particles,
                                     int k);

    void ComputeWeightsGaussian(vecPseudoJet const& particles,
                                vector<vector<double> >* weights,
                                int cluster_count,
                                vecPseudoJet const& mGMM_jets,
                                vector<MatTwo> const& mGMM_jets_params,
                                vector<double> const& mGMM_weights);

    vecPseudoJet UpdateJetsGaussian(vecPseudoJet const& particles,
                                    __attribute__((unused)) vecPseudoJet
                                    const& oldJets,
                                    vector<vector<double> > const& weights,
                                    int cluster_count,
                                    vector<MatTwo>* mGMM_jets_params,
                                    vector<double>* mGMM_weights);

    vecPseudoJet ClusterFuzzyGaussian(vecPseudoJet const& particles,
                                      vector<vector<double> >* weights,
                                      vector<MatTwo>* mGMM_jets_params_out,
                                      vector<double>* mGMM_weights_out,
                                      unsigned int &iter_count);

    vecPseudoJet UpdateJetsGaussianC(vecPseudoJet const& particles,
                                     __attribute__((unused)) vecPseudoJet
                                     const& oldJets,
                                     vector<vector<double> > const& weights,
                                     int cluster_count,
                                     vector<MatTwo>* mGMMc_jets_params,
                                     vector<double>* mGMMc_weights);

    vecPseudoJet ClusterFuzzyGaussianC(vecPseudoJet const& particles,
                                       vector<vector<double> >* weights,
                                       vector<MatTwo>* mGMMc_jets_params_out,
                                       vector<double>* mGMMc_weights_out,
                                       unsigned int &iter_count);

    void ComputeWeightsUniform(vecPseudoJet const& particles,
                               vector<vector<double> >* weights,
                               int cluster_count,
                               vecPseudoJet const& mUMM_jets,
                               vector<double> const& mUMM_weights);

    vecPseudoJet UpdateJetsUniform(vecPseudoJet const& particles,
                                   __attribute__((unused)) vecPseudoJet
                                   const& oldJets,
                                   vector<vector<double> > const& weights,
                                   int cluster_count,
                                   vector<double>* mUMM_weights);

    vecPseudoJet ClusterFuzzyUniform(vecPseudoJet const& particles,
                                     vector<vector<double> >* weights,
                                     vector<double>* mUMM_weights_out,
                                     unsigned int &iter_count);

    void ComputeWeightsTruncGaus(vecPseudoJet const& particles,
                                 vector<vector<double> >* weights,
                                 int cluster_count,
                                 vecPseudoJet const& mTGMM_jets,
                                 vector<MatTwo> const& mTGMM_jets_params,
                                 vector<double> const& mTGMM_weights);

    vecPseudoJet UpdateJetsTruncGaus(vecPseudoJet const& particles,
                                     __attribute__((unused)) vecPseudoJet
                                     const& oldJets,
                                     vector<vector<double> > const& weights,
                                     int cluster_count,
                                     vector<MatTwo>* mTGMM_jets_params,
                                     vector<double>* mTGMM_weights);

    vecPseudoJet ClusterFuzzyTruncGaus(vecPseudoJet const& particles,
                                       vector<vector<double> >* weights,
                                       vector<MatTwo>* mTGMM_jets_params_out,
                                       vector<double>* mTGMM_weights,
                                       unsigned int &iter_count);

    double LogLikelihoodGaussian(vecPseudoJet const& particles,
                                 vecPseudoJet const& mGMM_jets,
                                 vector<MatTwo> const& mGMM_jets_params,
                                 vector<double> const& mGMM_weights);

    double LogLikelihoodTruncGaus(vecPseudoJet const& particles,
                                  vecPseudoJet const& mTGMM_jets,
                                  vector<MatTwo> const& mTGMM_jets_params,
                                  vector<double> const& mTGMM_weights);

    double LogLikelihoodUniform(vecPseudoJet const& particles,
                                vecPseudoJet const& mUMM_jets,
                                vector<double> const& mUMM_weights);

    void SubsEventDisplay(vecPseudoJet const& particles,
                          vecPseudoJet const& mGMM_jets,
                          vector<vector<double> > const& weights,
                          int lead_mGMM_index,
                          vector<MatTwo> const& mGMM_jets_params,
                          std::string const& file_name,
                          std::string const& title);

    void EventJetDisplay(vecPseudoJet const& particles,
                                vecPseudoJet const& mGMM_jets,
                                vector<vector<double> > const& weights,
                                vector<MatTwo> const& mGMM_jets_params,
                                vector<double> const& mGMM_weights,
                                std::string const& label, int iter);

    void EventDisplay(vecPseudoJet const& particles,
                      vecPseudoJet const& ca_jets,
                      vecPseudoJet const& tops,
                      vecPseudoJet const& mGMM_jets,
                      vector<vector<double> > const& weights,
                      int which,
                      vector<MatTwo> const& mGMM_jets_params,
                      std::string const& out,
                      int iter);

    void ComparisonED(vecPseudoJet const& particles,
                      vecPseudoJet const& ca_jets,
                      vecPseudoJet const& tops,
                      vecPseudoJet const& mGMM_jets,
                      vector<vector<double> > const& weights,
                      int which,
                      vector<MatTwo> const& mGMM_jets_params,
                      vector<double> const& mGMM_weights,
                      std::string const& out,
                      int iter);

    void NewEventDisplay(vecPseudoJet const& particles,
                         vecPseudoJet const& ca_jets,
                         vecPseudoJet const& tops,
                         vecPseudoJet const& mGMM_jets,
                         vector<vector<double> > const& weights,
                         int which,
                         vector<MatTwo> const& mGMM_jets_params,
                         vector<double> const& mGMM_weights,
                         std::string const& out,
                         int iter);

    void NewEventDisplayPoster(vecPseudoJet const& particles,
                         vecPseudoJet const& ca_jets,
                         vecPseudoJet const& tops,
                         vecPseudoJet const& mGMM_jets,
                         vector<vector<double> > const& weights,
                         int which,
                         vector<MatTwo> const& mGMM_jets_params,
                         vector<double> const& mGMM_weights,
                         std::string const& out,
                         int iter);

    void NewEventDisplayUniform(vecPseudoJet const& particles,
                                vecPseudoJet const& ca_jets,
                                vecPseudoJet const& tops,
                                vecPseudoJet const& mUMM_jets,
                                vector<vector<double> > const& weights,
                                int which,
                                vector<double> const& mUMM_weights,
                                std::string const& out,
                                int iter);

    void JetContributionDisplay(vecPseudoJet particles,
                                vector<vector<double> > weights,
                                int which,
                                int m_type,
                                std::string out,
                                int iter);

    int belongs_idx(vector<vector<double> > const& weights,
                    int particle_index,
                    bool ignore_event_jet);

    double MLpT(vecPseudoJet particles,
                vector<vector<double> > weights,
                int jet_index,
                int k,
                int m_type,
                bool ignore_event_jet);

    double MLlpTGaussian(vecPseudoJet const& particles,
                         fastjet::PseudoJet const& jet,
                         MatTwo const& jet_params,
                         double jet_weight,
                         int m_type);

    double MLlpTUniform(vecPseudoJet const& particles,
                        fastjet::PseudoJet const& jet,
                        double jet_weight,
                        int m_type);

    double MLlpTTruncGaus(vecPseudoJet const& particles,
                          fastjet::PseudoJet const& jet,
                          MatTwo const& jet_params,
                          double jet_weight,
                          int m_type);

    void Qjetmass(vecPseudoJet particles,
                  vector<vector<double> > weights,
                  int which,
                  std::string out);

    double SoftpT(vecPseudoJet const& particles,
                  vector<vector<double> > const& weights,
                  int jet_index,
                  int m_type);

};

#endif
