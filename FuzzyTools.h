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

#include "TMatrix.h"

typedef std::vector<fastjet::PseudoJet> vecPseudoJet;


// declare helper functions to CentralMoments so they are accessible
// to any client code
double totalMass(vecPseudoJet particles, vector<unsigned int> indices);
double totalpT(vecPseudoJet particles, vector<unsigned int> indices);

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

    enum Mode {
        NOCLUSTERINGMODE,
        RECOMBINATION,
        FIXED
    };

 private:
    int m_test;
    double alpha;
    double R;
    int max_iters;

    string directory_prefix;

    bool learn_weights;
    bool learn_shape;

    double merge_dist;
    double min_weight;
    double min_sigma;

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

    void SetLearnWeights(bool b) {
        learn_weights = b;
    }

    void SetR(double r) {
        R = r;
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

    set<unsigned int> ClustersForRemovalGaussian(vecPseudoJet const& mGMM_jets,
                                                 vector<TMatrix> const& mGMM_jets_params,
                                                 vector<double> const& mGMM_weights);

    set<unsigned int> ClustersForRemovalUniform(vecPseudoJet const& mUMM_jets,
                                                vector<double> const& mUMM_weights);

    void SetSeeds(vecPseudoJet s) {
        seeds = s;
    }

    vector<vector<double> > InitWeights(vecPseudoJet const& particles,int k);

    double MDist(double x1, double x2, double mu1, double mu2, TMatrix const& sigma);

    double doGaus(double x1, double x2, double mu1, double mu2, TMatrix const& sigma);

    double doTruncGaus(double x1, double x2, double mu1, double mu2, TMatrix const& sigma);

    vector<TMatrix> Initializeparams(vecPseudoJet const& particles,
                                     int k);

    void ComputeWeightsGaussian(vecPseudoJet const& particles,
                                vector<vector<double> >* weights,
                                int k,
                                vecPseudoJet const& mGMM_jets,
                                vector<TMatrix> const& mGMM_jets_params,
                                vector<double> const& mGMM_weights);

    vecPseudoJet UpdateJetsGaussian(vecPseudoJet const& particles,
                                    vector<vector<double> > const& weights,
                                    int k,
                                    vector<TMatrix>* mGMM_jets_params,
                                    vector<double>* mGMM_weights);

    vecPseudoJet ClusterFuzzyGaussian(vecPseudoJet const& particles,
                                      vector<vector<double> >* weights,
                                      vector<TMatrix>* mGMM_jets_params_out,
                                      vector<double>* mGMM_weights_out);

    void ComputeWeightsUniform(vecPseudoJet const& particles,
                               vector<vector<double> >* weights,
                               int k,
                               vecPseudoJet const& mUMM_jets,
                               vector<double> const& mUMM_weights);

    vecPseudoJet UpdateJetsUniform(vecPseudoJet const& particles,
                                   vector<vector<double> > const& weights,
                                   int k,
                                   vector<double>* mUMM_weights);

    vecPseudoJet ClusterFuzzyUniform(vecPseudoJet const& particles,
                                     vector<vector<double> >* weights,
                                     vector<double>* mUMM_weights_out);

    void ComputeWeightsTruncGaus(vecPseudoJet const& particles,
                                 vector<vector<double> >* weights,
                                 int k,
                                 vecPseudoJet const& mTGMM_jets,
                                 vector<TMatrix> const& mTGMM_jets_params,
                                 vector<double> const& mTGMM_weights);

    vecPseudoJet UpdateJetsTruncGaus(vecPseudoJet const& particles,
                                     vector<vector<double> > const& weights,
                                     int k,
                                     vector<TMatrix>* mTGMM_jets_params,
                                     vector<double>* mTGMM_weights);

    vecPseudoJet ClusterFuzzyTruncGaus(vecPseudoJet const& particles,
                                       vector<vector<double> >* weights,
                                       vector<TMatrix>* mTGMM_jets_params_out,
                                       vector<double>* mTGMM_weights);

    void EventDisplay(vecPseudoJet const& particles,
                      vecPseudoJet const& ca_jets,
                      vecPseudoJet const& tops,
                      vecPseudoJet const& mGMM_jets,
                      vector<vector<double> > const& weights,
                      int which,
                      vector<TMatrix> const& mGMM_jets_params,
                      TString out);

    void NewEventDisplay(vecPseudoJet const& particles,
                         vecPseudoJet const& ca_jets,
                         vecPseudoJet const& tops,
                         vecPseudoJet const& mGMM_jets,
                         vector<vector<double> > const& weights,
                         int which,
                         vector<TMatrix> const& mGMM_jets_params,
                         vector<double> const& mGMM_weights,
                         TString const& out);

    void NewEventDisplayUniform(vecPseudoJet const& particles,
                                vecPseudoJet const& ca_jets,
                                vecPseudoJet const& tops,
                                vecPseudoJet const& mUMM_jets,
                                vector<vector<double> > const& weights,
                                int which,
                                vector<double> const& mUMM_weights,
                                TString const& out);

    void JetContributionDisplay(vecPseudoJet particles,
                                vector<vector<double> > weights,
                                int which,
                                int m_type,
                                TString out);

    double MLpT(vecPseudoJet particles,
                vector<vector<double> > weights,
                int jet_index,
                int k,
                int m_type);

    double MLlpTGaussian(vecPseudoJet const& particles,
                         fastjet::PseudoJet const& jet,
                         TMatrix const& jet_params,
                         double jet_weight,
                         int m_type);

    double MLlpTUniform(vecPseudoJet const& particles,
                        fastjet::PseudoJet const& jet,
                        double jet_weight,
                        int m_type);

    double MLlpTTruncGaus(vecPseudoJet const& particles,
                          fastjet::PseudoJet const& jet,
                          TMatrix const& jet_params,
                          double jet_weight,
                          int m_type);

    void Qjetmass(vecPseudoJet particles,
                  vector<vector<double> > weights,
                  int which,
                  TString out);

    double SoftpT(vecPseudoJet const& particles,
                  vector<vector<double> > const& weights,
                  int jet_index,
                  int m_type);

};

#endif
