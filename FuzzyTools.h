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
    int maxiters;

    string directoryPrefix;

    bool learnWeights;
    bool learnShape;

    double mergeDist;
    double minWeight;
    double minSigma;

    int clusterCount;
    Mode clusteringMode;
    Kernel kernelType;
    vecPseudoJet seeds;

 public:
    FuzzyTools();

    // methods
    vecPseudoJet Initialize(vecPseudoJet particles,
                            int k,
                            vecPseudoJet jets);

    void SetLearnWeights(bool b) {
        learnWeights = b;
    }

    void SetR(double r) {
        R = r;
    }

    void SetKernelType(Kernel k) {
        kernelType = k;
    }

    void SetMergeDistance(double d) {
        mergeDist = d;
    }

    void SetMinimumWeight(double m) {
        minWeight = m;
    }

    void SetLearnShape(bool s) {
        learnShape = s;
    }

    void SetMinimumSigma(double s) {
        minSigma = s;
    }

    void SetPrefix(string p) {
        directoryPrefix = p;
    }

    vector<double> CentralMoments(vecPseudoJet const& particles,
                                  vector<vector<double> > const& Weights,
                                  unsigned int clusterID,
                                  unsigned int momentCount,
                                  double (*f) (vecPseudoJet, vector<unsigned int>));

    void SetClusteringMode(Mode m) {
        clusteringMode = m;
    }

    void SetAlpha(double a) {
        alpha = a;
    }

    set<unsigned int> ClustersForRemovalGaussian(vecPseudoJet const& mGMMjets,
                                                 vector<TMatrix> const& mGMMjetsparams,
                                                 vector<double> const& mGMMweights);

    set<unsigned int> ClustersForRemovalUniform(vecPseudoJet const& mUMMjets,
                                                vector<double> const& mUMMweights);

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
                                vector<vector<double> >* Weights,
                                int k,
                                vecPseudoJet const& mGMMjets,
                                vector<TMatrix> const& mGMMjetsparams,
                                vector<double> const& mGMMweights);

    vecPseudoJet UpdateJetsGaussian(vecPseudoJet const& particles,
                                    vector<vector<double> > const& Weights,
                                    int k,
                                    vector<TMatrix>* mGMMjetsparams,
                                    vector<double>* mGMMweights);

    vecPseudoJet ClusterFuzzyGaussian(vecPseudoJet const& particles,
                                      vector<vector<double> >* Weights,
                                      vector<TMatrix>* mGMMjetsparamsout,
                                      vector<double>* mGMMweightsout);

    void ComputeWeightsUniform(vecPseudoJet const& particles,
                               vector<vector<double> >* Weights,
                               int k,
                               vecPseudoJet const& mUMMjets,
                               vector<double> const& mUMMweights);

    vecPseudoJet UpdateJetsUniform(vecPseudoJet const& particles,
                                   vector<vector<double> > const& Weights,
                                   int k,
                                   vector<double>* mUMMweights);

    vecPseudoJet ClusterFuzzyUniform(vecPseudoJet const& particles,
                                     vector<vector<double> >* Weights,
                                     vector<double>* mUMMweightsout);

    void ComputeWeightsTruncGaus(vecPseudoJet const& particles,
                                 vector<vector<double> >* Weights,
                                 int k,
                                 vecPseudoJet const& mTGMMjets,
                                 vector<TMatrix> const& mTGMMjetsparams,
                                 vector<double> const& mTGMMweights);

    vecPseudoJet UpdateJetsTruncGaus(vecPseudoJet const& particles,
                                     vector<vector<double> > const& Weights,
                                     int k,
                                     vector<TMatrix>* mTGMMjetsparams,
                                     vector<double>* mTGMMweights);

    vecPseudoJet ClusterFuzzyTruncGaus(vecPseudoJet const& particles,
                                       vector<vector<double> >* Weights,
                                       vector<TMatrix>* mTGMMjetsparamsout,
                                       vector<double>* mTGMMweights);

    void EventDisplay(vecPseudoJet const& particles,
                      vecPseudoJet const& CAjets,
                      vecPseudoJet const& tops,
                      vecPseudoJet const& mGMMjets,
                      vector<vector<double> > const& Weights,
                      int which,
                      vector<TMatrix> const& mGMMjetsparams,
                      TString out);

    void NewEventDisplay(vecPseudoJet const& particles,
                         vecPseudoJet const& CAjets,
                         vecPseudoJet const& tops,
                         vecPseudoJet const& mGMMjets,
                         vector<vector<double> > const& Weights,
                         int which,
                         vector<TMatrix> const& mGMMjetsparams,
                         vector<double> const& mGMMweights,
                         TString const& out);

    void NewEventDisplayUniform(vecPseudoJet const& particles,
                                vecPseudoJet const& CAjets,
                                vecPseudoJet const& tops,
                                vecPseudoJet const& mUMMjets,
                                vector<vector<double> > const& Weights,
                                int which,
                                vector<double> const& mUMMweights,
                                TString const& out);

    void JetContributionDisplay(vecPseudoJet particles,
                                vector<vector<double> > Weights,
                                int which,
                                int mtype,
                                TString out);

    double MLpT(vecPseudoJet particles,
                vector<vector<double> > Weights,
                int jetindex,
                int k,
                int mtype);

    double MLlpTGaussian(vecPseudoJet const& particles,
                         fastjet::PseudoJet const& jet,
                         TMatrix const& jetParams,
                         double jetWeight,
                         int mtype);

    double MLlpTUniform(vecPseudoJet const& particles,
                        fastjet::PseudoJet const& jet,
                        double jetWeight,
                        int mtype);

    double MLlpTTruncGaus(vecPseudoJet const& particles,
                          fastjet::PseudoJet const& jet,
                          TMatrix const& jetParams,
                          double jetWeight,
                          int mtype);

    void Qjetmass(vecPseudoJet particles,
                  vector<vector<double> > Weights,
                  int which,
                  TString out);

    double SoftpT(vecPseudoJet const& particles,
                  vector<vector<double> > const& Weights,
                  int jetindex,
                  int mtype);

};

#endif
