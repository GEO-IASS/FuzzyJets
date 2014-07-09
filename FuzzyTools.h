#ifndef FUZZYTOOLS_H
#define FUZZYTOOLS_H

#include <vector>
#include <set>
#include <math.h>
#include <string>
#include "TMatrix.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

#include "Pythia8/Pythia.h"

#include "myFastJetBase.h"

using namespace std;

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

    vector<double> CentralMoments(vecPseudoJet particles,
                                  vector<vector<double> >Weights,
                                  unsigned int clusterID,
                                  unsigned int momentCount,
                                  double (*f) (vecPseudoJet, vector<unsigned int>));

    void SetClusteringMode(Mode m) {
        clusteringMode = m;
    }

    void SetAlpha(double a) {
        alpha = a;
    }

    set<unsigned int> ClustersForRemovalGaussian(vecPseudoJet& mGMMjets,
                                                 vector<TMatrix>& mGMMjetsparams,
                                                 vector<double>& mGMMweights);

    set<unsigned int> ClustersForRemovalUniform(vecPseudoJet& mUMMjets,
                                                vector<double>& mUMMweights);

    void SetSeeds(vecPseudoJet s) {
        seeds = s;
    }

    vector<vector<double> > InitWeights(vecPseudoJet particles,int k);

    double MDist(double x1, double x2, double mu1, double mu2, TMatrix sigma);

    double doGaus(double x1, double x2, double mu1, double mu2, TMatrix sigma);

    double doTruncGaus(double x1, double x2, double mu1, double mu2, TMatrix sigma);

    vector<TMatrix> Initializeparams(vecPseudoJet particles,
                                     int k);

    void ComputeWeightsGaussian(vecPseudoJet particles,
                                vector<vector<double> >* Weights,
                                int k,
                                vecPseudoJet mGMMjets,
                                vector<TMatrix> mGMMjetsparams,
                                vector<double> mGMMweights);

    vecPseudoJet UpdateJetsGaussian(vecPseudoJet particles,
                                    vector<vector<double> > Weights,
                                    int k,
                                    vector<TMatrix>* mGMMjetsparams,
                                    vector<double>* mGMMweights);

    vecPseudoJet ClusterFuzzyGaussian(vecPseudoJet particles,
                                      vector<vector<double> >* Weights,
                                      vector<TMatrix>* mGMMjetsparamsout,
                                      vector<double>* mGMMweightsout);

    void ComputeWeightsUniform(vecPseudoJet particles,
                               vector<vector<double> >* Weights,
                               int k,
                               vecPseudoJet mUMMjets,
                               vector<double> mUMMweights);

    vecPseudoJet UpdateJetsUniform(vecPseudoJet particles,
                                   vector<vector<double> > Weights,
                                   int k,
                                   vector<double>* mUMMweights);

    vecPseudoJet ClusterFuzzyUniform(vecPseudoJet particles,
                                     vector<vector<double> >* Weights,
                                     vector<double>* mUMMweightsout);

    void ComputeWeightsTruncGaus(vecPseudoJet particles,
                                 vector<vector<double> >* Weights,
                                 int k,
                                 vecPseudoJet mTGMMjets,
                                 vector<TMatrix> mTGMMjetsparams,
                                 vector<double> mTGMMweights);

    vecPseudoJet UpdateJetsTruncGaus(vecPseudoJet particles,
                                     vector<vector<double> > Weights,
                                     int k,
                                     vector<TMatrix>* mTGMMjetsparams,
                                     vector<double>* mTGMMweights);

    vecPseudoJet ClusterFuzzyTruncGaus(vecPseudoJet particles,
                                       vector<vector<double> >* Weights,
                                       vector<TMatrix>* mTGMMjetsparamsout,
                                       vector<double>* mTGMMweights);

    void EventDisplay(vecPseudoJet particles,
                      vecPseudoJet CAjets,
                      vecPseudoJet tops,
                      vecPseudoJet mGMMjets,
                      vector<vector<double> > Weights,
                      int which,
                      vector<TMatrix> mGMMjetsparams,
                      TString out);

    void NewEventDisplay(vecPseudoJet particles,
                         vecPseudoJet CAjets,
                         vecPseudoJet tops,
                         vecPseudoJet mGMMjets,
                         vector<vector<double> > Weights,
                         int which,
                         vector<TMatrix> mGMMjetsparams,
                         vector<double> mGMMweights,
                         TString out);

    void NewEventDisplayUniform(vecPseudoJet particles,
                                vecPseudoJet CAjets,
                                vecPseudoJet tops,
                                vecPseudoJet mUMMjets,
                                vector<vector<double> > Weights,
                                int which,
                                vector<double> mUMMweights,
                                TString out);

    double MLpT(vecPseudoJet particles,
                vector<vector<double> > Weights,
                int jetindex,
                int k,
                int mtype);

    void Qjetmass(vecPseudoJet particles,
                  vector<vector<double> > Weights,
                  int which,
                  TString out);
};

#endif
