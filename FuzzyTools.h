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

    bool learnWeights;

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

    void SetMinimumSigma(double s) {
        minSigma = s;
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
                                                    vector<TMatrix>& mGMMjetsparams);

    set<unsigned int> ClustersForRemovalUniform(vecPseudoJet& mUMMjets,
                                                   vector<double>& mUMMweights);

    void SetSeeds(vecPseudoJet s) {
        seeds = s;
    }

    vector<vector<double> > InitWeights(vecPseudoJet particles,int k);

    double MDist(double x1, double x2, double mu1, double mu2, TMatrix sigma);

    double doGaus(double x1, double x2, double mu1, double mu2, TMatrix sigma);

    vector<TMatrix> Initializeparams(vecPseudoJet particles,
                                     int k);

    void ComputeWeightsGaussian(vecPseudoJet particles,
                                vector<vector<double> >* Weights,
                                int k,
                                vecPseudoJet mGMMjets,
                                vector<TMatrix> mGMMjetsparams);

    vecPseudoJet UpdateJetsGaussian(vecPseudoJet particles,
                                    vector<vector<double> > Weights,
                                    int k,
                                    vector<TMatrix>* mGMMjetsparams);

    vecPseudoJet ClusterFuzzyGaussian(vecPseudoJet particles,
                                      vector<vector<double> >* Weights,
                                      vector<TMatrix>* mGMMjetsparamsout);

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


    void EventDisplay(vecPseudoJet particles,
                      vecPseudoJet CAjets,
                      vecPseudoJet tops,
                      vecPseudoJet mGMMjets,
                      vector<vector<double> > Weights,
                      int which,
                      vector<TMatrix> mGMMjetsparams,
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
