#ifndef FUZZYTOOLS_H
#define FUZZYTOOLS_H

#include <vector>
#include <math.h>
#include <string>
#include "TMatrix.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

#include "Pythia8/Pythia.h"

#include "myFastJetBase.h"

using namespace std;

typedef std::vector<fastjet::PseudoJet> vecPseudoJet;


class FuzzyTools {
 private:
    int m_test;



 public:
    FuzzyTools();

    // methods
    vecPseudoJet Initialize(vecPseudoJet particles,
                            int k,
                            vecPseudoJet jets);

    vector<vector<double> > InitWeights(vecPseudoJet particles,int k);

    double doGaus(double x1, double x2, double mu1, double mu2, TMatrix sigma);

    vector<TMatrix> Initializeparams(vecPseudoJet particles,
                                     int k);

    void ComputeWeights(vecPseudoJet particles,
                        vector<vector<double> >* Weights,
                        int k,
                        vecPseudoJet mGMMjets,
                        vector<TMatrix> mGMMjetsparams);

    vecPseudoJet UpdateJets(vecPseudoJet particles,
                            vector<vector<double> > Weights,
                            int k,
                            vector<TMatrix>* mGMMjetsparams);

    vecPseudoJet ClusterFuzzy(vecPseudoJet particles,
                              vecPseudoJet inits,
                              vector<vector<double> >* Weights,
                              vector<TMatrix>* mGMMjetsparamsout);

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
