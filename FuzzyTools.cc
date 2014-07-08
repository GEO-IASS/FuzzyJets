#include <math.h>
#include <assert.h>
#include <vector>
#include <numeric>
#include <functional>
#include <string>
#include <sstream>
#include <ostream>
#include <set>
#include <time.h>

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

#include "FuzzyTools.h"
#include "myFastJetBase.h"

#include "TRandom3.h"
#include "TError.h"
#include "TVector3.h"
#include "TMath.h"
#include "TMatrix.h"
#include <map>
#include "TEllipse.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLine.h"

using namespace std;

// Makes debugging a little easier, pretty print vectors
template<typename T>
ostream& operator<< (ostream& out, const vector<T> v) {
    int last = v.size() - 1;
    if (last == -1) {
        out << "[EMPTYVEC]";
        return out;
    }
    out << "[";
    for(int i = 0; i < last; i++)
        out << v[i] << ", ";
    out << v[last] << "]";
    return out;
}


// Constructor
FuzzyTools::FuzzyTools(){
    alpha = 1.0;
    m_test = 0;
    maxiters = 100;
    clusteringMode = FuzzyTools::NOCLUSTERINGMODE;
    kernelType = FuzzyTools::NOKERNEL;

    learnWeights = true;
    learnShape = true;

    mergeDist = 0.01;
    minWeight = 0.0001;
    minSigma = 0.01;

    directoryPrefix = "";
}

vecPseudoJet
FuzzyTools::Initialize(__attribute__((unused)) vecPseudoJet particles,
                       int k,
                       vecPseudoJet jets){
    vecPseudoJet out;

    for (int i=0; i<k; i++){
        out.push_back(jets[i]);
    }
    return out;
}

vector<vector<double> >
FuzzyTools::InitWeights(vecPseudoJet particles,
                        int k){
    vector<vector<double> > out;
    for (unsigned int i=0; i<particles.size(); i++){
        vector<double> hold;
        hold.push_back(1.); //need to sum to 1!
        for (int j=1; j<k; j++){
            hold.push_back(0.);
        }
        out.push_back(hold);
    }
    return out;
}


//  compute N(x1, x2, Sigma, mu1, mu2)
//  the value of the Gaussian distribution with covariance Sigma
double
FuzzyTools::doGaus(double x1, double x2, double mu1, double mu2,
                   TMatrix sigma){
    TMatrix summT(1,2); // (x-mu) tranpose
    TMatrix summ(2,1);  // (x-mu)
    summT(0,0) = x1-mu1;
    summT(0,1) = x2-mu2;
    summ(0,0) = x1-mu1;
    summ(1,0) = x2-mu2;
    TMatrix sigmaInverse(2,2);
    Double_t det;

    // check for singularity in Sigma
    sigmaInverse(0,0)=sigma(0,0);
    sigmaInverse(0,1)=sigma(0,1);
    sigmaInverse(1,0)=sigma(1,0);
    sigmaInverse(1,1)=sigma(1,1);
    if (sigma(0,0)*sigma(1,1)-sigma(1,0)*sigma(0,1) < 0.001*0.001){
        sigmaInverse(0,0)=0.01;
        sigmaInverse(1,1)=0.01;
    }
    sigmaInverse.Invert(&det);

    // compute the value of the pdf at (x1, x2)
    TMatrix hold = summT*sigmaInverse*summ;
    double exparg = -0.5*hold(0, 0);
    return exp(exparg)/(sqrt(fabs(det))*2*TMath::Pi());
}

// return the *square* of the Mahalanobis distance
double
FuzzyTools::MDist(double x1, double x2, double mu1, double mu2,
                  TMatrix sigma) {
    TMatrix summT(1,2); // (x-mu) tranpose
    TMatrix summ(2,1);  // (x-mu)
    summT(0,0) = x1-mu1;
    summT(0,1) = x2-mu2;
    summ(0,0) = x1-mu1;
    summ(1,0) = x2-mu2;
    TMatrix sigmaInverse(2,2);
    Double_t det;

    // check for singularity in Sigma
    sigmaInverse(0,0)=sigma(0,0);
    sigmaInverse(0,1)=sigma(0,1);
    sigmaInverse(1,0)=sigma(1,0);
    sigmaInverse(1,1)=sigma(1,1);
    if (sigma(0,0)*sigma(1,1)-sigma(1,0)*sigma(0,1) < 0.001*0.001){
        sigmaInverse(0,0)=0.01;
        sigmaInverse(1,1)=0.01;
    }
    sigmaInverse.Invert(&det);
    TMatrix hold = summT*sigmaInverse*summ;
    return hold(0, 0);
}

double
FuzzyTools::doTruncGaus(double x1, double x2, double mu1, double mu2,
                        TMatrix sigma) {
    double dist = MDist(x1, x2, mu1, mu2, sigma);
    if (dist > R*R) {
        return 0;
    }
    double scale = (1-exp(-dist/2));
    return doGaus(x1, x2, mu1, mu2, sigma) / scale;
}

vector<TMatrix>
FuzzyTools::Initializeparams(__attribute__((unused)) vecPseudoJet particles,
                             int k){
    vector<TMatrix> outparams;
    for (int i=0; i<k;i++){
        TMatrix hold(2,2);
        hold(0,0) = 0.5;
        hold(1,1) = 0.5;
        hold(0,1) = 0.0;
        hold(1,0) = 0.0;
        outparams.push_back(hold);
    }
    return outparams;
}
// Compute new membership probabilities for each particle beloing to each
// cluster. This ammounts to computing the conditional probability that
// particle i belongs to jet j given a fixed set of cluster parameters.
// We divide by the total conditional probability (denom) in accordance with
// Bayes rule.

// This series of functions, which constitutes the expectation step of EM
// would be better suited to a functional style, and would be quite trivial
// with a higher order function. I don't want to litter the code with function
// pointers though. A way around this is to introduce a class
// ProbabilityDistribution so that subtype polymorphism can stand in for the
// higher function. Unfortunately, all this adds complexity.
void
FuzzyTools::ComputeWeightsGaussian(vecPseudoJet particles,
                                   vector<vector<double> >* Weights,
                                   __attribute__((unused)) int k,
                                   vecPseudoJet mGMMjets,
                                   vector<TMatrix> mGMMjetsparams,
                                   vector<double> mGMMweights){
    double li = 0;
    for (unsigned int i=0; i<particles.size(); i++){
        double denom=0.;
        for (unsigned int j=0; j<mGMMjets.size(); j++){
            denom+=doGaus(particles[i].rapidity(),
                          particles[i].phi(),
                          mGMMjets[j].rapidity(),
                          mGMMjets[j].phi(),
                          mGMMjetsparams[j])
                * mGMMweights[j];
        }
        li = 0;
        for (unsigned int j=0; j<mGMMjets.size(); j++){
            Weights->at(i)[j] = doGaus(particles[i].rapidity(),
                                       particles[i].phi(),
                                       mGMMjets[j].rapidity(),
                                       mGMMjets[j].phi(),
                                       mGMMjetsparams[j])
                * mGMMweights[j] / denom;
            li += Weights->at(i)[j];
        }

    }
}

void
FuzzyTools::ComputeWeightsTruncGaus(vecPseudoJet particles,
                                    vector<vector<double> >* Weights,
                                    __attribute__((unused)) int k,
                                    vecPseudoJet mTGMMjets,
                                    vector<TMatrix> mTGMMjetsparams,
                                    vector<double> mTGMMweights) {
    for (unsigned int i=0; i < particles.size(); i++) {
        double denom = 0;
        for (unsigned int j = 0; j < mTGMMjets.size(); j++) {
            denom += doTruncGaus(particles[i].rapidity(),
                                 particles[i].phi(),
                                 mTGMMjets[j].rapidity(),
                                 mTGMMjets[j].phi(),
                                 mTGMMjetsparams[j])
                * mTGMMweights[j];
        }
        for (unsigned int j = 0; j < mTGMMjets.size(); j++) {
            Weights->at(i)[j] = doTruncGaus(particles[i].rapidity(),
                                            particles[i].phi(),
                                            mTGMMjets[j].rapidity(),
                                            mTGMMjets[j].phi(),
                                            mTGMMjetsparams[j])
                * mTGMMweights[j] / denom;
        }
    }
}

void
FuzzyTools::ComputeWeightsUniform(vecPseudoJet particles,
                                  vector<vector<double> >* Weights,
                                  __attribute__((unused)) int k,
                                  vecPseudoJet mUMMjets,
                                  vector<double> mUMMweights) {
    for (unsigned int i=0; i<particles.size(); i++) {
        double denom=0.;
        double dist;
        double t;
        for (unsigned int j=0; j < mUMMjets.size(); j++) {
            dist = particles[i].delta_R(mUMMjets[j]);
            if (dist <= R) {
                // pT to scale by cluster weight
                denom += mUMMweights[j]/(TMath::Pi() * R*R);
            }
        }
        for (unsigned int j=0; j <  mUMMjets.size(); j++) {
            dist = particles[i].delta_R(mUMMjets[j]);
            t = (dist <= R) ? mUMMweights[j] : 0;
            Weights->at(i)[j] = t/((TMath::Pi() * R*R) * denom);
        }
    }
}

vecPseudoJet
FuzzyTools::UpdateJetsUniform(vecPseudoJet particles,
                              vector<vector<double> > Weights,
                              int clusterCount,
                              vector<double>* mUMMweights) {
    vecPseudoJet outjets;

    double totalParticlePt = 0;
    unsigned int particleCount = Weights.size();
    for (unsigned int particleIter = 0; particleIter < particleCount; particleIter++) {
        totalParticlePt += pow(particles[particleIter].pt(), alpha);
    }
    // iterate over clusters and update parameters and the covariance matrix
    for (int clusterIter=0; clusterIter<clusterCount; clusterIter++){
        double jety=0;
        double jetphi=0;
        double clusterWeightedPt = 0;
        double clusterPi = 0;

        // build the pT fraction belonging to cluster clusterIter
        for (unsigned int particleIter=0; particleIter<particleCount; particleIter++){
            clusterWeightedPt += pow(particles[particleIter].pt(),alpha)*Weights[particleIter][clusterIter];
        }

        // compute new cluster location on the basis of the EM update steps
        for (unsigned int particleIter=0; particleIter<particleCount; particleIter++){
            jety+=pow(particles[particleIter].pt(),alpha) * Weights[particleIter][clusterIter]
                * particles[particleIter].rapidity() / clusterWeightedPt;
            jetphi+=pow(particles[particleIter].pt(),alpha) * Weights[particleIter][clusterIter]
                * particles[particleIter].phi() / clusterWeightedPt;
        }
        clusterPi = clusterWeightedPt / totalParticlePt;
        if (learnWeights) {
            mUMMweights->at (clusterIter) = clusterPi;
        }

        // compute new cluster weight on the basis of the EM update steps
        if (!(clusterWeightedPt > 0)){
            jety = 0;
            jetphi = 0;
        }

        fastjet::PseudoJet myjet;
        myjet.reset_PtYPhiM(1,jety,jetphi,0.);
        outjets.push_back(myjet);
    }

    return outjets;
}

vecPseudoJet
FuzzyTools::UpdateJetsTruncGaus(vecPseudoJet particles,
                                vector<vector<double> > Weights,
                                int clusterCount,
                                vector<TMatrix>* mTGMMjetsparams,
                                vector<double>* mTGMMweights) {
    vecPseudoJet outjets;

    //For now, we fix the sigma at 1.
    //Means:

    double totalParticlePt = 0;
    unsigned int particleCount = Weights.size();
    for (unsigned int particleIter = 0; particleIter < particleCount; particleIter++) {
        totalParticlePt += pow(particles[particleIter].pt(), alpha);
    }
    // iterate over clusters and update parameters and the covariance matrix
    for (int clusterIter=0; clusterIter<clusterCount; clusterIter++){
        double jety=0;
        double jetphi=0;
        double clusterWeightedPt = 0;
        double clusterPi = 0;
        unsigned int particleCount = Weights.size();

        // build the pT fraction belonging to cluster clusterIter
        for (unsigned int particleIter=0; particleIter<particleCount; particleIter++){
            clusterWeightedPt += pow(particles[particleIter].pt(),alpha)*Weights[particleIter][clusterIter];
        }

        // compute new cluster location on the basis of the EM update steps
        for (unsigned int particleIter=0; particleIter<particleCount; particleIter++){
            jety+=pow(particles[particleIter].pt(),alpha) * Weights[particleIter][clusterIter]
                * particles[particleIter].rapidity() / clusterWeightedPt;
            jetphi+=pow(particles[particleIter].pt(),alpha) * Weights[particleIter][clusterIter]
                * particles[particleIter].phi() / clusterWeightedPt;
        }
        clusterPi = clusterWeightedPt / totalParticlePt;
        if (learnWeights) {
            mTGMMweights->at(clusterIter) = clusterPi;
        }
        if (!(clusterWeightedPt > 0)){
            jety = 0;
            jetphi = 0;
        }

        fastjet::PseudoJet myjet;
        myjet.reset_PtYPhiM(1.,jety,jetphi,0.);
        outjets.push_back(myjet);

        //now, we update sigma
        TMatrix sigmaupdate(2,2);
        for (unsigned int particleIter=0; particleIter<particleCount; particleIter++){
            TMatrix hold(2,2);

            // pt scaled particle weight
            double q_ji = pow(particles[particleIter].pt(),alpha) * Weights[particleIter][clusterIter];
            hold(0,0) = q_ji
                * (particles[particleIter].rapidity()-myjet.rapidity())
                * (particles[particleIter].rapidity()-myjet.rapidity()) / clusterWeightedPt;

            hold(0,1) = q_ji
                * (particles[particleIter].rapidity()-myjet.rapidity())
                * (particles[particleIter].phi()-myjet.phi()) / clusterWeightedPt;

            hold(1,0) = q_ji
                * (particles[particleIter].rapidity()-myjet.rapidity())
                * (particles[particleIter].phi()-myjet.phi()) / clusterWeightedPt;

            hold(1,1) = q_ji
                * (particles[particleIter].phi()-myjet.phi())
                * (particles[particleIter].phi()-myjet.phi()) / clusterWeightedPt;

            sigmaupdate+=hold;
        }

        // if the matrix is looking singular...
        if (sigmaupdate(0,0)+sigmaupdate(1,1)+sigmaupdate(0,1) < 0.01){
            sigmaupdate(0,0)=0.1*0.1;
            sigmaupdate(1,1)=0.1*0.1;
        }

        // updated sigma is junk if it had almost no contained pT
        if (!(clusterWeightedPt > 0)){
            sigmaupdate(0,0)=0.001*0.001;
            sigmaupdate(1,1)=0.001*0.001;
            sigmaupdate(0,1)=0.;
            sigmaupdate(1,0)=0.;
        }

        mTGMMjetsparams->at(clusterIter)=sigmaupdate;

    }

    return outjets;
}

vecPseudoJet
FuzzyTools::UpdateJetsGaussian(vecPseudoJet particles,
                               vector<vector<double> > Weights,
                               int clusterCount,
                               vector<TMatrix>* mGMMjetsparams,
                               vector<double>* mGMMweights){
    vecPseudoJet outjets;

    double totalParticlePt = 0;
    unsigned int particleCount = Weights.size();
    for (unsigned int particleIter = 0; particleIter < particleCount; particleIter++) {
        totalParticlePt += pow(particles[particleIter].pt(), alpha);
    }
    // iterate over clusters and update parameters and the covariance matrix
    for (int clusterIter=0; clusterIter<clusterCount; clusterIter++){
        double jety=0;
        double jetphi=0;
        double clusterWeightedPt = 0;
        double clusterPi;
        unsigned int particleCount = Weights.size();

        // build the pT fraction belonging to cluster clusterIter
        for (unsigned int particleIter=0; particleIter<particleCount; particleIter++){
            clusterWeightedPt += pow(particles[particleIter].pt(),alpha)*Weights[particleIter][clusterIter];
        }

        // compute new cluster location on the basis of the EM update steps
        for (unsigned int particleIter=0; particleIter<particleCount; particleIter++){
            jety+=pow(particles[particleIter].pt(),alpha) * Weights[particleIter][clusterIter]
                * particles[particleIter].rapidity() / clusterWeightedPt;
            jetphi+=pow(particles[particleIter].pt(),alpha) * Weights[particleIter][clusterIter]
                * particles[particleIter].phi() / clusterWeightedPt;
        }

        clusterPi = clusterWeightedPt / totalParticlePt;
        if (learnWeights) {
            mGMMweights->at(clusterIter) = clusterPi;
        }
        if (!(clusterWeightedPt > 0)){
            jety = 0;
            jetphi = 0;
        }

        fastjet::PseudoJet myjet;
        myjet.reset_PtYPhiM(1.,jety,jetphi,0.);
        outjets.push_back(myjet);

        //now, we update sigma
        TMatrix sigmaupdate(2,2);
        for (unsigned int particleIter=0; particleIter<particleCount; particleIter++){
            TMatrix hold(2,2);

            // pt scaled particle weight
            double q_ji = pow(particles[particleIter].pt(),alpha) * Weights[particleIter][clusterIter];
            hold(0,0) = q_ji
                * (particles[particleIter].rapidity()-myjet.rapidity())
                * (particles[particleIter].rapidity()-myjet.rapidity()) / clusterWeightedPt;

            hold(0,1) = q_ji
                * (particles[particleIter].rapidity()-myjet.rapidity())
                * (particles[particleIter].phi()-myjet.phi()) / clusterWeightedPt;

            hold(1,0) = q_ji
                * (particles[particleIter].rapidity()-myjet.rapidity())
                * (particles[particleIter].phi()-myjet.phi()) / clusterWeightedPt;

            hold(1,1) = q_ji
                * (particles[particleIter].phi()-myjet.phi())
                * (particles[particleIter].phi()-myjet.phi()) / clusterWeightedPt;

            sigmaupdate+=hold;
        }

        // if the matrix is looking singular...
        if (sigmaupdate(0,0)+sigmaupdate(1,1)+sigmaupdate(0,1) < 0.01){
            sigmaupdate(0,0)=0.1*0.1;
            sigmaupdate(1,1)=0.1*0.1;
        }

        // updated sigma is junk if it had almost no contained pT
        if (!(clusterWeightedPt > 0)){
            sigmaupdate(0,0)=0.001*0.001;
            sigmaupdate(1,1)=0.001*0.001;
            sigmaupdate(0,1)=0.;
            sigmaupdate(1,0)=0.;
        }
        if (learnShape) {
            mGMMjetsparams->at(clusterIter)=sigmaupdate;
        }
    }

    return outjets;
}



set<unsigned int>
FuzzyTools::ClustersForRemovalGaussian(vecPseudoJet& mGMMjets,
                                       vector<TMatrix>& mGMMjetsparams,
                                       vector<double>& mGMMweights) {
    set<unsigned int>removalIndices;
    // remove any jets which are candidates for mergers
    for (unsigned int j=0; j<mGMMjets.size(); j++){
        if (removalIndices.count(j)) continue; // skip flagged indices
        for (unsigned int k=j+1; k<mGMMjets.size(); k++){
            if (mGMMjets[j].delta_R(mGMMjets[k]) < mergeDist){
                if(mGMMweights[k] < mGMMweights[j]) {
                    removalIndices.insert(k);
                } else {
                    removalIndices.insert(j);
                    continue;
                }
            }
        }
    }

    //Also remove jets that are too small if the size is learned
    double epsilon = minSigma*minSigma;
    for (unsigned int j=0; j<mGMMjets.size(); j++){
        if (mGMMjetsparams[j](0,0) < epsilon || mGMMjetsparams[j](1,1) < epsilon){
            removalIndices.insert(j);
        }
    }

    for (unsigned int j=0; j < mGMMjets.size(); j++) {
        if (mGMMweights[j] < minWeight) {
            removalIndices.insert(j);
        }
    }
    cout << removalIndices.size() << endl;
    return removalIndices;
}

set<unsigned int>
FuzzyTools::ClustersForRemovalUniform(vecPseudoJet& mUMMjets,
                                      vector<double>& mUMMweights) {
    set<unsigned int>removalIndices;
    // remove any jets which are candidates for mergers
    for (unsigned int j=0; j<mUMMjets.size(); j++){
        if (removalIndices.count(j)) continue;
        for (unsigned int k=j+1; k<mUMMjets.size(); k++){
            if (mUMMjets[j].delta_R(mUMMjets[k])<mergeDist){
                if(mUMMweights[k] < mUMMweights[j]) {
                    removalIndices.insert(k);
                } else {
                    removalIndices.insert(j);
                    continue;
                }
            }
        }
    }

    //Also remove jets that are too small if the size is learned
    for (unsigned int j=0; j<mUMMjets.size(); j++){
        if (mUMMweights[j] < minWeight) {
            removalIndices.insert(j);
        }
    }
    cout << removalIndices.size() << endl;
    return removalIndices;
}

vecPseudoJet
FuzzyTools::ClusterFuzzyUniform(vecPseudoJet particles,
                                vector<vector<double> >* Weightsout,
                                vector<double>* mUMMweightsout) {
    assert(kernelType == FuzzyTools::UNIFORM);

    int clusterCount = seeds.size();

    vector<vector<double> > Weights = InitWeights(particles, clusterCount);
    vecPseudoJet mUMMjets = Initialize(particles, clusterCount, seeds);

    vector<double> mUMMweights;
    for (int i = 0; i < clusterCount; i++) {
        mUMMweights.push_back(1.0/clusterCount);
    }

    for(int iter = 0; iter < maxiters; iter++) {
        ComputeWeightsUniform(particles, &Weights, clusterCount, mUMMjets, mUMMweights);
        mUMMjets = UpdateJetsUniform(particles, Weights, clusterCount, &mUMMweights);

        if (clusteringMode == FuzzyTools::FIXED) continue;

        set<unsigned int>repeats = ClustersForRemovalUniform(mUMMjets,
                                                             mUMMweights);

        vector<vector<double> >Weights_hold;
        vecPseudoJet mUMMjets_hold;
        vector<double> mUMMweights_hold;

        for (unsigned int q = 0; q < particles.size(); q++) {
            vector<double> hhh;
            Weights_hold.push_back(hhh);

        }
        for (unsigned int j=0; j < mUMMjets.size(); j++) {
            if (repeats.count(j) == 0) {
                for (unsigned int q=0; q < particles.size(); q++) {
                    Weights_hold[q].push_back(Weights[q][j]);
                }
                mUMMjets_hold.push_back(mUMMjets[j]);
                mUMMweights_hold.push_back(mUMMweights[j]);
            }
        }
        Weights.clear();
        mUMMjets.clear();
        mUMMweights.clear();
        Weights = Weights_hold;
        mUMMjets = mUMMjets_hold;
        mUMMweights = mUMMweights_hold;
        clusterCount = mUMMjets.size();

        // rescale the weights vector so that they sum to one
        double totalWeight = std::accumulate(mUMMweights.begin(), mUMMweights.end(), 0.0);
        std::transform(mUMMweights.begin(), mUMMweights.end(), mUMMweights.begin(),
                       std::bind2nd(std::multiplies<double>(), 1.0/totalWeight));
    }
    Weightsout->clear();
    for (unsigned int i=0; i<Weights.size(); i++){
        Weightsout->push_back(Weights[i]);
    }
    mUMMweightsout->clear();
    for (unsigned int i=0; i<mUMMweights.size(); i++){
        mUMMweightsout->push_back(mUMMweights[i]);
    }

    return mUMMjets;
}

vecPseudoJet
FuzzyTools::ClusterFuzzyTruncGaus(vecPseudoJet particles,
                                  vector<vector<double> >* Weightsout,
                                  vector<TMatrix>* mTGMMjetsparamsout,
                                  vector<double>* mTGMMweightsout){
    assert(kernelType == FuzzyTools::TRUNCGAUSSIAN);

    int clusterCount = seeds.size();

    vector<double> mTGMMweights;
    for (int i = 0; i < clusterCount; i++) {
        mTGMMweights.push_back(1.0/clusterCount);
    }

    vector<vector<double> > Weights = InitWeights(particles,clusterCount);
    vecPseudoJet mTGMMjets = Initialize(particles,clusterCount,seeds);
    vector<TMatrix> mTGMMjetsparams = Initializeparams(particles,clusterCount);
    for (int iter=0; iter<maxiters; iter++){
        // EM algorithm update steps
        ComputeWeightsTruncGaus(particles,&Weights,clusterCount,mTGMMjets,mTGMMjetsparams, mTGMMweights);
        mTGMMjets = UpdateJetsTruncGaus(particles,Weights,clusterCount,&mTGMMjetsparams,&mTGMMweights);

        // do not flag clusters for deletion and remove them
        // if we are not in recombination mode
        if (clusteringMode == FuzzyTools::FIXED) continue;

        // determine which if any clusters should be removed
        set<unsigned int>repeats = ClustersForRemovalGaussian(mTGMMjets,
                                                              mTGMMjetsparams,
                                                              mTGMMweights);

        vector<vector<double> >Weights_hold;
        vecPseudoJet mTGMMjets_hold;
        vector<TMatrix> mTGMMjetsparams_hold;
        vector<double> mTGMMweights_hold;

        for (unsigned int q=0; q<particles.size(); q++){
            vector<double> hhh;
            Weights_hold.push_back(hhh);
        }
        for (unsigned int j=0; j<mTGMMjets.size(); j++){
            if (repeats.count(j) == 0){
                for (unsigned int q=0; q<particles.size(); q++){
                    Weights_hold[q].push_back(Weights[q][j]);
                }
                mTGMMjets_hold.push_back(mTGMMjets[j]);
                mTGMMjetsparams_hold.push_back(mTGMMjetsparams[j]);
                mTGMMweights_hold.push_back(mTGMMweights[j]);
            }
        }

        // now replace and update weights and parameters vectors
        Weights.clear();
        mTGMMjets.clear();
        mTGMMjetsparams.clear();
        mTGMMweights.clear();

        Weights = Weights_hold;
        mTGMMjets = mTGMMjets_hold;
        mTGMMjetsparams = mTGMMjetsparams_hold;
        mTGMMweights = mTGMMweights_hold;
        double totalWeight = std::accumulate(mTGMMweights.begin(), mTGMMweights.end(), 0.0);
        std::transform(mTGMMweights.begin(), mTGMMweights.end(), mTGMMweights.begin(),
                       std::bind2nd(std::multiplies<double>(), 1.0/totalWeight));
        clusterCount = mTGMMjets_hold.size();

    }
    Weightsout->clear();
    for (unsigned int i=0; i<Weights.size(); i++){
        Weightsout->push_back(Weights[i]);
    }
    mTGMMjetsparamsout->clear();
    mTGMMweightsout->clear();
    for (unsigned int i=0; i<mTGMMjetsparams.size(); i++){
        mTGMMjetsparamsout->push_back(mTGMMjetsparams[i]);
        mTGMMweightsout->push_back(mTGMMweights[i]);
    }

    return mTGMMjets;

}


vecPseudoJet
FuzzyTools::ClusterFuzzyGaussian(vecPseudoJet particles,
                                 vector<vector<double> >* Weightsout,
                                 vector<TMatrix>* mGMMjetsparamsout,
                                 vector<double>* mGMMweightsout){
    assert(kernelType == FuzzyTools::GAUSSIAN);

    int clusterCount = seeds.size();

    vector<double> mGMMweights;
    for (int i = 0; i < clusterCount; i++) {
        mGMMweights.push_back(1.0/clusterCount);
    }

    vector<vector<double> > Weights = InitWeights(particles,clusterCount);
    vecPseudoJet mGMMjets = Initialize(particles,clusterCount,seeds);
    vector<TMatrix> mGMMjetsparams = Initializeparams(particles,clusterCount);
    for (int iter=0; iter<maxiters; iter++){
        // EM algorithm update steps
        ComputeWeightsGaussian(particles,&Weights,clusterCount,mGMMjets,mGMMjetsparams, mGMMweights);
        mGMMjets = UpdateJetsGaussian(particles,Weights,clusterCount,&mGMMjetsparams,&mGMMweights);

        // do not flag clusters for deletion and remove them
        // if we are not in recombination mode
        if (clusteringMode == FuzzyTools::FIXED) continue;

        // determine which if any clusters should be removed
        set<unsigned int>repeats = ClustersForRemovalGaussian(mGMMjets,
                                                              mGMMjetsparams,
                                                              mGMMweights);

        vector<vector<double> >Weights_hold;
        vecPseudoJet mGMMjets_hold;
        vector<TMatrix> mGMMjetsparams_hold;
        vector<double> mGMMweights_hold;

        for (unsigned int q=0; q<particles.size(); q++){
            vector<double> hhh;
            Weights_hold.push_back(hhh);
        }
        for (unsigned int j=0; j<mGMMjets.size(); j++){
            if (repeats.count(j) == 0){
                for (unsigned int q=0; q<particles.size(); q++){
                    Weights_hold[q].push_back(Weights[q][j]);
                }
                mGMMjets_hold.push_back(mGMMjets[j]);
                mGMMjetsparams_hold.push_back(mGMMjetsparams[j]);
                mGMMweights_hold.push_back(mGMMweights[j]);
            }
        }

        // now replace and update weights and parameters vectors
        Weights.clear();
        mGMMjets.clear();
        mGMMjetsparams.clear();
        mGMMweights.clear();

        Weights = Weights_hold;
        mGMMjets = mGMMjets_hold;
        mGMMjetsparams = mGMMjetsparams_hold;
        mGMMweights = mGMMweights_hold;
        double totalWeight = std::accumulate(mGMMweights.begin(), mGMMweights.end(), 0.0);
        std::transform(mGMMweights.begin(), mGMMweights.end(), mGMMweights.begin(),
                       std::bind2nd(std::multiplies<double>(), 1.0/totalWeight));
        clusterCount = mGMMjets_hold.size();

    }
    Weightsout->clear();
    for (unsigned int i=0; i<Weights.size(); i++){
        Weightsout->push_back(Weights[i]);
    }
    mGMMjetsparamsout->clear();
    mGMMweightsout->clear();
    for (unsigned int i=0; i<mGMMjetsparams.size(); i++){
        mGMMjetsparamsout->push_back(mGMMjetsparams[i]);
        mGMMweightsout->push_back(mGMMweights[i]);
    }

    return mGMMjets;

}

void
FuzzyTools::NewEventDisplay(vecPseudoJet particles,
                            __attribute__((unused)) vecPseudoJet CAjets,
                            __attribute__((unused)) vecPseudoJet tops,
                            __attribute__((unused)) vecPseudoJet mGMMjets,
                            __attribute__((unused)) vector<vector<double> > Weights,
                            __attribute__((unused)) int which,
                            __attribute__((unused)) vector<TMatrix> mGMMjetsparams,
                            __attribute__((unused)) vector<double> mGMMweights,
                            TString out) {
    double mineta = -5;
    double maxeta = 5;
    TCanvas canv("", "", 1200, 600);
    TH2F hist("hist", "TEST", 30, mineta, maxeta, 28, 0, 7);

    double eta, phi, pT;
    for (unsigned int i = 0; i < particles.size(); i++) {
        eta = particles[i].eta();
        phi = particles[i].phi();
        pT = particles[i].pt();
        hist.Fill(eta, phi, pT);
    }

    vector<TEllipse> ellipses;
    for (unsigned int i=0; i < mGMMjetsparams.size(); i++) {
        double vareta = mGMMjetsparams[i](0,0);
        double varphi = mGMMjetsparams[i](1,1);
        double covar  = mGMMjetsparams[i](1,0);
        double tempa  = 0.5*(vareta + varphi);
        double tempb  = 0.5*sqrt((vareta-varphi)*(vareta-varphi) + 4*covar*covar);
        double lambdaeta = tempa + tempb;
        double lambdaphi = tempa - tempb;
        double theta = 0;
        if(covar > 0) theta=atan((lambdaeta - vareta)/covar);
        double loceta = mGMMjets[i].eta();
        double locphi = mGMMjets[i].phi();
        TEllipse currentEllipse(loceta, locphi,
                                sqrt(lambdaeta), sqrt(lambdaphi),
                                0, 360, theta*180/TMath::Pi());
        currentEllipse.SetFillStyle(0);
        if(learnWeights) {
            currentEllipse.SetLineWidth(2);
            currentEllipse.SetLineWidth(mGMMweights[i]);
        } else {
            currentEllipse.SetLineWidth(2);
        }
        if(loceta < maxeta && loceta > mineta)
            ellipses.push_back(currentEllipse);
    }


    canv.Divide(2, 1);

    canv.cd(1);
    hist.Draw("colz");
    for (unsigned int eli = 0; eli < ellipses.size(); eli++) {
        ellipses[eli].Draw();
    }

    canv.cd(2);
    hist.Draw("colz");
    for (unsigned int eli = 0; eli < ellipses.size(); eli++) {
        ellipses[eli].Draw();
    }
    gPad->SetLogz();
    canv.Update();


    canv.Print(directoryPrefix + "NEWEvent"+out+".root");
}

void
FuzzyTools::EventDisplay(vecPseudoJet particles,
                         vecPseudoJet CAjets,
                         vecPseudoJet tops,
                         vecPseudoJet mGMMjets,
                         vector<vector<double> > Weights,
                         int which,
                         vector<TMatrix> mGMMjetsparams,
                         TString out){

    gStyle->SetOptStat(0);
    //gROOT->Reset();
    //gROOT->SetStyle("ATLAS");
    //gROOT->ForceStyle();
    //gStyle->SetPadLeftMargin(0.15);
    //gStyle->SetPadTopMargin(0.15);

    TCanvas *c = new TCanvas("","",500,500);

    __attribute__((unused)) double maxpt=-1;

    const int n=particles.size();
    map<int, TGraph*> Vars;
    for (int i=0; i<n; i++) {
        //std::cout << i << " " << std::endl;
        double x[1];
        double y[1];
        x[0] = particles[i].rapidity();
        y[0] = particles[i].phi();

        //std::cout << i << " " << x[0] << " " << y[0] << std::endl;
        //std::cout << which << " " << Weights[i].size() << std::endl;
        //std::cout << Weights[i][which] << std::endl;

        Vars[i] = new TGraph (1, x, y);
        int mycolor = 19-floor(Weights[i][which]*10);
        if (mycolor < 12) mycolor =1;
        Vars[i]->SetMarkerColor(1);//mycolor);
    }

    std::cout << "here1 ? " << std::endl;

    const int n2=CAjets.size();
    double x2[n2];
    double y2[n2];
    for (int i=0; i<n2; i++) {
        x2[i] = CAjets[i].rapidity();
        y2[i] = CAjets[i].phi();
    }
    TGraph *gr2 = new TGraph (n2, x2, y2);

    const int n3=tops.size();
    double x3[n3];
    double y3[n3];
    for (int i=0; i<n3; i++) {
        x3[i] = tops[i].rapidity();
        y3[i] = tops[i].phi();
    }
    TGraph *gr3 = new TGraph (n3, x3, y3);

    std::cout << "here2 ? " << std::endl;

    const int n4=mGMMjets.size();
    double x4[n4];
    double y4[n4];
    for (int i=0; i<n4; i++) {
        x4[i] = mGMMjets[i].rapidity();
        y4[i] = mGMMjets[i].phi();
    }
    TGraph *gr4 = new TGraph (n4, x4, y4);

    TH1F * background = new TH1F("","",100,-4,8);
    background->GetXaxis()->SetTitle("Rapidity");
    background->GetYaxis()->SetTitleOffset(1.4);
    background->GetYaxis()->SetTitle("Azimuthal Angle [rad]");
    background->GetXaxis()->SetNdivisions(505);
    background->GetYaxis()->SetRangeUser(0,7);
    background->Draw();

    for (int i=0; i<n; i++) {
        Vars[i]->SetMarkerSize(1);
        Vars[i]->SetMarkerStyle(20);
        Vars[i]->Draw("samep");
    }

    std::cout << "here3 ? " << std::endl;

    gr2->SetMarkerSize(2);
    gr2->SetMarkerStyle(5);
    gr2->SetMarkerColor(2);
    gr2->Draw("psame");

    gr3->SetMarkerSize(3);
    gr3->SetMarkerStyle(29);
    gr3->SetMarkerColor(4);
    gr3->Draw("psame");

    gr4->SetMarkerSize(2);
    gr4->SetMarkerStyle(3);
    gr4->SetMarkerColor(6);
    gr4->Draw("psame");

    std::cout << "here4 ? " << std::endl;

    for (unsigned int i=0; i<mGMMjetsparams.size(); i++){
        double a = mGMMjetsparams[i](0,0);
        double b = mGMMjetsparams[i](1,1);
        double c = mGMMjetsparams[i](1,0);
        double lambda1 = 0.5*(a+b)+0.5*sqrt((a-b)*(a-b)+4.*c*c);
        double lambda2 = 0.5*(a+b)-0.5*sqrt((a-b)*(a-b)+4.*c*c);
        double theta = 0.;
        if (c>0) theta=atan((lambda1-a)/c);
        std::cout << "yo " << i << " " << sqrt(lambda1) << " " << sqrt(lambda2) << " " << theta << std::endl;
        TEllipse *el4 = new TEllipse(x4[i],y4[i],sqrt(lambda1),sqrt(lambda2),0,360,theta*180/TMath::Pi());
        el4->SetFillStyle(0);
        el4->Draw("same");
    }

    std::cout << "here5 ? " << std::endl;

    //myText(0.2,0.9,kBlack,"#scale[0.9]{#sqrt{s} = 8 TeV PYTHIA Z' #rightarrow t#bar{t}, m_{Z'}=1.5 TeV}");

    TLegend* leggaa = new TLegend(.7,.33,0.95,.67);
    leggaa->SetTextFont(42);
    //leggaa->AddEntry(gr1,"Particles","p");
    leggaa->AddEntry(gr3,"Top Quarks","p");
    leggaa->AddEntry(gr2,"C/A R=1.0 Jets","p");
    leggaa->AddEntry(gr4,"mGMM Jets","p");
    leggaa->SetFillStyle(0);
    leggaa->SetFillColor(0);
    leggaa->SetBorderSize(0);
    leggaa->Draw();

    c->Print(directoryPrefix + "Event"+out+".root");
}

void
FuzzyTools::Qjetmass(vecPseudoJet particles, vector<vector<double> > Weights, int which, TString out){

    TH1F* qjetmass = new TH1F("","",100,0,250);
    TRandom3 rand = TRandom3(1);

    for (int j=0; j<10000; j++){

        fastjet::PseudoJet qmass;
        for (unsigned int i=0; i<particles.size(); i++){
            double mythrow = rand.Uniform(0,1);
            if (mythrow < Weights[i][which]){
                qmass+=particles[i];
            }
        }
        qjetmass->Fill(qmass.m());
    }

    //gStyle->SetOptStat(0);
    //gROOT->Reset();
    //gROOT->SetStyle("ATLAS");
    //gROOT->ForceStyle();
    //gStyle->SetPadLeftMargin(0.15);
    //gStyle->SetPadTopMargin(0.15);

    TCanvas *c = new TCanvas("","",500,500);

    qjetmass->Scale(1./qjetmass->Integral());
    qjetmass->GetXaxis()->SetTitle("Single Jet Mass [GeV]");
    qjetmass->GetYaxis()->SetTitle("(1/N)dN/d(2.5 GeV)");
    qjetmass->GetXaxis()->SetTitleOffset(1.4);
    qjetmass->GetYaxis()->SetTitleOffset(1.4);
    qjetmass->Draw();
    //myText(0.2,0.9,kBlack,"#scale[0.9]{#sqrt{s} = 8 TeV PYTHIA Z' #rightarrow t#bar{t}, m_{Z'}=1.5 TeV}");
    c->Print(directoryPrefix + "Qjetmass"+out+".root");

}


// returns the vector of central moments for the distribution
// given by f(P) where P is the collection of particles,
// f is a function operating on all the particles,
// and the randomness in the distribution comes from the fact
// that particles are assigned fuzzily to the cluster
// note that instead of the first central moment (which is
// always zero) this function returns the mean instead.
//
// NOTE: At the moment, this function overwrites the value
// of the random number generator's seed so it can guarantee
// it will use the same samples when computing higher moments.
vector<double>
FuzzyTools::CentralMoments(vecPseudoJet particles,
                           vector<vector<double> >Weights,
                           unsigned int clusterID,
                           unsigned int momentCount,
                           double (*f)(vecPseudoJet, vector<unsigned int>)) {
    unsigned int nTrials = 10000;

    vector<unsigned int> particleIndices;
    particleIndices.reserve(particles.size());

    vector<double> moments;
    moments.reserve(momentCount);

    TRandom3 rand = TRandom3(1);
    UInt_t seed = rand.GetSeed();

    // compute throws and the mean
    moments.push_back(0);
    for (unsigned int trial = 0; trial < nTrials; trial++) {
        for (unsigned int particleIter=0; particleIter<particles.size(); particleIter++) {
            double mythrow = rand.Uniform(0, 1);
            if (mythrow < Weights[particleIter][clusterID]) {
                particleIndices.push_back(particleIter);
            }
        }
        moments[0] += (*f)(particles, particleIndices);
        particleIndices.clear();
    }

    moments[0] /= nTrials;
    rand.SetSeed(seed);
    particleIndices.clear();
    for (unsigned int i = 1; i < momentCount; i++) {
        moments.push_back(0);
    }

    for(unsigned int trial = 0; trial < nTrials; trial++) {
        for (unsigned int particleIter = 0; particleIter < particles.size(); particleIter++) {
            double mythrow = rand.Uniform(0, 1);
            if (mythrow < Weights[particleIter][clusterID]) {
                particleIndices.push_back(particleIter);
            }
        }
        double v = (*f)(particles, particleIndices);

        for (unsigned int i = 1; i < momentCount; i++) {
            moments[i] += pow(v - moments[0], i+1);
        }
        particleIndices.clear();
    }
    for (unsigned int i = 1; i < momentCount; i++) {
        moments[i] /= nTrials;
    }
    return moments;
}

double
FuzzyTools::MLpT(vecPseudoJet particles,
                 vector<vector<double> > Weights,
                 int jetindex,
                 int k,
                 int mtype){
    fastjet::PseudoJet myjet;
    for (unsigned int i=0; i<particles.size(); i++){
        double mymax = -1;
        double whichjet = -1;
        for (int j=0; j<k; j++){
            if (Weights[i][j]>mymax){
                mymax = Weights[i][j];
                whichjet = j;
            }
        }
        if (whichjet==jetindex){
            myjet+=particles[i];
        }
    }
    if (mtype==0) return myjet.pt();
    return myjet.m(); //mtype == 1
}
