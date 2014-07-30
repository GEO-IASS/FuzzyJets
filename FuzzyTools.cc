#include <vector>
#include <map>
#include <numeric>
#include <functional>
#include <string>
#include <sstream>
#include <ostream>
#include <set>
#include <time.h>
#include <math.h>
#include <assert.h>

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

#include "FuzzyTools.h"
#include "myFastJetBase.h"

#include "TTree.h"
#include "TRandom3.h"
#include "TError.h"
#include "TVector3.h"
#include "TMath.h"
#include "TMatrix.h"
#include "TEllipse.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLine.h"

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
    max_iters = 100;
    clustering_mode = FuzzyTools::NOCLUSTERINGMODE;
    kernel_type = FuzzyTools::NOKERNEL;

    learn_weights = true;
    learn_shape = true;

    merge_dist = 0.01;
    min_weight = 0.0001;
    min_sigma = 0.01;

    directory_prefix = "";
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
FuzzyTools::InitWeights(vecPseudoJet const& particles,
                        int k){
    vector<vector<double> > out;
    for (unsigned int p_iter=0; p_iter<particles.size(); p_iter++){
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
                   TMatrix const& sigma){
    TMatrix summ_t(1,2); // (x-mu) tranpose
    TMatrix summ(2,1);  // (x-mu)
    summ_t(0,0) = x1-mu1;
    summ_t(0,1) = x2-mu2;
    summ(0,0) = x1-mu1;
    summ(1,0) = x2-mu2;
    TMatrix sigma_inverse(2,2);
    Double_t det;

    // check for singularity in Sigma
    sigma_inverse(0,0)=sigma(0,0);
    sigma_inverse(0,1)=sigma(0,1);
    sigma_inverse(1,0)=sigma(1,0);
    sigma_inverse(1,1)=sigma(1,1);
    if (sigma(0,0)*sigma(1,1)-sigma(1,0)*sigma(0,1) < 0.001*0.001){
        sigma_inverse(0,0)=0.01;
        sigma_inverse(1,1)=0.01;
    }
    sigma_inverse.Invert(&det);

    // compute the value of the pdf at (x1, x2)
    TMatrix hold = summ_t*sigma_inverse*summ;
    double exp_arg = -0.5*hold(0, 0);
    return exp(exp_arg)/(sqrt(fabs(det))*2*TMath::Pi());
}

// return the *square* of the Mahalanobis distance
double
FuzzyTools::MDist(double x1, double x2, double mu1, double mu2,
                  TMatrix const& sigma) {
    TMatrix summ_t(1,2); // (x-mu) tranpose
    TMatrix summ(2,1);  // (x-mu)
    summ_t(0,0) = x1-mu1;
    summ_t(0,1) = x2-mu2;
    summ(0,0) = x1-mu1;
    summ(1,0) = x2-mu2;
    TMatrix sigma_inverse(2,2);
    Double_t det;

    // check for singularity in Sigma
    sigma_inverse(0,0)=sigma(0,0);
    sigma_inverse(0,1)=sigma(0,1);
    sigma_inverse(1,0)=sigma(1,0);
    sigma_inverse(1,1)=sigma(1,1);
    if (sigma(0,0)*sigma(1,1)-sigma(1,0)*sigma(0,1) < 0.001*0.001){
        sigma_inverse(0,0)=0.01;
        sigma_inverse(1,1)=0.01;
    }
    sigma_inverse.Invert(&det);
    TMatrix hold = summ_t*sigma_inverse*summ;
    return hold(0, 0);
}

double
FuzzyTools::doTruncGaus(double x1, double x2, double mu1, double mu2,
                        TMatrix const& sigma) {
    double dist = MDist(x1, x2, mu1, mu2, sigma);
    if (dist > R*R) {
        return 0;
    }
    double scale = (1.0-exp(-R*R/2.0));
    return doGaus(x1, x2, mu1, mu2, sigma) / scale;
}

vector<TMatrix>
FuzzyTools::Initializeparams(__attribute__((unused)) vecPseudoJet const& particles,
                             int k){
    vector<TMatrix> out_params;
    for (int i=0; i<k;i++){
        TMatrix hold(2,2);
        hold(0,0) = 0.5;
        hold(1,1) = 0.5;
        hold(0,1) = 0.0;
        hold(1,0) = 0.0;
        out_params.push_back(hold);
    }
    return out_params;
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
FuzzyTools::ComputeWeightsGaussian(vecPseudoJet const& particles,
                                   vector<vector<double> >* weights,
                                   __attribute__((unused)) int k,
                                   vecPseudoJet const& mGMM_jets,
                                   vector<TMatrix> const& mGMM_jets_params,
                                   vector<double> const& mGMM_weights){
    for (unsigned int i=0; i<particles.size(); i++){
        double denom=0.;
        for (unsigned int j=0; j<mGMM_jets.size(); j++){
            denom+=doGaus(particles[i].rapidity(),
                          particles[i].phi(),
                          mGMM_jets[j].rapidity(),
                          mGMM_jets[j].phi(),
                          mGMM_jets_params[j])
                * mGMM_weights[j];
        }
        for (unsigned int j=0; j<mGMM_jets.size(); j++){
            weights->at(i)[j] = doGaus(particles[i].rapidity(),
                                       particles[i].phi(),
                                       mGMM_jets[j].rapidity(),
                                       mGMM_jets[j].phi(),
                                       mGMM_jets_params[j])
                * mGMM_weights[j] / denom;
        }
    }
}

void
FuzzyTools::ComputeWeightsTruncGaus(vecPseudoJet const& particles,
                                    vector<vector<double> >* weights,
                                    __attribute__((unused)) int k,
                                    vecPseudoJet const& mTGMM_jets,
                                    vector<TMatrix> const& mTGMM_jets_params,
                                    vector<double> const& mTGMM_weights) {
    for (unsigned int i=0; i < particles.size(); i++) {
        double denom = 0;
        for (unsigned int j = 0; j < mTGMM_jets.size(); j++) {
            denom += doTruncGaus(particles[i].rapidity(),
                                 particles[i].phi(),
                                 mTGMM_jets[j].rapidity(),
                                 mTGMM_jets[j].phi(),
                                 mTGMM_jets_params[j])
                * mTGMM_weights[j];
        }
        for (unsigned int j = 0; j < mTGMM_jets.size(); j++) {
            double new_weight = doTruncGaus(particles[i].rapidity(),
                                           particles[i].phi(),
                                           mTGMM_jets[j].rapidity(),
                                           mTGMM_jets[j].phi(),
                                           mTGMM_jets_params[j])
                * mTGMM_weights[j] / denom;
            if(new_weight < 0 || new_weight > 1 || isnan(new_weight)) {
                new_weight = 0.;
            }
            weights->at(i)[j] = new_weight;
        }
    }
}

void
FuzzyTools::ComputeWeightsUniform(vecPseudoJet const& particles,
                                  vector<vector<double> >* weights,
                                  __attribute__((unused)) int k,
                                  vecPseudoJet const& mUMM_jets,
                                  vector<double> const& mUMM_weights) {
    for (unsigned int i=0; i<particles.size(); i++) {
        double denom=0.;
        double dist;
        double t;
        for (unsigned int j=0; j < mUMM_jets.size(); j++) {
            dist = particles[i].delta_R(mUMM_jets[j]);
            if (dist <= R) {
                // pT to scale by cluster weight
                denom += mUMM_weights[j]/(TMath::Pi() * R*R);
            }
        }
        for (unsigned int j=0; j <  mUMM_jets.size(); j++) {
            dist = particles[i].delta_R(mUMM_jets[j]);
            t = (dist <= R) ? mUMM_weights[j] : 0;
            double new_weight = t/((TMath::Pi() * R*R) * denom);
            if(new_weight < 0 || new_weight > 1 || isnan(new_weight)) {
                new_weight = 0.;
            }
            weights->at(i)[j] = new_weight;
        }
    }
}

vecPseudoJet
FuzzyTools::UpdateJetsUniform(vecPseudoJet const& particles,
                              vector<vector<double> > const& weights,
                              int cluster_count,
                              vector<double>* mUMM_weights) {
    vecPseudoJet out_jets;

    double total_particle_pt = 0;
    unsigned int particle_count = weights.size();
    for (unsigned int particle_iter = 0; particle_iter < particle_count; particle_iter++) {
        total_particle_pt += pow(particles[particle_iter].pt(), alpha);
    }
    // iterate over clusters and update parameters and the covariance matrix
    for (int cluster_iter=0; cluster_iter<cluster_count; cluster_iter++){
        double jet_y=0;
        double jet_phi=0;
        double cluster_weighted_pt = 0;
        double cluster_pi = 0;

        // build the pT fraction belonging to cluster cluster_iter
        for (unsigned int particle_iter=0; particle_iter<particle_count; particle_iter++) {
            cluster_weighted_pt += pow(particles[particle_iter].pt(),alpha)*weights[particle_iter][cluster_iter];
        }

        // compute new cluster location on the basis of the EM update steps
        for (unsigned int particle_iter=0; particle_iter<particle_count; particle_iter++){
            jet_y+=pow(particles[particle_iter].pt(),alpha) * weights[particle_iter][cluster_iter]
                * particles[particle_iter].rapidity() / cluster_weighted_pt;
            jet_phi+=pow(particles[particle_iter].pt(),alpha) * weights[particle_iter][cluster_iter]
                * particles[particle_iter].phi() / cluster_weighted_pt;
        }
        cluster_pi = cluster_weighted_pt / total_particle_pt;
        if (learn_weights) {
            mUMM_weights->at (cluster_iter) = cluster_pi;
        }

        // compute new cluster weight on the basis of the EM update steps
        if (!(cluster_weighted_pt > 0)){
            jet_y = 0;
            jet_phi = 0;
        }

        fastjet::PseudoJet my_jet;
        my_jet.reset_PtYPhiM(1,jet_y,jet_phi,0.);
        out_jets.push_back(my_jet);
    }

    return out_jets;
}

vecPseudoJet
FuzzyTools::UpdateJetsTruncGaus(vecPseudoJet const& particles,
                                vector<vector<double> > const& weights,
                                int cluster_count,
                                vector<TMatrix>* mTGMM_jets_params,
                                vector<double>* mTGMM_weights) {
    vecPseudoJet out_jets;

    double total_particle_pt = 0;
    unsigned int particle_count = weights.size();
    for (unsigned int particle_iter = 0; particle_iter < particle_count; particle_iter++) {
        total_particle_pt += pow(particles[particle_iter].pt(), alpha);
    }
    // iterate over clusters and update parameters and the covariance matrix
    for (int cluster_iter=0; cluster_iter<cluster_count; cluster_iter++){
        double jet_y=0;
        double jet_phi=0;
        double cluster_weighted_pt = 0;
        double cluster_pi = 0;
        unsigned int particle_count = weights.size();

        // build the pT fraction belonging to cluster cluster_iter
        for (unsigned int particle_iter=0; particle_iter<particle_count; particle_iter++){
            cluster_weighted_pt += pow(particles[particle_iter].pt(),alpha)*weights[particle_iter][cluster_iter];
        }

        // compute new cluster location on the basis of the EM update steps
        for (unsigned int particle_iter=0; particle_iter<particle_count; particle_iter++){
            jet_y+=pow(particles[particle_iter].pt(),alpha) * weights[particle_iter][cluster_iter]
                * particles[particle_iter].rapidity() / cluster_weighted_pt;
            jet_phi+=pow(particles[particle_iter].pt(),alpha) * weights[particle_iter][cluster_iter]
                * particles[particle_iter].phi() / cluster_weighted_pt;
        }
        cluster_pi = cluster_weighted_pt / total_particle_pt;
        if (learn_weights) {
            mTGMM_weights->at(cluster_iter) = cluster_pi;
        }
        if (!(cluster_weighted_pt > 0)){
            jet_y = 0;
            jet_phi = 0;
        }

        fastjet::PseudoJet my_jet;
        my_jet.reset_PtYPhiM(1.,jet_y,jet_phi,0.);
        out_jets.push_back(my_jet);

        //now, we update sigma
        TMatrix sigma_update(2,2);
        for (unsigned int particle_iter=0; particle_iter<particle_count; particle_iter++){
            TMatrix hold(2,2);

            // pt scaled particle weight
            double q_ji = pow(particles[particle_iter].pt(),alpha) * weights[particle_iter][cluster_iter];
            hold(0,0) = q_ji
                * (particles[particle_iter].rapidity()-my_jet.rapidity())
                * (particles[particle_iter].rapidity()-my_jet.rapidity()) / cluster_weighted_pt;

            hold(0,1) = q_ji
                * (particles[particle_iter].rapidity()-my_jet.rapidity())
                * (particles[particle_iter].phi()-my_jet.phi()) / cluster_weighted_pt;

            hold(1,0) = q_ji
                * (particles[particle_iter].rapidity()-my_jet.rapidity())
                * (particles[particle_iter].phi()-my_jet.phi()) / cluster_weighted_pt;

            hold(1,1) = q_ji
                * (particles[particle_iter].phi()-my_jet.phi())
                * (particles[particle_iter].phi()-my_jet.phi()) / cluster_weighted_pt;

            sigma_update+=hold;
        }

        // if the matrix is looking singular...
        if (sigma_update(0,0)+sigma_update(1,1)+sigma_update(0,1) < 0.01){
            sigma_update(0,0)=0.1*0.1;
            sigma_update(1,1)=0.1*0.1;
        }

        // updated sigma is junk if it had almost no contained pT
        if (!(cluster_weighted_pt > 0)){
            sigma_update(0,0)=0.001*0.001;
            sigma_update(1,1)=0.001*0.001;
            sigma_update(0,1)=0.;
            sigma_update(1,0)=0.;
        }
        if (learn_shape) {
            mTGMM_jets_params->at(cluster_iter) = sigma_update;
        }

    }

    return out_jets;
}

vecPseudoJet
FuzzyTools::UpdateJetsGaussian(vecPseudoJet const& particles,
                               vector<vector<double> > const& weights,
                               int cluster_count,
                               vector<TMatrix>* mGMM_jets_params,
                               vector<double>* mGMM_weights){
    vecPseudoJet out_jets;

    double total_particle_pt = 0;
    unsigned int particle_count = weights.size();
    for (unsigned int particle_iter = 0; particle_iter < particle_count; particle_iter++) {
        total_particle_pt += pow(particles[particle_iter].pt(), alpha);
    }
    // iterate over clusters and update parameters and the covariance matrix
    for (int cluster_iter=0; cluster_iter<cluster_count; cluster_iter++){
        double jet_y=0;
        double jet_phi=0;
        double cluster_weighted_pt = 0;
        double cluster_pi;
        unsigned int particle_count = weights.size();

        // build the pT fraction belonging to cluster cluster_iter
        for (unsigned int particle_iter=0; particle_iter<particle_count; particle_iter++){
            cluster_weighted_pt += pow(particles[particle_iter].pt(),alpha)*weights[particle_iter][cluster_iter];
        }

        // compute new cluster location on the basis of the EM update steps
        for (unsigned int particle_iter=0; particle_iter<particle_count; particle_iter++){
            jet_y+=pow(particles[particle_iter].pt(),alpha) * weights[particle_iter][cluster_iter]
                * particles[particle_iter].rapidity() / cluster_weighted_pt;
            jet_phi+=pow(particles[particle_iter].pt(),alpha) * weights[particle_iter][cluster_iter]
                * particles[particle_iter].phi() / cluster_weighted_pt;
        }

        cluster_pi = cluster_weighted_pt / total_particle_pt;
        if (learn_weights) {
            mGMM_weights->at(cluster_iter) = cluster_pi;
        }
        if (!(cluster_weighted_pt > 0)){
            jet_y = 0;
            jet_phi = 0;
        }

        fastjet::PseudoJet my_jet;
        my_jet.reset_PtYPhiM(1.,jet_y,jet_phi,0.);
        out_jets.push_back(my_jet);

        //now, we update sigma
        TMatrix sigma_update(2,2);
        for (unsigned int particle_iter=0; particle_iter<particle_count; particle_iter++){
            TMatrix hold(2,2);

            // pt scaled particle weight
            double q_ji = pow(particles[particle_iter].pt(),alpha) * weights[particle_iter][cluster_iter];
            hold(0,0) = q_ji
                * (particles[particle_iter].rapidity()-my_jet.rapidity())
                * (particles[particle_iter].rapidity()-my_jet.rapidity()) / cluster_weighted_pt;

            hold(0,1) = q_ji
                * (particles[particle_iter].rapidity()-my_jet.rapidity())
                * (particles[particle_iter].phi()-my_jet.phi()) / cluster_weighted_pt;

            hold(1,0) = q_ji
                * (particles[particle_iter].rapidity()-my_jet.rapidity())
                * (particles[particle_iter].phi()-my_jet.phi()) / cluster_weighted_pt;

            hold(1,1) = q_ji
                * (particles[particle_iter].phi()-my_jet.phi())
                * (particles[particle_iter].phi()-my_jet.phi()) / cluster_weighted_pt;

            sigma_update+=hold;
        }

        // if the matrix is looking singular...
        if (sigma_update(0,0)+sigma_update(1,1)+sigma_update(0,1) < 0.01){
            sigma_update(0,0)=0.1*0.1;
            sigma_update(1,1)=0.1*0.1;
        }

        // updated sigma is junk if it had almost no contained pT
        if (!(cluster_weighted_pt > 0)){
            sigma_update(0,0)=0.001*0.001;
            sigma_update(1,1)=0.001*0.001;
            sigma_update(0,1)=0.;
            sigma_update(1,0)=0.;
        }
        if (learn_shape) {
            mGMM_jets_params->at(cluster_iter)=sigma_update;
        }
    }

    return out_jets;
}



set<unsigned int>
FuzzyTools::ClustersForRemovalGaussian(vecPseudoJet const& mGMM_jets,
                                       vector<TMatrix> const& mGMM_jets_params,
                                       vector<double> const& mGMM_weights) {
    set<unsigned int>removal_indices;
    // remove any jets which are candidates for mergers
    for (unsigned int j=0; j<mGMM_jets.size(); j++){
        if (removal_indices.count(j)) continue; // skip flagged indices
        for (unsigned int k=j+1; k<mGMM_jets.size(); k++){
            if (mGMM_jets[j].delta_R(mGMM_jets[k]) < merge_dist){
                if(mGMM_weights[k] <= mGMM_weights[j]) {
                    removal_indices.insert(k);
                } else {
                    removal_indices.insert(j);
                    continue;
                }
            }
        }
    }

    //Also remove jets that are too small if the size is learned
    double epsilon = min_sigma*min_sigma;
    for (unsigned int j=0; j<mGMM_jets.size(); j++){
        if (mGMM_jets_params[j](0,0) < epsilon || mGMM_jets_params[j](1,1) < epsilon){
            removal_indices.insert(j);
        }
    }

    for (unsigned int j=0; j < mGMM_jets.size(); j++) {
        if (mGMM_weights[j] < min_weight) {
            removal_indices.insert(j);
        }
    }
    if (kernel_type == FuzzyTools::UNIFORM) {
        cout << removal_indices.size() << endl;
    }
    return removal_indices;
}

set<unsigned int>
FuzzyTools::ClustersForRemovalUniform(vecPseudoJet const& mUMM_jets,
                                      vector<double> const& mUMM_weights) {
    set<unsigned int> removal_indices;
    // remove any jets which are candidates for mergers
    for (unsigned int j=0; j<mUMM_jets.size(); j++){
        if (removal_indices.count(j)) continue;
        for (unsigned int k=j+1; k<mUMM_jets.size(); k++){
            if (mUMM_jets[j].delta_R(mUMM_jets[k])<merge_dist){
                if(mUMM_weights[k] <= mUMM_weights[j]) {
                    removal_indices.insert(k);
                } else {
                    removal_indices.insert(j);
                    continue;
                }
            }
        }
    }

    //Also remove jets that are too small if the size is learned
    for (unsigned int j=0; j<mUMM_jets.size(); j++){
        if (mUMM_weights[j] < min_weight) {
            removal_indices.insert(j);
        }
    }
    return removal_indices;
}

vecPseudoJet
FuzzyTools::ClusterFuzzyUniform(vecPseudoJet const& particles,
                                vector<vector<double> >* weights_out,
                                vector<double>* mUMM_weights_out) {
    assert(kernel_type == FuzzyTools::UNIFORM);

    int cluster_count = seeds.size();
    vector<vector<double> > weights = InitWeights(particles, cluster_count);
    vecPseudoJet mUMM_jets = Initialize(particles, cluster_count, seeds);

    vector<double> mUMM_weights;
    for (int i = 0; i < cluster_count; i++) {
        mUMM_weights.push_back(1.0/cluster_count);
    }

    for(int iter = 0; iter < max_iters; iter++) {
        ComputeWeightsUniform(particles, &weights, cluster_count, mUMM_jets, mUMM_weights);
        mUMM_jets = UpdateJetsUniform(particles, weights, cluster_count, &mUMM_weights);
        if (clustering_mode == FuzzyTools::FIXED) continue;

        set<unsigned int>repeats = ClustersForRemovalUniform(mUMM_jets,
                                                             mUMM_weights);

        vector<vector<double> >weights_hold;
        vecPseudoJet mUMM_jets_hold;
        vector<double> mUMM_weights_hold;

        for (unsigned int q = 0; q < particles.size(); q++) {
            vector<double> hhh;
            weights_hold.push_back(hhh);

        }
        for (unsigned int j=0; j < mUMM_jets.size(); j++) {
            if (repeats.count(j) == 0) {
                for (unsigned int q=0; q < particles.size(); q++) {
                    weights_hold[q].push_back(weights[q][j]);
                }
                mUMM_jets_hold.push_back(mUMM_jets[j]);
                mUMM_weights_hold.push_back(mUMM_weights[j]);
            }
        }
        weights.clear();
        mUMM_jets.clear();
        mUMM_weights.clear();
        weights = weights_hold;
        mUMM_jets = mUMM_jets_hold;
        mUMM_weights = mUMM_weights_hold;
        cluster_count = mUMM_jets.size();

        // rescale the weights vector so that they sum to one
        double total_weight = std::accumulate(mUMM_weights.begin(), mUMM_weights.end(), 0.0);
        std::transform(mUMM_weights.begin(), mUMM_weights.end(), mUMM_weights.begin(),
                       std::bind2nd(std::multiplies<double>(), 1.0/total_weight));
    }
    weights_out->clear();
    for (unsigned int i=0; i<weights.size(); i++){
        weights_out->push_back(weights[i]);
    }
    mUMM_weights_out->clear();
    for (unsigned int i=0; i<mUMM_weights.size(); i++){
        mUMM_weights_out->push_back(mUMM_weights[i]);
    }
    return mUMM_jets;
}

vecPseudoJet
FuzzyTools::ClusterFuzzyTruncGaus(vecPseudoJet const& particles,
                                  vector<vector<double> >* weights_out,
                                  vector<TMatrix>* mTGMM_jets_params_out,
                                  vector<double>* mTGMM_weights_out){
    assert(kernel_type == FuzzyTools::TRUNCGAUSSIAN);

    int cluster_count = seeds.size();

    vector<double> mTGMM_weights;
    for (int i = 0; i < cluster_count; i++) {
        mTGMM_weights.push_back(1.0/cluster_count);
    }

    vector<vector<double> > weights = InitWeights(particles,cluster_count);
    vecPseudoJet mTGMM_jets = Initialize(particles,cluster_count,seeds);
    vector<TMatrix> mTGMM_jets_params = Initializeparams(particles,cluster_count);
    for (int iter=0; iter<max_iters; iter++){
        // EM algorithm update steps
        ComputeWeightsTruncGaus(particles,&weights,cluster_count,mTGMM_jets,mTGMM_jets_params, mTGMM_weights);
        mTGMM_jets = UpdateJetsTruncGaus(particles,weights,cluster_count,&mTGMM_jets_params,&mTGMM_weights);

        // do not flag clusters for deletion and remove them
        // if we are not in recombination mode
        if (clustering_mode == FuzzyTools::FIXED) continue;

        // determine which if any clusters should be removed
        set<unsigned int>repeats = ClustersForRemovalGaussian(mTGMM_jets,
                                                              mTGMM_jets_params,
                                                              mTGMM_weights);

        vector<vector<double> >weights_hold;
        vecPseudoJet mTGMM_jets_hold;
        vector<TMatrix> mTGMM_jets_params_hold;
        vector<double> mTGMM_weights_hold;

        for (unsigned int q=0; q<particles.size(); q++){
            vector<double> hhh;
            weights_hold.push_back(hhh);
        }
        for (unsigned int j=0; j<mTGMM_jets.size(); j++){
            if (repeats.count(j) == 0){
                for (unsigned int q=0; q<particles.size(); q++){
                    weights_hold[q].push_back(weights[q][j]);
                }
                mTGMM_jets_hold.push_back(mTGMM_jets[j]);
                mTGMM_jets_params_hold.push_back(mTGMM_jets_params[j]);
                mTGMM_weights_hold.push_back(mTGMM_weights[j]);
            }
        }

        // now replace and update weights and parameters vectors
        weights.clear();
        mTGMM_jets.clear();
        mTGMM_jets_params.clear();
        mTGMM_weights.clear();

        weights = weights_hold;
        mTGMM_jets = mTGMM_jets_hold;
        mTGMM_jets_params = mTGMM_jets_params_hold;
        mTGMM_weights = mTGMM_weights_hold;
        double total_weight = std::accumulate(mTGMM_weights.begin(), mTGMM_weights.end(), 0.0);
        std::transform(mTGMM_weights.begin(), mTGMM_weights.end(), mTGMM_weights.begin(),
                       std::bind2nd(std::multiplies<double>(), 1.0/total_weight));
        cluster_count = mTGMM_jets_hold.size();

    }
    weights_out->clear();
    for (unsigned int i=0; i<weights.size(); i++){
        weights_out->push_back(weights[i]);
    }
    mTGMM_jets_params_out->clear();
    mTGMM_weights_out->clear();
    for (unsigned int i=0; i<mTGMM_jets_params.size(); i++){
        mTGMM_jets_params_out->push_back(mTGMM_jets_params[i]);
        mTGMM_weights_out->push_back(mTGMM_weights[i]);
    }

    return mTGMM_jets;

}


vecPseudoJet
FuzzyTools::ClusterFuzzyGaussian(vecPseudoJet const& particles,
                                 vector<vector<double> >* weights_out,
                                 vector<TMatrix>* mGMM_jets_params_out,
                                 vector<double>* mGMM_weights_out){
    assert(kernel_type == FuzzyTools::GAUSSIAN);

    int cluster_count = seeds.size();

    vector<double> mGMM_weights;
    for (int i = 0; i < cluster_count; i++) {
        mGMM_weights.push_back(1.0/cluster_count);
    }

    vector<vector<double> > weights = InitWeights(particles,cluster_count);
    vecPseudoJet mGMM_jets = Initialize(particles,cluster_count,seeds);
    vector<TMatrix> mGMM_jets_params = Initializeparams(particles,cluster_count);
    for (int iter=0; iter<max_iters; iter++){
        // EM algorithm update steps
        ComputeWeightsGaussian(particles,&weights,cluster_count,mGMM_jets,mGMM_jets_params, mGMM_weights);
        mGMM_jets = UpdateJetsGaussian(particles,weights,cluster_count,&mGMM_jets_params,&mGMM_weights);

        // do not flag clusters for deletion and remove them
        // if we are not in recombination mode
        if (clustering_mode == FuzzyTools::FIXED) continue;

        // determine which if any clusters should be removed
        set<unsigned int>repeats = ClustersForRemovalGaussian(mGMM_jets,
                                                              mGMM_jets_params,
                                                              mGMM_weights);

        vector<vector<double> >weights_hold;
        vecPseudoJet mGMM_jets_hold;
        vector<TMatrix> mGMM_jets_params_hold;
        vector<double> mGMM_weights_hold;

        for (unsigned int q=0; q<particles.size(); q++){
            vector<double> hhh;
            weights_hold.push_back(hhh);
        }
        for (unsigned int j=0; j<mGMM_jets.size(); j++){
            if (repeats.count(j) == 0){
                for (unsigned int q=0; q<particles.size(); q++){
                    weights_hold[q].push_back(weights[q][j]);
                }
                mGMM_jets_hold.push_back(mGMM_jets[j]);
                mGMM_jets_params_hold.push_back(mGMM_jets_params[j]);
                mGMM_weights_hold.push_back(mGMM_weights[j]);
            }
        }

        // now replace and update weights and parameters vectors
        weights.clear();
        mGMM_jets.clear();
        mGMM_jets_params.clear();
        mGMM_weights.clear();

        weights = weights_hold;
        mGMM_jets = mGMM_jets_hold;
        mGMM_jets_params = mGMM_jets_params_hold;
        mGMM_weights = mGMM_weights_hold;
        double total_weight = std::accumulate(mGMM_weights.begin(), mGMM_weights.end(), 0.0);
        std::transform(mGMM_weights.begin(), mGMM_weights.end(), mGMM_weights.begin(),
                       std::bind2nd(std::multiplies<double>(), 1.0/total_weight));
        cluster_count = mGMM_jets_hold.size();

    }
    weights_out->clear();
    for (unsigned int i=0; i<weights.size(); i++){
        weights_out->push_back(weights[i]);
    }
    mGMM_jets_params_out->clear();
    mGMM_weights_out->clear();
    for (unsigned int i=0; i<mGMM_jets_params.size(); i++){
        mGMM_jets_params_out->push_back(mGMM_jets_params[i]);
        mGMM_weights_out->push_back(mGMM_weights[i]);
    }

    return mGMM_jets;

}

void
FuzzyTools::NewEventDisplay(vecPseudoJet const& particles,
                            __attribute__((unused)) vecPseudoJet const& ca_jets,
                            __attribute__((unused)) vecPseudoJet const& tops,
                            __attribute__((unused)) vecPseudoJet const& mGMM_jets,
                            __attribute__((unused)) vector<vector<double> > const& weights,
                            __attribute__((unused)) int which,
                            __attribute__((unused)) vector<TMatrix> const& mGMM_jets_params,
                            __attribute__((unused)) vector<double> const& mGMM_weights,
                            TString const& out) {
    double min_eta = -5;
    double max_eta = 5;
    TCanvas canv("NEVC" + out, "", 1200, 600);
    TH2F hist("NEVH" + out, "", 30, min_eta, max_eta, 28, 0, 7);
    TTree aux("NEVT" + out, "");

    double loc_eta;
    double loc_phi;
    double theta;
    double x_r;
    double y_r;
    double w;

    aux.Branch("loc_eta", &loc_eta, "loc_eta/D");
    aux.Branch("loc_phi", &loc_phi, "loc_phi/D");
    aux.Branch("theta", &theta, "theta/D");
    aux.Branch("x_r", &x_r, "x_r/D");
    aux.Branch("y_r", &y_r, "y_r/D");
    aux.Branch("w", &w, "w/D");

    double eta, phi, pT;
    for (unsigned int i = 0; i < particles.size(); i++) {
        eta = particles[i].eta();
        phi = particles[i].phi();
        pT = particles[i].pt();
        hist.Fill(eta, phi, pT);
    }

    vector<TEllipse> ellipses;
    for (unsigned int i=0; i < mGMM_jets_params.size(); i++) {
        double var_eta = mGMM_jets_params[i](0,0);
        double var_phi = mGMM_jets_params[i](1,1);
        double covar  = mGMM_jets_params[i](1,0);
        double temp_a  = 0.5*(var_eta + var_phi);
        double temp_b  = 0.5*sqrt((var_eta-var_phi)*(var_eta-var_phi) + 4*covar*covar);
        double lambda_eta = temp_a + temp_b;
        double lambda_phi = temp_a - temp_b;
        theta = 0;
        if(covar > 0) theta=atan((lambda_eta - var_eta)/covar);
        loc_eta = mGMM_jets[i].eta();
        loc_phi = mGMM_jets[i].phi();
        x_r = sqrt(lambda_eta);
        y_r = sqrt(lambda_phi);
        w = mGMM_weights[i];
        theta = theta * 180 / TMath::Pi();
        aux.Fill();
        TEllipse current_ellipse(loc_eta, loc_phi,
                                x_r, y_r,
                                0, 360, theta);
        current_ellipse.SetFillStyle(0);
        if(learn_weights) {
            current_ellipse.SetLineWidth(2);
            current_ellipse.SetLineWidth(w);
        } else {
            current_ellipse.SetLineWidth(2);
        }
        if(loc_eta < max_eta && loc_eta > min_eta)
            ellipses.push_back(current_ellipse);
    }

    hist.Write();
    aux.Write();
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
}

void
FuzzyTools::NewEventDisplayUniform(vecPseudoJet const& particles,
                                   __attribute__((unused)) vecPseudoJet const& ca_jets,
                                   __attribute__((unused)) vecPseudoJet const& tops,
                                   __attribute__((unused)) vecPseudoJet const& mUMM_jets,
                                   __attribute__((unused)) vector<vector<double> > const& weights,
                                   __attribute__((unused)) int which,
                                   __attribute__((unused)) vector<double> const& mUMM_weights,
                                   TString const& out) {
    double min_eta = -5;
    double max_eta = 5;
    TCanvas canv("NEVC" + out, "", 1200, 600);
    TTree aux("NEVT" + out, "");

    double loc_eta;
    double loc_phi;
    double theta;
    double x_r;
    double y_r;
    __attribute__((unused)) double w;

    aux.Branch("loc_eta", &loc_eta, "loc_eta/D");
    aux.Branch("loc_phi", &loc_phi, "loc_phi/D");
    aux.Branch("theta", &theta, "theta/D");
    aux.Branch("x_r", &x_r, "x_r/D");
    aux.Branch("y_r", &y_r, "y_r/D");
    aux.Branch("w", &w, "w/D");


    TH2F hist("NEVH" + out, "TEST", 30, min_eta, max_eta, 28, 0, 7);

    double eta, phi, pT;
    for (unsigned int i = 0; i < particles.size(); i++) {
        eta = particles[i].eta();
        phi = particles[i].phi();
        pT = particles[i].pt();
        hist.Fill(eta, phi, pT);
    }

    vector<TEllipse> ellipses;
    for (unsigned int i=0; i < mUMM_jets.size(); i++) {
        loc_eta = mUMM_jets[i].eta();
        loc_phi = mUMM_jets[i].phi();
        x_r = R;
        y_r = R;
        w = mUMM_weights[i];
        theta = 0;
        aux.Fill();
        TEllipse current_ellipse(loc_eta, loc_phi,
                                x_r, y_r,
                                0, 360, theta);
        current_ellipse.SetFillStyle(0);
        if(learn_weights) {
            current_ellipse.SetLineWidth(2);
            current_ellipse.SetLineWidth(mUMM_weights[i]);
        } else {
            current_ellipse.SetLineWidth(2);
        }
        if(loc_eta < max_eta && loc_eta > min_eta)
            ellipses.push_back(current_ellipse);
    }

    hist.Write();
    aux.Write();

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

    //canv.Print(directory_prefix + "NEWEventUniform"+out+".root");
}

void
FuzzyTools::EventDisplay(vecPseudoJet const& particles,
                         vecPseudoJet const& ca_jets,
                         vecPseudoJet const& tops,
                         vecPseudoJet const& mGMM_jets,
                         vector<vector<double> > const& weights,
                         int which,
                         vector<TMatrix> const& mGMM_jets_params,
                         TString out){

    gStyle->SetOptStat(0);
    //gROOT->Reset();
    //gROOT->SetStyle("ATLAS");
    //gROOT->ForceStyle();
    //gStyle->SetPadLeftMargin(0.15);
    //gStyle->SetPadTopMargin(0.15);

    TCanvas *c = new TCanvas("EventDisplayOld" + out,"",500,500);

    __attribute__((unused)) double max_pt=-1;

    const int n=particles.size();
    map<int, TGraph*> vars;
    for (int i=0; i<n; i++) {
        //std::cout << i << " " << std::endl;
        double x[1];
        double y[1];
        x[0] = particles[i].rapidity();
        y[0] = particles[i].phi();

        //std::cout << i << " " << x[0] << " " << y[0] << std::endl;
        //std::cout << which << " " << weights[i].size() << std::endl;
        //std::cout << weights[i][which] << std::endl;

        vars[i] = new TGraph (1, x, y);
        int mycolor = 19-floor(weights[i][which]*10);
        if (mycolor < 12) mycolor =1;
        vars[i]->SetMarkerColor(1);//mycolor);
    }

    // std::cout << "here1 ? " << std::endl;

    const int n2=ca_jets.size();
    double x2[n2];
    double y2[n2];
    for (int i=0; i<n2; i++) {
        x2[i] = ca_jets[i].rapidity();
        y2[i] = ca_jets[i].phi();
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

    // std::cout << "here2 ? " << std::endl;

    const int n4=mGMM_jets.size();
    double x4[n4];
    double y4[n4];
    for (int i=0; i<n4; i++) {
        x4[i] = mGMM_jets[i].rapidity();
        y4[i] = mGMM_jets[i].phi();
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
        vars[i]->SetMarkerSize(1);
        vars[i]->SetMarkerStyle(20);
        vars[i]->Draw("samep");
    }

    // std::cout << "here3 ? " << std::endl;

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

    // std::cout << "here4 ? " << std::endl;

    for (unsigned int i=0; i<mGMM_jets_params.size(); i++){
        double a = mGMM_jets_params[i](0,0);
        double b = mGMM_jets_params[i](1,1);
        double c = mGMM_jets_params[i](1,0);
        double lambda1 = 0.5*(a+b)+0.5*sqrt((a-b)*(a-b)+4.*c*c);
        double lambda2 = 0.5*(a+b)-0.5*sqrt((a-b)*(a-b)+4.*c*c);
        double theta = 0.;
        if (c>0) theta=atan((lambda1-a)/c);
        // std::cout << "yo " << i << " " << sqrt(lambda1) << " " << sqrt(lambda2) << " " << theta << std::endl;
        TEllipse *el4 = new TEllipse(x4[i],y4[i],sqrt(lambda1),sqrt(lambda2),0,360,theta*180/TMath::Pi());
        el4->SetFillStyle(0);
        el4->Draw("same");
    }

    // std::cout << "here5 ? " << std::endl;

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

    c->Print(directory_prefix + "Event"+out+".root");
    delete c;
}

void
FuzzyTools::Qjetmass(vecPseudoJet particles, vector<vector<double> > weights, int which, TString out){

    TH1F q_jet_mass("qjetmass" + out,"",100,0,250);
    TRandom3 rand = TRandom3(1);

    for (int j=0; j<10000; j++){

        fastjet::PseudoJet qmass;
        for (unsigned int i=0; i<particles.size(); i++){
            double my_throw = rand.Uniform(0,1);
            if (my_throw < weights[i][which]){
                qmass+=particles[i];
            }
        }
        q_jet_mass.Fill(qmass.m());
    }

    q_jet_mass.Scale(1./q_jet_mass.Integral());
    q_jet_mass.GetXaxis()->SetTitle("Single Jet Mass [GeV]");
    q_jet_mass.GetYaxis()->SetTitle("(1/N)dN/d(2.5 GeV)");
    q_jet_mass.GetXaxis()->SetTitleOffset(1.4);
    q_jet_mass.GetYaxis()->SetTitleOffset(1.4);
    q_jet_mass.Write();
}

void
FuzzyTools::JetContributionDisplay(vecPseudoJet particles,
                                   vector<vector<double> > weights,
                                   int which,
                                   int m_type,
                                   TString out) {
    double min_eta = -5;
    double max_eta = 5;
    unsigned int k = weights[0].size();
    TH2F hist("JCH_hard" + out, "", 30, min_eta, max_eta, 28, 0, 7);
    for (unsigned int i = 0; i < particles.size(); i++) {
        int m_idx = -1;
        double mweight = -1;
        for (unsigned int j = 0; j < k; j++) {
            if (weights[i][j] > mweight) {
                mweight = weights[i][j];
                m_idx = j;
            }
        }
        if (m_idx == which) {
            double fill_val = m_type == 1 ? particles[i].m() : particles[i].pt();
            hist.Fill(particles[i].eta(), particles[i].phi(), fill_val);
        }
    }
    hist.Write();

    TH2F hist2("JCH_soft" + out, "", 30, min_eta, max_eta, 28, 0, 7);
    for (unsigned int i=0; i < particles.size(); i++) {
        double fill_val = m_type == 1 ? particles[i].m() : particles[i].pt();
        fill_val *= weights[i][which];
        hist2.Fill(particles[i].eta(), particles[i].phi(), fill_val);
    }

    hist2.Write();
}

double totalMass(vecPseudoJet particles, vector<unsigned int> indices) {
    double mass = 0;
    for (unsigned int i = 0; i < indices.size(); i++) {
        mass += particles[indices[i]].m();
    }
    return mass;
}

double totalpT(vecPseudoJet particles, vector<unsigned int> indices) {
    double pT = 0;
    for (unsigned int i = 0; i < indices.size(); i++) {
        pT += particles[indices[i]].pt();
    }
    return pT;
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
FuzzyTools::CentralMoments(vecPseudoJet const& particles,
                           vector<vector<double> > const& weights,
                           unsigned int cluster_id,
                           unsigned int moment_count,
                           double (*f)(vecPseudoJet, vector<unsigned int>)) {
    unsigned int n_trials = 1000;

    vector<unsigned int> particle_indices;
    particle_indices.reserve(particles.size());

    vector<double> moments;
    moments.reserve(moment_count);

    TRandom3 rand = TRandom3(1);
    UInt_t seed = rand.GetSeed();

    // compute throws and the mean
    moments.push_back(0);
    for (unsigned int trial = 0; trial < n_trials; trial++) {
        for (unsigned int particle_iter=0; particle_iter<particles.size(); particle_iter++) {
            double my_throw = rand.Uniform(0, 1);
            if (my_throw < weights[particle_iter][cluster_id]) {
                particle_indices.push_back(particle_iter);
            }
        }
        moments[0] += (*f)(particles, particle_indices);
        particle_indices.clear();
    }

    moments[0] /= n_trials;
    rand.SetSeed(seed);
    particle_indices.clear();
    for (unsigned int i = 1; i < moment_count; i++) {
        moments.push_back(0);
    }

    for(unsigned int trial = 0; trial < n_trials; trial++) {
        for (unsigned int particle_iter = 0; particle_iter < particles.size(); particle_iter++) {
            double my_throw = rand.Uniform(0, 1);
            if (my_throw < weights[particle_iter][cluster_id]) {
                particle_indices.push_back(particle_iter);
            }
        }
        double v = (*f)(particles, particle_indices);

        for (unsigned int i = 1; i < moment_count; i++) {
            moments[i] += pow(v - moments[0], i+1);
        }
        particle_indices.clear();
    }
    for (unsigned int i = 1; i < moment_count; i++) {
        moments[i] /= n_trials;
    }
    return moments;
}

double
FuzzyTools::MLpT(vecPseudoJet particles,
                 vector<vector<double> > weights,
                 int jet_index,
                 int k,
                 int m_type){
    fastjet::PseudoJet my_jet;
    for (unsigned int i=0; i<particles.size(); i++){
        double my_max = -1;
        double which_jet = -1;
        for (int j=0; j<k; j++){
            if (weights[i][j]>my_max){
                my_max = weights[i][j];
                which_jet = j;
            }
        }
        if (which_jet==jet_index){
            my_jet+=particles[i];
        }
    }
    if (m_type==0) return my_jet.pt();
    return my_jet.m(); //m_type == 1
}

double
FuzzyTools::SoftpT(vecPseudoJet const& particles,
                   vector<vector<double> > const& weights,
                   int jet_index,
                   int m_type) {
    fastjet::PseudoJet my_jet;
    for (unsigned int i=0; i < particles.size(); i++) {
        my_jet += particles[i]*weights[i][jet_index];
    }
    if (m_type == 0) return my_jet.pt();
    return my_jet.m();
}

// has potentially very bad properties, might not conserve pT or mass!
double
FuzzyTools::MLlpTGaussian(vecPseudoJet const& particles,
                          fastjet::PseudoJet const& jet,
                          TMatrix const& jet_params,
                          double jetWeight, int m_type) {
    fastjet::PseudoJet my_jet;
    for(unsigned int i = 0; i < particles.size(); i++) {
        const fastjet::PseudoJet p = particles[i];
        double l = doGaus(p.eta(), p.phi(), jet.eta(), jet.phi(), jet_params);
        l /= doGaus(jet.eta(), jet.phi(), jet.eta(), jet.phi(), jet_params);
        my_jet += l * p;
    }
    if (m_type == 0) return my_jet.pt()*jetWeight;
    return my_jet.m() * jetWeight;
}

double
FuzzyTools::MLlpTUniform(vecPseudoJet const& particles,
                         fastjet::PseudoJet const& jet,
                         double jetWeight, int m_type) {
    fastjet::PseudoJet my_jet;
    for(unsigned int i = 0; i < particles.size(); i++) {
        const fastjet::PseudoJet p = particles[i];
        double l = p.delta_R(jet) <= R ? 1 : 0;
        my_jet += l * p;
    }
    if (m_type == 0) return my_jet.pt()*jetWeight;
    return my_jet.m() * jetWeight;
}

double
FuzzyTools::MLlpTTruncGaus(vecPseudoJet const& particles,
                           fastjet::PseudoJet const& jet,
                           TMatrix const& jet_params,
                           double jetWeight, int m_type) {
    fastjet::PseudoJet my_jet;
    for(unsigned int i = 0; i < particles.size(); i++) {
        const fastjet::PseudoJet p = particles[i];
        double l = doTruncGaus(p.eta(), p.phi(), p.eta(), p.phi(), jet_params);
        l /= doTruncGaus(jet.eta(), jet.phi(), jet.eta(), jet.phi(), jet_params);
        my_jet += l * p;
    }
    if (m_type == 0) return my_jet.pt()*jetWeight;
    return my_jet.m() * jetWeight;
}
