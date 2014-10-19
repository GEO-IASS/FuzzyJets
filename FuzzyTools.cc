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
#include "ROOTConf.h"

#ifdef WITHROOT
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
#endif

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
FuzzyTools::FuzzyTools()
    : default_sigma(MatTwo(0.5, 0, 0, 0.5)) {
    alpha = 1.0;
    m_test = 0;
    max_iters = 100;
    clustering_mode = FuzzyTools::NOCLUSTERINGMODE;
    kernel_type = FuzzyTools::NOKERNEL;

    learn_weights = true;
    learn_shape = true;

    log_log_likelihood_limit = -10;

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

#ifdef NEW
void
GaussianKernel::AdditionalUpdate(vecPseudoJet const& particles,
                                 FuzzyTools const& tool,
                                 double cluster_weighted_pt) {
    _sigma.xx = 0;
    _sigma.yy = 0;
    _sigma.xy = 0;
    _sigma.yx = 0;
    for (unsigned int particle_iter=0; particle_iter<particles.size(); particle_iter++){
        // pt scaled particle weight
        double q_ji = pow(particles[particle_iter].pt(), tool.GetAlpha()) * _weights[particle_iter];
        fastjet::PseudoJet const& particle = particles[particle_iter];
        double x_diff = particle.rapidity() - _mu_x;
        double y_diff = particle.rapidity() - _mu_y;
        _sigma.xx += q_ji * x_diff * x_diff / cluster_weighted_pt;
        _sigma.xy += q_ji * x_diff * y_diff / cluster_weighted_pt;
        _sigma.yx += q_ji * x_diff * y_diff / cluster_weighted_pt;
        _sigma.yy += q_ji * y_diff * y_diff / cluster_weighted_pt;
    }

    // if the matrix is looking singular...
    if (_sigma.xx+_sigma.yy+_sigma.xy < 0.01){
        _sigma.xx=0.5*0.5;
        _sigma.yy=0.5*0.5;
    }

    // updated sigma is junk if it had almost no contained pT
    if (!(cluster_weighted_pt > 0)){
        _sigma.xx=0.001*0.001;
        _sigma.yy=0.001*0.001;
        _sigma.xy=0.;
        _sigma.yx=0.;
    }
}

double
GaussianKernel::PDF(double x, double y) {
    const double limit = 0.001*0.001;
    double det = _sigma.determinant();

    double invdet;
    double ixx, iyx, iyy;

    if (fabs(det) < limit) {
        const double sigmasubsx = 0.01;
        const double sigmasubsy = 0.01;

        det = sigmasubsx * sigmasubsy - _sigma.xy * _sigma.yx;
        invdet = 1.0/det;

        ixx = sigmasubsy * invdet;
        iyy = sigmasubsx * invdet;
        iyx = -_sigma.yx * invdet;
    } else {
        invdet = 1.0/det;
        ixx = _sigma.yy * invdet;
        iyy = _sigma.xx * invdet;
        iyx = -_sigma.yx * invdet;
    }

    const double mxmu_x = x - _mu_x;
    const double mxmu_y = y - _mu_y;
    const double expval = exp(-0.5 * (ixx * mxmu_x * mxmu_x + iyy * mxmu_y * mxmu_y) - iyx * mxmu_x * mxmu_y);

    return expval / (2 * M_PI * sqrt(fabs(det)));
}
#endif

//  compute N(x1, x2, Sigma, mu1, mu2)
//  the value of the Gaussian distribution with covariance Sigma
double
FuzzyTools::doGaus(double x1, double x2, double mu1, double mu2,
                   MatTwo const& sigma){
    /*#ifdef FIXED_MOD_PI
    if (x2 - mu2 > M_PI) {
        mu2 += 2*M_PI
    }
    #endif*/

    const double limit = 0.001 * 0.001;
    double det = sigma.determinant();

    double invdet;
    double ixx, iyx, iyy;

    if (fabs(det) < limit) {
      const double sigmasubsx = 0.01;
      const double sigmasubsy = 0.01;
      det = sigmasubsx * sigmasubsy - sigma.xy * sigma.yx;
      invdet = 1.0/det;
      ixx = sigmasubsy * invdet;
      iyy = sigmasubsx * invdet;
      iyx = -sigma.yx * invdet;
    } else {
      invdet = 1.0/det;
      ixx = sigma.yy * invdet;
      iyy = sigma.xx * invdet;
      iyx = -sigma.yx * invdet;
    }

    const double mxmu1 = x1 - mu1;
    double mxmu2 = x2 - mu2;
    const double expvalA = exp(-0.5 * (ixx * mxmu1 * mxmu1 + iyy * mxmu2 * mxmu2) - iyx * mxmu1 * mxmu2);
#ifdef FIXED_MOD_PI
    mxmu2 = x2 - mu2 + 2 * M_PI * (x2 > mu2 ? -1 : 1);
    const double expvalB = exp(-0.5 * (ixx * mxmu1 * mxmu1 + iyy * mxmu2 * mxmu2) - iyx * mxmu1 * mxmu2);
    return max(expvalA, expvalB) / (2 * M_PI * sqrt(fabs(det)));
#else
    return expvalA / (2 * M_PI * sqrt(fabs(det)));
#endif
}

// return the *square* of the Mahalanobis distance
double
FuzzyTools::MDist(double x1, double x2, double mu1, double mu2,
                  MatTwo const& sigma) {
    const double limit = 0.001 * 0.001;
    double det = sigma.determinant();

    double invdet;
    double ixx, iyx, iyy;

    if (fabs(det) < limit) {
      const double sigmasubsx = 0.01;
      const double sigmasubsy = 0.01;
      det = sigmasubsx * sigmasubsy - sigma.xy * sigma.yx;
      invdet = 1.0/det;
      ixx = sigmasubsy * invdet;
      iyy = sigmasubsx * invdet;
      iyx = -sigma.yx * invdet;
    } else {
      invdet = 1.0/det;
      ixx = sigma.yy * invdet;
      iyy = sigma.xx * invdet;
      iyx = -sigma.yx * invdet;
    }

    const double mxmu1 = x1 - mu1;
    double mxmu2 = x2 - mu2;
    const double retA = (ixx * mxmu1 * mxmu1 + iyy * mxmu2 * mxmu2 + 2 * iyx * mxmu1 * mxmu2);
#ifdef FIXED_MOD_PI
    mxmu2 = x2 - mu2 + 2 * M_PI * (x2 > mu2 ? -1 : 1);
    const double retB = (ixx * mxmu1 * mxmu1 + iyy * mxmu2 * mxmu2 + 2 * iyx * mxmu1 * mxmu2);
    return min(retA, retB);
#else
    return retA;
#endif
}

double
FuzzyTools::doTruncGaus(double x1, double x2, double mu1, double mu2,
                        MatTwo const& sigma) {
    double dist = MDist(x1, x2, mu1, mu2, sigma);

    if (dist > R*R) {
        return 0;
    }
    double scale = (1.0-exp(-R*R/2.0));
    return doGaus(x1, x2, mu1, mu2, sigma) / scale;
}

vector<MatTwo>
FuzzyTools::Initializeparams(__attribute__((unused)) vecPseudoJet const& particles,
                             int k){
    vector<MatTwo> out_params;

    for (int i=0; i<k;i++){
        out_params.push_back(default_sigma);
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
                                   __attribute__((unused)) int cluster_count,
                                   vecPseudoJet const& mGMM_jets,
                                   vector<MatTwo> const& mGMM_jets_params,
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
            weights->at(i).at(j) = doGaus(particles[i].rapidity(),
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
                                    __attribute__((unused)) int cluster_count,
                                    vecPseudoJet const& mTGMM_jets,
                                    vector<MatTwo> const& mTGMM_jets_params,
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
            if(new_weight < 0 || new_weight > 1 || std::isnan(new_weight)) {
                new_weight = 0.;
            }
            weights->at(i).at(j) = new_weight;
        }
    }
}

void
FuzzyTools::ComputeWeightsUniform(vecPseudoJet const& particles,
                                  vector<vector<double> >* weights,
                                  __attribute__((unused)) int cluster_count,
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
                denom += mUMM_weights[j]/(M_PI * R*R);
            }
        }
        for (unsigned int j=0; j <  mUMM_jets.size(); j++) {
            dist = particles[i].delta_R(mUMM_jets[j]);
            t = (dist <= R) ? mUMM_weights[j] : 0;
            double new_weight = t/((M_PI * R*R) * denom);
            if(new_weight < 0 || new_weight > 1 || std::isnan(new_weight)) {
                new_weight = 0.;
            }
            weights->at(i).at(j) = new_weight;
        }
    }
}

vecPseudoJet
FuzzyTools::UpdateJetsUniform(vecPseudoJet const& particles,
                              __attribute__((unused)) vecPseudoJet const& oldJets,
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
            double delta_jet_y = pow(particles[particle_iter].pt(),alpha)
                * weights[particle_iter][cluster_iter] / cluster_weighted_pt;
            double delta_jet_phi = delta_jet_y;
            delta_jet_y *= particles[particle_iter].rapidity();

#ifdef FIXED_MOD_PI
            double k = particles[particle_iter].phi();
            if (fabs(oldJets[cluster_iter].phi() - k) > M_PI) {
                k += 2*M_PI*(oldJets[cluster_iter].phi() > k ? 1 : -1);
            }
            delta_jet_phi *= k;
#else
            delta_jet_phi *= particles[particle_iter].phi();
#endif

            jet_y += delta_jet_y;
            jet_phi += delta_jet_phi;
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
                                __attribute__((unused))
                                vecPseudoJet const& oldJets,
                                vector<vector<double> > const& weights,
                                int cluster_count,
                                vector<MatTwo>* mTGMM_jets_params,
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
            double delta_jet_y = pow(particles[particle_iter].pt(),alpha)
                * weights[particle_iter][cluster_iter] / cluster_weighted_pt;
            double delta_jet_phi = delta_jet_y;
            delta_jet_y *= particles[particle_iter].rapidity();

#ifdef FIXED_MOD_PI
            double k = particles[particle_iter].phi();
            if (fabs(oldJets[cluster_iter].phi() - k) > M_PI) {
                k += 2*M_PI*(oldJets[cluster_iter].phi() > k ? 1 : -1);
            }
            delta_jet_phi *= k;
#else
            delta_jet_phi *= particles[particle_iter].phi();
#endif

            jet_y += delta_jet_y;
            jet_phi += delta_jet_phi;
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

        if (learn_shape) {
            //now, we update sigma
            MatTwo sigma_update(0, 0, 0, 0);
            for (unsigned int particle_iter=0; particle_iter<particle_count; particle_iter++){
                // pt scaled particle weight
                double q_ji = pow(particles[particle_iter].pt(),alpha) * weights[particle_iter][cluster_iter];
                double eff_delta_phi = particles[particle_iter].phi()
                    -my_jet.phi();
#ifdef FIXED_MOD_PI
                if (fabs(eff_delta_phi) > M_PI) {
                    eff_delta_phi += 2*M_PI*(eff_delta_phi > 0 ? -1 : 1);
                }
#endif
                sigma_update.xx += q_ji
                    * (particles[particle_iter].rapidity()-my_jet.rapidity())
                    * (particles[particle_iter].rapidity()-my_jet.rapidity()) / cluster_weighted_pt;

                sigma_update.xy += q_ji
                    * (particles[particle_iter].rapidity()-my_jet.rapidity())
                    * (eff_delta_phi) / cluster_weighted_pt;

                sigma_update.yx += q_ji
                    * (particles[particle_iter].rapidity()-my_jet.rapidity())
                    * (eff_delta_phi) / cluster_weighted_pt;

                sigma_update.yy += q_ji
                    * (eff_delta_phi)
                    * (eff_delta_phi) / cluster_weighted_pt;
            }

            // if the matrix is looking singular...
            if (sigma_update.xx+sigma_update.yy+sigma_update.xy < 0.01){
                sigma_update.xx=0.5*0.5;
                sigma_update.yy=0.5*0.5;
            }

             // updated sigma is junk if it had almost no contained pT
            if (!(cluster_weighted_pt > 0)){
                sigma_update.xx=0.001*0.001;
                sigma_update.yy=0.001*0.001;
                sigma_update.xy=0.;
                sigma_update.yx=0.;
            }
            mTGMM_jets_params->at(cluster_iter) = sigma_update;
        }
    }

    return out_jets;
}

vecPseudoJet
FuzzyTools::UpdateJetsGaussianC(vecPseudoJet const& particles,
                                __attribute__((unused))
                                vecPseudoJet const& oldJets,
                                vector<vector<double> > const& weights,
                                int cluster_count,
                                vector<MatTwo>* mGMMc_jets_params,
                                vector<double>* mGMMc_weights) {
    vecPseudoJet out_jets;

    double total_particle_pt = 0;
    unsigned int particle_count = weights.size();
    for (unsigned int particle_iter = 0; particle_iter < particle_count; particle_iter++) {
        total_particle_pt += pow(particles[particle_iter].pt(), alpha);
    }

    for (int cluster_iter = 0; cluster_iter < cluster_count; cluster_iter++) {
        double jet_y = 0;
        double jet_phi = 0;
        double cluster_weighted_pt = 0;
        double cluster_pi = 0;

        for (unsigned int particle_iter = 0; particle_iter < particle_count; particle_iter++) {
            cluster_weighted_pt += pow(particles[particle_iter].pt(), alpha) * weights[particle_iter][cluster_iter];
        }

        for (unsigned int particle_iter = 0; particle_iter < particle_count; particle_iter++) {
            double delta_jet_y = pow(particles[particle_iter].pt(),alpha)
                * weights[particle_iter][cluster_iter] / cluster_weighted_pt;
            double delta_jet_phi = delta_jet_y;
            delta_jet_y *= particles[particle_iter].rapidity();

#ifdef FIXED_MOD_PI
            double k = particles[particle_iter].phi();
            if (fabs(oldJets[cluster_iter].phi() - k) > M_PI) {
                k += 2*M_PI*(oldJets[cluster_iter].phi() > k ? 1 : -1);
            }
            delta_jet_phi *= k;
#else
            delta_jet_phi *= particles[particle_iter].phi();
#endif

            jet_y += delta_jet_y;
            jet_phi += delta_jet_phi;
        }

        cluster_pi = cluster_weighted_pt / total_particle_pt;
        if (learn_weights) {
            mGMMc_weights->at(cluster_iter) = cluster_pi;
        }

        if(!(cluster_weighted_pt > 0)) {
            // setup jet to be pruned
            jet_y = 0;
            jet_phi = 0;
        }

        fastjet::PseudoJet my_jet;
        my_jet.reset_PtYPhiM(1.0, jet_y, jet_phi, 0.);
        out_jets.push_back(my_jet);

        if(learn_shape) {
            MatTwo sigma_update(0, 0, 0, 0);
            double sum = 0;
            for (unsigned int particle_iter = 0; particle_iter < particle_count; particle_iter++) {
                double q_ji = pow(particles[particle_iter].pt(), alpha) * weights[particle_iter][cluster_iter]
                    / cluster_weighted_pt;
                double v_off_x = particles[particle_iter].rapidity() - my_jet.rapidity();
                double v_off_y = particles[particle_iter].phi() - my_jet.phi();
#ifdef FIXED_MOD_PI
                if (fabs(v_off_y) > M_PI) {
                    v_off_y += 2*M_PI*(v_off_y > 0 ? -1 : 1);
                }
#endif
                sum += q_ji * (v_off_x * v_off_x + v_off_y * v_off_y);
            }
            double little_sigma = sqrt(sum / 2);
            if (2 * little_sigma < 0.001) {
                little_sigma = 0.01;
            }
            if (!(cluster_weighted_pt > 0)) {
                little_sigma = 0.00001;
            }
            sigma_update.xx = little_sigma;
            sigma_update.yy = little_sigma;
            mGMMc_jets_params->at(cluster_iter) = sigma_update;
        }
    }

    return out_jets;
}

vecPseudoJet
FuzzyTools::UpdateJetsGaussian(vecPseudoJet const& particles,
                               __attribute__((unused))
                               vecPseudoJet const& oldJets,
                               vector<vector<double> > const& weights,
                               int cluster_count,
                               vector<MatTwo>* mGMM_jets_params,
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
            double delta_jet_y = pow(particles[particle_iter].pt(),alpha)
                * weights[particle_iter][cluster_iter] / cluster_weighted_pt;
            double delta_jet_phi = delta_jet_y;
            delta_jet_y *= particles[particle_iter].rapidity();

#ifdef FIXED_MOD_PI
            double k = particles[particle_iter].phi();
            if (fabs(oldJets[cluster_iter].phi() - k) > M_PI) {
                k += 2*M_PI*(oldJets[cluster_iter].phi() > k ? 1 : -1);
            }
            delta_jet_phi *= k;
#else
            delta_jet_phi *= particles[particle_iter].phi();
#endif

            jet_y += delta_jet_y;
            jet_phi += delta_jet_phi;
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

        if (learn_shape) {
            //now, we update sigma
            MatTwo sigma_update(0, 0, 0, 0);
            for (unsigned int particle_iter=0; particle_iter<particle_count; particle_iter++){
                // pt scaled particle weight
                double q_ji = pow(particles[particle_iter].pt(),alpha) * weights[particle_iter][cluster_iter];
                double eff_delta_phi = particles[particle_iter].phi()
                    -my_jet.phi();
#ifdef FIXED_MOD_PI
                if (fabs(eff_delta_phi) > M_PI) {
                    eff_delta_phi += 2*M_PI*(eff_delta_phi > 0 ? -1 : 1);
                }
#endif
                sigma_update.xx += q_ji
                    * (particles[particle_iter].rapidity()-my_jet.rapidity())
                    * (particles[particle_iter].rapidity()-my_jet.rapidity()) / cluster_weighted_pt;

                sigma_update.xy += q_ji
                    * (particles[particle_iter].rapidity()-my_jet.rapidity())
                    * (eff_delta_phi) / cluster_weighted_pt;

                sigma_update.yx += q_ji
                    * (particles[particle_iter].rapidity()-my_jet.rapidity())
                    * (eff_delta_phi) / cluster_weighted_pt;

                sigma_update.yy += q_ji
                    * (eff_delta_phi)
                    * (eff_delta_phi) / cluster_weighted_pt;
            }

            // if the matrix is looking singular...
            if (sigma_update.xx+sigma_update.yy+sigma_update.xy < 0.01){
                sigma_update.xx=0.5*0.5;
                sigma_update.yy=0.5*0.5;
            }

            // updated sigma is junk if it had almost no contained pT
            if (!(cluster_weighted_pt > 0)){
                sigma_update.xx=0.001*0.001;
                sigma_update.yy=0.001*0.001;
                sigma_update.xy=0.;
                sigma_update.yx=0.;
            }
            mGMM_jets_params->at(cluster_iter) = sigma_update;
        }
    }

    return out_jets;
}



set<unsigned int>
FuzzyTools::ClustersForRemovalGaussian(vecPseudoJet const& mGMM_jets,
                                       vector<MatTwo> const& mGMM_jets_params,
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
        if (mGMM_jets_params[j].xx < epsilon || mGMM_jets_params[j].yy < epsilon){
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
                                vector<double>* mUMM_weights_out,
                                unsigned int &iter_count) {
    assert(kernel_type == FuzzyTools::UNIFORM);

    int cluster_count = seeds.size();
    vector<vector<double> > weights = InitWeights(particles, cluster_count);
    vecPseudoJet mUMM_jets = Initialize(particles, cluster_count, seeds);

    vector<double> mUMM_weights;
    for (int i = 0; i < cluster_count; i++) {
        mUMM_weights.push_back(1.0/cluster_count);
    }
    double log_likelihood_norm_last = 0;
    int iter = 0;
    for(; iter < max_iters; iter++) {
        ComputeWeightsUniform(particles, &weights, cluster_count, mUMM_jets, mUMM_weights);
        mUMM_jets = UpdateJetsUniform(particles, mUMM_jets, weights, cluster_count, &mUMM_weights);
        double log_likelihood_norm = LogLikelihoodUniform(particles, mUMM_jets, mUMM_weights) / particles.size();
        if (log(fabs(log_likelihood_norm_last - log_likelihood_norm)) < log_log_likelihood_limit) break;
        log_likelihood_norm_last = log_likelihood_norm;

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
    iter_count = iter;
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
                                  vector<MatTwo>* mTGMM_jets_params_out,
                                  vector<double>* mTGMM_weights_out,
                                  unsigned int &iter_count){
    assert(kernel_type == FuzzyTools::TRUNCGAUSSIAN);

    int cluster_count = seeds.size();

    vector<double> mTGMM_weights;
    for (int i = 0; i < cluster_count; i++) {
        mTGMM_weights.push_back(1.0/cluster_count);
    }

    vector<vector<double> > weights = InitWeights(particles,cluster_count);
    vecPseudoJet mTGMM_jets = Initialize(particles,cluster_count,seeds);
    vector<MatTwo> mTGMM_jets_params = Initializeparams(particles,cluster_count);
    double log_likelihood_norm_last = 0;
    int iter = 0;
    for (; iter<max_iters; iter++){
        // EM algorithm update steps
        ComputeWeightsTruncGaus(particles,&weights,cluster_count,mTGMM_jets,mTGMM_jets_params, mTGMM_weights);
        mTGMM_jets = UpdateJetsTruncGaus(particles,mTGMM_jets,weights,cluster_count,&mTGMM_jets_params,&mTGMM_weights);
        double log_likelihood_norm = LogLikelihoodTruncGaus(particles, mTGMM_jets, mTGMM_jets_params, mTGMM_weights) / particles.size();
        if (log(fabs(log_likelihood_norm_last - log_likelihood_norm)) < log_log_likelihood_limit) break;
        log_likelihood_norm_last = log_likelihood_norm;
        // do not flag clusters for deletion and remove them
        // if we are not in recombination mode
        if (clustering_mode == FuzzyTools::FIXED) continue;

        // determine which if any clusters should be removed
        set<unsigned int>repeats = ClustersForRemovalGaussian(mTGMM_jets,
                                                              mTGMM_jets_params,
                                                              mTGMM_weights);

        vector<vector<double> >weights_hold;
        vecPseudoJet mTGMM_jets_hold;
        vector<MatTwo> mTGMM_jets_params_hold;
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
    iter_count = iter;
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
FuzzyTools::ClusterFuzzyGaussianC(vecPseudoJet const& particles,
                                 vector<vector<double> >* weights_out,
                                 vector<MatTwo>* mGMMc_jets_params_out,
                                 vector<double>* mGMMc_weights_out,
                                 unsigned int &iter_count){
    assert(kernel_type == FuzzyTools::GAUSSIAN);

    int cluster_count = seeds.size();

    vector<double> mGMMc_weights;
    for (int i = 0; i < cluster_count; i++) {
        mGMMc_weights.push_back(1.0/cluster_count);
    }

    vector<vector<double> > weights = InitWeights(particles,cluster_count);
    vecPseudoJet mGMMc_jets = Initialize(particles,cluster_count,seeds);
    vector<MatTwo> mGMMc_jets_params = Initializeparams(particles,cluster_count);
    double log_likelihood_norm_last = 0;
    int iter=0;
    for (; iter<max_iters; iter++){
        // EM algorithm update steps
        ComputeWeightsGaussian(particles,&weights,cluster_count,mGMMc_jets,mGMMc_jets_params, mGMMc_weights);
        mGMMc_jets = UpdateJetsGaussianC(particles,mGMMc_jets,weights,cluster_count,&mGMMc_jets_params,&mGMMc_weights);

        double log_likelihood_norm = LogLikelihoodGaussian(particles, mGMMc_jets, mGMMc_jets_params, mGMMc_weights) / particles.size();
        if (log(fabs(log_likelihood_norm_last - log_likelihood_norm)) < log_log_likelihood_limit) break;
        log_likelihood_norm_last = log_likelihood_norm;
        // do not flag clusters for deletion and remove them
        // if we are not in recombination mode
        if (clustering_mode == FuzzyTools::FIXED) continue;

        // determine which if any clusters should be removed
        set<unsigned int>repeats = ClustersForRemovalGaussian(mGMMc_jets,
                                                              mGMMc_jets_params,
                                                              mGMMc_weights);

        vector<vector<double> >weights_hold;
        vecPseudoJet mGMMc_jets_hold;
        vector<MatTwo> mGMMc_jets_params_hold;
        vector<double> mGMMc_weights_hold;

        for (unsigned int q=0; q<particles.size(); q++){
            vector<double> hhh;
            weights_hold.push_back(hhh);
        }
        for (unsigned int j=0; j<mGMMc_jets.size(); j++){
            if (repeats.count(j) == 0){
                for (unsigned int q=0; q<particles.size(); q++){
                    weights_hold[q].push_back(weights[q][j]);
                }
                mGMMc_jets_hold.push_back(mGMMc_jets[j]);
                mGMMc_jets_params_hold.push_back(mGMMc_jets_params[j]);
                mGMMc_weights_hold.push_back(mGMMc_weights[j]);
            }
        }

        // now replace and update weights and parameters vectors
        weights.clear();
        mGMMc_jets.clear();
        mGMMc_jets_params.clear();
        mGMMc_weights.clear();

        weights = weights_hold;
        mGMMc_jets = mGMMc_jets_hold;
        mGMMc_jets_params = mGMMc_jets_params_hold;
        mGMMc_weights = mGMMc_weights_hold;
        double total_weight = std::accumulate(mGMMc_weights.begin(), mGMMc_weights.end(), 0.0);
        std::transform(mGMMc_weights.begin(), mGMMc_weights.end(), mGMMc_weights.begin(),
                       std::bind2nd(std::multiplies<double>(), 1.0/total_weight));
        cluster_count = mGMMc_jets_hold.size();

    }
    iter_count= iter;
    weights_out->clear();
    for (unsigned int i=0; i<weights.size(); i++){
        weights_out->push_back(weights[i]);
    }
    mGMMc_jets_params_out->clear();
    mGMMc_weights_out->clear();
    for (unsigned int i=0; i<mGMMc_jets_params.size(); i++){
        mGMMc_jets_params_out->push_back(mGMMc_jets_params[i]);
        mGMMc_weights_out->push_back(mGMMc_weights[i]);
    }

    return mGMMc_jets;

}

#ifdef NEW
set<unsigned int> ClustersForRemoval(__attribute__((unused)) FuzzyTools const& tool,
                                     vector<AbstractKernel *>& jets) {
    assert(jets.size());
    AbstractKernel *first = jets[0];
    if (dynamic_cast<GaussianKernel *> (first)) {
        // we really have a whole bunch of Gaussian kernels
        //cout << "GAUSSIAN" << endl;
    }
    return set<unsigned int>();
}

vector<GaussianKernel *>
MakeGaussianKernels(FuzzyTools const& tool) {
    vector<GaussianKernel *> out;
    vecPseudoJet const& seeds = tool.GetSeeds();
    double initial_weight = 1.0/seeds.size();
    for (unsigned int iter = 0; iter < seeds.size(); iter++) {
        out.push_back(new GaussianKernel(initial_weight, seeds[iter].rapidity(), seeds[iter].phi(), tool.GetDefaultSigma()));
    }
    return out;
}


void ClusterFuzzy(vecPseudoJet const& particles,
                  FuzzyTools const& tool,
                  vector<AbstractKernel *> & jets) {
    double total_particle_pt = 0;
    for (fastjet::PseudoJet const& particle : particles) {
        total_particle_pt += pow(particle.pt(), tool.GetAlpha());
    }

    for (int iter = 0; iter < tool.GetMaxIters(); iter++) {
        // E step
        for (unsigned int particle_iter = 0; particle_iter < particles.size(); particle_iter++) {
            double denom = 0;
            fastjet::PseudoJet particle = particles[particle_iter];
            for (AbstractKernel *p_jet : jets) {
                denom += p_jet->PDF(particle.rapidity(), particle.phi()) * p_jet->Weight();
            }
            for (AbstractKernel *p_jet : jets) {
                double new_weight = p_jet->PDF(particle.rapidity(), particle.phi()) * p_jet->Weight() / denom;
                p_jet->SetParticleWeight(particle_iter, new_weight);
            }
        }

        // M step
        for (AbstractKernel *p_jet : jets) {
            double cluster_weighted_pt = 0;
            vector<double> const& particle_weights = p_jet->ParticleWeights();

            for (unsigned int particle_iter = 0; particle_iter < particles.size(); particle_iter++) {
                cluster_weighted_pt += pow(particles[particle_iter].pt(), tool.GetAlpha()) * particle_weights[particle_iter];
            }

            p_jet->AdditionalUpdate(particles, tool, cluster_weighted_pt);

            double jet_y = 0;
            double jet_phi = 0;
            if (cluster_weighted_pt > 0) {
                for (unsigned int particle_iter = 0; particle_iter < particles.size(); particle_iter++) {
                    fastjet::PseudoJet const& particle = particles[particle_iter];
                    jet_y += pow(particle.pt(), tool.GetAlpha()) * particle_weights[particle_iter]
                        * particle.rapidity() / cluster_weighted_pt;
                    jet_phi += pow(particle.pt(), tool.GetAlpha()) * particle_weights[particle_iter]
                        * particle.phi() / cluster_weighted_pt;
                }
            }
            p_jet->SetLocation(jet_y, jet_phi);

            if (tool.GetLearnWeights()) {
                p_jet->SetWeight(cluster_weighted_pt / total_particle_pt);
            }
        }

        // remove clusters if necessary
        if (tool.GetClusteringMode() == FuzzyTools::FIXED) continue;

        set<unsigned int> repeats = ClustersForRemoval(tool, jets);

        vector<AbstractKernel *> new_jets;
        for (unsigned int jet_iter = 0; jet_iter < jets.size(); jet_iter++) {
            if (repeats.count(jet_iter) != 0) {
                delete jets[jet_iter];
            } else {
                new_jets.push_back(jets[jet_iter]);
            }
        }
        jets = new_jets;
    }
}
#endif

vecPseudoJet
FuzzyTools::ClusterFuzzyGaussian(vecPseudoJet const& particles,
                                 vector<vector<double> >* weights_out,
                                 vector<MatTwo>* mGMM_jets_params_out,
                                 vector<double>* mGMM_weights_out,
                                 unsigned int &iter_count){
    assert(kernel_type == FuzzyTools::GAUSSIAN);

    int cluster_count = seeds.size();

    vector<double> mGMM_weights;
    for (int i = 0; i < cluster_count; i++) {
        mGMM_weights.push_back(1.0/cluster_count);
    }

    vector<vector<double> > weights = InitWeights(particles,cluster_count);
    vecPseudoJet mGMM_jets = Initialize(particles,cluster_count,seeds);
    vector<MatTwo> mGMM_jets_params = Initializeparams(particles,cluster_count);
    double log_likelihood_norm_last = 0;
    int iter=0;
    for (; iter<max_iters; iter++){
        // EM algorithm update steps
        ComputeWeightsGaussian(particles,&weights,cluster_count,mGMM_jets,mGMM_jets_params, mGMM_weights);
        mGMM_jets = UpdateJetsGaussian(particles,mGMM_jets,weights,cluster_count,&mGMM_jets_params,&mGMM_weights);
        double log_likelihood_norm = LogLikelihoodGaussian(particles, mGMM_jets, mGMM_jets_params, mGMM_weights) / particles.size();
        if (log(fabs(log_likelihood_norm_last - log_likelihood_norm)) < log_log_likelihood_limit) break;
        log_likelihood_norm_last = log_likelihood_norm;
        // do not flag clusters for deletion and remove them
        // if we are not in recombination mode
        if (clustering_mode == FuzzyTools::FIXED) continue;

        // determine which if any clusters should be removed
        set<unsigned int>repeats = ClustersForRemovalGaussian(mGMM_jets,
                                                              mGMM_jets_params,
                                                              mGMM_weights);

        vector<vector<double> >weights_hold;
        vecPseudoJet mGMM_jets_hold;
        vector<MatTwo> mGMM_jets_params_hold;
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
    iter_count = iter;
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

double FuzzyTools::LogLikelihoodGaussian(vecPseudoJet const& particles,
                                         vecPseudoJet const& mGMM_jets,
                                         vector<MatTwo> const& mGMM_jets_params,
                                         vector<double> const& mGMM_weights) {
    double log_likelihood = 0;
    for (unsigned int particle_iter = 0; particle_iter < particles.size(); particle_iter++) {
        fastjet::PseudoJet p = particles[particle_iter];
        double likelihood_summand = 0;
        for (unsigned int cluster_iter = 0; cluster_iter < mGMM_weights.size(); cluster_iter++) {
            fastjet::PseudoJet cluster = mGMM_jets[cluster_iter];
            likelihood_summand += mGMM_weights[cluster_iter] *
                doGaus(p.rapidity(), p.phi(), cluster.rapidity(), cluster.phi(), mGMM_jets_params[cluster_iter]);

        }
        log_likelihood += log(likelihood_summand);
    }
    return log_likelihood;
}

double FuzzyTools::LogLikelihoodTruncGaus(vecPseudoJet const& particles,
                                          vecPseudoJet const& mTGMM_jets,
                                          vector<MatTwo> const& mTGMM_jets_params,
                                          vector<double> const& mTGMM_weights) {
    double log_likelihood = 0;
    for (unsigned int particle_iter = 0; particle_iter < particles.size(); particle_iter++) {
        fastjet::PseudoJet p = particles[particle_iter];
        double likelihood_summand = 0;
        for (unsigned int cluster_iter = 0; cluster_iter < mTGMM_weights.size(); cluster_iter++) {
            fastjet::PseudoJet cluster = mTGMM_jets[cluster_iter];
            likelihood_summand += mTGMM_weights[cluster_iter] *
                doTruncGaus(p.rapidity(), p.phi(), cluster.rapidity(), cluster.phi(), mTGMM_jets_params[cluster_iter]);
        }
        log_likelihood += log(likelihood_summand);
    }
    return log_likelihood;
}

double FuzzyTools::LogLikelihoodUniform(vecPseudoJet const& particles,
                                        vecPseudoJet const& mUMM_jets,
                                        vector<double> const& mUMM_weights) {
    double log_likelihood = 0;
    for (unsigned int particle_iter = 0; particle_iter < particles.size(); particle_iter++) {
        fastjet::PseudoJet p = particles[particle_iter];
        double likelihood_summand = 0;
        for (unsigned int cluster_iter = 0; cluster_iter < mUMM_weights.size(); cluster_iter++) {
            fastjet::PseudoJet cluster = mUMM_jets[cluster_iter];
            if (p.delta_R(cluster) > R) continue;
            likelihood_summand += mUMM_weights[cluster_iter] * M_PI / (R*R);
        }
        log_likelihood += log(likelihood_summand);
    }
    return log_likelihood;
}

void
FuzzyTools::NewEventDisplay(__attribute__((unused)) vecPseudoJet const& particles,
                            __attribute__((unused)) vecPseudoJet const& ca_jets,
                            __attribute__((unused)) vecPseudoJet const& tops,
                            __attribute__((unused)) vecPseudoJet const& mGMM_jets,
                            __attribute__((unused)) vector<vector<double> > const& weights,
                            __attribute__((unused)) int which,
                            __attribute__((unused)) vector<MatTwo> const& mGMM_jets_params,
                            __attribute__((unused)) vector<double> const& mGMM_weights,
                            __attribute__((unused)) std::string const& out,
                            __attribute__((unused)) int iter) {
    #ifdef WITHROOT
    double min_eta = -5;
    double max_eta = 5;
    TCanvas canv(TString::Format("NEVC_%s_%d", out.c_str(), iter), "", 1200, 600);
    TH2F hist(TString::Format("NEVH_%s_%d", out.c_str(), iter), "", 30, min_eta, max_eta, 28, 0, 7);
    TTree aux(TString::Format("NEVT_%s_%d", out.c_str(), iter), "");

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
        double var_eta = mGMM_jets_params[i].xx;
        double var_phi = mGMM_jets_params[i].yy;
        double covar  = mGMM_jets_params[i].yx;
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
    canv.Write(TString::Format("%sNewEventDisplay_%s_%d",
                               directory_prefix.c_str(),
                               out.c_str(), iter));
    canv.Print(TString::Format("%sNewEventDisplay_%s_%d.pdf",
                               directory_prefix.c_str(),
                               out.c_str(), iter));
    #endif
}

void
FuzzyTools::NewEventDisplayPoster(__attribute__((unused)) vecPseudoJet const& particles,
                            __attribute__((unused)) vecPseudoJet const& ca_jets,
                            __attribute__((unused)) vecPseudoJet const& tops,
                            __attribute__((unused)) vecPseudoJet const& mGMM_jets,
                            __attribute__((unused)) vector<vector<double> > const& weights,
                            __attribute__((unused)) int which,
                            __attribute__((unused)) vector<MatTwo> const& mGMM_jets_params,
                            __attribute__((unused)) vector<double> const& mGMM_weights,
                            __attribute__((unused)) std::string const& out,
                            __attribute__((unused)) int iter) {
    #ifdef WITHROOT
    double min_eta = -5;
    double max_eta = 5;
    TCanvas canv(TString::Format("NEVC_%s_%d", out.c_str(), iter), "", 1000, 1000);
    TH2F hist(TString::Format("NEVH_%s_%d", out.c_str(), iter), "", 50, min_eta, max_eta, 35, 0, 7);
    TTree aux(TString::Format("NEVT_%s_%d", out.c_str(), iter), "");

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
        double var_eta = mGMM_jets_params[i].xx;
        double var_phi = mGMM_jets_params[i].yy;
        double covar  = mGMM_jets_params[i].yx;
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

    gStyle->SetPadRightMargin(0.15);
    hist.SetXTitle("#eta (Rapidity)");
    hist.SetYTitle("#phi");
    hist.Draw("colz");
    for (unsigned int eli = 0; eli < ellipses.size(); eli++) {
        ellipses[eli].Draw();
    }
    gPad->SetLogz();
    canv.Update();

    canv.Print(TString::Format("%sNewEventDisplayPoster_%s_%d.pdf",
                               directory_prefix.c_str(),
                               out.c_str(), iter));
    #endif
}

void
FuzzyTools::NewEventDisplayUniform(__attribute__((unused)) vecPseudoJet const& particles,
                                   __attribute__((unused)) vecPseudoJet const& ca_jets,
                                   __attribute__((unused)) vecPseudoJet const& tops,
                                   __attribute__((unused)) vecPseudoJet const& mUMM_jets,
                                   __attribute__((unused)) vector<vector<double> > const& weights,
                                   __attribute__((unused)) int which,
                                   __attribute__((unused)) vector<double> const& mUMM_weights,
                                   __attribute__((unused)) std::string const& out,
                                   __attribute__((unused)) int iter) {
    #ifdef WITHROOT
    double min_eta = -5;
    double max_eta = 5;
    TCanvas canv(TString::Format("NEVC_%s_%d", out.c_str(), iter), "", 1200, 600);
    TTree aux(TString::Format("NEVT_%s_%d", out.c_str(), iter), "");

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


    TH2F hist(TString::Format("NEVH_%s_%d", out.c_str(), iter), "TEST", 30, min_eta, max_eta, 28, 0, 7);

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

    canv.Write(TString::Format("%sNewEventDisplay_%s_%d",
                               directory_prefix.c_str(),
                               out.c_str(), iter));
    #endif
}

void
FuzzyTools::SubsEventDisplay(vecPseudoJet const& particles,
                             vecPseudoJet const& mGMM_jets,
                             vector<vector<double> > const& weights,
                             int lead_mGMM_index,
                             vector<MatTwo> const& mGMM_jets_params,
                             std::string const& file_name,
                             std::string const& title) {
    #ifdef WITHROOT
    gStyle->SetOptStat(0);
    TCanvas * c = new TCanvas("c", title.c_str(), 500, 500);

    const int particle_count = particles.size();

    map<int, TGraph*> vars;
    for (int particle_iter = 0; particle_iter < particle_count; particle_iter++) {
        double x = particles[particle_iter].rapidity();
        double y = particles[particle_iter].phi();

        vars[particle_iter] = new TGraph(1, &x, &y);
        int mycolor = 17 - floor(weights[particle_iter][lead_mGMM_index] * 8);
        if (mycolor < 12) mycolor = 1;
        vars[particle_iter]->SetMarkerColor(mycolor);
    }

    const int mGMM_jet_count = mGMM_jets.size();
    double x_jet[mGMM_jet_count];
    double y_jet[mGMM_jet_count];
    for (int jet_iter = 0; jet_iter < mGMM_jet_count; jet_iter++) {
        x_jet[jet_iter] = mGMM_jets[jet_iter].rapidity();
        y_jet[jet_iter] = mGMM_jets[jet_iter].phi();
    }
    TGraph *jet_graph = new TGraph(mGMM_jet_count, x_jet, y_jet);

    TH1F * background = new TH1F("", "", 100, -4, 8);
    background->GetXaxis()->SetTitle("Rapidity");
    background->GetYaxis()->SetTitleOffset(1.4);
    background->GetYaxis()->SetTitle("Azimuthal Angle [rad]");
    background->GetXaxis()->SetNdivisions(505);
    background->GetYaxis()->SetRangeUser(0,7);
    background->Draw();

    for (int i=0; i<particle_count; i++) {
        vars[i]->SetMarkerSize(1);
        vars[i]->SetMarkerStyle(20);
        vars[i]->Draw("samep");
    }

    jet_graph->SetMarkerSize(2);
    jet_graph->SetMarkerStyle(3);
    jet_graph->SetMarkerColor(6);
    jet_graph->Draw("psame");

    for (unsigned int i=0; i<mGMM_jets_params.size(); i++){
        double a = mGMM_jets_params[i].xx;
        double b = mGMM_jets_params[i].yy;
        double c = mGMM_jets_params[i].yx;
        double lambda1 = 0.5*(a+b)+0.5*sqrt((a-b)*(a-b)+4.*c*c);
        double lambda2 = 0.5*(a+b)-0.5*sqrt((a-b)*(a-b)+4.*c*c);
        double theta = 0.;
        if (c>0) theta=atan((lambda1-a)/c);
        // std::cout << "yo " << i << " " << sqrt(lambda1) << " " << sqrt(lambda2) << " " << theta << std::endl;
        TEllipse *el4 = new TEllipse(x_jet[i],y_jet[i],sqrt(lambda1),sqrt(lambda2),0,360,theta*180/TMath::Pi());
        el4->SetFillStyle(0);
        el4->Draw("same");
    }

    TLegend* leggaa = new TLegend(.7,.33,0.95,.67);
    leggaa->SetTextFont(42);

    leggaa->AddEntry(jet_graph,"Fuzzy Jets","p");
    leggaa->SetFillStyle(0);
    leggaa->SetFillColor(0);
    leggaa->SetBorderSize(0);
    leggaa->Draw();


    c->Print(file_name.c_str(), "pdf");


    delete c;
    #endif
}

void
FuzzyTools::EventDisplay(__attribute__((unused)) vecPseudoJet const& particles,
                         __attribute__((unused)) vecPseudoJet const& ca_jets,
                         __attribute__((unused)) vecPseudoJet const& tops,
                         __attribute__((unused)) vecPseudoJet const& mGMM_jets,
                         __attribute__((unused)) vector<vector<double> > const& weights,
                         __attribute__((unused)) int which,
                         __attribute__((unused)) vector<MatTwo> const& mGMM_jets_params,
                         __attribute__((unused)) std::string const& out,
                         __attribute__((unused)) int iter) {
    #ifdef WITHROOT
    gStyle->SetOptStat(0);
    //gROOT->Reset();
    //gROOT->SetStyle("ATLAS");
    //gROOT->ForceStyle();
    //gStyle->SetPadLeftMargin(0.15);
    //gStyle->SetPadTopMargin(0.15);

    TCanvas *c = new TCanvas(TString::Format("EventDisplayOld_%s_%d", out.c_str(), iter),"",500,500);

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
        double a = mGMM_jets_params[i].xx;
        double b = mGMM_jets_params[i].yy;
        double c = mGMM_jets_params[i].yx;
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

    TString file_name = TString::Format("%sEvent%s_%d.root", directory_prefix.c_str(), out.c_str(), iter);
    c->Print(file_name);
    delete c;
    #endif
}

void
FuzzyTools::Qjetmass(__attribute__((unused)) vecPseudoJet particles,
                     __attribute__((unused)) vector<vector<double> > weights,
                     __attribute__((unused)) int which,
                     __attribute__((unused)) std::string out){
    #ifdef WITHROOT
    TH1F q_jet_mass(TString::Format("qjetmass%s", out.c_str()),"",100,0,250);
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
    #endif
}

void
FuzzyTools::JetContributionDisplay(__attribute__((unused)) vecPseudoJet particles,
                                   __attribute__((unused)) vector<vector<double> > weights,
                                   __attribute__((unused)) int which,
                                   __attribute__((unused)) int m_type,
                                   __attribute__((unused)) std::string out,
                                   __attribute__((unused)) int iter) {
    #ifdef WTIHROOT
    double min_eta = -5;
    double max_eta = 5;
    unsigned int k = weights[0].size();
    TH2F hist(TString::Format("JCH_hard_%s_%d", out, iter), "", 30, min_eta, max_eta, 28, 0, 7);
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

    TH2F hist2(TString::Format("JCH_soft_%s_%d", out, iter), "", 30, min_eta, max_eta, 28, 0, 7);
    for (unsigned int i=0; i < particles.size(); i++) {
        double fill_val = m_type == 1 ? particles[i].m() : particles[i].pt();
        fill_val *= weights[i][which];
        hist2.Fill(particles[i].eta(), particles[i].phi(), fill_val);
    }

    hist2.Write();
    #endif
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

    #ifdef WITHROOT
    TRandom3 rand = TRandom3(1);
    UInt_t seed = rand.GetSeed();
    #else
    time_t seed = time(NULL);
    srand(seed);
    #endif
    // compute throws and the mean
    moments.push_back(0);
    for (unsigned int trial = 0; trial < n_trials; trial++) {
        for (unsigned int particle_iter=0; particle_iter<particles.size(); particle_iter++) {
            #ifdef WITHROOT
            double my_throw = rand.Uniform(0, 1);
            #else
            double my_throw = ((double) rand() / (RAND_MAX));
            #endif
            if (my_throw < weights[particle_iter][cluster_id]) {
                particle_indices.push_back(particle_iter);
            }
        }
        moments[0] += (*f)(particles, particle_indices);
        particle_indices.clear();
    }

    moments[0] /= n_trials;
    #ifdef WITHROOT
    rand.SetSeed(seed);
    #else
    srand(seed);
    #endif
    particle_indices.clear();
    for (unsigned int i = 1; i < moment_count; i++) {
        moments.push_back(0);
    }

    for(unsigned int trial = 0; trial < n_trials; trial++) {
        for (unsigned int particle_iter = 0; particle_iter < particles.size(); particle_iter++) {
            #ifdef WITHROOT
            double my_throw = rand.Uniform(0, 1);
            #else
            double my_throw = ((double) rand() / (RAND_MAX));
            #endif
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
                          MatTwo const& jet_params,
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
                           MatTwo const& jet_params,
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
