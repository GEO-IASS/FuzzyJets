#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <set>


#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TMath.h"
#include "TVector3.h"
#include "TMatrix.h"

#include "FuzzyAnalysis.h"
#include "FuzzyTools.h"

#include "myFastJetBase.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/Selector.hh"

#include "Pythia8/Pythia.h"

using namespace std;


// Privately separate the logic of different analysis modes from using them
namespace {
    void WeightDistribution(vector<vector<double> > const& Weights,
                            int which,
                            TString out) {
        TTree aux("WT" + out, "");
        double w;
        int pCount = Weights.size();
        aux.Branch("w", &w, "w/D");
        for (int pidx = 0; pidx < pCount; pidx++) {
            w = Weights[pidx][which];
            aux.Fill();
        }
        aux.Write();
    }

    // is the reference particle given by index pidx clustered at all?
    bool isClustered(vector<vector<double> > const& Weights,
                     int pidx) {
        bool clustered = false;
        int nClusters = Weights[0].size();
        for (int j = 0; j < nClusters; j++) {
            double w = Weights[pidx][j];
            if (!isnan(w) && w > 0) {
                clustered = true;
            }
        }
        return clustered;
    }

    // what is the index of the fuzzy jet that the particle belongs to?
    int belongsTo(vector<vector<double> > const& Weights,
                  int pidx) {
        double wmax = -1;
        int bestidx = -1;
        const int nClusters = Weights[0].size();
        for (int j = 0; j < nClusters; j++) {
            double w = Weights[pidx][j];
            if(!isnan(w) && w > 0 && w > wmax) {
                wmax = w;
                bestidx = j;
            }
        }
        return bestidx;
    }

    // Compute the mass of a Fuzzy Jet which comes from pileup
    // Please note that this is somewhat ill defined, due to how particles add
    double JetPuMassHard(vecPseudoJet const& particles,
                             vector<vector<double> > const& Weights,
                             int jetidx) {
        fastjet::PseudoJet pujet;
        const int nParticles = particles.size();
        for (int i = 0; i < nParticles; i++) {
            if (belongsTo(Weights, i) == jetidx &&
                particles[i].user_info<MyUserInfo>().isPU()) {
                pujet += particles[i];
            }
        }
        return pujet.m();
    }

    // Compute the soft mass due to pileup
    // Please note that this is somewhat ill defined, due to how particles add
    double JetPuMassSoft(vecPseudoJet const& particles,
                             vector<vector<double> > const& Weights,
                             int jetidx) {
        fastjet::PseudoJet pujet;
        const int nParticles = particles.size();
        for (int i = 0; i < nParticles; i++) {
            if (particles[i].user_info<MyUserInfo>().isPU()) {
                pujet += particles[i] * Weights[i][jetidx];
            }
        }
        return pujet.m();
    }

    // Compute the mass of a fastjet jet which is due to pileup
    // Please note that this is somewhat ill defined, due to how particles add
    double JetPuMassFastjet(fastjet::PseudoJet const& jet) {
        const vecPseudoJet c = jet.constituents();
        const int nParticles = c.size();
        fastjet::PseudoJet myjet;
        for (int i = 0; i < nParticles; i++) {
            if(c[i].user_info<MyUserInfo>().isPU()) {
                myjet += c[i];
            }
        }
        return myjet.m();
    }

    // Compute the fraction of particles in a jet which are pileup
    double JetPuFracHard(vecPseudoJet const& particles,
                         vector<vector<double> > const& Weights,
                         int jetidx) {
        const int nParticles = particles.size();
        int clusteredParticles = 0;
        int clusteredPu = 0;
        for (int i = 0; i < nParticles; i++) {
            if (belongsTo(Weights, i) == jetidx) {
                clusteredParticles++;
                if(particles[i].user_info<MyUserInfo>().isPU()) {
                    clusteredPu++;
                }
            }
        }
        return (double)clusteredPu / clusteredParticles;
    }

    // Compute the fraction of particles in a jet which are pileup by soft assignment
    double JetPuFracSoft(vecPseudoJet const& particles,
                         vector<vector<double> > const& Weights,
                         int jetidx) {
        const int nParticles = particles.size();
        double clusteredParticles = 0;
        double clusteredPu = 0;
        for (int i = 0; i < nParticles; i++) {
            if (isnan(Weights[i][jetidx])) continue;
            clusteredParticles += Weights[i][jetidx];
            if (particles[i].user_info<MyUserInfo>().isPU()) {
                clusteredPu += Weights[i][jetidx];
            }
        }
        return clusteredPu / clusteredParticles;
    }

    // Compute the fraction of particles in a jet which are pileup for standard jets
    double JetPuFracFastjet(fastjet::PseudoJet const& jet) {
        const vecPseudoJet c = jet.constituents();
        const int nParticles = c.size();
        int pu = 0;
        for (int i = 0; i < nParticles; i++) {
            if (c[i].user_info<MyUserInfo>().isPU()) {
                pu++;
            }
        }
        return (double)pu / nParticles;
    }

    // Compute (#clustered pileup particles) / (#pileup particles)
    double ClusteredPileupFrac(vecPseudoJet const& particles,
                               vector<vector<double> > const& Weights) {
        int pileup = 0;
        int clusteredPileup = 0;
        for(unsigned int i = 0; i < particles.size(); i++) {
            fastjet::PseudoJet p = particles[i];
            if (p.user_info<MyUserInfo>().isPU()) {
                pileup++;
                if(isClustered(Weights, i))
                    clusteredPileup++;
            }
        }
        return 1.0*clusteredPileup / pileup;
    }

    // Compute (#unclustered pileup particles) / (#unclustered particles)
    double UnclusteredPileupComposition(vecPseudoJet particles,
                                        vector<vector<double> > const& Weights) {
        int unclusteredParticles = 0;
        int unclusteredPileup = 0;
        for(unsigned int i = 0; i < particles.size(); i++) {
            if(!isClustered(Weights, i)) {
                unclusteredParticles++;
                fastjet::PseudoJet p = particles[i];
                if(p.user_info<MyUserInfo>().isPU()) {
                    unclusteredPileup++;
                }
            }
        }
        return (double)unclusteredPileup / unclusteredParticles;
    }

    void DoFastJetFinding(vecPseudoJet particles,
                          fastjet::JetDefinition *pJetDef,
                          double pTmin,
                          vecPseudoJet& jets) {
        fastjet::ClusterSequence csLargeR_ca(particles, *pJetDef);
        jets = fastjet::sorted_by_pt(csLargeR_ca.inclusive_jets(pTmin));
    }

    void FindLeadingJet(vecPseudoJet& particles,
                        vecPseudoJet& jets,
                        vector<vector<double> >& particleWeights,
                        FuzzyTools *tool,
                        int& index,
                        double& pT) {
        pT = -1;
        index = -1;
        int clusterCount = particleWeights[0].size();
        for (unsigned int i=0; i < jets.size(); i++) {
            double holdpT = tool->MLpT(particles, particleWeights, i,
                                       clusterCount, 0);
            if (holdpT > pT) {
                pT = holdpT;
                index = i;
            }
        }
    }

    void DoMUMMJetFinding(vecPseudoJet& particles,
                          bool learnWeights,
                          double size,
                          int& leadIndex,
                          double& leadpT,
                          FuzzyTools *tool,
                          vecPseudoJet& jets,
                          vector<vector<double> >& particleWeights,
                          vector<double>& jetWeights) {
        vecPseudoJet parts = particles;
        tool->SetKernelType(FuzzyTools::UNIFORM);
        tool->SetSeeds(parts);
        tool->SetLearnWeights(learnWeights);
        tool->SetClusteringMode(FuzzyTools::RECOMBINATION);
        tool->SetR(size);
        jets = tool->ClusterFuzzyUniform(particles,
                                         &particleWeights,
                                         &jetWeights);
        FindLeadingJet(particles, jets, particleWeights, tool, leadIndex, leadpT);
    }

    void DoMTGMMJetFinding(vecPseudoJet& particles,
                           vecPseudoJet& seeds,
                           bool learnWeights,
                           bool learnShape,
                           bool size,
                           int& leadIndex,
                           double& leadpT,
                           FuzzyTools *tool,
                           vecPseudoJet& jets,
                           vector<vector<double> >& particleWeights,
                           vector<TMatrix>& parameters,
                           vector<double>& jetWeights) {
        tool->SetKernelType(FuzzyTools::TRUNCGAUSSIAN);
        tool->SetSeeds(seeds);
        tool->SetLearnWeights(learnWeights);
        if(learnShape) {
            tool->SetClusteringMode(FuzzyTools::FIXED);
        } else {
            tool->SetClusteringMode(FuzzyTools::RECOMBINATION);
        }
        tool->SetLearnShape(learnShape);
        tool->SetR(size);
        jets = tool->ClusterFuzzyTruncGaus(particles,
                                           &particleWeights,
                                           &parameters,
                                           &jetWeights);
        FindLeadingJet(particles, jets, particleWeights, tool, leadIndex, leadpT);
    }

    void DoMGMMJetFinding(vecPseudoJet& particles,
                          vecPseudoJet& seeds,
                          bool learnWeights,
                          bool learnShape,
                          int& leadIndex,
                          double& leadpT,
                          FuzzyTools *tool,
                          vecPseudoJet& jets,
                          vector<vector<double> >& particleWeights,
                          vector<TMatrix>& parameters,
                          vector<double>& jetWeights) {
        tool->SetKernelType(FuzzyTools::GAUSSIAN);
        tool->SetSeeds(seeds);
        tool->SetLearnWeights(learnWeights);
        if(learnShape) {
            tool->SetClusteringMode(FuzzyTools::FIXED);
        } else {
            tool->SetClusteringMode(FuzzyTools::RECOMBINATION);
        }
        tool->SetLearnShape(learnShape);
        jets = tool->ClusterFuzzyGaussian(particles,
                                          &particleWeights,
                                          &parameters,
                                          &jetWeights);
        FindLeadingJet(particles, jets, particleWeights, tool, leadIndex, leadpT);
    }
}

// Constructor
FuzzyAnalysis::FuzzyAnalysis(){
    fDebug = false;
    if(fDebug) cout << "FuzzyAnalysis::FuzzyAnalysis Start " << endl;
    ftest = 0;
    fOutName = "test.root";
    directoryPrefix = "results/";

    batched = false;

    tool = new FuzzyTools();

    tool->SetClusteringMode(FuzzyTools::RECOMBINATION);


    // jet def
    m_jet_def                = new fastjet::JetDefinition(fastjet::antikt_algorithm, 0.4);
    m_jet_def_largeR_antikt  = new fastjet::JetDefinition(fastjet::antikt_algorithm, 1.0);
    m_jet_def_largeR_ca      = new fastjet::JetDefinition(fastjet::cambridge_algorithm, 1.0);

    if(fDebug) cout << "FuzzyAnalysis::FuzzyAnalysis End " << endl;
}

// Destructor
FuzzyAnalysis::~FuzzyAnalysis(){
    delete tool;
    delete m_jet_def;
    delete m_jet_def_largeR_antikt;
    delete m_jet_def_largeR_ca;
}

// Begin method
void FuzzyAnalysis::Begin(){
    // Declare TTree
    string fullName = directoryPrefix + fOutName;
    tF = new TFile(fullName.c_str(), "RECREATE");
    tT = new TTree("EventTree", "Event Tree for Fuzzy");

    DeclareBranches();
    ResetBranches();

    tool->SetPrefix(directoryPrefix);

    return;
}

// End
void FuzzyAnalysis::End(){

    tT->Write();
    tF->Close();
    return;
}

// Analyze
void FuzzyAnalysis::AnalyzeEvent(int ievt, Pythia8::Pythia* pythia8, Pythia8::Pythia* pythia_MB, int NPV){
    if(fDebug) cout << "FuzzyAnalysis::AnalyzeEvent Begin " << endl;

    // generate a new event
    if (!pythia8->next()) return;
    if (!pythia_MB->next()) return;
    if(fDebug) cout << "FuzzyAnalysis::AnalyzeEvent Event Number " << ievt << endl;

    // reset branches
    ResetBranches();

    // new event-----------------------
    fTEventNumber = ievt;
    std::vector <fastjet::PseudoJet>           particlesForJets;
    vector<fastjet::PseudoJet> tops;
    fastjet::PseudoJet dl;
    tops.push_back(dl);
    fastjet::PseudoJet dl2;
    tops.push_back(dl);

    fTNPV = NPV;
    // Pileup loop -------------------------------------------------------------
    double px, py, pz, e;
    for (int iPU = 0; iPU <= NPV; ++iPU) {
        for (unsigned int ip = 0; ip < (unsigned) pythia_MB->event.size(); ++ip) {
            if (!pythia_MB->event[ip].isFinal()) continue;
            if (fabs(pythia_MB->event[ip].id())==12) continue;
            if (fabs(pythia_MB->event[ip].id())==14) continue;
            if (fabs(pythia_MB->event[ip].id())==13) continue;
            if (fabs(pythia_MB->event[ip].id())==16) continue;
            //if (pythia_MB->event[ip].pT() < 0.5)     continue;
            px = pythia_MB->event[ip].px();
            py = pythia_MB->event[ip].py();
            pz = pythia_MB->event[ip].pz();
            e  = pythia_MB->event[ip].e();

            fastjet::PseudoJet p(px, py, pz, e);
            p.reset_PtYPhiM(p.pt(), p.rapidity(), p.phi(), 0.);
            // note that we don't really keep the particle number! only store ip, should really store ip and iPU
            p.set_user_info(new MyUserInfo(pythia_MB->event[ip].id(),ip,pythia_MB->event[ip].charge(),true));
            particlesForJets.push_back(p);
        }
        if (!pythia_MB->next()) continue;
    }

    // Particle loop -----------------------------------------------------------
    // The Pythia event listing contains a lot more than we want to process,
    // we prune out certain particles (muons / neutrinos) and only add final
    // state particles
    for (unsigned int ip=0; ip < (unsigned) pythia8->event.size(); ++ip){
        px = pythia8->event[ip].px();
        py = pythia8->event[ip].py();
        pz = pythia8->event[ip].pz();
        e  = pythia8->event[ip].e();
        fastjet::PseudoJet p(px, py, pz, e);
        p.set_user_info(new MyUserInfo(pythia8->event[ip].id(),ip,pythia8->event[ip].charge(),false));

        // In reality we should be more careful about finding tops,
        // but this will do for now. In the future consider refactoring
        // and tracing to find a top quark with no daughters
        if (pythia8->event[ip].id()  ==6) tops[0]=p;
        if (pythia8->event[ip].id()  ==-6) tops[1]=p;

        // prune uninteresting particles
        if (!pythia8->event[ip].isFinal() )      continue; // only final state
        if (fabs(pythia8->event[ip].id())  ==12) continue; // prune nu-e
        if (fabs(pythia8->event[ip].id())  ==13) continue; // ...   mu
        if (fabs(pythia8->event[ip].id())  ==14) continue; // ...   nu-mu
        if (fabs(pythia8->event[ip].id())  ==16) continue; // ...   nu-tau
        if (pythia8->event[ip].pT()       < 0.5) continue; // ...   low pT

        particlesForJets.push_back(p);

    } // end particle loop -----------------------------------------------

    // large-R jets: C/A --------------------
    double pTmin = 50;
    fastjet::ClusterSequence csLargeR_ca(particlesForJets, *m_jet_def_largeR_ca);
    vecPseudoJet myJetsLargeR_ca = fastjet::sorted_by_pt(csLargeR_ca.inclusive_jets(pTmin));

    fTCA_m = myJetsLargeR_ca[0].m();
    fTCA_pt = myJetsLargeR_ca[0].pt();
    fTCA_pufrac = JetPuFracFastjet(myJetsLargeR_ca[0]);
    fTCA_m_pu = JetPuMassFastjet(myJetsLargeR_ca[0]);

    // anti-kt R:1.0 trimmed ----------------
    fastjet::Filter filter_two(0.2, fastjet::SelectorPtFractionMin(0.05));
    fastjet::Filter filter_three(0.3, fastjet::SelectorPtFractionMin(0.05));

    fastjet::ClusterSequence csLargeR_antikt(particlesForJets, *m_jet_def_largeR_antikt);
    vecPseudoJet myJetsLargeR_antikt = fastjet::sorted_by_pt(csLargeR_antikt.inclusive_jets(pTmin));

    fastjet::PseudoJet leadAkt = myJetsLargeR_antikt[0];
    const fastjet::PseudoJet leadAkt_filter_two = filter_two(leadAkt);
    const fastjet::PseudoJet leadAkt_filter_three = filter_three(leadAkt);
    fTantikt_m = leadAkt.m();
    fTantikt_pt = leadAkt.pt();
    fTantikt_m_trimmed_two = leadAkt_filter_two.m();
    fTantikt_pt_trimmed_two = leadAkt_filter_two.pt();
    fTantikt_m_trimmed_three = leadAkt_filter_three.m();
    fTantikt_pt_trimmed_three = leadAkt_filter_three.pt();
    fTantikt_pufrac_trimmed_two = JetPuFracFastjet(leadAkt_filter_two);
    fTantikt_pufrac_trimmed_three = JetPuFracFastjet(leadAkt_filter_three);
    fTantikt_m_pu = JetPuMassFastjet(leadAkt);
    fTantikt_m_pu_trimmed_two = JetPuMassFastjet(leadAkt_filter_two);
    fTantikt_m_pu_trimmed_three = JetPuMassFastjet(leadAkt_filter_three);

    // ======================================
    // Various mixture models ---------------
    // ======================================
    tool->SetMergeDistance(0.05);

    // which jets to run
    bool mUMMon = true;
    bool mGMMon = true;
    __attribute__((unused)) bool mGMMson = true;
    bool mTGMMon = true;
    __attribute__((unused)) bool mTGMMson = true;


    // Fuzzy Jets: mGMMs --------------------
    vector<vector<double> > mGMMsparticleWeights;
    vector<TMatrix> mGMMsjetsparams;
    vector<double> mGMMsweights;
    vecPseudoJet mGMMsjets;
    int leadmGMMsindex;
    double maxpTmGMMs;
    if(mGMMson) {
        DoMGMMJetFinding(particlesForJets, myJetsLargeR_ca,
                         fLearnWeights, true,
                         leadmGMMsindex, maxpTmGMMs,
                         tool, mGMMsjets, mGMMsparticleWeights,
                         mGMMsjetsparams, mGMMsweights);
    }

    // Fuzzy Jets: mTGMMs -------------------
    vector<vector<double > > mTGMMsparticleWeights;
    vector<TMatrix> mTGMMsjetsparams;
    vector<double> mTGMMsweights;
    vecPseudoJet mTGMMsjets;
    int leadmTGMMsindex;
    double maxpTmTGMMs;
    if(mTGMMson) {
        DoMTGMMJetFinding(particlesForJets, myJetsLargeR_ca,
                          fLearnWeights, true,
                          fSize, leadmTGMMsindex, maxpTmTGMMs,
                          tool, mTGMMsjets, mTGMMsparticleWeights,
                          mTGMMsjetsparams, mTGMMsweights);
    }

    // Fuzzy Jets: mGMM ---------------------
    vector<vector<double> > Weights;
    vector<TMatrix> mGMMjetsparams;
    vector<double> mGMMweights;
    vector<fastjet::PseudoJet> mGMMjets;
    int leadmGMMindex;
    double maxpTmGMM;
    if(mGMMon) {
        DoMGMMJetFinding(particlesForJets, particlesForJets,
                         fLearnWeights, false,
                         leadmGMMindex, maxpTmGMM,
                         tool, mGMMjets, Weights,
                         mGMMjetsparams, mGMMweights);
    }

    // Fuzzy Jets: mUMM ---------------------
    vector<vector<double> > mUMMparticleWeights;
    vector<double> mUMMweights;
    vecPseudoJet mUMMjets;

    int leadmUMMindex;
    double maxpTmUMM;
    if(mUMMon) {
        DoMUMMJetFinding(particlesForJets, fLearnWeights, fSize,
                         leadmUMMindex, maxpTmUMM, tool, mUMMjets,
                         mUMMparticleWeights, mUMMweights);
    }

    // Fuzzy Jets: mTGMM --------------------
    vector<vector<double> > mTGMMparticleWeights;
    vector<double> mTGMMweights;
    vecPseudoJet mTGMMjets;
    vector<TMatrix> mTGMMjetsparams;
    int leadmTGMMindex;
    double maxpTmTGMM;
    if(mTGMMon) {
        DoMTGMMJetFinding(particlesForJets, particlesForJets,
                          fLearnWeights, false, fSize,
                          leadmTGMMindex, maxpTmTGMM, tool, mTGMMjets,
                          mTGMMparticleWeights, mTGMMjetsparams, mTGMMweights);
    }

    // Having found jets, now do a bit of data logging and analysis
    bool doWeightDistributions = true;
    if (doWeightDistributions && ievt < 10 && !batched) {

        if(mGMMon) {
            WeightDistribution(Weights, leadmGMMindex,
                               TString::Format("_mGMM_%d", ievt));
        }
        if(mTGMMon) {
            WeightDistribution(mTGMMparticleWeights, leadmTGMMindex,
                               TString::Format("_mTGMM_%d", ievt));
        }
        if(mUMMon) {
            WeightDistribution(mUMMparticleWeights, leadmUMMindex,
                               TString::Format("_mUMM_%d", ievt));
        }
        if(mGMMson) {
            WeightDistribution(mGMMsparticleWeights, leadmGMMsindex,
                               TString::Format("_mGMMs_%d", ievt));
        }
        if(mTGMMson) {
            WeightDistribution(mTGMMsparticleWeights, leadmTGMMsindex,
                               TString::Format("_mTGMMs_%d", ievt));
        }
    }

    bool doEventDisplays = false;
    if (doEventDisplays && !batched) {
        if(ievt < 10 && mGMMon) {
            tool->EventDisplay(particlesForJets,
                               myJetsLargeR_ca,tops,
                               mGMMjets,
                               Weights,
                               leadmGMMindex,
                               mGMMjetsparams,
                               TString::Format("%i",ievt));
        }
        if(mGMMon) {
            tool->NewEventDisplay(particlesForJets,
                                  myJetsLargeR_ca,tops,
                                  mGMMjets,
                                  Weights,
                                  leadmGMMindex,
                                  mGMMjetsparams,
                                  mGMMweights,
                                  TString::Format("_mGMM_%i",ievt));
            tool->JetContributionDisplay(particlesForJets,
                                         Weights,
                                         leadmGMMindex,
                                         1, TString::Format("_m_mGMM_%i", ievt));
            tool->JetContributionDisplay(particlesForJets,
                                         Weights,
                                         leadmGMMindex,
                                         0, TString::Format("_pt_mGMM_%i", ievt));

        }
        if(mTGMMon) {
            tool->NewEventDisplay(particlesForJets,
                                  myJetsLargeR_ca,tops,
                                  mTGMMjets,
                                  mTGMMparticleWeights,
                                  leadmTGMMindex,
                                  mTGMMjetsparams,
                                  mTGMMweights,
                                  TString::Format("_mTGMM_%i",ievt));
        }
        if(mUMMon) {
            tool->NewEventDisplayUniform(particlesForJets,
                                         myJetsLargeR_ca,tops,
                                         mUMMjets,
                                         mUMMparticleWeights,
                                         leadmUMMindex,
                                         mUMMweights,
                                         TString::Format("_mUMM_%i",ievt));
        }
        if(mGMMson) {
            tool->NewEventDisplay(particlesForJets,
                                  myJetsLargeR_ca, tops,
                                  mGMMsjets,
                                  mGMMsparticleWeights,
                                  leadmGMMsindex,
                                  mGMMsjetsparams,
                                  mGMMsweights,
                                  TString::Format("_mGMMs_%i", ievt));
        }
        if(mTGMMson) {
            tool->NewEventDisplay(particlesForJets,
                                  myJetsLargeR_ca,tops,
                                  mTGMMsjets,
                                  mTGMMsparticleWeights,
                                  leadmTGMMsindex,
                                  mTGMMsjetsparams,
                                  mTGMMsweights,
                                  TString::Format("_mTGMMs_%i", ievt));
        }
    }

    if(mGMMon) {
        fTmGMM_m = tool->MLpT(particlesForJets,Weights,
                              leadmGMMindex,Weights[0].size(),1);
        fTmGMM_pt = tool->MLpT(particlesForJets,Weights,
                               leadmGMMindex,Weights[0].size(),0);
        fTmGMM_ml = tool->MLlpTGaussian(particlesForJets,mGMMjets[leadmGMMindex],
                                        mGMMjetsparams[leadmGMMindex], mGMMweights[leadmGMMindex], 1);
        fTmGMM_ptl = tool->MLlpTGaussian(particlesForJets,mGMMjets[leadmGMMindex],
                                         mGMMjetsparams[leadmGMMindex], mGMMweights[leadmGMMindex], 0);
        fTmGMM_pufrac_soft = JetPuFracSoft(particlesForJets, Weights, leadmGMMindex);
        fTmGMM_pufrac_hard = JetPuFracHard(particlesForJets, Weights, leadmGMMindex);
        fTmGMM_m_pu_soft = JetPuMassSoft(particlesForJets, Weights, leadmGMMindex);
        fTmGMM_m_pu_hard = JetPuMassHard(particlesForJets, Weights, leadmGMMindex);
}
    if(mUMMon) {
        fTmUMM_m = tool->MLpT(particlesForJets, mUMMparticleWeights,
                              leadmUMMindex, mUMMparticleWeights[0].size(), 1);
        fTmUMM_pt = maxpTmUMM;
        fTmUMM_ml = tool->MLlpTUniform(particlesForJets,mUMMjets[leadmUMMindex],
                                       mUMMweights[leadmUMMindex], 1);
        fTmUMM_ptl = tool->MLlpTUniform(particlesForJets,mUMMjets[leadmUMMindex],
                                        mUMMweights[leadmUMMindex], 0);
        fTmUMM_pufrac_soft = JetPuFracSoft(particlesForJets, mUMMparticleWeights, leadmUMMindex);
        fTmUMM_pufrac_hard = JetPuFracHard(particlesForJets, mUMMparticleWeights, leadmUMMindex);
        fTmUMM_m_pu_soft = JetPuMassSoft(particlesForJets, mUMMparticleWeights, leadmUMMindex);
        fTmUMM_m_pu_hard = JetPuMassHard(particlesForJets, mUMMparticleWeights, leadmUMMindex);
    }
    if(mTGMMon) {
        fTmTGMM_m = tool->MLpT(particlesForJets, mTGMMparticleWeights,
                               leadmTGMMindex, mTGMMparticleWeights[0].size(), 1);
        fTmTGMM_pt = maxpTmTGMM;
        fTmTGMM_ml = tool->MLlpTTruncGaus(particlesForJets,mTGMMjets[leadmTGMMindex],
                                          mTGMMjetsparams[leadmTGMMindex], mTGMMweights[leadmTGMMindex],1);
        fTmTGMM_ptl = tool->MLlpTTruncGaus(particlesForJets,mTGMMjets[leadmTGMMindex],
                                           mTGMMjetsparams[leadmTGMMindex], mTGMMweights[leadmTGMMindex],0);
        fTmTGMM_pufrac_soft = JetPuFracSoft(particlesForJets, mTGMMparticleWeights, leadmTGMMindex);
        fTmTGMM_pufrac_hard = JetPuFracHard(particlesForJets, mTGMMparticleWeights, leadmTGMMindex);
        fTmTGMM_m_pu_soft = JetPuMassSoft(particlesForJets, mTGMMparticleWeights, leadmTGMMindex);
        fTmTGMM_m_pu_hard = JetPuMassHard(particlesForJets, mTGMMparticleWeights, leadmTGMMindex);
    }
    if(mTGMMson) {
        fTmTGMMs_m = tool->MLpT(particlesForJets, mTGMMsparticleWeights,
                                leadmTGMMsindex, mTGMMsparticleWeights[0].size(), 1);
        fTmTGMMs_pt = maxpTmTGMMs;
        fTmTGMMs_ml = tool->MLlpTTruncGaus(particlesForJets,mTGMMsjets[leadmTGMMsindex],
                                           mTGMMsjetsparams[leadmTGMMsindex], mTGMMsweights[leadmTGMMsindex],1);
        fTmTGMMs_ptl = tool->MLlpTTruncGaus(particlesForJets,mTGMMsjets[leadmTGMMsindex],
                                            mTGMMsjetsparams[leadmTGMMsindex], mTGMMsweights[leadmTGMMsindex],0);
        fTmTGMMs_pufrac_soft = JetPuFracSoft(particlesForJets, mTGMMsparticleWeights, leadmTGMMsindex);
        fTmTGMMs_pufrac_hard = JetPuFracHard(particlesForJets, mTGMMsparticleWeights, leadmTGMMsindex);
        fTmTGMMs_m_pu_soft = JetPuMassSoft(particlesForJets, mTGMMsparticleWeights, leadmTGMMsindex);
        fTmTGMMs_m_pu_hard = JetPuMassHard(particlesForJets, mTGMMsparticleWeights, leadmTGMMsindex);
    }
    if(mGMMson) {
        fTmGMMs_m = tool->MLpT(particlesForJets, mGMMsparticleWeights,
                               leadmGMMsindex, mGMMsparticleWeights[0].size(), 1);
        fTmGMMs_pt = maxpTmGMMs;
        fTmGMMs_ml = tool->MLlpTGaussian(particlesForJets,mGMMsjets[leadmGMMsindex],
                                         mGMMsjetsparams[leadmGMMsindex], mGMMsweights[leadmGMMsindex],1);
        fTmGMMs_ptl = tool->MLlpTGaussian(particlesForJets,mGMMsjets[leadmGMMsindex],
                                          mGMMsjetsparams[leadmGMMsindex], mGMMsweights[leadmGMMsindex],0);
        fTmGMMs_pufrac_soft = JetPuFracSoft(particlesForJets, mGMMsparticleWeights, leadmGMMsindex);
        fTmGMMs_pufrac_hard = JetPuFracHard(particlesForJets, mGMMsparticleWeights, leadmGMMsindex);
        fTmGMMs_m_pu_soft = JetPuMassSoft(particlesForJets, mGMMsparticleWeights, leadmGMMsindex);
        fTmGMMs_m_pu_hard = JetPuMassHard(particlesForJets, mGMMsparticleWeights, leadmGMMsindex);
    }

    int mytop=0;
    if (tops[1].pt()> tops[0].pt()){
        mytop=1;
    }

    fTtoppt = tops[0].pt();

    if(mGMMon) {
        fTdeltatop_mGMM = tops[0].delta_R(mGMMjets[leadmGMMindex]);

        if (tops[1].delta_R(mGMMjets[leadmGMMindex]) < fTdeltatop_mGMM){
            fTdeltatop_mGMM = tops[1].delta_R(mGMMjets[leadmGMMindex]);
            fTtoppt = tops[1].pt();
        }
        fTmGMM_m_soft = tool->SoftpT(particlesForJets,
                                     Weights,
                                     leadmGMMindex,
                                     1);
        fTmGMM_pt_soft = tool->SoftpT(particlesForJets,
                                      Weights,
                                      leadmGMMindex,
                                      0);
        fTmGMM_ucpu = UnclusteredPileupComposition(particlesForJets,
                                                   Weights);
        fTmGMM_clpu = ClusteredPileupFrac(particlesForJets,
                                          Weights);
    }
    if(mUMMon) {
        fTdeltatop_mUMM = tops[0].delta_R(mUMMjets[leadmUMMindex]);
        if (tops[1].delta_R(mUMMjets[leadmUMMindex]) <  fTdeltatop_mUMM) {
            fTdeltatop_mUMM = tops[1].delta_R(mUMMjets[leadmUMMindex]);
        }
        fTmUMM_m_soft = tool->SoftpT(particlesForJets,
                                     mUMMparticleWeights,
                                     leadmUMMindex,
                                     1);
        fTmUMM_pt_soft = tool->SoftpT(particlesForJets,
                                      mUMMparticleWeights,
                                      leadmUMMindex,
                                      0);
        fTmUMM_ucpu = UnclusteredPileupComposition(particlesForJets,
                                                   mUMMparticleWeights);
        fTmUMM_clpu = ClusteredPileupFrac(particlesForJets,
                                          mUMMparticleWeights);
    }
    if(mTGMMon) {
        fTdeltatop_mTGMM = tops[0].delta_R(mTGMMjets[leadmTGMMindex]);
        if (tops[1].delta_R(mTGMMjets[leadmTGMMindex]) < fTdeltatop_mTGMM) {
            fTdeltatop_mTGMM = tops[1].delta_R(mTGMMjets[leadmTGMMindex]);
        }
        fTmTGMM_m_soft = tool->SoftpT(particlesForJets,
                                      mTGMMparticleWeights,
                                      leadmTGMMindex,
                                      1);
        fTmTGMM_pt_soft = tool->SoftpT(particlesForJets,
                                       mTGMMparticleWeights,
                                       leadmTGMMindex,
                                       0);
        fTmTGMM_ucpu = UnclusteredPileupComposition(particlesForJets,
                                                    mTGMMparticleWeights);
        fTmTGMM_clpu = ClusteredPileupFrac(particlesForJets,
                                           mTGMMparticleWeights);
    }
    if(mTGMMson) {
        fTdeltatop_mTGMMs = tops[0].delta_R(mTGMMsjets[leadmTGMMsindex]);
        if (tops[1].delta_R(mTGMMsjets[leadmTGMMsindex]) < fTdeltatop_mTGMMs) {
            fTdeltatop_mTGMMs = tops[1].delta_R(mTGMMsjets[leadmTGMMsindex]);
        }
        fTmTGMMs_m_soft = tool->SoftpT(particlesForJets,
                                       mTGMMsparticleWeights,
                                       leadmTGMMsindex,
                                       1);
        fTmTGMMs_pt_soft = tool->SoftpT(particlesForJets,
                                        mTGMMsparticleWeights,
                                        leadmTGMMsindex,
                                        0);
        fTmTGMMs_ucpu = UnclusteredPileupComposition(particlesForJets,
                                                     mTGMMsparticleWeights);
        fTmTGMMs_clpu = ClusteredPileupFrac(particlesForJets,
                                            mTGMMsparticleWeights);
    }
    if(mGMMson) {
        fTdeltatop_mGMMs = tops[0].delta_R(mGMMsjets[leadmGMMsindex]);
        if (tops[1].delta_R(mGMMsjets[leadmGMMsindex]) < fTdeltatop_mGMMs) {
            fTdeltatop_mGMMs = tops[1].delta_R(mGMMsjets[leadmGMMsindex]);
        }
        fTmGMMs_m_soft = tool->SoftpT(particlesForJets,
                                      mGMMsparticleWeights,
                                      leadmGMMsindex,
                                      1);
        fTmGMMs_pt_soft = tool->SoftpT(particlesForJets,
                                       mGMMsparticleWeights,
                                       leadmGMMsindex,
                                       0);
        fTmGMMs_ucpu = UnclusteredPileupComposition(particlesForJets,
                                                    mGMMsparticleWeights);
        fTmGMMs_clpu = ClusteredPileupFrac(particlesForJets,
                                           mGMMsparticleWeights);
    }

    // Moments
    vector<double> moments_m;
    vector<double> moments_pt;
    double sig;
    if(mGMMon) {
        moments_m = tool->CentralMoments(particlesForJets, Weights,
                                         leadmGMMindex, 3, &totalMass);

        moments_pt = tool->CentralMoments(particlesForJets, Weights,
                                          leadmGMMindex, 3, &totalpT);

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

    if(mUMMon) {
        moments_m = tool->CentralMoments(particlesForJets, mUMMparticleWeights,
                                         leadmUMMindex, 3, &totalMass);

        moments_pt = tool->CentralMoments(particlesForJets, mUMMparticleWeights,
                                          leadmUMMindex, 3, &totalpT);

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

    if(mTGMMon) {
        moments_m = tool->CentralMoments(particlesForJets, mTGMMparticleWeights,
                                         leadmTGMMindex, 3, &totalMass);

        moments_pt = tool->CentralMoments(particlesForJets, mTGMMparticleWeights,
                                          leadmTGMMindex, 3, &totalpT);

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

    if(mTGMMson) {
        moments_m = tool->CentralMoments(particlesForJets, mTGMMsparticleWeights,
                                         leadmTGMMsindex, 3, &totalMass);

        moments_pt = tool->CentralMoments(particlesForJets, mTGMMsparticleWeights,
                                          leadmTGMMsindex, 3, &totalpT);

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

    if(mGMMson) {
        moments_m = tool->CentralMoments(particlesForJets, mGMMsparticleWeights,
                                         leadmGMMsindex, 3, &totalMass);

        moments_pt = tool->CentralMoments(particlesForJets, mGMMsparticleWeights,
                                          leadmGMMsindex, 3, &totalpT);

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

    tT->Fill();

    if(fDebug) cout << "FuzzyAnalysis::AnalyzeEvent End " << endl;
    return;
}



// declate branches
void FuzzyAnalysis::DeclareBranches(){

    // Event Properties
    tT->Branch("EventNumber",    &fTEventNumber,    "EventNumber/I");
    tT->Branch("NPV",            &fTNPV,            "NPV/I");
    tT->Branch("CA_m",           &fTCA_m,           "CA_m/F");
    tT->Branch("CA_pt",          &fTCA_pt,          "CA_pt/F");
    tT->Branch("CA_pufrac",      &fTCA_pufrac,      "CA_pufrac/F");
    tT->Branch("CA_m_pu",        &fTCA_m_pu,        "CA_m_pu/F");

    tT->Branch("toppt",          &fTtoppt,          "toppt/F");

    tT->Branch("antikt_m",       &fTantikt_m,       "antikt_m/F");
    tT->Branch("antikt_pt",      &fTantikt_pt,      "antikt_pt/F");
    tT->Branch("antikt_m_trimmed_two", &fTantikt_m_trimmed_two, "antikt_m_trimmed_two/F");
    tT->Branch("antikt_pt_trimmed_two", &fTantikt_pt_trimmed_two, "antikt_pt_trimmed_two/F");
    tT->Branch("antikt_m_trimmed_three", &fTantikt_m_trimmed_three, "antikt_m_trimmed_three/F");
    tT->Branch("antikt_pt_trimmed_three", &fTantikt_pt_trimmed_three, "antikt_pt_trimmed_three/F");
    tT->Branch("antikt_pufrac_trimmed_two", &fTantikt_pufrac_trimmed_two, "antikt_pufrac_trimmed_two/F");
    tT->Branch("antikt_pufrac_trimmed_three", &fTantikt_pufrac_trimmed_three, "antikt_pufrac_trimmed_three/F");
    tT->Branch("antikt_m_pu", &fTantikt_m_pu, "antikt_m_pu/F");
    tT->Branch("antikt_m_pu_trimmed_two", &fTantikt_m_pu_trimmed_two, "antikt_m_pu_trimmed_two/F");
    tT->Branch("antikt_m_pu_trimmed_three", &fTantikt_m_pu_trimmed_three, "antikt_m_pu_trimmed_three/F");

    tT->Branch("mGMMs_m",        &fTmGMMs_m,        "mGMMs_m/F");
    tT->Branch("mGMMs_pt",       &fTmGMMs_pt,       "mGMMs_pt/F");
    tT->Branch("deltatop_mGMMs", &fTdeltatop_mGMMs, "deltatop_mGMMs/F");
    tT->Branch("mGMMs_m_mean",   &fTmGMMs_m_mean,   "mGMMs_m_mean/F");
    tT->Branch("mGMMs_m_var",    &fTmGMMs_m_var,    "mGMMs_m_var/F");
    tT->Branch("mGMMs_m_skew",   &fTmGMMs_m_skew,   "mGMMs_m_skew/F");
    tT->Branch("mGMMs_pt_mean",  &fTmGMMs_pt_mean,  "mGMMs_pt_mean/F");
    tT->Branch("mGMMs_pt_var",   &fTmGMMs_pt_var,   "mGMMs_pt_var/F");
    tT->Branch("mGMMs_pt_skew",  &fTmGMMs_pt_skew,  "mGMMs_pt_skew/F");
    tT->Branch("mGMMs_ptl",      &fTmGMMs_ptl,      "mGMMs_ptl/F");
    tT->Branch("mGMMs_ml",       &fTmGMMs_ml,       "mGMMs_ml/F");
    tT->Branch("mGMMs_m_soft",    &fTmGMMs_m_soft,    "mGMMs_m_soft/F");
    tT->Branch("mGMMs_pt_soft",   &fTmGMMs_pt_soft,   "mGMMs_pt_soft/F");
    tT->Branch("mGMMs_ucpu",     &fTmGMMs_ucpu,      "mGMMs_ucpu/F");
    tT->Branch("mGMMs_clpu",     &fTmGMMs_clpu,      "mGMMs_clpu/F");
    tT->Branch("mGMMs_pufrac_soft", &fTmGMMs_pufrac_soft, "mGMMs_pufrac_soft/F");
    tT->Branch("mGMMs_pufrac_hard", &fTmGMMs_pufrac_hard, "mGMMs_pufrac_hard/F");
    tT->Branch("mGMMs_m_pu_soft", &fTmGMMs_m_pu_soft, "mGMMs_m_pu_soft/F");
    tT->Branch("mGMMs_m_pu_hard", &fTmGMMs_m_pu_hard, "mGMMs_m_pu_hard/F");

    tT->Branch("mTGMMs_m",       &fTmTGMMs_m,       "mTGMMs_m/F");
    tT->Branch("mTGMMs_pt",      &fTmTGMMs_pt,      "mTGMMs_pt/F");
    tT->Branch("deltatop_mTGMMs",&fTdeltatop_mTGMMs,"deltatop_mTGMMs/F");
    tT->Branch("mTGMMs_m_mean",  &fTmTGMMs_m_mean,  "mTGMMs_m_mean/F");
    tT->Branch("mTGMMs_m_var",   &fTmTGMMs_m_var,   "mTGMMs_m_var/F");
    tT->Branch("mTGMMs_m_skew",  &fTmTGMMs_m_skew,  "mTGMMs_m_skew/F");
    tT->Branch("mTGMMs_pt_mean", &fTmTGMMs_pt_mean, "mTGMMs_pt_mean/F");
    tT->Branch("mTGMMs_pt_var",  &fTmTGMMs_pt_var,  "mTGMMs_pt_var/F");
    tT->Branch("mTGMMs_pt_skew", &fTmTGMMs_pt_skew, "mTGMMs_pt_skew/F");
    tT->Branch("mTGMMs_ptl",     &fTmTGMMs_ptl,     "mTGMMs_ptl/F");
    tT->Branch("mTGMMs_ml",      &fTmTGMMs_ml,      "mTGMMs_ml/F");
    tT->Branch("mTGMMs_m_soft",    &fTmTGMMs_m_soft,    "mTGMMs_m_soft/F");
    tT->Branch("mTGMMs_pt_soft",   &fTmTGMMs_pt_soft,   "mTGMMs_pt_soft/F");
    tT->Branch("mTGMMs_ucpu",     &fTmTGMMs_ucpu,      "mTGMMs_ucpu/F");
    tT->Branch("mTGMMs_clpu",     &fTmTGMMs_clpu,      "mTGMMs_clpu/F");
    tT->Branch("mTGMMs_pufrac_soft", &fTmTGMMs_pufrac_soft, "mTGMMs_pufrac_soft/F");
    tT->Branch("mTGMMs_pufrac_hard", &fTmTGMMs_pufrac_hard, "mTGMMs_pufrac_hard/F");
    tT->Branch("mTGMMs_m_pu_soft", &fTmTGMMs_m_pu_soft, "mTGMMs_m_pu_soft/F");
    tT->Branch("mTGMMs_m_pu_hard", &fTmTGMMs_m_pu_hard, "mTGMMs_m_pu_hard/F");

    tT->Branch("mUMM_m",         &fTmUMM_m,         "mUMM_m/F");
    tT->Branch("mUMM_pt",        &fTmUMM_pt,        "mUMM_pt/F");
    tT->Branch("deltatop_mUMM",  &fTdeltatop_mUMM,  "deltatop_mUMM/F");
    tT->Branch("mUMM_m_mean",    &fTmUMM_m_mean,    "mUMM_m_mean/F");
    tT->Branch("mUMM_m_var",     &fTmUMM_m_var,     "mUMM_m_var/F");
    tT->Branch("mUMM_m_skew",    &fTmUMM_m_skew,    "mUMM_m_skew/F");
    tT->Branch("mUMM_pt_mean",   &fTmUMM_pt_mean,   "mUMM_pt_mean/F");
    tT->Branch("mUMM_pt_var",    &fTmUMM_pt_var,    "mUMM_pt_var/F");
    tT->Branch("mUMM_pt_skew",   &fTmUMM_pt_skew,   "mUMM_pt_skew/F");
    tT->Branch("mUMM_ptl",       &fTmUMM_ptl,       "mUMM_ptl/F");
    tT->Branch("mUMM_ml",        &fTmUMM_ml,        "mUMM_ml/F");
    tT->Branch("mUMM_m_soft",    &fTmUMM_m_soft,    "mUMM_m_soft/F");
    tT->Branch("mUMM_pt_soft",   &fTmUMM_pt_soft,   "mUMM_pt_soft/F");
    tT->Branch("mUMM_ucpu",     &fTmUMM_ucpu,      "mUMM_ucpu/F");
    tT->Branch("mUMM_clpu",     &fTmUMM_clpu,      "mUMM_clpu/F");
    tT->Branch("mUMM_pufrac_soft", &fTmUMM_pufrac_soft, "mUMM_pufrac_soft/F");
    tT->Branch("mUMM_pufrac_hard", &fTmUMM_pufrac_hard, "mUMM_pufrac_hard/F");
    tT->Branch("mUMM_m_pu_soft", &fTmUMM_m_pu_soft, "mUMM_m_pu_soft/F");
    tT->Branch("mUMM_m_pu_hard", &fTmUMM_m_pu_hard, "mUMM_m_pu_hard/F");

    tT->Branch("mGMM_m",         &fTmGMM_m,         "mGMM_m/F");
    tT->Branch("mGMM_pt",        &fTmGMM_pt,        "mGMM_pt/F");
    tT->Branch("deltatop_mGMM",  &fTdeltatop_mGMM,  "deltatop_mGMM/F");
    tT->Branch("mGMM_m_mean",    &fTmGMM_m_mean,    "mGMM_m_mean/F");
    tT->Branch("mGMM_m_var",     &fTmGMM_m_var,     "mGMM_m_var/F");
    tT->Branch("mGMM_m_skew",    &fTmGMM_m_skew,    "mGMM_m_skew/F");
    tT->Branch("mGMM_pt_mean",   &fTmGMM_pt_mean,   "mGMM_pt_mean/F");
    tT->Branch("mGMM_pt_var",    &fTmGMM_pt_var,    "mGMM_pt_var/F");
    tT->Branch("mGMM_pt_skew",   &fTmGMM_pt_skew,   "mGMM_pt_skew/F");
    tT->Branch("mGMM_ptl",       &fTmGMM_ptl,       "mGMM_ptl/F");
    tT->Branch("mGMM_ml",        &fTmGMM_ml,        "mGMM_ml/F");
    tT->Branch("mGMM_m_soft",    &fTmUMM_m_soft,    "mGMM_m_soft/F");
    tT->Branch("mGMM_pt_soft",   &fTmUMM_pt_soft,   "mGMM_pt_soft/F");
    tT->Branch("mGMM_ucpu",     &fTmGMM_ucpu,      "mGMM_ucpu/F");
    tT->Branch("mGMM_clpu",     &fTmGMM_clpu,      "mGMM_clpu/F");
    tT->Branch("mGMM_pufrac_soft", &fTmGMM_pufrac_soft, "mGMM_pufrac_soft/F");
    tT->Branch("mGMM_pufrac_hard", &fTmGMM_pufrac_hard, "mGMM_pufrac_hard/F");
    tT->Branch("mGMM_m_pu_soft", &fTmGMM_m_pu_soft, "mGMM_m_pu_soft/F");
    tT->Branch("mGMM_m_pu_hard", &fTmGMM_m_pu_hard, "mGMM_m_pu_hard/F");

    tT->Branch("mTGMM_m",        &fTmTGMM_m,        "mTGMM_m/F");
    tT->Branch("mTGMM_pt",       &fTmTGMM_pt,       "mTGMM_pt/F");
    tT->Branch("deltatop_mTGMM", &fTdeltatop_mTGMM, "deltatop_mTGMM/F");
    tT->Branch("mTGMM_m_mean",   &fTmTGMM_m_mean,   "mTGMM_m_mean/F");
    tT->Branch("mTGMM_m_var",    &fTmTGMM_m_var,    "mTGMM_m_var/F");
    tT->Branch("mTGMM_m_skew",   &fTmTGMM_m_skew,   "mTGMM_m_skew/F");
    tT->Branch("mTGMM_pt_mean",  &fTmTGMM_pt_mean,  "mTGMM_pt_mean/F");
    tT->Branch("mTGMM_pt_var",   &fTmTGMM_pt_var,   "mTGMM_pt_var/F");
    tT->Branch("mTGMM_pt_skew",  &fTmTGMM_pt_skew,  "mTGMM_pt_skew/F");
    tT->Branch("mTGMM_ptl",      &fTmTGMM_ptl,      "mTGMM_ptl/F");
    tT->Branch("mTGMM_ml",       &fTmTGMM_ml,       "mTGMM_ml/F");
    tT->Branch("mTGMM_m_soft",    &fTmTGMM_m_soft,    "mTGMM_m_soft/F");
    tT->Branch("mTGMM_pt_soft",   &fTmTGMM_pt_soft,   "mTGMM_pt_soft/F");
    tT->Branch("mTGMM_ucpu",     &fTmTGMM_ucpu,      "mTGMM_ucpu/F");
    tT->Branch("mTGMM_clpu",     &fTmTGMM_clpu,      "mTGMM_clpu/F");
    tT->Branch("mTGMM_pufrac_soft", &fTmTGMM_pufrac_soft, "mTGMM_pufrac_soft/F");
    tT->Branch("mTGMM_pufrac_hard", &fTmTGMM_pufrac_hard, "mTGMM_pufrac_hard/F");
    tT->Branch("mTGMM_m_pu_soft", &fTmTGMM_m_pu_soft, "mTGMM_m_pu_soft/F");
    tT->Branch("mTGMM_m_pu_hard", &fTmTGMM_m_pu_hard, "mTGMM_m_pu_hard/F");


    tT->GetListOfBranches()->ls();

    return;
}


// resets vars
void FuzzyAnalysis::ResetBranches(){
    // reset branches
    fTEventNumber     = -999;
    fTNPV = -999;
    fTCA_m            = -1.;
    fTCA_pt           = -1.;
    fTCA_pufrac = -1;
    fTCA_m_pu = -1;

    fTtoppt           = -1.;

    fTantikt_m = -1;
    fTantikt_pt = -1;
    fTantikt_m_trimmed_two = -1;
    fTantikt_pt_trimmed_two = -1;
    fTantikt_m_trimmed_three = -1;
    fTantikt_pt_trimmed_three = -1;
    fTantikt_pufrac_trimmed_two = -1;
    fTantikt_pufrac_trimmed_three = -1;
    fTantikt_m_pu = -1;
    fTantikt_m_pu_trimmed_two = -1;
    fTantikt_m_pu_trimmed_three = -1;

    fTmGMMs_m         = -1.;
    fTmGMMs_pt        = -1.;
    fTdeltatop_mGMMs  = -1.;
    fTmGMMs_m_mean    = -1.;
    fTmGMMs_m_var     = -1.;
    fTmGMMs_m_skew    = -99999;
    fTmGMMs_pt_mean   = -1.;
    fTmGMMs_pt_var    = -1.;
    fTmGMMs_pt_skew   = -99999;
    fTmGMMs_ptl       = -1.;
    fTmGMMs_ml        = -1.;
    fTmGMMs_m_soft = -1;
    fTmGMMs_pt_soft = -1;
    fTmGMMs_ucpu = -1;
    fTmGMMs_clpu = -1;
    fTmGMMs_pufrac_soft = -1;
    fTmGMMs_pufrac_hard = -1;
    fTmGMMs_m_pu_soft = -1;
    fTmGMMs_m_pu_hard = -1;

    fTmTGMMs_m        = -1.;
    fTmTGMMs_pt       = -1.;
    fTdeltatop_mTGMMs = -1.;
    fTmTGMMs_m_mean   = -1.;
    fTmTGMMs_m_var    = -1.;
    fTmTGMMs_m_skew   = -99999;
    fTmTGMMs_pt_mean  = -1.;
    fTmTGMMs_pt_var   = -1.;
    fTmTGMMs_pt_skew  = -99999;
    fTmTGMMs_ptl      = -1.;
    fTmTGMMs_ml = -1.;
    fTmTGMMs_m_soft = -1;
    fTmTGMMs_pt_soft = -1;
    fTmTGMMs_ucpu = -1;
    fTmTGMMs_clpu = -1;
    fTmTGMMs_pufrac_soft = -1;
    fTmTGMMs_pufrac_hard = -1;
    fTmTGMMs_m_pu_soft = -1;
    fTmTGMMs_m_pu_hard = -1;

    fTmGMM_m         = -1.;
    fTmGMM_pt        = -1.;
    fTdeltatop_mGMM  = -1.;
    fTmGMM_m_mean    = -1.;
    fTmGMM_m_var     = -1.;
    fTmGMM_m_skew    = -99999;
    fTmGMM_pt_mean    = -1.;
    fTmGMM_pt_var     = -1.;
    fTmGMM_pt_skew    = -99999;
    fTmGMM_ptl = -1.;
    fTmGMM_ml = -1.;
    fTmGMM_m_soft = -1;
    fTmGMM_pt_soft = -1;
    fTmGMM_ucpu = -1;
    fTmGMM_clpu = -1;
    fTmGMM_pufrac_soft = -1;
    fTmGMM_pufrac_hard = -1;
    fTmGMM_m_pu_soft = -1;
    fTmGMM_m_pu_hard = -1;

    fTmUMM_m         = -1.;
    fTmUMM_pt        = -1.;
    fTdeltatop_mUMM  = -1.;
    fTmUMM_m_mean    = -1.;
    fTmUMM_m_var     = -1.;
    fTmUMM_m_skew    = -99999;
    fTmUMM_pt_mean    = -1.;
    fTmUMM_pt_var     = -1.;
    fTmUMM_pt_skew    = -99999;
    fTmUMM_ptl = -1.;
    fTmUMM_ml = -1.;
    fTmUMM_m_soft = -1;
    fTmUMM_pt_soft = -1;
    fTmUMM_ucpu = -1;
    fTmUMM_clpu = -1;
    fTmUMM_pufrac_soft = -1;
    fTmUMM_pufrac_hard = -1;
    fTmUMM_m_pu_soft = -1;
    fTmUMM_m_pu_hard = -1;

    fTmTGMM_m        = -1.;
    fTmTGMM_pt       = -1.;
    fTdeltatop_mTGMM = -1.;
    fTmTGMM_m_mean    = -1.;
    fTmTGMM_m_var     = -1.;
    fTmTGMM_m_skew    = -99999;
    fTmTGMM_pt_mean    = -1.;
    fTmTGMM_pt_var     = -1.;
    fTmTGMM_pt_skew    = -99999;
    fTmTGMM_ptl = -1.;
    fTmTGMM_ml = -1.;
    fTmTGMM_m_soft = -1;
    fTmTGMM_pt_soft = -1;
    fTmTGMM_ucpu = -1;
    fTmTGMM_clpu = -1;
    fTmTGMM_pufrac_soft = -1;
    fTmTGMM_pufrac_hard = -1;
    fTmTGMM_m_pu_soft = -1;
    fTmTGMM_m_pu_hard = -1;
}
