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
    void DoFastJetFinding(vecPseudoJet particles,
                          fastjet::JetDefinition *pJetDef,
                          double size,
                          vecPseudoJet& jets) {
        fastjet::ClusterSequence csLargeR_ca(particles, *pJetDef);
        jets = fastjet::sorted_by_pt(csLargeR_ca.inclusive_jets(size));
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
void FuzzyAnalysis::AnalyzeEvent(int ievt, Pythia8::Pythia* pythia8){
    if(fDebug) cout << "FuzzyAnalysis::AnalyzeEvent Begin " << endl;

    // generate a new event
    if (!pythia8->next()) return;
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


    // Particle loop -----------------------------------------------------------
    // The Pythia event listing contains a lot more than we want to process,
    // we prune out certain particles (muons / neutrinos) and only add final
    // state particles

    double px, py, pz, e;
    for (unsigned int ip=0; ip < (unsigned) pythia8->event.size(); ++ip){
        px = pythia8->event[ip].px();
        py = pythia8->event[ip].py();
        pz = pythia8->event[ip].pz();
        e  = pythia8->event[ip].e();
        fastjet::PseudoJet p(px, py, pz, e);
        p.set_user_info(new MyUserInfo(pythia8->event[ip].id(),ip,pythia8->event[ip].charge()));

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
    vecPseudoJet myJetsLargeR_ca;
    DoFastJetFinding(particlesForJets, m_jet_def_largeR_ca, 1.0, myJetsLargeR_ca);
    fTCA_m = myJetsLargeR_ca[0].m();
    fTCA_pt = myJetsLargeR_ca[0].pt();

    // Various mixture models ---------------
    // ======================================
    tool->SetMergeDistance(0.05);
    double size = 0.8;
    bool learnWeights = true;

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
                         learnWeights, true,
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
                          learnWeights, true,
                          size, leadmTGMMsindex, maxpTmTGMMs,
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
                         learnWeights, false,
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
        DoMUMMJetFinding(particlesForJets, learnWeights, size,
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
                          learnWeights, false, size,
                          leadmTGMMindex, maxpTmTGMM, tool, mTGMMjets,
                          mTGMMparticleWeights, mTGMMjetsparams, mTGMMweights);
    }

    // Having found jets, now do a bit of data logging and analysis
    bool doEventDisplays = true;
    if (doEventDisplays) {
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
    }
    if(mUMMon) {
        fTmUMM_m = tool->MLpT(particlesForJets, mUMMparticleWeights,
                              leadmUMMindex, mUMMparticleWeights[0].size(), 1);
        fTmUMM_pt = maxpTmUMM;
        fTmUMM_ml = tool->MLlpTUniform(particlesForJets,mUMMjets[leadmUMMindex],
                                 mUMMweights[leadmUMMindex], 1);
        fTmUMM_ptl = tool->MLlpTUniform(particlesForJets,mUMMjets[leadmUMMindex],
                                   mUMMweights[leadmUMMindex], 0);
    }
    if(mTGMMon) {
        fTmTGMM_m = tool->MLpT(particlesForJets, mTGMMparticleWeights,
                               leadmTGMMindex, mTGMMparticleWeights[0].size(), 1);
        fTmTGMM_pt = maxpTmTGMM;
        fTmTGMM_ml = tool->MLlpTTruncGaus(particlesForJets,mTGMMjets[leadmTGMMindex],
                                  mTGMMjetsparams[leadmTGMMindex], mTGMMweights[leadmTGMMindex],1);
        fTmTGMM_ptl = tool->MLlpTTruncGaus(particlesForJets,mTGMMjets[leadmTGMMindex],
                                   mTGMMjetsparams[leadmTGMMindex], mTGMMweights[leadmTGMMindex],0);
    }
    if(mTGMMson) {
        fTmTGMMs_m = tool->MLpT(particlesForJets, mTGMMsparticleWeights,
                                leadmTGMMsindex, mTGMMsparticleWeights[0].size(), 1);
        fTmTGMMs_pt = maxpTmTGMMs;
        fTmTGMM_ml = tool->MLlpTTruncGaus(particlesForJets,mTGMMjets[leadmTGMMindex],
                                  mTGMMjetsparams[leadmTGMMindex], mTGMMweights[leadmTGMMindex],1);
        fTmTGMM_ptl = tool->MLlpTTruncGaus(particlesForJets,mTGMMjets[leadmTGMMindex],
                                   mTGMMjetsparams[leadmTGMMindex], mTGMMweights[leadmTGMMindex],0);
    }
    if(mGMMson) {
        fTmGMMs_m = tool->MLpT(particlesForJets, mGMMsparticleWeights,
                               leadmGMMsindex, mGMMsparticleWeights[0].size(), 1);
        fTmGMMs_pt = maxpTmGMMs;
        fTmGMMs_ml = tool->MLlpTGaussian(particlesForJets,mGMMsjets[leadmGMMsindex],
                                  mGMMsjetsparams[leadmGMMsindex], mGMMsweights[leadmGMMsindex],1);
        fTmGMMs_ptl = tool->MLlpTGaussian(particlesForJets,mGMMsjets[leadmGMMsindex],
                                   mGMMsjetsparams[leadmGMMsindex], mGMMsweights[leadmGMMsindex],0);

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
    }
    if(mUMMon) {
        fTdeltatop_mUMM = tops[0].delta_R(mUMMjets[leadmUMMindex]);
        if (tops[1].delta_R(mUMMjets[leadmUMMindex]) <  fTdeltatop_mUMM) {
            fTdeltatop_mUMM = tops[1].delta_R(mUMMjets[leadmUMMindex]);
        }
    }
    if(mTGMMon) {
        fTdeltatop_mTGMM = tops[0].delta_R(mTGMMjets[leadmTGMMindex]);
        if (tops[1].delta_R(mTGMMjets[leadmTGMMindex]) < fTdeltatop_mTGMM) {
            fTdeltatop_mTGMM = tops[1].delta_R(mTGMMjets[leadmTGMMindex]);
        }
    }
    if(mTGMMson) {
        fTdeltatop_mTGMMs = tops[0].delta_R(mTGMMsjets[leadmTGMMsindex]);
        if (tops[1].delta_R(mTGMMsjets[leadmTGMMsindex]) < fTdeltatop_mTGMMs) {
            fTdeltatop_mTGMMs = tops[1].delta_R(mTGMMsjets[leadmTGMMsindex]);
        }
    }
    if(mGMMson) {
        fTdeltatop_mGMMs = tops[0].delta_R(mGMMsjets[leadmGMMsindex]);
        if (tops[1].delta_R(mGMMsjets[leadmGMMsindex]) < fTdeltatop_mGMMs) {
            fTdeltatop_mGMMs = tops[1].delta_R(mGMMsjets[leadmGMMsindex]);
        }
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
    tT->Branch("CA_m",           &fTCA_m,           "CA_m/F");
    tT->Branch("CA_pt",          &fTCA_pt,          "CA_pt/F");
    tT->Branch("toppt",          &fTtoppt,          "toppt/F");

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


    tT->Branch("mTGMMs_m",       &fTmTGMMs_m,       "mTGMMs_m/F");
    tT->Branch("mTGMMs_pt",      &fTmTGMMs_pt,      "mTGMMs_pt/F");
    tT->Branch("deltatop_mTGMMs",&fTdeltatop_mTGMMs,"deltatop_mTGMMs/F");
    tT->Branch("mTGMMs_m_mean",  &fTmTGMMs_m_mean,  "mTGMMs_m_mean/F");
    tT->Branch("mTGMMs_m_var",   &fTmTGMMs_m_var,   "mTGMMs_m_var/F");
    tT->Branch("mTGMMs_m_skew",  &fTmTGMMs_m_skew,  "mTGMMs_m_skew/F");
    tT->Branch("mTGMMs_pt_mean", &fTmTGMMs_pt_mean, "mTGMMs_pt_mean/F");
    tT->Branch("mTGMMs_pt_var",  &fTmTGMMs_pt_var,  "mTGMMs_pt_var/F");
    tT->Branch("mTGMMs_pt_skew", &fTmTGMMs_pt_skew, "mTGMMs_pt_skew/F");
    tT->Branch("mTGMMs_ptl",     &fTmTGMMs_ptl,     "mTGMMs_ml/F");


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


    tT->GetListOfBranches()->ls();

    return;
}


// resets vars
void FuzzyAnalysis::ResetBranches(){
    // reset branches
    fTEventNumber     = -999;
    fTCA_m            = -1.;
    fTCA_pt           = -1.;
    fTtoppt           = -1.;

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
}
