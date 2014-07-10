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
        tool->SetR(size);
        jets = tool->ClusterFuzzyUniform(parts,
                                         &particleWeights,
                                         &jetWeights);
        FindLeadingJet(particles, jets, particleWeights, tool, leadIndex, leadpT);
    }

    void DoMTGMMJetFinding(vecPseudoJet& particles,
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
        vecPseudoJet parts = particles;
        tool->SetKernelType(FuzzyTools::TRUNCGAUSSIAN);
        tool->SetSeeds(parts);
        tool->SetLearnWeights(learnWeights);
        tool->SetLearnShape(learnShape);
        tool->SetR(size);
        jets = tool->ClusterFuzzyTruncGaus(parts,
                                           &particleWeights,
                                           &parameters,
                                           &jetWeights);
        FindLeadingJet(particles, jets, particleWeights, tool, leadIndex, leadpT);
    }

    void DoMGMMJetFinding(vecPseudoJet& particles,
                          bool learnWeights,
                          bool learnShape,
                          int& leadIndex,
                          double& leadpT,
                          FuzzyTools *tool,
                          vecPseudoJet& jets,
                          vector<vector<double> >& particleWeights,
                          vector<TMatrix>& parameters,
                          vector<double>& jetWeights) {
        vecPseudoJet parts = particles;
        tool->SetKernelType(FuzzyTools::GAUSSIAN);
        tool->SetSeeds(parts);
        tool->SetLearnWeights(learnWeights);
        tool->SetLearnShape(learnShape);
        jets = tool->ClusterFuzzyGaussian(parts,
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
    tool->SetMergeDistance(0.05);
    // Fuzzy Jets: mGMM ---------------------
    vector<vector<double> > Weights;
    vector<TMatrix> mGMMjetsparams;
    vector<double> mGMMweights;
    vector<fastjet::PseudoJet> mGMMjets;
    int leadmGMMindex;
    double maxpTmGMM;
    bool learnWeights = false;
    bool learnShape = false;
    DoMGMMJetFinding(particlesForJets, learnWeights, learnShape,
                     leadmGMMindex, maxpTmGMM,
                     tool, mGMMjets, Weights,
                     mGMMjetsparams, mGMMweights);


    int learnedClusterCount = Weights[0].size();
    if (learnedClusterCount == 0) {
        cout << "At event " << ievt << " returned no clusters. Continuing." << endl;
        return;
    }

    // Fuzzy Jets: mUMM ---------------------
    vector<vector<double> > mUMMparticleWeights;
    vector<double> mUMMweights;
    vecPseudoJet mUMMjets;

    int leadmUMMindex;
    double maxpTmUMM;
    double size = 1.0;
    DoMUMMJetFinding(particlesForJets, learnWeights, size,
                     leadmUMMindex, maxpTmUMM, tool, mUMMjets,
                     mUMMparticleWeights, mUMMweights);

    // Fuzzy Jets: mTGMM --------------------
    vector<vector<double> > mTGMMparticleWeights;
    vector<double> mTGMMweights;
    vecPseudoJet mTGMMjets;
    vector<TMatrix> mTGMMjetsparams;
    int leadmTGMMindex;
    double maxpTmTGMM;
    DoMTGMMJetFinding(particlesForJets, learnWeights, learnShape, size,
                      leadmTGMMindex, maxpTmTGMM, tool, mTGMMjets,
                      mTGMMparticleWeights, mTGMMjetsparams, mTGMMweights);


    // Having found jets, now do a bit of data logging and analysis
    bool doEventDisplays = true;
    if (doEventDisplays) {
        if(ievt < 10) {
            tool->EventDisplay(particlesForJets,
                               myJetsLargeR_ca,tops,
                               mGMMjets,
                               Weights,
                               leadmGMMindex,
                               mGMMjetsparams,
                               TString::Format("%i",ievt));
        }
        tool->NewEventDisplay(particlesForJets,
                              myJetsLargeR_ca,tops,
                              mGMMjets,
                              Weights,
                              leadmGMMindex,
                              mGMMjetsparams,
                              mGMMweights,
                              TString::Format("_mGMM_%i",ievt));

        tool->NewEventDisplay(particlesForJets,
                              myJetsLargeR_ca,tops,
                              mTGMMjets,
                              mTGMMparticleWeights,
                              leadmTGMMindex,
                              mTGMMjetsparams,
                              mTGMMweights,
                              TString::Format("_mTGMM_%i",ievt));

        tool->NewEventDisplayUniform(particlesForJets,
                                     myJetsLargeR_ca,tops,
                                     mUMMjets,
                                     mUMMparticleWeights,
                                     leadmUMMindex,
                                     mUMMweights,
                                     TString::Format("_mUMM_%i",ievt));
    }


    fTmGMM_m = tool->MLpT(particlesForJets,Weights,
                          leadmGMMindex,Weights[0].size(),1);
    fTmGMM_pt = tool->MLpT(particlesForJets,Weights,
                           leadmGMMindex,Weights[0].size(),0);

    fTmUMM_m = tool->MLpT(particlesForJets, mUMMparticleWeights,
                          leadmUMMindex, mUMMparticleWeights[0].size(), 1);
    fTmUMM_pt = maxpTmUMM;

    fTmTGMM_m = tool->MLpT(particlesForJets, mTGMMparticleWeights,
                           leadmTGMMindex, mTGMMparticleWeights[0].size(), 1);
    fTmTGMM_pt = maxpTmTGMM;

    int mytop=0;
    if (tops[1].pt()> tops[0].pt()){
        mytop=1;
    }

    fTtoppt = tops[0].pt();
    fTdeltatop_mGMM = tops[0].delta_R(mGMMjets[leadmGMMindex]);

    if (tops[1].delta_R(mGMMjets[leadmGMMindex]) < fTdeltatop_mGMM){
        fTdeltatop_mGMM = tops[1].delta_R(mGMMjets[leadmGMMindex]);
        fTtoppt = tops[1].pt();
    }

    fTdeltatop_mUMM = tops[0].delta_R(mUMMjets[leadmUMMindex]);
    if (tops[1].delta_R(mUMMjets[leadmUMMindex]) <  fTdeltatop_mUMM) {
        fTdeltatop_mUMM = tops[1].delta_R(mUMMjets[leadmUMMindex]);
    }

    fTdeltatop_mTGMM = tops[0].delta_R(mTGMMjets[leadmTGMMindex]);
    if (tops[1].delta_R(mTGMMjets[leadmTGMMindex]) < fTdeltatop_mTGMM) {
        fTdeltatop_mTGMM = tops[1].delta_R(mTGMMjets[leadmTGMMindex]);
    }

    // Moments
    vector<double> moments_m;
    vector<double> moments_pt;
    moments_m = tool->CentralMoments(particlesForJets, Weights,
                                     leadmGMMindex, 3, &totalMass);

    moments_pt = tool->CentralMoments(particlesForJets, Weights,
                                      leadmGMMindex, 3, &totalpT);

    double sig;
    fTmGMM_m_mean = moments_m[0];
    fTmGMM_m_var  = moments_m[1];
    sig = sqrt(fTmGMM_m_var);
    fTmGMM_m_skew = moments_m[2];
    fTmGMM_pt_mean = moments_pt[0];
    fTmGMM_pt_var = moments_pt[1];
    sig = sqrt(fTmGMM_m_var);
    fTmGMM_pt_skew = moments_pt[2];

    moments_m.clear();
    moments_pt.clear();

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

    moments_m.clear();
    moments_pt.clear();

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

    tT->Branch("mTGMMs_m",       &fTmTGMMs_m,       "mTGMMs_m/F");
    tT->Branch("mTGMMs_pt",      &fTmTGMMs_pt,      "mTGMMs_pt/F");
    tT->Branch("deltatop_mTGMMs",&fTdeltatop_mTGMMs,"deltatop_mTGMMs/F");
    tT->Branch("mTGMMs_m_mean",  &fTmTGMMs_m_mean,  "mTGMMs_m_mean/F");
    tT->Branch("mTGMMs_m_var",   &fTmTGMMs_m_var,   "mTGMMs_m_var/F");
    tT->Branch("mTGMMs_m_skew",  &fTmTGMMs_m_skew,  "mTGMMs_m_skew/F");
    tT->Branch("mTGMMs_pt_mean", &fTmTGMMs_pt_mean, "mTGMMs_pt_mean/F");
    tT->Branch("mTGMMs_pt_var",  &fTmTGMMs_pt_var,  "mTGMMs_pt_var/F");
    tT->Branch("mTGMMs_pt_skew", &fTmTGMMs_pt_skew, "mTGMMs_pt_skew/F");


    tT->Branch("mUMM_m",         &fTmUMM_m,         "mUMM_m/F");
    tT->Branch("mUMM_pt",        &fTmUMM_pt,        "mUMM_pt/F");
    tT->Branch("deltatop_mUMM",  &fTdeltatop_mUMM,  "deltatop_mUMM/F");
    tT->Branch("mUMM_m_mean",    &fTmUMM_m_mean,    "mUMM_m_mean/F");
    tT->Branch("mUMM_m_var",     &fTmUMM_m_var,     "mUMM_m_var/F");
    tT->Branch("mUMM_m_skew",    &fTmUMM_m_skew,    "mUMM_m_skew/F");
    tT->Branch("mUMM_pt_mean",   &fTmUMM_pt_mean,   "mUMM_pt_mean/F");
    tT->Branch("mUMM_pt_var",    &fTmUMM_pt_var,    "mUMM_pt_var/F");
    tT->Branch("mUMM_pt_skew",   &fTmUMM_pt_skew,   "mUMM_pt_skew/F");


    tT->Branch("mGMM_m",         &fTmGMM_m,         "mGMM_m/F");
    tT->Branch("mGMM_pt",        &fTmGMM_pt,        "mGMM_pt/F");
    tT->Branch("deltatop_mGMM",  &fTdeltatop_mGMM,  "deltatop_mGMM/F");
    tT->Branch("mGMM_m_mean",    &fTmGMM_m_mean,    "mGMM_m_mean/F");
    tT->Branch("mGMM_m_var",     &fTmGMM_m_var,     "mGMM_m_var/F");
    tT->Branch("mGMM_m_skew",    &fTmGMM_m_skew,    "mGMM_m_skew/F");
    tT->Branch("mGMM_pt_mean",   &fTmGMM_pt_mean,   "mGMM_pt_mean/F");
    tT->Branch("mGMM_pt_var",    &fTmGMM_pt_var,    "mGMM_pt_var/F");
    tT->Branch("mGMM_pt_skew",   &fTmGMM_pt_skew,   "mGMM_pt_skew/F");


    tT->Branch("mTGMM_m",        &fTmTGMM_m,        "mTGMM_m/F");
    tT->Branch("mTGMM_pt",       &fTmTGMM_pt,       "mTGMM_pt/F");
    tT->Branch("deltatop_mTGMM", &fTdeltatop_mTGMM, "deltatop_mTGMM/F");
    tT->Branch("mTGMM_m_mean",   &fTmTGMM_m_mean,   "mTGMM_m_mean/F");
    tT->Branch("mTGMM_m_var",    &fTmTGMM_m_var,    "mTGMM_m_var/F");
    tT->Branch("mTGMM_m_skew",   &fTmTGMM_m_skew,   "mTGMM_m_skew/F");
    tT->Branch("mTGMM_pt_mean",  &fTmTGMM_pt_mean,  "mTGMM_pt_mean/F");
    tT->Branch("mTGMM_pt_var",   &fTmTGMM_pt_var,   "mTGMM_pt_var/F");
    tT->Branch("mTGMM_pt_skew",  &fTmTGMM_pt_skew,  "mTGMM_pt_skew/F");


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

    fTmTGMMs_m        = -1.;
    fTmTGMMs_pt       = -1.;
    fTdeltatop_mTGMMs = -1.;
    fTmTGMMs_m_mean   = -1.;
    fTmTGMMs_m_var    = -1.;
    fTmTGMMs_m_skew   = -99999;
    fTmTGMMs_pt_mean  = -1.;
    fTmTGMMs_pt_var   = -1.;
    fTmTGMMs_pt_skew  = -99999;

    fTmGMM_m         = -1.;
    fTmGMM_pt        = -1.;
    fTdeltatop_mGMM  = -1.;
    fTmGMM_m_mean    = -1.;
    fTmGMM_m_var     = -1.;
    fTmGMM_m_skew    = -99999;
    fTmGMM_pt_mean    = -1.;
    fTmGMM_pt_var     = -1.;
    fTmGMM_pt_skew    = -99999;

    fTmUMM_m         = -1.;
    fTmUMM_pt        = -1.;
    fTdeltatop_mUMM  = -1.;
    fTmUMM_m_mean    = -1.;
    fTmUMM_m_var     = -1.;
    fTmUMM_m_skew    = -99999;
    fTmUMM_pt_mean    = -1.;
    fTmUMM_pt_var     = -1.;
    fTmUMM_pt_skew    = -99999;

    fTmTGMM_m        = -1.;
    fTmTGMM_pt       = -1.;
    fTdeltatop_mTGMM = -1.;
    fTmTGMM_m_mean    = -1.;
    fTmTGMM_m_var     = -1.;
    fTmTGMM_m_skew    = -99999;
    fTmTGMM_pt_mean    = -1.;
    fTmTGMM_pt_var     = -1.;
    fTmTGMM_pt_skew    = -99999;
}
