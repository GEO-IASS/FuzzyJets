#include "myFastJetBase.h"

#include "FuzzyUtil.h"

void generate_particles(vecPseudoJet &out,
                        vecPseudoJet &tops_out,
                        float &max_rap,
                        bool is_pileup,
                        Pythia8::Pythia *pythia_inst) {
    assert(pythia_inst->next());
    while(tops_out.size() < 2) {
        tops_out.push_back(fastjet::PseudoJet());
    }

    // Particle loop -----------------------------------------------------------
    // The Pythia event listing contains a lot more than we want to process,
    // we prune out certain particles (muons / neutrinos) and only add final
    // state particles
    double px, py, pz, e;
        for (unsigned int particle_idx = 0; particle_idx < (unsigned) pythia_inst->event.size(); ++particle_idx){
            px = pythia_inst->event[particle_idx].px();
            py = pythia_inst->event[particle_idx].py();
            pz = pythia_inst->event[particle_idx].pz();
            e  = pythia_inst->event[particle_idx].e();
            fastjet::PseudoJet p(px, py, pz, e);
            p.reset_PtYPhiM(p.pt(), p.rapidity(), p.phi(), 0.);
            p.set_user_info(new MyUserInfo(pythia_inst->event[particle_idx].id(),particle_idx,pythia_inst->event[particle_idx].charge(),is_pileup,false));

            // In reality we should be more careful about finding tops,
            // but this will do for now. In the future consider refactoring
            // and tracing to find a top quark with no daughters
            if (pythia_inst->event[particle_idx].id()  ==6) tops_out.at(0)=p;
            if (pythia_inst->event[particle_idx].id()  ==-6) tops_out.at(1)=p;

            // prune uninteresting particles
            if (!pythia_inst->event[particle_idx].isFinal() )      continue; // only final state
            if (fabs(pythia_inst->event[particle_idx].id())  ==11) continue; // ...   electron
            if (fabs(pythia_inst->event[particle_idx].id())  ==12) continue; // prune nu-e
            if (fabs(pythia_inst->event[particle_idx].id())  ==13) continue; // ...   mu
            if (fabs(pythia_inst->event[particle_idx].id())  ==14) continue; // ...   nu-mu
            if (fabs(pythia_inst->event[particle_idx].id())  ==16) continue; // ...   nu-tau
            if (pythia_inst->event[particle_idx].pT()       < 0.5) continue; // ...   low pT

            max_rap = max_rap < fabs(p.rap()) ? fabs(p.rap()) : max_rap;
            out.push_back(p);
        } // end particle loop -----------------------------------------------
}

void generate_particles_multi(vecPseudoJet &out,
                              vecPseudoJet &tops_out,
                              float &max_rap,
                              bool is_pileup,
                              Pythia8::Pythia *pythia_inst,
                              size_t n_events) {
    while(tops_out.size() < 2) {
        tops_out.push_back(fastjet::PseudoJet());
    }
    for (size_t i = 0; i < n_events; i++) {
        generate_particles(out, tops_out, max_rap, is_pileup, pythia_inst);
    }
}
