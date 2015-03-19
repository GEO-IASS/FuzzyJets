#ifndef FUZZYUTIL_H
#define FUZZYUTIL_H

#include "Pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"

typedef std::vector<fastjet::PseudoJet> vecPseudoJet;

void generate_particles(vecPseudoJet &out,
                        vecPseudoJet &tops_out,
                        float &max_rap,
                        bool is_pileup,
                        Pythia8::Pythia *pythia_inst);

void generate_particles_multi(vecPseudoJet &out,
                              vecPseudoJet &tops_out,
                              float &max_rap,
                              bool is_pileup,
                              Pythia8::Pythia *pythia_inst,
                              size_t n_events);

#endif
