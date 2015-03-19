#include "ComputationManager.h"

#include <assert.h>
#include <iostream>
#include <unordered_set>

#include "boost/foreach.hpp"

void ComputationManager::compile() {
    assert(!_compiled);
    // compile providers, results, and final state computations
    BOOST_FOREACH(ComputationBase *c, _all_computations) {
        BOOST_FOREACH(label l, c->provided()) {
            if (_providers.find(l) == _providers.end()) {
                // can safely install the provided computation
                _providers[l] = c;
            } else {
                std::cout << "CONFLICT: " << c->name() << " and "
                          << _providers[l]->name() << " both provide "
                          << l << "." << std::endl;
                assert(0);
            }
        }
    }

    std::unordered_set<label> actually_required;
    std::unordered_set<ComputationBase *>
        inserted_computations(_final_state_computations.begin(),
                              _final_state_computations.end());

    BOOST_FOREACH(ComputationBase *c, _final_state_computations) {
        BOOST_FOREACH(label l, c->provided()) {
            actually_required.insert(l);
            _results[l] = c->how_to_provide(l);
        }
    }

    for (size_t iter = 0; iter < _final_state_computations.size(); iter++) {
        ComputationBase *c = _final_state_computations.at(iter);

        // record all the provided computations
        BOOST_FOREACH(label l, c->provided()) {
            if(_results.find(l) == _results.end()) {
                _results[l] = c->how_to_provide(l);
            }
        }

        BOOST_FOREACH(label l, c->required()) {
            if (_providers.find(l) == _providers.end()) {
                std::cout << "No computation provides " << l << "." << std::endl;
                assert(0);
            }

            actually_required.insert(l);
            if (inserted_computations.find(_providers[l]) == inserted_computations.end()) {
                inserted_computations.insert(_providers[l]);
                _final_state_computations.push_back(_providers[l]);
            }
        }
    }

    _toggle = true;
    _compiled = true;

    BOOST_FOREACH(ComputationBase *c, _final_state_computations) {
        if (c->toggled() != _toggle) c->toggle();
        BOOST_FOREACH(label l, c->provided()) {
            if (actually_required.find(l) == actually_required.end()) {
                c->may_disable(l);
            }
        }
    }
}

void ComputationManager::compute_step() {
    _toggle = !_toggle;

    assert(_compiled);
    BOOST_FOREACH(ComputationBase *c, _final_state_computations) {
        if (c->toggled() == _toggle) continue;

        c->toggle();
        c->compute(this);
    }
}

void ComputationManager::publish() {
    _toggle = !_toggle;

    _publishing = true;
    BOOST_FOREACH(ComputationBase *c, _final_state_computations) {
        if (c->toggled() == _toggle) continue;
        c->toggle();
        c->publish(this);
    }
}
