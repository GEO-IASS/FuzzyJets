#ifndef COMPUTATIONMANAGER_H
#define COMPUTATIONMANAGER_H

#include <list>
#include <iostream>
#include <cstddef>
#include <string>
#include <vector>
#include <map>
#include <assert.h>

typedef std::string label;

class ComputationManager;

class ComputationBase {
protected:
    // essentially a binary clock
    bool _toggle;

    std::string _name;
    std::vector<label> _provides;
    std::vector<label> _requires;

public:
    ComputationBase() {
        _name = "";
        _provides = {};
        _requires = {};
        _toggle = false;
    }

    void toggle() {
        _toggle = !_toggle;
    }

    bool toggled() {
        return _toggle;
    }

    virtual void *how_to_provide(label l) = 0;

    std::vector<label> const& provided() {
        return _provides;
    }

    std::vector<label> const& required() {
        return _requires;
    }

    virtual void compute(__attribute__((unused)) ComputationManager *manager) {};

    virtual void publish(__attribute__((unused)) ComputationManager *manager) {};

    // allows disabling branches in trees, for instance
    virtual void may_disable(__attribute__((unused)) label l) {};

    std::string name() {
        return _name;
    }
};

class ComputationManager {
protected:
    // master clock to check whether computations have proceeded
    bool _toggle;
    bool _pushing_to_final_state;
    bool _compiled;
    bool _publishing;

    std::map<label, ComputationBase *> _providers;
    std::map<label, void*> _results;

    std::vector<ComputationBase *> _final_state_computations;
    std::vector<ComputationBase *> _all_computations;

public:
    ComputationManager() {
        _compiled = false;
        _publishing = false;
        _pushing_to_final_state = false;
        _toggle = false;
    };

    void compile();

    void publish();

    void compute_step();

    void compute(size_t steps) {
        for (size_t iter = 0; iter < steps; iter++) {
            compute_step();
        }
    }

    void push_to_final_state() {
        _pushing_to_final_state = true;
    }

    void install_computation(ComputationBase *c) {
        _all_computations.push_back(c);
        if (_pushing_to_final_state)
            _final_state_computations.push_back(c);
    }

    // quite dangerous, yay type systems :P
    template <typename T> const T Get(label l) {
        if (_results.find(l) == _results.end()) {
            assert(0);
        }
        if (_providers[l]->toggled() != _toggle) {
            _providers[l]->toggle();
            if (_publishing) {
                _providers[l]->publish(this);
            } else {
                _providers[l]->compute(this);
            }
        }
        return static_cast<T>(_results.find(l)->second);
    }
};

#endif
