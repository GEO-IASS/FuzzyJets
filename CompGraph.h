#ifndef COMPGRAPH_H
#define COMPGRAPH_H

#include <string>
#include <assert.h>

#define NOTCOMPUTED 0
#define COMPUTED 1

// supports float values, due to C++ restrictions on
// type heterogeneity in STL containers. This could be circumvented
// with some Pointer Magic, but I don't want to bother at the moment.
template <typename T>
struct ComputerWrapper {
    typedef T (*_computer) (ComputationManager const& manager);
};

class ComputationBase {
protected:
    unsigned int _status;
    unsigned int _dependency;
    unsigned int _internal;
    std::string _label;
    
public:
    ComputationBase() {};
    
    virtual void doComputation(ComputationManager const& manager);
};

template <typename T>
class Computation : public ComputationBase {
protected:
    T _value;
    ComputerWrapper _computer_wrapper;
    
    
public:
    Computation(std::string label, ComputerWrapper wrapper)
        : _computer_wrapper(wrapper)
        , _label(label) {
        _status = NOTCOMPUTED;
    }

    T provide() {
        assert(_status = COMPUTED);
        return _value;
    }
    
    void doComputation(ComputationManager const& manager) {
        _value = *(_computer_wrapper.computer) (manager);
    }
};


#endif
