#!/bin/bash

setup_PYTHIA() {
    export PYTHIA8LOCATION=/u/at/pnef/Work/Code/pythia8183/
    export PYTHIA8DATA=${PYTHIA8LOCATION}xmldoc/
    export LD_LIBRARY_PATH=${PYTHIA8LOCATION}lib/:$LD_LIBRARY_PATH
}

setup_ROOT() {
    source /u/at/pnef/Work/Code/root_v5.34.17/bin/thisroot.sh
}

setup_fastjet() {
    export FASTJETLOCATION=/u/at/chstan/nfs/src/fastjet-install/
    export LD_LIBRARY_PATH=${FASTJETPATH}lib/:$LD_LIBRARY_PATH
}

setup_boost() {
    export BOOSTINCDIR=/nfs/slac/g/atlas/u01/users/chstan/src/boost_1_56_0
    export BOOSTLIBLOCATION=/nfs/slac/g/atlas/u01/users/chstan/src/boost_1_56_0/stage/lib
    export LD_LIBRARY_PATH=${BOOSTLIBLOCATION}:$LD_LIBRARY_PATH

}

setup_ROOT
setup_PYTHIA
setup_fastjet
setup_boost
