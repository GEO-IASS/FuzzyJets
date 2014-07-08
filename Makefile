# --------------------------------------------- #
# Makefile for FuzzyJets code                   #
# Ben N, May 3, 2014                            #
# Based on code by Pascal Nef, March 6th 2014   #
#                                               #
# Note: source setup.sh before make             #
# --------------------------------------------- #

CXXFLAGS =   -O2 -Wall -Wextra -equal -Wshadow -Werror
CXX = g++

.PHONY: clean debug all clear

all: Fuzzy

Fuzzy:  Fuzzy.so FuzzyTools.so FuzzyAnalysis.so
	$(CXX) Fuzzy.so FuzzyTools.so FuzzyAnalysis.so -o $@ \
	$(CXXFLAGS) -Wno-shadow  \
	`root-config --glibs`  \
	-L$(FASTJETLOCATION)/lib `$(FASTJETLOCATION)/bin/fastjet-config --libs --plugins ` \
	-L$(PYTHIA8LOCATION)/lib -lpythia8 -llhapdfdummy \
	-L$(BOOSTLIBLOCATION) -lboost_program_options

Fuzzy.so: Fuzzy.C    FuzzyTools.so FuzzyAnalysis.so
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -Wno-shadow -fPIC -shared \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins` \
	-isystem $(PYTHIA8LOCATION)/include \
	-isystem $(BOOSTINCDIR) \
	`root-config --cflags`

FuzzyTools.so : FuzzyTools.cc FuzzyTools.h
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -Wno-shadow -fPIC -shared \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins` \
	-isystem $(PYTHIA8LOCATION)/include \
	`root-config --cflags`

FuzzyAnalysis.so : FuzzyAnalysis.cc FuzzyAnalysis.h
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -Wno-shadow -fPIC -shared \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins`  \
	-isystem $(PYTHIA8LOCATION)/include \
	`root-config --cflags`

clean:
	rm -rf *.exe
	rm -rf *.o
	rm -rf *.so
	rm -f *~

clear:
	rm ./results/tmp/*.root
