# --------------------------------------------- #
# Makefile for FuzzyJets code                   #
# Ben N, May 3, 2014                            #
# Based on code by Pascal Nef, March 6th 2014   #
#                                               #
# Note: source setup.sh before make             #
# --------------------------------------------- #

CXXFLAGS =   -O2 -Wall -Wextra -equal -Wshadow -Werror -Wno-shadow
CXX = g++

.PHONY: clean debug all clear

all: Fuzzy Histogrammer

Histogrammer: Histogrammer.so AnalyzeFuzzyTools.so AtlasUtils.so
	$(CXX) Histogrammer.so AnalyzeFuzzyTools.so AtlasUtils.so -o $@ \
	$(CXXFLAGS) `root-config --glibs` \
	-L$(BOOSTLIBLOCATION) -lboost_program_options

Histogrammer.so: Histogrammer.cc AnalyzeFuzzyTools.so AtlasUtils.so
	$(CXX) -o $@ -c $< \
	$(CXXFLAGS) -fPIC -shared `root-config --cflags` \
	-isystem $(BOOSTINCDIR)

AnalyzeFuzzyTools.so: AnalyzeFuzzyTools.cc AnalyzeFuzzyTools.h
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -fPIC -shared `root-config --cflags`

AtlasUtils.so: AtlasUtils.cc AtlasUtils.h
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -fPIC -shared `root-config --cflags`

Fuzzy:  Fuzzy.so FuzzyTools.so FuzzyAnalysis.so
	$(CXX) Fuzzy.so FuzzyTools.so FuzzyAnalysis.so -o $@ \
	$(CXXFLAGS) `root-config --glibs`  \
	-L$(FASTJETLOCATION)/lib \
	`$(FASTJETLOCATION)/bin/fastjet-config --libs --plugins ` \
	-L$(PYTHIA8LOCATION)/lib -lpythia8 -llhapdfdummy \
	-L$(BOOSTLIBLOCATION) -lboost_program_options

Fuzzy.so: Fuzzy.C    FuzzyTools.so FuzzyAnalysis.so
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -fPIC -shared \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins` \
	-isystem $(PYTHIA8LOCATION)/include \
	-isystem $(BOOSTINCDIR) \
	`root-config --cflags`

FuzzyTools.so : FuzzyTools.cc FuzzyTools.h
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -fPIC -shared \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins` \
	-isystem $(PYTHIA8LOCATION)/include \
	`root-config --cflags`

FuzzyAnalysis.so : FuzzyAnalysis.cc FuzzyAnalysis.h
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -fPIC -shared \
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
