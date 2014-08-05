# --------------------------------------------- #
# Makefile for FuzzyJets code                   #
# Ben N, May 3, 2014                            #
# Based on code by Pascal Nef, March 6th 2014   #
#                                               #
# Note: source setup.sh before make             #
# --------------------------------------------- #

CXXFLAGS =   -O2 -Wall -Wextra -equal -Wshadow -Werror -Wno-shadow
CXX = g++

# Are we profiling?
PROFILER = -pg
#PROFILER =

# Use these for ROOTless compiles
#ROOTLIBS =
#ROOTFLAGS = 

# Use these for actually logging data
ROOTLIBS = `root-config --glibs`
ROOTFLAGS = `root-config --cflags`

# Use these always, of course
ROOTLIBSREAL = `root-config --glibs`
ROOTFLAGSREAL = `root-config --cflags`
FASTJETLIBS = `$(FASTJETLOCATION)/bin/fastjet-config --libs --plugins`
FASTJETFLAGS = `$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins`

.PHONY: clean debug all clear

all: Fuzzy Histogrammer

Histogrammer: Histogrammer.so AnalyzeFuzzyTools.so AtlasUtils.so
	$(CXX) Histogrammer.so AnalyzeFuzzyTools.so AtlasUtils.so -o $@ \
	$(CXXFLAGS) $(ROOTLIBSREAL) \
	-L$(BOOSTLIBLOCATION) -lboost_program_options

Histogrammer.so: Histogrammer.cc AnalyzeFuzzyTools.so AtlasUtils.so
	$(CXX) -o $@ -c $< \
	$(CXXFLAGS) -fPIC -shared $(ROOTFLAGSREAL) \
	-isystem $(BOOSTINCDIR)

AnalyzeFuzzyTools.so: AnalyzeFuzzyTools.cc AnalyzeFuzzyTools.h
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -fPIC -shared $(ROOTFLAGSREAL) \

AtlasUtils.so: AtlasUtils.cc AtlasUtils.h
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -fPIC -shared $(ROOTFLAGSREAL)

Fuzzy:  Fuzzy.so FuzzyTools.so FuzzyAnalysis.so
	$(CXX) Fuzzy.so FuzzyTools.so FuzzyAnalysis.so -o $@ $(PROFILER) \
	$(CXXFLAGS) $(ROOTLIBS) \
	-L$(FASTJETLOCATION)/lib $(FASTJETLIBS) \
	-L$(PYTHIA8LOCATION)/lib -lpythia8 -llhapdfdummy \
	-L$(BOOSTLIBLOCATION) -lboost_program_options

Fuzzy.so: Fuzzy.C    FuzzyTools.so FuzzyAnalysis.so
	$(CXX) -o $@ -c $< $(PROFILER)  \
	$(CXXFLAGS) -fPIC -shared $(FASTJETFLAGS) \
	-isystem $(PYTHIA8LOCATION)/include \
	-isystem $(BOOSTINCDIR) $(ROOTFLAGS)

FuzzyTools.so : FuzzyTools.cc FuzzyTools.h
	$(CXX) -o $@ -c $< $(PROFILER) \
	$(CXXFLAGS) -fPIC -shared $(FASTJETFLAGS) \
	-isystem $(PYTHIA8LOCATION)/include $(ROOTFLAGS)

FuzzyAnalysis.so : FuzzyAnalysis.cc FuzzyAnalysis.h
	$(CXX) -o $@ -c $< $(PROFILER) \
	$(CXXFLAGS) -fPIC -shared $(FASTJETFLAGS) \
	-isystem $(PYTHIA8LOCATION)/include $(ROOTFLAGS)

clean:
	rm -rf *.exe
	rm -rf *.o
	rm -rf *.so
	rm -f *~

clear:
	rm ./results/tmp/*.root
