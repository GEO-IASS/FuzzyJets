# --------------------------------------------- #
# Makefile for FuzzyJets code                   #
# Ben N, May 3, 2014                            #
# Based on code by Pascal Nef, March 6th 2014   #
#                                               #
# Note: source setup.sh before make             #
# --------------------------------------------- #
OPTIMIZATION = -O3
CXXFLAGS = -Wall -Wextra -equal -Werror -Wno-shadow -fstack-protector-all -Wno-unused-but-set-variable
CXXFLAGSEXTRA = -g -std=c++0x
CXX = g++

# Are we profiling?
#PROFILER = -pg
PROFILER =

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

Histogrammer: Histogrammer.so AnalyzeFuzzyTools.so AtlasUtils.so Event.so Histogram.so
	$(CXX) Histogrammer.so AnalyzeFuzzyTools.so AtlasUtils.so Event.so Histogram.so -o $@ \
	$(CXXFLAGS) $(CXXFLAGSEXTRA) $(ROOTLIBSREAL) \
	-L$(BOOSTLIBLOCATION) -lboost_program_options

Histogrammer.so: Histogrammer.cc AnalyzeFuzzyTools.so AtlasUtils.so Event.so Histogram.so
	$(CXX) -o $@ -c $< \
	$(CXXFLAGS) $(CXXFLAGSEXTRA) -fPIC -shared $(ROOTFLAGSREAL) \
	-isystem $(BOOSTINCDIR)

AnalyzeFuzzyTools.so: AnalyzeFuzzyTools.cc AnalyzeFuzzyTools.h
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) $(CXXFLAGSEXTRA) -fPIC -shared $(ROOTFLAGSREAL) \

Event.so: Event.cc Event.h
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) $(CXXFLAGSEXTRA) -fPIC -shared $(ROOTFLAGSREAL) \

Histogram.so: Histogram.cc Histogram.h
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) $(CXXFLAGSEXTRA) -fPIC -shared $(ROOTFLAGSREAL) \

AtlasUtils.so: AtlasUtils.cc AtlasUtils.h
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) $(CXXFLAGSEXTRA) -fPIC -shared $(ROOTFLAGSREAL)

Fuzzy:  Fuzzy.so FuzzyTools.so FuzzyAnalysis.so
	$(CXX) Fuzzy.so FuzzyTools.so FuzzyAnalysis.so -o $@ $(PROFILER) \
	$(CXXFLAGS) $(OPTIMIZATION) $(ROOTLIBS) \
	-L$(FASTJETLOCATION)/lib $(FASTJETLIBS) \
	-L$(PYTHIA8LOCATION)/lib -lpythia8 -llhapdfdummy \
	-L$(BOOSTLIBLOCATION) -lboost_program_options -lNsubjettiness


Fuzzy.so: Fuzzy.C    FuzzyTools.so FuzzyAnalysis.so
	$(CXX) -o $@ -c $< $(PROFILER)  \
	$(CXXFLAGS) $(OPTIMIZATION) -fPIC -shared $(FASTJETFLAGS) \
	-isystem $(PYTHIA8LOCATION)/include \
	-isystem $(BOOSTINCDIR) $(ROOTFLAGS)

FuzzyTools.so : FuzzyTools.cc FuzzyTools.h
	$(CXX) -o $@ -c $< $(PROFILER) \
	$(CXXFLAGS) $(OPTIMIZATION) -fPIC -shared $(FASTJETFLAGS) \
	-isystem $(PYTHIA8LOCATION)/include $(ROOTFLAGS)

FuzzyAnalysis.so : FuzzyAnalysis.cc FuzzyAnalysis.h
	$(CXX) -o $@ -c $< $(PROFILER) \
	$(CXXFLAGS) $(OPTIMIZATION) -fPIC -shared $(FASTJETFLAGS) \
	-isystem $(PYTHIA8LOCATION)/include $(ROOTFLAGS) -INsubjettiness

clean:
	rm -rf *.exe
	rm -rf *.o
	rm -rf *.so
	rm -f *~

clear:
	rm ./results/tmp/*.root
