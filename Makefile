.PHONY: clean

ROOTLIBS = $(shell root-config --glibs)
ROOTCXXF = $(shell root-config --cflags)
LDFLAGS  := -Wl,--no-as-needed $(LDFLAGS) $(ROOTLIBS)
CXXFLAGS := $(CXXFLAGS) -std=c++11 -O3 -mavx $(ROOTCXXF) -Iinclude/

TGT =	Delineate	\
	SVD		\
	Plotter
	#GeneratePlot	\
	RandomScan	\
	Test		\
	OneGeneration

all: $(TGT)

old:
	g++ mass.c -lgsl -lgomp -lgslcblas
