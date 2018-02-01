.PHONY: clean

#ROOTLIBS = $(shell root-config --glibs)
#ROOTCXXF = $(shell root-config --cflags)
LDFLAGS  := $(LDFLAGS) $(ROOTLIBS)
CXXFLAGS := $(CXXFLAGS) -std=c++11 $(ROOTCXXF) -Iinclude/

TGT =	RandomScan
#TGT =	Test
#TGT =	OneGeneration

all: $(TGT)

old:
	g++ mass.c -lgsl -lgomp -lgslcblas
