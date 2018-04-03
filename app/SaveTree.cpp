#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <random>
#include <chrono>
#include <vector>
#include <algorithm>
#include <getopt.h>

#include "TTree.h"
#include "TFile.h"
#include "TChain.h"

#include "Tools.h"

int main(int argc, char** argv)
{
	TChain Ch("Eigen");

	for (unsigned int i = 2; i < argc; ++i)
		Ch.Add(argv[i]);

	Ch.Merge(argv[1]);

	return 0;
}
