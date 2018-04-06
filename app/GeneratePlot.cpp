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

#include "Tools.h"

int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;

	std::ofstream OutFile;
	unsigned int nMAX = 10000, nS = 3, nR = 3;
	TFile *RootFile;
	bool TwoPerturbations = false;

	while((iarg = getopt_long(argc,argv, "I:R:S:yo:r:h", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'I':
				nMAX = strtol(optarg, NULL, 10);
				break;
			case 'y':
				TwoPerturbations = true;
				break;
			case 'R':
				nR = strtol(optarg, NULL, 10);
				break;
			case 'S':
				nS = strtol(optarg, NULL, 10);
				break;
			case 'r':
				RootFile = new TFile(optarg, "RECREATE");
				break;
			case 'o':
				OutFile.open(optarg);
				break;
			case 'h':
				return 1;
			default:
				break;
		}
	}
	std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;
	std::ofstream OutFile(argv[1]);
	std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	TFile *FileTree;
	TTree *Scan;
	TBranch *b_Dim, *b_MM, *b_VE, *b_VM, *b_VT, *b_BB0, *b_NSI, *b_MEG;   //!
	int Dim;
	double *MM = new double[Dim];
	double *VE = new double[Dim];
	double *VM = new double[Dim];
	double *VT = new double[Dim];
	//double OSC[6], EWS[6], MM[Dim], VE[Dim], VM[Dim], VT[Dim];
	bool NH, EXP, BB0, MEG, NSI;
	
	for (unsigned int i = 2; i < argc; ++i)
	{
		FileTree = new TFile(argv[i], "OPEN");
		FileTree->GetObject("Eigen", Scan);

		Scan->SetBranchAddress("Dim", &Dim, &b_Dim);
		Scan->SetBranchAddress("MM", MM, &b_MM);
		Scan->SetBranchAddress("VE", VE, &b_VE);
		Scan->SetBranchAddress("VM", VM, &b_VM);
		Scan->SetBranchAddress("VT", VT, &b_VT);
		Scan->SetBranchAddress("BB0", &BB0, &b_BB0);
		Scan->SetBranchAddress("NSI", &NSI, &b_NSI);
		Scan->SetBranchAddress("MEG", &MEG, &b_MEG);

		for (unsigned int j = 0; j < Scan->GetEntriesFast(); ++j)
		{
			Scan->GetEntry(j);
			if (MM[4])
			bool Pass = BB0 && NSI && MEG;
			if (Pass)
				Out << MM[4]*1e-9 << "\t" << pow(VE[4], 2) << "\t" << pow(VM[4], 2) << "\t" << pow(VT[4], 2) << std::endl;
				//Out << MM[3]*1e-9 << "\t" << pow(VE[3], 2) << "\t" << pow(VM[3], 2) << "\t" << pow(VT[3], 2) << std::endl;
		}
		FileTree->Close();
	}

	return 0;
}
