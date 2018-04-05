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
		std::cout << "H0" << std::endl;
		FileTree = new TFile(argv[i], "OPEN");
		std::cout << "H1" << std::endl;
		FileTree->GetObject("Eigen", Scan);

		std::cout << "H2" << std::endl;
		Scan->SetBranchAddress("Dim", &Dim, &b_Dim);
		Scan->SetBranchAddress("MM", MM, &b_MM);
		Scan->SetBranchAddress("VE", VE, &b_VE);
		Scan->SetBranchAddress("VM", VM, &b_VM);
		Scan->SetBranchAddress("VT", VT, &b_VT);
		Scan->SetBranchAddress("BB0", &BB0, &b_BB0);
		Scan->SetBranchAddress("NSI", &NSI, &b_NSI);
		Scan->SetBranchAddress("MEG", &MEG, &b_MEG);

		std::cout << "H3 " << Scan->GetEntriesFast() << std::endl;
		for (unsigned int j = 480; j < Scan->GetEntriesFast(); ++j)
		{
			std::cout << "H3 " << Scan->GetEntriesFast() << std::endl;
			std::cout << "F0 " << j << std::endl;
			std::cout << Scan->GetEntry(j) << std::endl;
			std::cout << "F1" << std::endl;
			std::cout << "BB0 " << BB0 << "\tNSI " << NSI << "\tMEG " << MEG << std::endl;
			std::cout << "MM " << MM << "\tMM[0] " << MM[0] << std::endl;
			if (BB0 && NSI && MEG)
				Out << MM[3] << "\t" << VE[3] << "\t" << VM[3] << "\t" << VT[3] << std::endl;
			std::cout << "F2" << std::endl;
		}
		std::cout << "H4" << std::endl;
		FileTree->Close();
	}

	return 0;
}
