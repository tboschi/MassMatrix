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
	TBranch *b_MM, *b_VE, *b_VM, *b_VT;   //!
	double MM[8], VE[8], VM[8], VT[8];   //[iDim]
	
	for (unsigned int i = 2; i < argc; ++i)
	{
		FileTree = new TFile(argv[i], "OPEN");
		FileTree->GetObject("Eigen", Scan);

		Scan->SetBranchAddress("MM", MM, &b_MM);
		Scan->SetBranchAddress("VE", VE, &b_VE);
		Scan->SetBranchAddress("VM", VM, &b_VM);
		Scan->SetBranchAddress("VT", VT, &b_VT);
		for (unsigned int j = 0; j < Scan->GetEntriesFast(); ++j)
		{
			Scan->GetEntry(j);
			Out << MM[3] << "\t" << VE[3] << "\t" << VM[3] << "\t" << VT[3] << std::endl;
		}
		FileTree->Close();
	}

	return 0;
}
