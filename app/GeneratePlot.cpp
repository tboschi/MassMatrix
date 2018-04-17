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

bool TestLimit(const double M0, const double V0, std::vector<double> &vM, std::vector<double> &vU, unsigned int &cc);
int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"elimit", 	required_argument,	0, 'E'},
		{"mlimit", 	required_argument,	0, 'M'},
		{"tlimit", 	required_argument,	0, 'T'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;

	std::ofstream OutFile;
	std::ifstream FileE, FileM, FileT;
	int xI;
	TFile *TreeFile;

	while((iarg = getopt_long(argc,argv, "i:o:I:E:M:T:h", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'i':
				TreeFile = new TFile(optarg, "OPEN");
				break;
			case 'o':
				OutFile.open(optarg);
				break;
			case 'I':
				xI = std::strtol(optarg, NULL, 10);
				break;
			case 'E':
				FileE.open(optarg);
				break;
			case 'M':
				FileM.open(optarg);
				break;
			case 'T':
				FileT.open(optarg);
				break;
			case 'h':
				return 1;
			default:
				break;
		}
	}
	std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	TTree *Scan;
	TBranch *b_Dim, *b_MM, *b_VE, *b_VM, *b_VT, *b_BB0, *b_NSI, *b_MEG;   //!
	int Dim;
	double *MM = new double[Dim];
	double *VE = new double[Dim];
	double *VM = new double[Dim];
	double *VT = new double[Dim];
	//double OSC[6], EWS[6], MM[Dim], VE[Dim], VM[Dim], VT[Dim];
	bool NH, EXP, BB0, MEG, NSI;
	
	std::vector<double> vME, vMM, vMT;
	std::vector<double> vVE, vVM, vVT;

	std::string Line;
	std::stringstream ssL;
	double mm, uu;

	if (FileE.is_open())
		while (getline(FileE, Line))
		{
			if (Line[0] == '#')
				continue;
			ssL.str("");
			ssL.clear();

			ssL << Line;
			ssL >> mm >> uu;

			vME.push_back(mm);
			vVE.push_back(uu);
		}
	if (FileM.is_open())
		while (getline(FileM, Line))
		{
			if (Line[0] == '#')
				continue;
			ssL.str("");
			ssL.clear();

			ssL << Line;
			ssL >> mm >> uu;

			vMM.push_back(mm);
			vVM.push_back(uu);
		}
	if (FileT.is_open())
		while (getline(FileT, Line))
		{
			if (Line[0] == '#')
				continue;
			ssL.str("");
			ssL.clear();

			ssL << Line;
			ssL >> mm >> uu;

			vMT.push_back(mm);
			vVT.push_back(uu);
		}


	TreeFile->GetObject("Eigen", Scan);

	Scan->SetBranchAddress("Dim", &Dim, &b_Dim);
	Scan->SetBranchAddress("MM", MM, &b_MM);
	Scan->SetBranchAddress("VE", VE, &b_VE);
	Scan->SetBranchAddress("VM", VM, &b_VM);
	Scan->SetBranchAddress("VT", VT, &b_VT);
	Scan->SetBranchAddress("BB0", &BB0, &b_BB0);
	Scan->SetBranchAddress("NSI", &NSI, &b_NSI);
	Scan->SetBranchAddress("MEG", &MEG, &b_MEG);

	unsigned int ec = 0, mc = 0, tc = 0;
	for (unsigned int j = 0; j < Scan->GetEntriesFast(); ++j)
	{
		Scan->GetEntry(j);

		bool Pass = BB0 && NSI && MEG;
		if (Pass)
		{
			double m0 = MM[xI] * 1e-9;
			double ve = VE[xI] * VE[xI];
			double vm = VM[xI] * VM[xI];
			double vt = VT[xI] * VT[xI];

			bool ELimit = TestLimit(m0, ve, vME, vVE, ec); 
			bool MLimit = TestLimit(m0, vm, vMM, vVM, mc); 
			bool TLimit = TestLimit(m0, vt, vMT, vVT, tc); 

			//Out << m0 << "\t" << (ELimit ? ve : -1) << "\t" << (MLimit ? vm : -1) << "\t" << (TLimit ? vt : -1) << std::endl;
			Out << m0 << "\t" << ve << "\t" << vm << "\t" << vt << std::endl;
		}
	}

	TreeFile->Close();

	return 0;
}

bool TestLimit(const double M0, const double V0, std::vector<double> &vM, std::vector<double> &vU, unsigned int &cc)
{
	if (M0 < vM.front() || M0 > vM.back())
		return true;

	double m1, m2, v1, v2;
	if (M0 > vM.at(cc))
	{
		while (M0 > vM.at(cc) && cc < vM.size()-1)
			++(cc);
		m1 = vM.at(cc-1);
		m2 = vM.at(cc);
		v1 = vU.at(cc-1);
		v2 = vU.at(cc);
	}
	else
	{
		while (M0 < vM.at(cc) && cc >= 0)
			--cc;
		m1 = vM.at(cc);
		m2 = vM.at(cc+1);
		v1 = vU.at(cc);
		v2 = vU.at(cc+1);
	}

	//double Vx = (v2-v1)/(m2-m1) * (M0 - m1) + v1;
	double lVx = log10(v2/v1)/(m2-m1) * (M0 - m1) + log10(v1);
	if (log10(V0) < lVx)
		return true;
	else
		return false;
}
