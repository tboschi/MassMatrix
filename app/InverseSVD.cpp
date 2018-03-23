#include <iostream>
#include <algorithm>
#include <string>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <getopt.h>

#include "TFile.h"
#include "TTree.h"
#include "TH2I.h"
#include "Tools.h"
#include "InverseMatrix.h"

void FillHistogram(TH2I *&hM, std::vector<int> vMag, unsigned int &k, unsigned int ix, unsigned int iy, unsigned int Dim)
{
	hM->SetBinContent(ix, iy, vMag.at(k));
	hM->SetBinContent(Dim-iy+1, Dim-ix+1, vMag.at(k));
	++k;
}

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

	InverseMatrix *ISS = new InverseMatrix(nR, nS);

	unsigned int Cap = 1e7;
	unsigned int Dim = ISS->nM();
	unsigned int Var = 3 * nR + (nR + nS) * (nR + nS + 1);
	int Mag[Var];
	double Mbb, bMG;
	double OSC[6], EWS[6], MM[Dim], VE[Dim], VM[Dim], VT[Dim];
	bool NH, EXP, BB0, MEG, NSI;

	TTree *tEigen = new TTree("Eigen", "eigen");

	tEigen->Branch("Dim", &Dim, "iDim/I");
	//tEigen->Branch("Var", &Var, "iVar/I");
	tEigen->Branch("NH",  &NH,  "bNH/O");
	//tEigen->Branch("EXP", &EXP, "bEXP/O");
	tEigen->Branch("BB0", &BB0, "bBB0/O");	//bool neutrinoless doublebeta
	tEigen->Branch("Mbb", &Mbb, "fMbb/D");	//effective doubel beta mass
	tEigen->Branch("MEG", &MEG, "bMEG/O");	//bool mu to e gamma
	tEigen->Branch("bMG", &bMG, "fbMG/D");	//mu to e gamma branch
	tEigen->Branch("NSI", &NSI, "bNSI/O");	//bool if nonunitarity is satisfied
	tEigen->Branch("OSC", &OSC, "fOSC[6]/D");	//deviation from osc
	tEigen->Branch("EWS", &EWS, "fEWS[6]/D");	//deviation from EW precision
	//tEigen->Branch("Mag", Mag,  "fMag[iVar]/I");
	tEigen->Branch("MM",  MM,   "fMM[iDim]/D");
	tEigen->Branch("VE",  VE,   "fVE[iDim]/D");
	tEigen->Branch("VM",  VM,   "fVM[iDim]/D");
	tEigen->Branch("VT",  VT,   "fVT[iDim]/D");

	Out << "Realisation ISS(" << ISS->nR() << "," << ISS->nS() << ")" << std::endl;
	Out << "Saving to file every " << Cap << " entries" << std::endl;
	Out << "Total number of savings expected is " << nMAX/Cap << std::endl;
	
	ISS->Clean(Block::Full);

	ISS->Set(Block::Mr, 4, 6);
	ISS->Set(Block::Ms, 6, 10);
	ISS->Set(Block::Ur, -4, 3);
	ISS->Set(Block::Us, -4, 3);

	ISS->Clean(Block::Mr);
	ISS->Clean(Block::Ms);
	ISS->Clean(Block::Ur);
	ISS->Clean(Block::Us);

	std::vector<double> vVal, vOsc, vEws;
	//std::vector<int> vMag;

	for (unsigned int i = 0; i < nMAX; ++i)
	{
		ISS->Clean(Block::Full);
		//vMag = ISS->Populate(Block::Full);
		ISS->Populate(Block::Full);

		//ISS->Clean(Block::Ur);
		//ISS->Clean(Block::Us);
		//ISS->Show(Block::Full, 0);

		vVal.clear();
		Eigen::MatrixXcd VV = ISS->MassMatrixSVD(vVal);

		if (ISS->FindDeltaM2(vVal, NH))
		{
			vOsc.clear();
			vEws.clear();
			EXP = ISS->FindMass(vVal, 1e6, 2e9);
			BB0 = ISS->BB0(vVal, VV, Mbb);
			MEG = ISS->MEG(vVal, VV, bMG);
			NSI = ISS->NSI(vVal, VV, vOsc, vEws);

			/*
			for (unsigned int i = 0; i < vMag.size(); ++i)
				Mag[i] = vMag.at(i);
			*/

			for (unsigned int i = 0; i < 6; ++i)
			{
				OSC[i] = vOsc.at(i);
				EWS[i] = vEws.at(i);
			}

			for (unsigned int i = 0; i < Dim; ++i)
			{
				MM[i] = vVal.at(i);
				VE[i] = std::abs(VV(0, i));
				VM[i] = std::abs(VV(1, i));
				VT[i] = std::abs(VV(2, i));
			}

			Out << "Filling " << i << " at " <<  tEigen->GetEntries() << std::endl;
			tEigen->Fill();
		}

		if (i % Cap == Cap-1)
		{
			/*
			unsigned int k = 0;
			for (unsigned int ix = 4; ix < 4+nR; ++ix)
				for (unsigned int iy = Dim; iy > Dim-3; --iy)
					FillHistogram(hMatrix, vMag, k, ix, iy, Dim);
			for (unsigned int ix = 4+nR; ix < Dim+1; ++ix)
				for (unsigned int iy = Dim-3; iy > Dim-3-nR; --iy)
					FillHistogram(hMatrix, vMag, k, ix, iy, Dim);
			for (unsigned int ix = 4; ix < 4+nR; ++ix)
				for (unsigned int iy = Dim-3; iy > Dim-4-(ix-4); --iy)
					FillHistogram(hMatrix, vMag, k, ix, iy, Dim);
			for (unsigned int ix = 4+nR; ix < Dim+1; ++ix)
				for (unsigned int iy = Dim-3-nR; iy > Dim-4-nR-(ix-4-nR); --iy)
					FillHistogram(hMatrix, vMag, k, ix, iy, Dim);
			*/

			Out << "Saving " << i/Cap << std::endl;
			tEigen->Write();

			//std::stringstream ssl;
			//ssl << "matrix" << int((i+1)/Cap); 
			//hMatrix->Write(ssl.str().c_str());
		}
	}

	tEigen->Write();

	return 0;
}
