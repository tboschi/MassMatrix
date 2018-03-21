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

	unsigned int Cap = 1000000;
	unsigned int Dim = ISS->nM();
	unsigned int Var = 3 * nR + (nR + nS) * (nR + nS + 1);
	int Mag[Var];
	double MM[Dim], VE[Dim], VM[Dim], VT[Dim];
	bool NH, EXP, BB0, MEG, NSI;

	TTree *tEigen = new TTree("Eigen", "eigen");

	tEigen->Branch("Dim", &Dim, "iDim/I");
	tEigen->Branch("Var", &Var, "iVar/I");
	tEigen->Branch("NH",  &NH,  "bNH/O");
	tEigen->Branch("EXP", &EXP, "bEXP/O");
	tEigen->Branch("BB0", &BB0, "bBB0/O");
	tEigen->Branch("MEG", &MEG, "bMEG/O");
	tEigen->Branch("NSI", &NSI, "bNSI/O");
	tEigen->Branch("Mag", Mag,  "fMag[iVar]/I");
	tEigen->Branch("MM",  MM,   "fMM[iDim]/D");
	tEigen->Branch("VE",  VE,   "fVE[iDim]/D");
	tEigen->Branch("VM",  VM,   "fVM[iDim]/D");
	tEigen->Branch("VT",  VT,   "fVT[iDim]/D");

	std::cout << "Realisation ISS(" << ISS->nR() << "," << ISS->nS() << ")" << std::endl;
	std::cout << "Saving to file every " << Cap << " entries" << std::endl;
	std::cout << "Total number of savings expected is " << nMAX/Cap+1 << std::endl;
	
	ISS->Clean(Block::Full);
	ISS->Set(Block::Mr,  4,  6);
	ISS->Set(Block::Ms,  6, 15);
	ISS->Set(Block::Ur, -5,  4);
	ISS->Set(Block::Us, -5,  4);

	ISS->Clean(Block::Mr);
	ISS->Clean(Block::Ms);
	ISS->Clean(Block::Ur);
	ISS->Clean(Block::Us);

	std::vector<double> vVal;
	std::vector<int> vMag;

	for (unsigned int i = 0; i < nMAX; ++i)
	{
		ISS->Clean(Block::Full);
		vMag = ISS->Populate(Block::Full);

		//ISS->Show(Block::Full, 0);

		Eigen::MatrixXcd VV = ISS->MassMatrixSVD(vVal);
		if (ISS->FindDeltaM2(vVal, NH))
		{
			EXP = ISS->FindMass(vVal, 1e6, 2e9);
			BB0 = ISS->BB0(vVal, VV);
			MEG = ISS->MEG(vVal, VV);
			NSI = ISS->NSI(vVal, VV);

			for (unsigned int i = 0; i < vMag.size(); ++i)
				Mag[i] = vMag.at(i);

			for (unsigned int i = 0; i < Dim; ++i)
			{
				MM[i] = vVal.at(i);
				VE[i] = std::abs(VV(0, i));
				VM[i] = std::abs(VV(1, i));
				VT[i] = std::abs(VV(2, i));
			}

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

			std::cout << "Saving " << i/Cap << std::endl;
			tEigen->Write();

			//std::stringstream ssl;
			//ssl << "matrix" << int((i+1)/Cap); 
			//hMatrix->Write(ssl.str().c_str());
		}
	}

	tEigen->Write();

	/*
	std::cout <<   "Mr\t";
	for (unsigned int i = 0; i < ISS->Get(Block::Mr).size(); ++i)
		std::cout << std::abs(*(ISS->Get(Block::Mr).data() + i)) << "\t";
	std::cout << "\nlogMr\t";
	for (auto p : vMr)
		std::cout << p << "\t";
	std::cout << std::endl;

	std::cout <<   "Ms\t";
	for (unsigned int i = 0; i < ISS->Get(Block::Ms).size(); ++i)
		std::cout << std::abs(*(ISS->Get(Block::Ms).data() + i)) << "\t";
	std::cout << "\nlogMs\t";
	for (auto p : vMs)
		std::cout << p << "\t";
	std::cout << std::endl << std::endl;

	Eigen::MatrixXcd VV = ISS->MassMatrixSVD(vValues);
	std::cout << "0th order" << "\t";
	for (unsigned int i = 0; i < ISS->n0(); ++i)
		std::cout << vValues.at(2*i) << "\t";
	std::cout << "and " << 3+ISS->nR()+ISS->nS() - 2*ISS->n0() << " massless";
	std::cout << std::endl << std::endl << VV.cwiseAbs() << std::endl << std::endl;
	std::cout << "____________________________________________________________" << std::endl;
	std::complex<double> dUrij = ISS->RandomAs(Block::Ur);
	std::complex<double> dUsij = ISS->RandomAs(Block::Us);
	std::complex<double> dUr00 = dUrij;
	std::complex<double> dUs00 = dUsij;
	std::cout << "Perturbation " << std::abs(dUrij) << "\t" << std::abs(dUsij) << std::endl;
	std::cout << "____________________________________________________________" << std::endl;

	for (unsigned int i = 0; i < nMAX; ++i)
	{
		for (unsigned int j = i; j < ISS->nS(); ++j)
		{
			ISS->Clean(Block::Ur);
			ISS->Clean(Block::Us);
			ISS->Manual(Block::Us, dUsij, i, j);
			ISS->Manual(Block::Us, dUsij, j, i);
			if (TwoPerturbations)
			{
				ISS->Manual(Block::Us, dUs00, 0, 1);
				ISS->Manual(Block::Us, dUs00, 1, 0);
			}
			std::cout << "Perturbed at (" << i << "," << j << ")" << std::endl;
			std::cout << ISS->Get(Block::Us).cwiseAbs() << std::endl << std::endl;

			VV = ISS->MassMatrixSVD(vValues);
			for (auto p : vValues)
				std::cout << p << "\t";
			std::cout << std::endl;
			std::cout << std::endl << std::endl << VV.cwiseAbs() << std::endl << std::endl;
			//for (unsigned int i = 0; i < ISS->n0(); ++i)
			//	std::cout << (vValues.at(2*i)+vValues.at(2*i+1))*0.5 << "+/-" << vValues.at(2*i)-vValues.at(2*i+1) << "\t";
			//for (unsigned int i = 2*ISS->n0(); i < vValues.size(); ++i)
			//	std::cout << vValues.at(i) << "\t";
			//std::cout << std::endl << std::endl << VV.cwiseAbs() << std::endl << std::endl;
			//std::cout << "____________________________________________________________" << std::endl;
		}
	}
	*/

	return 0;
}
