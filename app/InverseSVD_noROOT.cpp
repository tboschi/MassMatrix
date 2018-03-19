#include <iostream>
#include <algorithm>
#include <string>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <getopt.h>

#include "Tools.h"
#include "InverseMatrix.h"

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
	bool NH, BB0, MEG, NSI;

	ISS->Clean(Block::Full);
	ISS->Set(Block::Mr,  4, 15);
	ISS->Set(Block::Ms,  4, 15);
	ISS->Set(Block::Ur, -4,  4);
	ISS->Set(Block::Us, -4,  4);

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

		Eigen::MatrixXcd VV = ISS->MassMatrixSVD(vVal);
		if (ISS->FindDeltaM2(vVal, NH))
		{
			if (ISS->FindMass(vVal, 1e6, 2e9))
			{
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
	
				//Place here necessary std::cout or Out (streamer defined with -o option)
				//std::cout << MM[3] << "\t" << VE[3] << std::endl;
				//Out << MM[3] << "\t" << VE[3] << std::endl;
			}
		}

	}

	return 0;
}
