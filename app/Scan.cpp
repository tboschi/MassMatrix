#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <complex>
#include <random>
#include <chrono>
#include <vector>
#include <algorithm>
#include <getopt.h>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/SVD>

#include "TTree.h"
#include "TFile.h"

#include "Tools.h"

template <typename MM>
void Print(const Eigen::MatrixBase<MM> &S)
{
	std::cout << std::fixed << std::setprecision(5);
	for (unsigned int i = 0; i < S.rows(); ++i)
	{
		for (unsigned int j = 0; j < S.cols(); ++j)
			std::cout << "(" << std::abs(S(i, j)) << "," << Const::fDeg * std::arg(S(i,j)) << ")\t";
		std::cout << std::endl;
	}
}

template <typename MM>
void MajoranaPhase(Eigen::MatrixBase<MM> &MAJ, Eigen::MatrixBase<MM> &Diag)
{
	double phase;
	std::complex<double> maj;
	for (unsigned int i = 0; i < Diag.rows(); ++i)
	{
		MAJ(i, i).real( cos( -0.5 * std::arg( Diag(i,i) ) ) );
		MAJ(i, i).imag( sin( -0.5 * std::arg( Diag(i,i) ) ) );
	}
}

void Pol2Cart(std::complex<double> &w, double lMod, double Phs)
{
	std::cout << "Module " << lMod << "\tPhase " << Phs*Const::fDeg << std::endl;
	w.real(pow(10, lMod)*cos(Phs));
	w.imag(pow(10, lMod)*sin(Phs));
}

int main(int argc, char** argv)
{
	const double Pi = 3.14159;
	const double deg = Pi/180;
	//TRandom3 *Ran = new TRandom3();

	const struct option longopts[] = 
	{
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;

	std::ofstream OutFile;
	TFile *RootFile;
	int dd = 9, nn = 12, mm = 0, uu = 6;
	unsigned int nMAX = 10000;
	unsigned long Seed = std::chrono::system_clock::now().time_since_epoch()/std::chrono::milliseconds(1);
	double CL = 0.90;	//90% C.L.

	while((iarg = getopt_long(argc,argv, "I:S:o:r:h", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'I':
				nMAX = strtol(optarg, NULL, 10);
				break;
			case 'S':
				Seed = strtol(optarg, NULL, 10);
				break;
			case 'o':
				OutFile.open(optarg);
				break;
			case 'r':
				RootFile = new TFile(optarg, "RECREATE");
				break;
			case 'h':
				return 1;
			default:
				break;
		}
	}
	std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	const unsigned int nD = 8;
	int Dim;
	TTree *tEigen = new TTree("Eigen", "eigen");

	int N, D, U, M = 0;
	double MM[nD], VE[nD], VM[nD], VT[nD];

	tEigen->Branch("Dim", &Dim, "iDim/I");

	tEigen->Branch("N", &N, "fN/I");
	tEigen->Branch("D", &D, "fD/I");
	tEigen->Branch("U", &U, "fU/I");
	tEigen->Branch("M", &M, "fM/I");
	
	tEigen->Branch("MM", MM, "fMM[iDim]/D");

	tEigen->Branch("VE", VE, "fVE[iDim]/D");
	tEigen->Branch("VM", VM, "fVM[iDim]/D");
	tEigen->Branch("VT", VT, "fVT[iDim]/D");

	std::mt19937 MT(Seed);
	
	std::uniform_real_distribution<double> nPow(3, 15);
	std::uniform_real_distribution<double> dPow(0,  9);
	//std::uniform_real_distribution<double> uPow(4,  9);

	//std::uniform_real_distribution<double> Val(0, 1);
	std::uniform_real_distribution<double> Phase(0, 2*Pi);

	std::complex<double> n11, n12, n13, n21, n22, n23;
	std::complex<double> d11, d12, d21, d22, d31, d32;
	std::complex<double> u11, u12, u13, u22, u23, u33;
	std::complex<double> m11, m12, m22;

	unsigned int nIter = 0;
	while (nIter < nMAX)
	{
		Dim = nD;

		Pol2Cart(n11, mPow(MT)+Val(MT), Phase(MT));
		Pol2Cart(n12, mPow(MT)+Val(MT), Phase(MT));
		Pol2Cart(n13, mPow(MT)+Val(MT), Phase(MT));
		Pol2Cart(n21, mPow(MT)+Val(MT), Phase(MT));
		Pol2Cart(n22, mPow(MT)+Val(MT), Phase(MT));
		Pol2Cart(n23, mPow(MT)+Val(MT), Phase(MT));

		Pol2Cart(d11, mPow(MT)+Val(MT), Phase(MT));
		Pol2Cart(d12, mPow(MT)+Val(MT), Phase(MT));
		Pol2Cart(d21, mPow(MT)+Val(MT), Phase(MT));
		Pol2Cart(d22, mPow(MT)+Val(MT), Phase(MT));
		Pol2Cart(d31, mPow(MT)+Val(MT), Phase(MT));
		Pol2Cart(d32, mPow(MT)+Val(MT), Phase(MT));

		Pol2Cart(u11, mPow(MT)+Val(MT), Phase(MT));
		Pol2Cart(u12, mPow(MT)+Val(MT), Phase(MT));
		Pol2Cart(u13, mPow(MT)+Val(MT), Phase(MT));
		Pol2Cart(u22, mPow(MT)+Val(MT), Phase(MT));
		Pol2Cart(u23, mPow(MT)+Val(MT), Phase(MT));
		Pol2Cart(u33, mPow(MT)+Val(MT), Phase(MT));

		Pol2Cart(m11, mPow(MT)+Val(MT), Phase(MT));
		Pol2Cart(m12, mPow(MT)+Val(MT), Phase(MT));
		Pol2Cart(m22, mPow(MT)+Val(MT), Phase(MT));

		Eigen::MatrixXcd M(nD, nD);
		M <<	0,	0,	0,	d11,	d12,	0,	0,	0,
		  	0,	0,	0,	d21,	d22,	0,	0,	0,
			0,	0,	0,	d31,	d32,	0,	0,	0,
			d11,	d21,	d31,	m11,	m12,	n11,	n12,	n13,
			d12,	d22,	d32,	m12,	m22,	n21,	n22,	n23,
			0,	0,	0,	n11,	n21,	u11,	u12,	u13,
			0,	0,	0,	n12,	n22,	u12,	u22,	u23,
			0,	0,	0,	n13,	n23,	u13,	u23,	u33;
			
		//V is PMNS matrix
		Eigen::JacobiSVD<Eigen::MatrixXcd> svd(M, Eigen::ComputeFullV | Eigen::ComputeFullU);

		Eigen::MatrixXcd SVDmass(nD, nD);
		Eigen::MatrixXcd SVDdiag(nD, nD);
		Eigen::MatrixXcd UMUdiag(nD, nD);
		Eigen::MatrixXcd MAJphase(nD, nD);
		Eigen::MatrixXcd PMNSreal(nD, nD);
		Print(PMNSreal);
		Eigen::MatrixXcd CorrectD(nD, nD);

		SVDmass = svd.singularValues().asDiagonal();
		SVDdiag = svd.matrixU().adjoint() * M * svd.matrixV();
		UMUdiag = svd.matrixV().transpose() * M * svd.matrixV();

		MajoranaPhase(MAJphase, UMUdiag);
		PMNSreal = svd.matrixV() * MAJphase;
		CorrectD = PMNSreal.transpose() * M * PMNSreal;

		std::cout << " M \n";
		Print(M);
		std::cout << "\n SVDmass \n";
		Print(SVDmass);
		std::cout << "\n SVDdiag \n";
		Print(SVDdiag);
		std::cout << "\n UMUdiag \n";
		Print(UMUdiag);
		std::cout << "\n MAJphase \n";
		Print(MAJphase);
		std::cout << "\n PMNSreal \n";
		Print(PMNSreal);
		std::cout << "\n Correct diag \n";
		Print(CorrectD);
		std::cout << "_________________________" << std::endl;

		Eigen::MatrixXcd Test(nD, nD);
		Test << 1;
		Test << 2;
		Test << 3;
		Test << 4;
		Test << 5;

		/*
		if (Pass(svd))
		{
			for (unsigned int i = 0; i < nD; ++i)
			{
				MM[i] = std::abs(svd.singularValues()[nD-1-i]);
				VE[i] = std::norm(svd.matrixV()(0, nD-1-i));
				VM[i] = std::norm(svd.matrixV()(1, nD-1-i));
				VT[i] = std::norm(svd.matrixV()(2, nD-1-i));
			}

			tEigen->Fill();
		}
		*/

		++nIter;
	}

	//tEigen->Write();

	return 0;
}
