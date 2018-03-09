#include <iostream>
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

double Gloop(double x)
{
	return - 0.25 * (2*x*x*x + 5*x*x - x) / pow((1-x), 3) -
	       - 1.50 * (x*x*x) / pow(1-x, 4) * log(x);
}

template <typename Derived>
bool PassFIT(const Eigen::MatrixBase<Derived> &V, double *mm)
{
	bool M21 = (mm[1]-mm[0] > 6.8e-5	&& mm[1]-mm[0] < 8.02e-5);
	bool M31 = (mm[2]-mm[0] > 2.399e-3	&& mm[2]-mm[0] < 2.593e-3);
	bool Ms4 = (sqrt(mm[3]) > 1e5		&& sqrt(mm[3]) < 5e9);

	/*
	Eigen::Matrix3d PMNSmin, PMNSmax, Vpmns;
	PMNSmin <<	0.76,	0.50,	0.13,
			0.21,	0.42,	0.61,
			0.18,	0.38,	0.40;
	PMNSmax <<	0.85,	0.60,	0.16,
			0.54,	0.70,	0.79,
			0.58,	0.72,	0.78;

	bool PMNS = true;
	for (unsigned int c = 0; c < PMNSmin.cols(); ++c)
		for (unsigned int r = 0; r < PMNSmin.rows(); ++r)
			PMNS *= V.cwiseAbs()(r, c) > PMNSmin(r, c) && V.cwiseAbs()(r, c) < PMNSmax(r, c);
	*/

	return M21 && M31 && Ms4;
}

template <typename Derived>
bool PassMEG(const Eigen::MatrixBase<Derived> &V, double *mm)
{
	std::complex<double> MEGamp;
	for (unsigned int i = 0; i < V.cols(); ++i)
		MEGamp += std::conj(V(0, i)) * V(1, i) * Gloop( mm[i] / pow(Const::fMW, 2) );
	double MEGbranch = 3 * Const::fAem * std::norm(MEGamp) / (32 * Const::fPi);

	//bool MEG = MEGbranch < 4.2e-13;		//present
	bool MEG = MEGbranch < 5e-14;		//future

	return MEG;
}

template <typename Derived>
bool PassBB0(const Eigen::MatrixBase<Derived> &V, double *mm)
{
	double p2 = -pow(125e6, 2);
	std::complex<double> BBeff;
	for (unsigned int i = 0; i < V.cols(); ++i)
		BBeff += V(0, i) * V(0, i) * p2 * sqrt(mm[i]) / (p2 - mm[i]);

	//bool BB0 = std::abs(BBeff) < 150e-3;	//present
	bool BB0 = std::abs(BBeff) < 20e-3;	//future

	return BB0;
}

template <typename Derived>
bool PassNSI(const Eigen::MatrixBase<Derived> &V, double *mm)
{
	Eigen::Matrix3d NSIabove, NSIbelow;
	Eigen::Matrix3cd Kab;
	NSIabove <<	4.0e-3,	1.2e-4,	3.2e-3,
			1.2e-4,	1.6e-3,	2.1e-3,
			3.2e-3,	2.1e-3,	5.3e-3;
	NSIbelow <<	4.0e-3,	1.8e-3,	3.2e-3,
		 	1.8e-3,	1.6e-3,	2.1e-3,
			3.2e-3,	2.1e-3,	5.3e-3;
	for (unsigned int i = 4; i < V.cols(); ++i)
	{
		Kab(0,0) += V(0, i) * std::conj(V(0, i));
		Kab(0,1) += V(0, i) * std::conj(V(1, i));
		Kab(0,2) += V(0, i) * std::conj(V(2, i));
		Kab(1,0) += V(1, i) * std::conj(V(0, i));
		Kab(1,1) += V(1, i) * std::conj(V(1, i));
		Kab(1,2) += V(1, i) * std::conj(V(2, i));
		Kab(2,0) += V(2, i) * std::conj(V(0, i));
		Kab(2,1) += V(2, i) * std::conj(V(1, i));
		Kab(2,2) += V(2, i) * std::conj(V(2, i));
	}

	bool NSI = true;
	bool EWscale = mm[4] > 246e9 && mm[5] > 249e9 && mm[6] > 249e9 && mm[7] > 249e9;
	for (unsigned int c = 0; c < Kab.cols(); ++c)
		for (unsigned int r = 0; r < Kab.rows(); ++r)
			NSI *= EWscale ? Kab.cwiseAbs()(r, c) < NSIabove(r, c) : Kab.cwiseAbs()(r, c) < NSIbelow(r, c);
	return true;
}

template <typename Derived>
bool Pass(const Eigen::JacobiSVD<Derived> &S, bool &vbb0)	//SVD are in descending order
{
	const unsigned int nS = S.singularValues().size();
	double mm[nS];
	Eigen::MatrixXcd V(nS, nS);
	for (unsigned int i = 0; i < nS; ++i)
	{
		mm[i] = pow(S.singularValues()[nS-1-i], 2);
		V.col(i) = S.matrixV().col(nS-1-i);
	}	

	bool temp = PassFIT(V, mm) && PassNSI(V, mm) &&	PassMEG(V, mm);
	if (temp)
		vbb0 = PassBB0(V, mm);
	return temp;

	/*
	return PassFIT(V, mm) &&
	       PassBB0(V, mm) &&	
	       PassMEG(V, mm) &&	
	       PassNSI(V, mm) ;
	*/
}

void Pol2Cart(std::complex<double> &w, double lMod, double Phs)
{
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

	int N[6], D[6], U[6], M[6];
	double MM[nD], VE[nD], VM[nD], VT[nD];
	bool BB0;

	tEigen->Branch("Dim", &Dim, "iDim/I");
	tEigen->Branch("BB0", &BB0, "bDim/O");

	tEigen->Branch("N", N, "fN[6]/I");
	tEigen->Branch("D", D, "fD[6]/I");
	tEigen->Branch("U", U, "fU[6]/I");
	tEigen->Branch("M", M, "fM[6]/I");
	
	tEigen->Branch("MM", MM, "fMM[iDim]/D");

	tEigen->Branch("VE", VE, "fVE[iDim]/D");
	tEigen->Branch("VM", VM, "fVM[iDim]/D");
	tEigen->Branch("VT", VT, "fVT[iDim]/D");

	std::mt19937 MT(Seed);
	
	std::uniform_int_distribution<int> nPow( 3, 15);
	std::uniform_int_distribution<int> dPow( 3, 15);
	std::uniform_int_distribution<int> uPow(-9,  9);
	std::uniform_int_distribution<int> mPow(-9,  9);

	std::uniform_real_distribution<double> Val(0, 1);
	std::uniform_real_distribution<double> Phase(0, 2*Pi);

	std::complex<double> n11, n12, n13, n21, n22, n23;
	std::complex<double> d11, d12, d21, d22, d31, d32;
	std::complex<double> u11, u12, u13, u22, u23, u33;
	std::complex<double> m11, m12, m22;

	unsigned int nIter = 0;
	while (nIter < nMAX)
	{
		Dim = nD;

		//naturalness conditions
		//U < N and (D > N-6 & D < N)
		//if ((U > N) || (N - D < 3) || (N - D > 6))
		//	continue;

		Pol2Cart(n11, nPow(MT)+Val(MT), Phase(MT));
		Pol2Cart(n12, nPow(MT)+Val(MT), Phase(MT));
		Pol2Cart(n13, nPow(MT)+Val(MT), Phase(MT));
		Pol2Cart(n21, nPow(MT)+Val(MT), Phase(MT));
		Pol2Cart(n22, nPow(MT)+Val(MT), Phase(MT));
		Pol2Cart(n23, nPow(MT)+Val(MT), Phase(MT));

		Pol2Cart(d11, dPow(MT)+Val(MT), Phase(MT));
		Pol2Cart(d12, dPow(MT)+Val(MT), Phase(MT));
		Pol2Cart(d21, dPow(MT)+Val(MT), Phase(MT));
		Pol2Cart(d22, dPow(MT)+Val(MT), Phase(MT));
		Pol2Cart(d31, dPow(MT)+Val(MT), Phase(MT));
		Pol2Cart(d32, dPow(MT)+Val(MT), Phase(MT));

		Pol2Cart(u11, uPow(MT)+Val(MT), Phase(MT));
		Pol2Cart(u12, uPow(MT)+Val(MT), Phase(MT));
		Pol2Cart(u13, uPow(MT)+Val(MT), Phase(MT));
		Pol2Cart(u22, uPow(MT)+Val(MT), Phase(MT));
		Pol2Cart(u23, uPow(MT)+Val(MT), Phase(MT));
		Pol2Cart(u33, uPow(MT)+Val(MT), Phase(MT));

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
		Eigen::JacobiSVD<Eigen::MatrixXcd> svd(M, Eigen::ComputeFullV);

		D[0] = log10(std::abs(d11));
		D[1] = log10(std::abs(d12));
		D[2] = log10(std::abs(d21));
		D[3] = log10(std::abs(d22));
		D[4] = log10(std::abs(d31));
		D[5] = log10(std::abs(d32));

		N[0] = log10(std::abs(n11));
		N[1] = log10(std::abs(n12));
		N[2] = log10(std::abs(n13));
		N[3] = log10(std::abs(n21));
		N[4] = log10(std::abs(n22));
		N[5] = log10(std::abs(n23));

		U[0] = log10(std::abs(u11));
		U[1] = log10(std::abs(u12));
		U[2] = log10(std::abs(u13));
		U[3] = log10(std::abs(u22));
		U[4] = log10(std::abs(u23));
		U[5] = log10(std::abs(u33));

		M[0] = log10(std::abs(m11));
		M[1] = log10(std::abs(m12));
		M[2] = log10(std::abs(m22));
		M[3] = 0;
		M[4] = 0;
		M[5] = 0;

		if (Pass(svd, BB0))
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

		++nIter;
	}

	tEigen->Write();

	return 0;
}
