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

//#include "TTree.h"

//from nufit 18
template <typename Derived>
bool Pass(const Eigen::JacobiSVD<Derived> &S)	//SVD are in descending order
{
	unsigned int nS = S.singularValues().size();
	double m1 = S.singularValues()[nS-1];
	double m2 = S.singularValues()[nS-2];
	double m3 = S.singularValues()[nS-3];
	double m4 = S.singularValues()[nS-4];

	double mm1 = m1*m1;
	double mm2 = m2*m2;
	double mm3 = m3*m3;

	//should be ok!
	double uue1 = std::abs(S.matrixV()(0, nS-1));
	double uue2 = std::abs(S.matrixV()(0, nS-2));
	double uue3 = std::abs(S.matrixV()(0, nS-3));
	double uum1 = std::abs(S.matrixV()(1, nS-1));
	double uum2 = std::abs(S.matrixV()(1, nS-2));
	double uum3 = std::abs(S.matrixV()(1, nS-3));
	double uut1 = std::abs(S.matrixV()(2, nS-1));
	double uut2 = std::abs(S.matrixV()(2, nS-2));
	double uut3 = std::abs(S.matrixV()(2, nS-3));

	
	Eigen::Matrix3d PMNSmin, PMNSmax, Vpmns;
	PMNSmin <<	0.76,	0.50,	0.13,
			0.21,	0.42,	0.61,
			0.18,	0.38,	0.40;
	PMNSmax <<	0.85,	0.60,	0.16,
			0.54,	0.70,	0.79,
			0.58,	0.72,	0.78;
	Vpmns <<	uue1, 	uue2, 	uue3,
	      		uum1, 	uum2, 	uum3,
	      		uut1, 	uut2, 	uut3;

	std::cout << "pmns check" << std::endl;
	std::cout << PMNSmin << std::endl;
	std::cout << std::endl;
	std::cout << PMNSmax << std::endl;
	std::cout << std::endl;
	std::cout << Vpmns << std::endl;
	std::cout << "____________________________________\n" << std::endl;

	bool Ue1 = (uue1 > PMNSmin(0,0) && uue1 < PMNSmax(0,0));
	bool Ue2 = (uue2 > PMNSmin(0,1) && uue2 < PMNSmax(0,1));
	bool Ue3 = (uue3 > PMNSmin(0,2) && uue3 < PMNSmax(0,2));
	bool Um1 = (uum1 > PMNSmin(1,0) && uum1 < PMNSmax(1,0));
	bool Um2 = (uum2 > PMNSmin(1,1) && uum2 < PMNSmax(1,1));
	bool Um3 = (uum3 > PMNSmin(1,2) && uum3 < PMNSmax(1,2));
	bool Ut1 = (uut1 > PMNSmin(2,0) && uut1 < PMNSmax(2,0));
	bool Ut2 = (uut2 > PMNSmin(2,1) && uut2 < PMNSmax(2,1));
	bool Ut3 = (uut3 > PMNSmin(2,2) && uut3 < PMNSmax(2,2));
	
	bool M21 = (mm2-mm1 > 6.8e-5	&& mm2-mm1 < 8.02e-5);
	bool M31 = (mm3-mm1 > 2.399e-3	&& mm3-mm1 < 2.593e-3);
	bool Ms4 = (m4 > 1e5		&& m4 < 5e9);
	/*
	if (M21 && M31 && Ms4)
	{
		Eigen::Matrix3d Vpmns;
		std::cout << std::endl;
		std::cout << S.matrixV().cwiseAbs() << std::endl;
		std::cout << std::endl;
		std::cout << Vpmns << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
	}
	*/
	
	return  M21 &&	M31 &&	Ms4;	//&&	
	//	Ue1 &&	Ue2 &&	Ue3 && 
	//	Um1 &&	Um2 &&	Um3 && 
	//	Ut1 &&	Ut2 &&	Ut3 ;
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
	int dd = 9, nn = 12, mm = 0, uu = 6;
	unsigned int nMAX = 10000;
	unsigned long Seed = std::chrono::system_clock::now().time_since_epoch()/std::chrono::milliseconds(1);
	double CL = 0.90;	//90% C.L.

	while((iarg = getopt_long(argc,argv, "I:d:n:m:u:S:o:h", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'I':
				nMAX = strtol(optarg, NULL, 10);
				break;
			case 'd':
				dd = strtol(optarg, NULL, 10);
				break;
			case 'n':
				nn = strtol(optarg, NULL, 10);
				break;
			case 'u':
				uu = strtol(optarg, NULL, 10);
				break;
			case 'm':
				mm = strtol(optarg, NULL, 10);
				break;
			case 'S':
				Seed = strtol(optarg, NULL, 10);
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

	/*
	double dmm21 = 7.4e-5;		//eV2
	double dmm31 = 2.494e-3;		//eV2

	double t12 = 33.62 * deg;
	double t23 = 47.2  * deg;
	double t13 = 8.54  * deg;
	double cp  = 234   * deg;
	std::complex<double> dcp(cos(cp), -sin(cp));
	Eigen::Matrix3cd U1, U2, U3, PMNS;

	U1 <<	1, 		0, 		0,
	   	0, 		cos(t23),	sin(t23),
		0, 		-sin(t23), 	cos(t23);

	U2 <<	cos(t13),	0, 		sin(t23)*dcp,
	   	0, 		1,		0,
		-sin(t13)*dcp, 	0,		cos(t13);

	U3 <<	cos(t12), 	sin(t12),	0,
	   	-sin(t12), 	cos(t12),	0,
		0,		0,		1;

	PMNS = U1 * U2 * U3;
	//std::cout << "PMNS matrix is " << std::endl;
	//std::cout << PMNS << std::endl;
	//std::cout << PMNS.cwiseAbs() << std::endl;

	double mm1 = 0;
	double mm2 = mm1 + dmm21;
	double mm3 = mm1 + dmm31;

	double mme = std::abs(PMNS(1,1)) * mm1 +
		     std::abs(PMNS(1,2)) * mm2 +
		     std::abs(PMNS(1,2)) * mm2 ;

	//std::cout << "light masses : " << mm1 << "\t" << mm2 << "\t" << mm3 << std::endl;
	//if (sqrt(mme) < 2.05)
	//	ok;
	*/

	const unsigned int nD = 8;
	//typedef Eigen::Matrix<std::complex<double>, nD, nD> LDMatrixXcd;

	std::mt19937 MT(Seed);
	
	//U < D < N	is a good condition
	std::uniform_int_distribution<int> nPow( 7, 15);
	std::uniform_int_distribution<int> dPow( 3, 11);
	std::uniform_int_distribution<int> uPow( 4, 10);
	std::uniform_int_distribution<int> mPow(-3,  3);

	std::uniform_real_distribution<double> Val(0, 1);
	std::uniform_real_distribution<double> Phase(0, 2*Pi);

	std::complex<double> d11, d12, d21, d22, d31, d32;
	std::complex<double> n11, n12, n13, n21, n22, n23;
	std::complex<double> u11, u12, u13, u22, u23, u33;
	std::complex<double> m11, m12, m22;

	unsigned int nIter = 0;
	while (nIter < nMAX)
	{
		double dBase = dPow(MT);
		double nBase = nPow(MT);
		double uBase = uPow(MT);
		double mBase = mPow(MT);

		//naturalness conditions
		if (dBase > nBase)
			continue;
		if (uBase > nBase)
			continue;
		if (nBase - dBase > 6)
			continue;

		Pol2Cart(d11, dBase+Val(MT), Phase(MT));
		Pol2Cart(d12, dBase+Val(MT), Phase(MT));
		Pol2Cart(d21, dBase+Val(MT), Phase(MT));
		Pol2Cart(d22, dBase+Val(MT), Phase(MT));
		Pol2Cart(d31, dBase+Val(MT), Phase(MT));
		Pol2Cart(d32, dBase+Val(MT), Phase(MT));
		Pol2Cart(n11, nBase+Val(MT), Phase(MT));
		Pol2Cart(n12, nBase+Val(MT), Phase(MT));
		Pol2Cart(n13, nBase+Val(MT), Phase(MT));
		Pol2Cart(n21, nBase+Val(MT), Phase(MT));
		Pol2Cart(n22, nBase+Val(MT), Phase(MT));
		Pol2Cart(n23, nBase+Val(MT), Phase(MT));
		Pol2Cart(u11, uBase+Val(MT), Phase(MT));
		Pol2Cart(u12, uBase+Val(MT), Phase(MT));
		Pol2Cart(u13, uBase+Val(MT), Phase(MT));
		Pol2Cart(u22, uBase+Val(MT), Phase(MT));
		Pol2Cart(u23, uBase+Val(MT), Phase(MT));
		Pol2Cart(u33, uBase+Val(MT), Phase(MT));
		Pol2Cart(m11, mBase+Val(MT), Phase(MT));
		Pol2Cart(m12, mBase+Val(MT), Phase(MT));
		Pol2Cart(m22, mBase+Val(MT), Phase(MT));

		//Eigen::MatrixXcd M0(nD,nD), dM(nD,nD);
		Eigen::MatrixXcd M(nD, nD), M2(nD, nD);
		M <<	0,	0,	0,	d11,	d12,	0,	0,	0,
		  	0,	0,	0,	d21,	d22,	0,	0,	0,
			0,	0,	0,	d31,	d32,	0,	0,	0,
			d11,	d21,	d31,	m11,	m12,	n11,	n12,	n13,
			d12,	d22,	d32,	m12,	m22,	n21,	n22,	n23,
			0,	0,	0,	n11,	n21,	u11,	u12,	u13,
			0,	0,	0,	n12,	n22,	u12,	u22,	u23,
			0,	0,	0,	n13,	n23,	u13,	u23,	u33;
			
		M2 = M.adjoint() * M;
	
		//Eigen::MatrixXcd M2 = M0.adjoint() * M0 + 
		//		      dM.adjoint()*M0 + M0.adjoint()*dM ;//+
				      //dM.adjoint()*dM;

		Eigen::JacobiSVD<Eigen::MatrixXcd> svd(M, Eigen::ComputeFullV);		//V is PMNS matrix..?
	
		if (Pass(svd))
		{
			Out << dBase << "\t"; 
			Out << nBase << "\t"; 
			Out << uBase << "\t" ;
			Out << mBase << "\t";
			for (unsigned int i = 0; i < nD; ++i)
				Out << std::abs(svd.singularValues()[i]) << "\t"; 
			Out << std::abs(svd.matrixV()(0, nD-4)) << "\t";	//e4
			Out << std::abs(svd.matrixV()(1, nD-4)) << "\t";	//m4
			Out << std::abs(svd.matrixV()(2, nD-4)) << "\t";	//t4
			Out << std::endl;
		}

		++nIter;
	}

	return 0;
}
