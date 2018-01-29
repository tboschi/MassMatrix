#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <random>
#include <chrono>
#include <getopt.h>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

//#include "TRandom3.h"

//from nufit 18
template <typename Derived>
//bool Pass(const Eigen::MatrixBase<Derived1> &U, const Eigen::ArrayBase<Derived2> &L)
bool Pass(const Eigen::SelfAdjointEigenSolver<Derived> &C)
{
	double mm1 = C.eigenvalues()[0];
	double mm2 = C.eigenvalues()[1];
	double mm3 = C.eigenvalues()[2];

	if (mm1 < 0 || mm2 < 0 || mm3 < 0)
		return false;
	else
	{


	/*
	double uue1 = std::abs(C.eigenvectors()(0,0));
	double uue2 = std::abs(C.eigenvectors()(0,1));
	double uue3 = std::abs(C.eigenvectors()(0,2));
	double uum1 = std::abs(C.eigenvectors()(1,0));
	double uum2 = std::abs(C.eigenvectors()(1,1));
	double uum3 = std::abs(C.eigenvectors()(1,2));
	double uut1 = std::abs(C.eigenvectors()(2,0));
	double uut2 = std::abs(C.eigenvectors()(2,1));
	double uut3 = std::abs(C.eigenvectors()(2,2));
	*/

	
	/*
	Eigen::Matrix3d PMNSmin, PMNSmax;
	PMNSmin <<	0.76,	0.50,	0.13,
			0.21,	0.42,	0.61,
			0.18,	0.38,	0.40;
	PMNSmax <<	0.85,	0.60,	0.16,
			0.54,	0.70,	0.79,
			0.58,	0.72,	0.78;
	bool Ue1 = (uue1 > PMNSmin(0,0) && uue1 < PMNSmax(0,0));
	bool Ue2 = (uue2 > PMNSmin(0,1) && uue2 < PMNSmax(0,1));
	bool Ue3 = (uue3 > PMNSmin(0,2) && uue3 < PMNSmax(0,2));
	bool Um1 = (uum1 > PMNSmin(1,0) && uum1 < PMNSmax(1,0));
	bool Um2 = (uum2 > PMNSmin(1,1) && uum2 < PMNSmax(1,1));
	bool Um3 = (uum3 > PMNSmin(1,2) && uum3 < PMNSmax(1,2));
	bool Ut1 = (uut1 > PMNSmin(2,0) && uut1 < PMNSmax(2,0));
	bool Ut2 = (uut2 > PMNSmin(2,1) && uut2 < PMNSmax(2,1));
	bool Ut3 = (uut3 > PMNSmin(2,2) && uut3 < PMNSmax(2,2));
	*/
	
		bool M21 = (mm2-mm1 > 6.8e-5	&& mm2-mm1 < 8.02e-5);
		bool M31 = (mm3-mm1 > 2.399e-3	&& mm3-mm1 < 2.593e-3);
		return  M21 &&	M31 ;
		//	Ue1 &&	Ue2 &&	Ue3 && 
		//	Um1 &&	Um2 &&	Um3 && 
		//	Ut1 &&	Ut2 &&	Ut3 ;
	}
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
	double CL = 0.90;	//90% C.L.

	while((iarg = getopt_long(argc,argv, "I:d:n:m:u:o:h", longopts, &index)) != -1)
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

	unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::mt19937 MT(seed);
	
	//std::uniform_int_distribution<int> dExp(3, 9);
	//std::uniform_int_distribution<int> nExp(6, 15);
	//std::uniform_int_distribution<int> uExp(3, 9);
	//std::uniform_int_distribution<int> mExp(-5, 0);

	std::uniform_int_distribution<int> Pow(-5, 15);
	std::uniform_real_distribution<double> Val(0, 1);
	std::uniform_real_distribution<double> Phase(0, 2*Pi);

	std::complex<double> d11, d12, d21, d22, d31, d32;
	std::complex<double> n11, n12, n13, n21, n22, n23;
	std::complex<double> u11, u12, u13, u22, u23, u33;
	std::complex<double> m11, m12, m22;

	unsigned int nD = 8;
	unsigned int nIter = 0;
	while (nIter < nMAX)
	{
		double dPow = Pow(MT);
		Pol2Cart(d11, dPow+Val(MT), Phase(MT));
		Pol2Cart(d12, dPow+Val(MT), Phase(MT));
		Pol2Cart(d21, dPow+Val(MT), Phase(MT));
		Pol2Cart(d22, dPow+Val(MT), Phase(MT));
		Pol2Cart(d31, dPow+Val(MT), Phase(MT));
		Pol2Cart(d32, dPow+Val(MT), Phase(MT));
		double nPow = Pow(MT);
		Pol2Cart(n11, nPow+Val(MT), Phase(MT));
		Pol2Cart(n12, nPow+Val(MT), Phase(MT));
		Pol2Cart(n13, nPow+Val(MT), Phase(MT));
		Pol2Cart(n21, nPow+Val(MT), Phase(MT));
		Pol2Cart(n22, nPow+Val(MT), Phase(MT));
		Pol2Cart(n23, nPow+Val(MT), Phase(MT));
		double uPow = Pow(MT);
		Pol2Cart(u11, uPow+Val(MT), Phase(MT));
		Pol2Cart(u12, uPow+Val(MT), Phase(MT));
		Pol2Cart(u13, uPow+Val(MT), Phase(MT));
		Pol2Cart(u22, uPow+Val(MT), Phase(MT));
		Pol2Cart(u23, uPow+Val(MT), Phase(MT));
		Pol2Cart(u33, uPow+Val(MT), Phase(MT));
		double mPow = Pow(MT);
		Pol2Cart(m11, mPow+Val(MT), Phase(MT));
		Pol2Cart(m12, mPow+Val(MT), Phase(MT));
		Pol2Cart(m22, mPow+Val(MT), Phase(MT));

		Eigen::MatrixXcd M(nD,nD);
		M <<	0,	0,	0,	d11,	d12,	0,	0,	0,
		  	0,	0,	0,	d21,	d22,	0,	0,	0,
			0,	0,	0,	d31,	d32,	0,	0,	0,
			d11,	d21,	d31,	m11,	m12,	n11,	n12,	n13,
			d12,	d22,	d32,	m12,	m22,	n21,	n22,	n23,
			0,	0,	0,	n11,	n21,	u11,	u12,	u13,
			0,	0,	0,	n12,	n22,	u12,	u22,	u23,
			0,	0,	0,	n13,	n23,	u13,	u23,	u33;
			
	
		//Diagonalise;
		//std::cout << M << std::endl;
		//std::cout << M.cwiseAbs() << std::endl << std::endl;
		//std::cout << std::endl;
	
		Eigen::MatrixXcd M2 = M * M.adjoint();
		//std::cout << M2 << std::endl;
		//std::cout << M2.cwiseAbs() << std::endl;
		//std::cout << std::endl;

		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> Ces(M2);
	
		//std::cout << Ces.eigenvalues() << std::endl;
		//std::cout << Ces.eigenvectors().real() << std::endl;
	
		//Out << nPow << "\t" << dPow << "\t" << uPow << "\t" << mPow << "\t";
		//for (unsigned int i = 0 ;  i < nD; ++i)
		//	Out << sqrt(std::abs(Ces.eigenvalues()[i])) << "\t";
		//std::cout << std::abs(Ces.eigenvalues()[2]) << "\t" << std::abs(Ces.eigenvalues()[1]) << std::endl;

		//double asd21 = std::abs(Ces.eigenvalues()[1]) - std::abs(Ces.eigenvalues()[0]);
		//double asd31 = std::abs(Ces.eigenvalues()[2]) - std::abs(Ces.eigenvalues()[0]);
		//std::cout << "diff masses1 : " << dmm21 << "\t" << dmm31 << std::endl;
		//std::cout << "diff masses2 : " << asd21 << "\t" << asd31 << std::endl;
		//std::cout << Ces.eigenvalues()[0] << "\t"; 
		//std::cout << Ces.eigenvalues()[3] << "\t"; 
		//std::cout << std::abs(Ces.eigenvectors()(0,3)) << "\t"; 
		//std::cout << std::abs(Ces.eigenvectors()(3,0)) << std::endl;
		if (Pass(Ces))
			Out << dPow << "\t" << nPow << "\t" << uPow << "\t" << mPow << std::endl;

		++nIter;
	}

	return 0;
}

/*
void Populate(std::complex<double> &w, TRandom3 *Gen, double Min, double Max)	//radius in log scale
{
	double x, y, r = pow(10, Gen->Uniform(Min, Max));
	Gen->Circle(x, y, r);

	std::cout << "min " << Min << "\t" << r << std::endl;

	w.real(x);
	w.imag(y);
}
*/
