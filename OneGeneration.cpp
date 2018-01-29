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

void Pol2Cart(std::complex<double> &w, double lMod, double Phs);

//void Populate(std::complex<double> &w, TRandom3 *Gen, double Min, double Max);

int main(int argc, char** argv)
{
	const double Pi = 3.14159;
	const double deg = Pi/180;
	//TRandom3 *Ran = new TRandom3();

	const struct option longopts[] = 
	{
		{"output", 	required_argument,	0, 'o'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;

	std::ofstream OutFile;
	int dd = 9, nn = 12, mm = 0, uu = 6;
	double CL = 0.90;	//90% C.L.

	while((iarg = getopt_long(argc,argv, "d:n:m:u:o:h", longopts, &index)) != -1)
	{
		switch(iarg)
		{
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
	std::uniform_real_distribution<double> Phase(0, 2*Pi);
	std::uniform_real_distribution<double> dExp(dd, dd+2);
	std::uniform_real_distribution<double> nExp(nn, nn+2);
	std::uniform_real_distribution<double> uExp(uu, uu+2);
	std::uniform_real_distribution<double> mExp(mm, mm+2);

	unsigned int nD = 3;
	unsigned int nIter = 0;
	while (nIter < 10000)
	{
		std::complex<double> ddd;
		std::complex<double> nnn;
		std::complex<double> uuu;
		std::complex<double> mmm;
		Pol2Cart(ddd, dExp(MT), Phase(MT));
		Pol2Cart(nnn, nExp(MT), Phase(MT));
		Pol2Cart(uuu, uExp(MT), Phase(MT));
		Pol2Cart(mmm, mExp(MT), Phase(MT));

		Out << std::abs(ddd) << "\t" << std::abs(nnn) << "\t";
		Out << std::abs(uuu) << "\t" << std::abs(mmm) << "\t";
		Eigen::MatrixXcd M(nD,nD);
		M <<	0,	ddd,	0,
			ddd,	mmm,	nnn,
			0,	nnn,	uuu;
			
	
		//Diagonalise;
		std::cout << M << std::endl;
		//std::cout << M.cwiseAbs() << std::endl << std::endl;
		std::cout << std::endl;
	
		Eigen::MatrixXcd M2 = M * M.adjoint();
		std::cout << M2 << std::endl;
		//std::cout << M2.cwiseAbs() << std::endl;
		std::cout << std::endl;

		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> Ces(M2);
	
		//std::cout << Ces.eigenvalues() << std::endl;
		//std::cout << Ces.eigenvectors().real() << std::endl;
	
		for (unsigned int i = 0 ;  i < nD; ++i)
			Out << sqrt(std::abs(Ces.eigenvalues()[i])) << "\t";
		Out << std::endl;
		//std::cout << std::abs(Ces.eigenvalues()[2]) << "\t" << std::abs(Ces.eigenvalues()[1]) << std::endl;

		//double asd21 = std::abs(Ces.eigenvalues()[1]) - std::abs(Ces.eigenvalues()[0]);
		//double asd31 = std::abs(Ces.eigenvalues()[2]) - std::abs(Ces.eigenvalues()[0]);
		//std::cout << "diff masses1 : " << dmm21 << "\t" << dmm31 << std::endl;
		//std::cout << "diff masses2 : " << asd21 << "\t" << asd31 << std::endl;
		//std::cout << Ces.eigenvalues()[0] << "\t"; 
		//std::cout << Ces.eigenvalues()[3] << "\t"; 
		//std::cout << std::abs(Ces.eigenvectors()(0,3)) << "\t"; 
		//std::cout << std::abs(Ces.eigenvectors()(3,0)) << std::endl;

		++nIter;
	}

	return 0;
}

void Pol2Cart(std::complex<double> &w, double lMod, double Phs)
{
	w.real(pow(10, lMod)*cos(Phs));
	w.imag(pow(10, lMod)*sin(Phs));
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
