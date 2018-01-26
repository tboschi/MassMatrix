#include <iostream>
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
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;

	int dd = 9, nn = 12, mm = 0, uu = 6;
	double CL = 0.90;	//90% C.L.

	while((iarg = getopt_long(argc,argv, "d:n:m:u:h", longopts, &index)) != -1)
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
			case 'h':
				return 1;
			default:
				break;
		}
	}

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
	std::uniform_real_distribution<double> dExp(6.0,7.0);
	std::uniform_real_distribution<double> nExp(7.0,8.0);
	std::uniform_real_distribution<double> uExp(-1.0,1.0);
	std::uniform_real_distribution<double> mExp(-1.0,1.0);

	unsigned int nD = 8;
	unsigned int nIter = 0;
	while (nIter < 10000)
	{
	std::complex<double> d11(pow(10,dExp(MT)), pow(10,dExp(MT))),
			     d12(pow(10,dExp(MT)), pow(10,dExp(MT))),
			     d21(pow(10,dExp(MT)), pow(10,dExp(MT))),
			     d22(pow(10,dExp(MT)), pow(10,dExp(MT))),
			     d31(pow(10,dExp(MT)), pow(10,dExp(MT))),
			     d32(pow(10,dExp(MT)), pow(10,dExp(MT)));
	std::complex<double> n11(pow(10,nExp(MT)), pow(10,nExp(MT))),
			     n12(pow(10,nExp(MT)), pow(10,nExp(MT))),
			     n13(pow(10,nExp(MT)), pow(10,nExp(MT))),
			     n21(pow(10,nExp(MT)), pow(10,nExp(MT))),
			     n22(pow(10,nExp(MT)), pow(10,nExp(MT))),
			     n23(pow(10,nExp(MT)), pow(10,nExp(MT)));
	std::complex<double> u11(pow(10,uExp(MT)), pow(10,uExp(MT))),
			     u12(pow(10,uExp(MT)), pow(10,uExp(MT))),
			     u13(pow(10,uExp(MT)), pow(10,uExp(MT))),
			     u22(pow(10,uExp(MT)), pow(10,uExp(MT))),
			     u23(pow(10,uExp(MT)), pow(10,uExp(MT))),
			     u33(pow(10,uExp(MT)), pow(10,uExp(MT)));
	std::complex<double> m11(pow(10,mExp(MT)), pow(10,mExp(MT))),
			     m12(pow(10,mExp(MT)), pow(10,mExp(MT))),
			     m22(pow(10,mExp(MT)), pow(10,mExp(MT)));
		/*
		Populate(d11, Ran, 11, 12);
		Populate(d12, Ran, 11, 12);
		Populate(d21, Ran, 11, 12);
		Populate(d22, Ran, 11, 12);
		Populate(d31, Ran, 11, 12);
		Populate(d32, Ran, 11, 12);

		Populate(n11, Ran, 11, 12);
		Populate(n12, Ran, 11, 12);
		//Populate(n13, Ran, 11, 12);
		Populate(n21, Ran, 11, 12);
		Populate(n22, Ran, 11, 12);
		//Populate(n23, Ran, 11, 12);

		Populate(u11, Ran, 6, 7);
		Populate(u12, Ran, 6, 7);
		//Populate(u13, Ran, 6, 7);
		Populate(u22, Ran, 6, 7);
		//Populate(u23, Ran, 6, 7);
		//Populate(u33, Ran, 6, 7);
	
		Populate(m11, Ran, -1, 0);
		Populate(m12, Ran, -1, 0);
		Populate(m22, Ran, -1, 0);
		*/

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
	
		//for (unsigned int i = 0 ;  i < nD; ++i)
		//	std::cout << Ces.eigenvalues()[i] << "\t";
		//std::cout << std::endl;
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
