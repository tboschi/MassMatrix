#include <iostream>
#include <cmath>
#include <complex>
#include <random>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "TRandom3.h"

void Pol2Cart(std::complex<double> &w, double lMod, double Phs);

void Populate(std::complex<double> &w, TRandom3 *Gen, double Min, double Max);

int main(int argc, char** argv)
{
	const double Pi = 3.14159;
	const double deg = Pi/180;
	TRandom3 *Ran = new TRandom3();

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
	std::cout << "PMNS matrix is " << std::endl;
	std::cout << PMNS << std::endl;
	std::cout << PMNS.cwiseAbs() << std::endl;

	double mm1 = 0;
	double mm2 = mm1 + dmm21;
	double mm3 = mm1 + dmm31;

	double mme = std::abs(PMNS(1,1)) * mm1 +
		     std::abs(PMNS(1,2)) * mm2 +
		     std::abs(PMNS(1,2)) * mm2 ;

	std::cout << "light masses : " << mm1 << "\t" << mm2 << "\t" << mm3 << std::endl;
	//if (sqrt(mme) < 2.05)
	//	ok;

	std::complex<double> z;			//this is zero (0, 0)
	std::complex<double> d11(1e6,0), d12(2e6,0), d21(3e6,0), d22(4e6,0), d31(5e6,0), d32(6e6,0);
	std::complex<double> n11(1e9,0), n12(3e9,0), n13(5e9,0), n21(6e9,0), n22(2e9,0), n23(8e9,0);
	std::complex<double> u11(4e2,0), u12(3e2,0), u13(5e2,0), u22(2e2,0), u23(6e2,0), u33(1e2,0);
	std::complex<double> m11, m12, m22;	//do not populate, so they are zero

	unsigned int nD = 8;
	unsigned int nIter = 0;
	while (nIter < 1)
	{
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
		std::cout << M.cwiseAbs() << std::endl;
		std::cout << std::endl;
	
		Eigen::MatrixXcd M2 = M.adjoint() * M;
		std::cout << M2 << std::endl;
		std::cout << M2.cwiseAbs() << std::endl;
		std::cout << std::endl;

		Eigen::ComplexEigenSolver<Eigen::MatrixXcd> Ces;
		Ces.compute(M2);
	
		//std::cout << Ces.eigenvalues() << std::endl;
		//std::cout << Ces.eigenvectors().real() << std::endl;
	
		for (unsigned int i = 0 ;  i < nD; ++i)
			std::cout << sqrt(std::abs(Ces.eigenvalues()[i])) << "\t";

		std::cout << std::endl;

		double asd21 = std::abs(Ces.eigenvalues()[1]) - std::abs(Ces.eigenvalues()[0]);
		double asd31 = std::abs(Ces.eigenvalues()[2]) - std::abs(Ces.eigenvalues()[0]);
		std::cout << "diff masses1 : " << dmm21 << "\t" << dmm31 << std::endl;
		std::cout << "diff masses2 : " << asd21 << "\t" << asd31 << std::endl;

		++nIter;
	}

	return 0;
}

void Populate(std::complex<double> &w, TRandom3 *Gen, double Min, double Max)	//radius in log scale
{
	double x, y, r = pow(10, Gen->Uniform(Min, Max));
	Gen->Circle(x, y, r);

	std::cout << "min " << Min << "\t" << r << std::endl;

	w.real(x);
	w.imag(y);
}
