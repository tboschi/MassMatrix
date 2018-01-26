#include <iostream>
#include <cmath>
#include <complex>
#include <random>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

std::complex<double> Pol2Cart(double lMod, double Phs);
int main(int argc, char** argv)
{
	const double Pi = 3.14159;
	const double deg = Pi/180;
	std::random_device Ran;
	std::mt19937_64 Gen(Ran());
	std::uniform_real_distribution<double> Phase(0, 2*Pi);
	std::uniform_real_distribution<double> lMod_d(6, 7);	//log eV
	std::uniform_real_distribution<double> lMod_n(7, 8);	//log eV
	std::uniform_real_distribution<double> lMod_m(-1, 1);	//log eV
	std::uniform_real_distribution<double> lMod_u(-1, 1);	//log eV
	unsigned int nIter = 0;

	double dmm21 = 7.4e-5;		//eV2
	double dmm31 = 2.494e-3;		//eV2

	double t12 = 33.62 * deg;
	double t23 = 47.2  * deg;
	double t13 = 8.54  * deg;
	double cp  = 234   * deg;
	std::complex<double> dcp(cos(cp), -sin(cp));

	std::complex<double> U11(cos(t12)*cos(t13), 0);
	std::complex<double> U12(sin(t12)*cos(t13), 0);
	std::complex<double> U13(sin(t13)*dcp.real(), sin(t13)*dcp.imag());
	std::complex<double> U21(-sin(t12)*cos(23)-cos(t12)*cos(t12)*cos(t13), 0);
	std::complex<double> U22(cos(t12)*cos(t13), 0);
	std::complex<double> U23(cos(t12)*cos(t13), 0);
	std::complex<double> U31(cos(t12)*cos(t13), 0);
	std::complex<double> U32(cos(t12)*cos(t13), 0);
	std::complex<double> U33(cos(t12)*cos(t13), 0);
	std::complex<double> z;					//this is zero
	std::complex<double> p1, p2, p3; 			//theta1, 2, 3
	std::complex<double> d11, d12, d21, d22, d31, d32;
	std::complex<double> m11, m12, m22;
	std::complex<double> n11, n12, n13, n21, n22, n23;
	std::complex<double> MR1, MR2, MR3;
	std::complex<double> u11, u12, u13, u22, u23, u33;	//muX

	Eigen::ComplexEigenSolver<Eigen::MatrixXcd> CesX;
	Eigen::ComplexEigenSolver<Eigen::Matrix3cd> Ces3;

	Eigen::Matrix3cd U1, U2, U3, PMNS;
	Eigen::Matrix3cd R, M, V, md, UX;
	Eigen::Vector3cd MD, MR, mn;

	U1 <<	1, 		z, 		z,
	   	z, 		cos(t23),	sin(t23),
		z, 		-sin(t23), 	cos(t23);

	U2 <<	cons(t13),	z, 		sin(t23)*dcp,
	   	z, 		1,		z,
		-sin(t13)*dcp, 	z,		cos(t13);

	U3 <<	cos(t12), 	sin(t12),	z,
	   	-sin(t12), 	cos(t12),	z,
		z,		z,		1;

	PMNS = U1 * (U2 * U3);
	std::cout << "PMNS matrix is " << PMNS << std::endl;

	double mm1;
	double mm2 = mm1 + dmm21;
	double mm3 = mm1 + dmm31;

	double mme = std::abs(PMNS(1,1)) * mm1 +
		     std::abs(PMNS(1,2)) * mm2 +
		     std::abs(PMNS(1,2)) * mm2 ;

	if (sqrt(mme) < 2.05)
		ok;

	//Vector3
	MR <<	MR1, 	MR2,	MR3;
	mn <<	sqrt(mm1),	sqrt(mm2),	sqrt(mm3);

	M = M.asDiagonal() * UX.inverse() * M.asDiagonal();
	Ces.compute(M);
	MD = Ces.eigenvalues();
	V  = Ces.eigenvectors();

	md = V.adjoint() * MD.asDiagonal().sqrt() * R * mn.asDiagonal().sqrt() * PMNS.adjoint();

	//Mass Matrix is this one!
	Eigen::MatrixXcd M9(9, 9);
	M9 <<	Z3, 	md.transpose(), 	Z3,
		md,	Z3,			MR.asDiagonal(),
		Z3,	MR.asDiagonal(),	UX;

	while (nIter < 100)
	{

		std::complex<double> d11(Pol2Cart(lMod_d(Gen), Phase(Gen)));
		std::complex<double> d12(Pol2Cart(lMod_d(Gen), Phase(Gen)));
		std::complex<double> d21(Pol2Cart(lMod_d(Gen), Phase(Gen)));
		std::complex<double> d22(Pol2Cart(lMod_d(Gen), Phase(Gen)));
		std::complex<double> d31(Pol2Cart(lMod_d(Gen), Phase(Gen)));
		std::complex<double> d32(Pol2Cart(lMod_d(Gen), Phase(Gen)));

		std::complex<double> n11(Pol2Cart(lMod_n(Gen), Phase(Gen)));
		std::complex<double> n12(Pol2Cart(lMod_n(Gen), Phase(Gen)));
		std::complex<double> n13(Pol2Cart(lMod_n(Gen), Phase(Gen)));
		std::complex<double> n21(Pol2Cart(lMod_n(Gen), Phase(Gen)));
		std::complex<double> n22(Pol2Cart(lMod_n(Gen), Phase(Gen)));
		std::complex<double> n23(Pol2Cart(lMod_n(Gen), Phase(Gen)));

		std::complex<double> m11(Pol2Cart(lMod_m(Gen), Phase(Gen)));
		std::complex<double> m12(Pol2Cart(lMod_m(Gen), Phase(Gen)));
		std::complex<double> m22(Pol2Cart(lMod_m(Gen), Phase(Gen)));

		std::complex<double> u11(Pol2Cart(lMod_u(Gen), Phase(Gen)));
		std::complex<double> u12(Pol2Cart(lMod_u(Gen), Phase(Gen)));
		std::complex<double> u13(Pol2Cart(lMod_u(Gen), Phase(Gen)));
		std::complex<double> u22(Pol2Cart(lMod_u(Gen), Phase(Gen)));
		std::complex<double> u23(Pol2Cart(lMod_u(Gen), Phase(Gen)));
		std::complex<double> u33(Pol2Cart(lMod_u(Gen), Phase(Gen)));
	
		Eigen::MatrixXcd M(8,8);
		M <<	z,	z,	z,	d11,	d12,	z,	z,	z,
		  	z,	z,	z,	d21,	d22,	z,	z,	z,
			z,	z,	z,	d31,	d32,	z,	z,	z,
			d11,	d21,	d31,	m11,	m12,	n11,	n12,	n13,
			d12,	d22,	d32,	m12,	m22,	n21,	n22,	n23,
			z,	z,	z,	n11,	n21,	u11,	u12,	u13,
			z,	z,	z,	n12,	n22,	u12,	u22,	u23,
			z,	z,	z,	n13,	n23,	u13,	u23,	u33;
			
	
		//Diagonalise;
//		std::cout << M << std::endl;
	
		Eigen::MatrixXcd M2 = M.adjoint() * M;
//		std::cout << M2.real() << std::endl;
	
		Ces.compute(M2);
	
		//std::cout << Ces.eigenvalues() << std::endl;
		//std::cout << Ces.eigenvectors().real() << std::endl;
	
		for (unsigned int i = 0 ;  i < 8; ++i)
			std::cout << sqrt(std::abs(Ces.eigenvalues()[i])) << "\t";

		std::cout << std::endl;

		++nIter;
	}

	return 0;
}

std::complex<double> Pol2Cart(double lMod, double Phs)
{
	double Re = pow(10, lMod) * cos(Phs);
	double Im = pow(10, lMod) * sin(Phs);
	std::complex<double> w(Re, Im);

	return w;
}

//bool Oscillation();
