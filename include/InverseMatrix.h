///////////////////////////////////////////////////////////////////////////////
// a Inverse Matrix is in the form of 
//
// /	0	Mr	0	\		/	0	0	Mr	\
// |	MrT	Ur	Ms	|	~	|	0	Us	MsT	| 
// \	0	MsT	Us	/		\	MrT	Ms	Ur	/
//
// after construction each elemenet must be populated
///////////////////////////////////////////////////////////////////////////////

#ifndef INVERSEMATRIX_H
#define INVERSEMATRIX_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <random>
#include <cmath>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/SVD>

#include "Tools.h"

enum class Block
{
	Full,
	Mr,
	Ms,
	Ur,
	Us,
};

class InverseMatrix
{
	public:
		InverseMatrix(unsigned int inR, unsigned int inS);
		Eigen::MatrixXcd MassMatrix();
		void Show(Block BN, bool Abs = false);
		void Print(const Eigen::MatrixXcd &S);
		void Set(Block BN, int min, int max);
		void Clean(Block BN);

		std::complex<double> RandomAs(Block BN);
		int RandomEntry(Block BN, unsigned int i, unsigned int j);
		int Manual(Block BN, std::complex<double> z, unsigned int i, unsigned int j);
		std::vector<int> Populate(Block BN, bool Fixed = false);
		std::vector<int> PopPerRow(Block BN);
		std::vector<int> PopPerCol(Block BN);
		std::vector<int> PerRowAll(Eigen::Ref<Eigen::MatrixXcd> mBlock, std::uniform_int_distribution<int> &Pow);
		std::vector<int> PerColAll(Eigen::Ref<Eigen::MatrixXcd> mBlock, std::uniform_int_distribution<int> &Pow);
		std::vector<int> PopulateAll(Eigen::Ref<Eigen::MatrixXcd> mBlock, std::uniform_int_distribution<int> &Pow, bool Fixed = false);
		std::vector<int> PopulateSup(Eigen::Ref<Eigen::MatrixXcd> mBlock, std::uniform_int_distribution<int> &Pow, bool Fixed = false);
		std::vector<int> PopulateReduceMr(std::uniform_int_distribution<int> &Pow, bool Fixed = false);
		std::vector<int> PopulateReduceMs(std::uniform_int_distribution<int> &Pow, bool Fixed = false);
		std::vector<int> PopulateReduceUs(std::uniform_int_distribution<int> &Pow, bool Fixed = false);

		Eigen::MatrixXcd Get(Block BN);
		Eigen::MatrixXcd ColumnBlock();	
		Eigen::MatrixXcd ColumnBlockSVD(std::vector<double> &SingularValues);
		Eigen::MatrixXcd MassMatrixSVD(std::vector<double> &SingularValues);

		unsigned int nR();
		unsigned int nS();
		unsigned int nM();
		unsigned int n0();

		bool FindDeltaM2(std::vector<double> &vM, bool &Hierarchy);
		bool FindMass(std::vector<double> &vM, double Min, double Max);
		bool BB0(std::vector<double> &vM, Eigen::MatrixXcd &VA, double &Mbb);
		bool MEG(std::vector<double> &vM, Eigen::MatrixXcd &VA, double &MEGbranch);
		bool NSI(std::vector<double> &vM, Eigen::MatrixXcd &VA, std::vector<double> &vOsc, std::vector<double> &vEws);

	private:
		unsigned int inR, inS;
		bool ReduceParameters;
	
		std::mt19937 MT;

		std::uniform_int_distribution<int> Pow_Mr,
						   Pow_Ms,
						   Pow_Ur,
						   Pow_Us;

		std::uniform_real_distribution<double> Valor, Phase;

		Eigen::MatrixXcd Block_Mr, 
				 Block_Ms, 
				 Block_Ur, 
				 Block_Us;
};

#endif
