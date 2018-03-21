#include "InverseMatrix.h"

InverseMatrix::InverseMatrix(unsigned int nR_, unsigned int nS_) :
	inR(nR_),
	inS(nS_),
	MT(std::random_device()()),
	Valor(std::uniform_real_distribution<double>(0, 1)),
	Phase(std::uniform_real_distribution<double>(0, 2*Const::fPi))
{
	Block_Mr.resize(  3, inR);
	Block_Ms.resize(inR, inS);
	Block_Ur.resize(inR, inR);
	Block_Us.resize(inS, inS);

	if ( nS()-nR() == pow(nS()-nR(), 2) )
		ReduceParameters = true;
	else
		ReduceParameters = false;

	//if ReduceParameters is true then some dof can be removed easily
	//it will speed up computation
	//Mr can have 3 entries which are real
	//Ms is diagonal and real -> great for fixing LO mass states
	//Ur does not change, but it is useless anyway
	//Us has a real diagonal
}

Eigen::MatrixXcd InverseMatrix::MassMatrix()
{
	Eigen::MatrixXcd MM(nM(), nM());
	MM << Eigen::MatrixXcd::Zero(3,3),	Block_Mr,		Eigen::MatrixXcd::Zero(3,nS()),
	      Block_Mr.transpose(),		Block_Ur,		Block_Ms,
	      Eigen::MatrixXcd::Zero(nS(),3),	Block_Ms.transpose(),	Block_Us;

	return MM;
}

//Print blocks
void InverseMatrix::Show(Block BN, bool Abs)
{
	switch (BN)
	{
		case Block::Full:
			if (Abs)
				std::cout << std::endl << MassMatrix().cwiseAbs() << std::endl << std::endl; 
			else
				std::cout << std::endl << MassMatrix() << std::endl << std::endl; 
			break;
		case Block::Mr:
			if (Abs)
				std::cout << std::endl << Block_Mr.cwiseAbs() << std::endl << std::endl; 
			else
				std::cout << std::endl << Block_Mr << std::endl << std::endl; 
			break;
		case Block::Ms:
			if (Abs)
				std::cout << std::endl << Block_Ms.cwiseAbs() << std::endl << std::endl; 
			else
				std::cout << std::endl << Block_Ms << std::endl << std::endl; 
			break;
		case Block::Ur:
			if (Abs)
				std::cout << std::endl << Block_Ur.cwiseAbs() << std::endl << std::endl; 
			else
				std::cout << std::endl << Block_Ur << std::endl << std::endl; 
			break;
		case Block::Us:
			if (Abs)
				std::cout << std::endl << Block_Us.cwiseAbs() << std::endl << std::endl; 
			else
				std::cout << std::endl << Block_Us << std::endl << std::endl; 
			break;
		default:
			std::cerr << "Wrong block" << std::endl;
			break;
	}
}

void InverseMatrix::Set(Block BN, int min, int max)
{
	std::uniform_int_distribution<int> Pow(min, max);

	switch (BN)
	{
		case Block::Full:
			Pow_Mr = Pow; 
			Pow_Ms = Pow; 
			Pow_Us = Pow; 
			Pow_Ur = Pow; 
			break;
		case Block::Mr:
			Pow_Mr = Pow; 
			break;
		case Block::Ms:
			Pow_Ms = Pow; 
			break;
		case Block::Ur:
			Pow_Ur = Pow; 
			break;
		case Block::Us:
			Pow_Us = Pow; 
			break;
		default:
			std::cerr << "Wrong block" << std::endl;
			break;
	}
}

//Set all block entries to zero
void InverseMatrix::Clean(Block BN)
{
	switch (BN)
	{
		case Block::Full:
			Block_Mr.setZero();
			Block_Ms.setZero();
			Block_Ur.setZero();
			Block_Us.setZero();
			break;
		case Block::Mr:
			Block_Mr.setZero();
			break;
		case Block::Ms:
			Block_Ms.setZero();
			break;
		case Block::Ur:
			Block_Ur.setZero();
			break;
		case Block::Us:
			Block_Us.setZero();
			break;
		default:
			std::cerr << "Wrong block" << std::endl;
			break;
	}
}

//MM values to entries of block. Fixed decides if the overall magnitued is the same for each block and return that value
std::vector<int> InverseMatrix::Populate(Block BN, bool Fixed)
{
	std::vector<int> vRet;
	std::vector<int> vMr, vMs, vUr, vUs;
	switch (BN)
	{
		case Block::Full:
			vMr = Populate(Block::Mr, Fixed);
			vMs = Populate(Block::Ms, Fixed);
			vUr = Populate(Block::Ur, Fixed);
			vUs = Populate(Block::Us, Fixed);
			vRet.insert(vRet.end(), vMr.begin(), vMr.end());
			vRet.insert(vRet.end(), vMs.begin(), vMs.end());
			vRet.insert(vRet.end(), vUr.begin(), vUr.end());
			vRet.insert(vRet.end(), vUs.begin(), vUs.end());
			break;
		case Block::Mr:
			if (ReduceParameters)
				vRet = PopulateReduceMr(Pow_Mr, Fixed);
			else 
				vRet = PopulateAll(Block_Mr, Pow_Mr, Fixed);
			break;
		case Block::Ms:
			if (ReduceParameters)
				vRet = PopulateReduceMs(Pow_Ms, Fixed);
			else
				vRet = PopulateAll(Block_Ms, Pow_Ms, Fixed);
			break;
		case Block::Ur:
			vRet = PopulateAll(Block_Ur, Pow_Ur, Fixed);
			break;
		case Block::Us:
			if (ReduceParameters)
				vRet = PopulateReduceUs(Pow_Us, Fixed);
			else
				vRet = PopulateAll(Block_Us, Pow_Us, Fixed);
			break;
		default:
			std::cerr << "Wrong block" << std::endl;
			break;
	}

	return vRet;
}

//Assign values to entries of block, varying magnitude in each row. No full option
std::vector<int> InverseMatrix::PopPerRow(Block BN)
{
	switch (BN)
	{
		case Block::Mr:
			return PerRowAll(Block_Mr, Pow_Mr);
		case Block::Ms:
			return PerRowAll(Block_Ms, Pow_Ms);
		default:
			std::cerr << "Wrong block" << std::endl;
			//return NULL;
	}
}

//Assign values to entries of block, varying magnitude in each column. No full option
std::vector<int> InverseMatrix::PopPerCol(Block BN)
{
	switch (BN)
	{
		case Block::Mr:
			return PerColAll(Block_Mr, Pow_Mr);
		case Block::Ms:
			return PerColAll(Block_Ms, Pow_Ms);
		default:
			std::cerr << "Wrong block" << std::endl;
			//return 0;
	}
}

//Assign randomly block entires, one by one
std::complex<double> InverseMatrix::RandomAs(Block BN)
{
	switch (BN)
	{
		case Block::Mr:
			return Const::Pol2Cart(Pow_Mr(MT)+Valor(MT), Phase(MT));
		case Block::Ms:
			return Const::Pol2Cart(Pow_Ms(MT)+Valor(MT), Phase(MT));
		case Block::Ur:
			return Const::Pol2Cart(Pow_Ur(MT)+Valor(MT), Phase(MT));
		case Block::Us:
			return Const::Pol2Cart(Pow_Us(MT)+Valor(MT), Phase(MT));
		default:
			std::cerr << "Wrong block" << std::endl;
			break;
	}
}

//Assign randomly block entires, one by one
int InverseMatrix::RandomEntry(Block BN, unsigned int i, unsigned int j)
{
	return Manual(BN, RandomAs(BN), i, j);
}

//Assign manually block entires, one by one
int InverseMatrix::Manual(Block BN, std::complex<double> z, unsigned int i, unsigned int j)
{
	switch (BN)
	{
		case Block::Mr:
			Block_Mr(i, j) = z;
			break;
		case Block::Ms:
			Block_Ms(i, j) = z;
			break;
		case Block::Ur:
			Block_Ur(i, j) = z;
			break;
		case Block::Us:
			Block_Us(i, j) = z;
			break;
		default:
			std::cerr << "Wrong block" << std::endl;
			break;
	}
	return int( floor( log10( std::abs(z) ) ) );	//magnitude of |z|
}

//Populates all entries of a block, keeping same magnitude each row
std::vector<int> InverseMatrix::PerRowAll(Eigen::Ref<Eigen::MatrixXcd> mBlock, std::uniform_int_distribution<int> &Pow)
{
	std::vector<int> RowMod;
	for (unsigned int i = 0; i < mBlock.rows(); ++i)
	{
		int GlobMod = Pow(MT);
		RowMod.push_back(GlobMod);
		for (unsigned int j = 0; j < mBlock.cols(); ++j)
			mBlock(i,j) = Const::Pol2Cart(GlobMod+Valor(MT), Phase(MT));
	}
	return RowMod;
}

//Populates all entries of a block, keeping same magnitude each col
std::vector<int> InverseMatrix::PerColAll(Eigen::Ref<Eigen::MatrixXcd> mBlock, std::uniform_int_distribution<int> &Pow)
{
	std::vector<int> ColMod;
	for (unsigned int j = 0; j < mBlock.cols(); ++j)
	{
		int GlobMod = Pow(MT);
		ColMod.push_back(GlobMod);
		for (unsigned int i = 0; i < mBlock.rows(); ++i)
			mBlock(i,j) = Const::Pol2Cart(GlobMod+Valor(MT), Phase(MT));
	}
	return ColMod;
}

//Populates all entries of a block, used if there is no symmetry
std::vector<int> InverseMatrix::PopulateAll(Eigen::Ref<Eigen::MatrixXcd> mBlock, std::uniform_int_distribution<int> &Pow, bool Fixed)
{
	std::vector<int> vRet;
	int GlobMod;
	if (Fixed)
		GlobMod = Pow(MT);

	for (unsigned int j = 0; j < mBlock.cols(); ++j)
	{
		for (unsigned int i = 0; i < mBlock.rows(); ++i)
		{
			if (!Fixed)
				GlobMod = Pow(MT);
			vRet.push_back(GlobMod);

			mBlock(i,j) = Const::Pol2Cart(GlobMod + Valor(MT), Phase(MT));
		}
	}

	return vRet;
}

//Populates superior triangle of block and the copy for inferior, used for symmetric matrices
std::vector<int> InverseMatrix::PopulateSup(Eigen::Ref<Eigen::MatrixXcd> mBlock, std::uniform_int_distribution<int> &Pow, bool Fixed)
{
	std::vector<int> vRet;
	int GlobMod;
	if (Fixed)
		GlobMod = Pow(MT);

	for (unsigned int j = 0; j < mBlock.cols(); ++j)
	{
		for (unsigned int i = 0; i < j+1; ++i)
		{
			if (!Fixed)
				GlobMod = Pow(MT);
			vRet.push_back(GlobMod);

			mBlock(i, j) = Const::Pol2Cart(GlobMod + Valor(MT), Phase(MT));
			if (j > i)
				mBlock(j, i) = mBlock(i, j);
		}
	}
	return vRet;
}

//Populates all entries of blocks Mr, Ms, Us if in ReducedParameters mode, will speed uyp things
//Mr : the first column is set to be real (j == 0 has no phase)
std::vector<int> InverseMatrix::PopulateReduceMr(std::uniform_int_distribution<int> &Pow, bool Fixed)
{
	std::vector<int> vRet;
	int GlobMod;
	if (Fixed)
		GlobMod = Pow(MT);

	for (unsigned int j = 0; j < Block_Mr.cols(); ++j)
	{
		for (unsigned int i = 0; i < Block_Mr.rows(); ++i)
		{
			if (!Fixed)
				GlobMod = Pow(MT);
			vRet.push_back(GlobMod);

			Block_Mr(i,j) = Const::Pol2Cart(GlobMod + Valor(MT), j == 0 ? 0 : Phase(MT));
		}
	}

	return vRet;
}

//Ms : is set to be diagonal and real
std::vector<int> InverseMatrix::PopulateReduceMs(std::uniform_int_distribution<int> &Pow, bool Fixed)
{
	std::vector<int> vRet;
	int GlobMod;
	if (Fixed)
		GlobMod = Pow(MT);

	for (unsigned int i = 0; i < Block_Ms.rows(); ++i)
	{
		if (!Fixed)
			GlobMod = Pow(MT);
		vRet.push_back(GlobMod);

		Block_Ms(i,i) = Const::Pol2Cart(GlobMod + Valor(MT), 0);
	}

	return vRet;
}

//Us : is set to have a real diagonal
std::vector<int> InverseMatrix::PopulateReduceUs(std::uniform_int_distribution<int> &Pow, bool Fixed)
{
	std::vector<int> vRet;
	int GlobMod;
	if (Fixed)
		GlobMod = Pow(MT);

	for (unsigned int j = 0; j < Block_Us.cols(); ++j)
	{
		for (unsigned int i = 0; i < j+1; ++i)
		{
			if (!Fixed)
				GlobMod = Pow(MT);
			vRet.push_back(GlobMod);

			Block_Us(i, j) = Const::Pol2Cart(GlobMod + Valor(MT), i == j ? 0 : Phase(MT));
			if (j > i)
				Block_Us(j, i) = Block_Us(i, j);
		}
	}

	return vRet;
}

//Return a block
Eigen::MatrixXcd InverseMatrix::Get(Block BN)
{
	switch (BN)
	{
		case Block::Mr:
			return Block_Mr;
		case Block::Ms:
			return Block_Ms;
		case Block::Ur:
			return Block_Ur;
		case Block::Us:
			return Block_Us;
		default:
			std::cerr << "Wrong block" << std::endl;
			//return 0;
	}
}

//Return the column 	MCB^T = (Mr^T, Ms), used in SVD
Eigen::MatrixXcd InverseMatrix::ColumnBlock()
{
	Eigen::MatrixXcd MCB(3+nS(), nR());
	MCB << Block_Mr, Block_Ms.transpose();
	return MCB;
}

//Compute singular values of ColumnBlock
Eigen::MatrixXcd InverseMatrix::ColumnBlockSVD(std::vector<double> &SingularValues)
{
	SingularValues.clear();

	Eigen::JacobiSVD<Eigen::MatrixXcd> ColumnSVD(ColumnBlock(), Eigen::ComputeFullV | Eigen::ComputeFullU);
	//Eigen::JacobiSVD<Eigen::MatrixXcd> ColumnSVD(ColumnBlock());

	for (unsigned int i = 0; i < ColumnSVD.singularValues().size(); ++i)
		SingularValues.push_back(ColumnSVD.singularValues()[i]);

	unsigned int nU = ColumnSVD.matrixU().cols();
	unsigned int nV = ColumnSVD.matrixV().cols();
	Eigen::MatrixXcd UU(nU+nV, nU+nV);
	UU << ColumnSVD.matrixU(),	Eigen::MatrixXcd::Zero(nU, nV),
	      Eigen::MatrixXcd::Zero(nV, nU), 	ColumnSVD.matrixV().conjugate();

	return UU;
}

//SVD decomposition of the full mass matrix
Eigen::MatrixXcd InverseMatrix::MassMatrixSVD(std::vector<double> &SingularValues)
{
	SingularValues.clear();

	Eigen::JacobiSVD<Eigen::MatrixXcd> MassSVD(MassMatrix(), Eigen::ComputeFullV);

	for (unsigned int i = 0; i < MassSVD.singularValues().size(); ++i)
		SingularValues.push_back(MassSVD.singularValues()[i]);
	std::reverse(SingularValues.begin(), SingularValues.end());

	Eigen::MatrixXcd Diag = MassSVD.matrixV().transpose() * MassMatrix() * MassSVD.matrixV();
	Eigen::MatrixXcd Majorana  = Eigen::MatrixXcd::Zero(nM(), nM());
	Eigen::MatrixXcd Projector = Eigen::MatrixXcd::Zero(nM(), nM());
	for (unsigned int i = 0; i < Diag.rows(); ++i)
	{
		Majorana( i, i) = std::exp(- Const::i * std::arg(Diag(i, i)) / 2.0);
		Projector(i, nM()-1-i) = std::complex<double>(1,0);
	}

	return MassSVD.matrixV() * Projector * Majorana;
}

//number of right-handed neutrinos
unsigned int InverseMatrix::nR()
{
	return inR;
}

//number of sterile fermions
unsigned int InverseMatrix::nS()
{
	return inS;
}

unsigned int InverseMatrix::nM()
{
	return 3+inR+inS;
}

//number of non-zero SVD for the leading order matrix
unsigned int InverseMatrix::n0()
{
	return nR() < 3 ? nR() : std::min(nR(), nS());
}

//return true if the Dm2 from nufit oscillation data is found
//Hierarchy is true if is normal ordering, false if inverted
bool InverseMatrix::FindDeltaM2(std::vector<double> &vMass, bool &Hierarchy)
{
	bool N21 = (vMass.at(1)-vMass.at(0) > 6.8e-5	&& vMass.at(1)-vMass.at(0) < 8.02e-5);
	bool N31 = (vMass.at(2)-vMass.at(0) > 2.399e-3	&& vMass.at(2)-vMass.at(0) < 2.593e-3);

	bool I21 = (vMass.at(2)-vMass.at(1) > 6.8e-5	&& vMass.at(2)-vMass.at(1) < 8.02e-5);
	bool I32 = (vMass.at(2)-vMass.at(0) > 2.369e-3	&& vMass.at(2)-vMass.at(0) < 2.562e-3);

	if (N21 && N31)
		Hierarchy = true;
	if (I21 && I32)
		Hierarchy = false;

	return (N21 && N31) || (I21 && I32);
}

//returns true if there is a mass state in the range Min-Max
bool InverseMatrix::FindMass(std::vector<double> &vMass, double Min, double Max)
{
	bool M4 = false; 
	for (auto p : vMass)
		if (!M4 && p > Min && p < Max)
			M4 = true;

	return M4;
}

//return true if satisfies GERDA
bool InverseMatrix::BB0(std::vector<double> &vMass, Eigen::MatrixXcd &VA)
{
	double p2 = -pow(125e6, 2);
	std::complex<double> BBeff;
	for (unsigned int i = 0; i < nM(); ++i)
		BBeff += VA(0, i) * VA(0, i) * p2 * vMass.at(i) / (p2 - vMass.at(i));

	//bool BB0 = std::abs(BBeff) < 150e-3;	//present
	bool BB0 = std::abs(BBeff) < 20e-3;	//future

	return BB0;
}

//return true if satisfies MEG
bool InverseMatrix::MEG(std::vector<double> &vMass, Eigen::MatrixXcd &VA)
{
	std::complex<double>  MEGamp;
	for (unsigned int i = 0; i < nM(); ++i)
		MEGamp += std::conj(VA(0, i)) * VA(1, i) * Const::LoopG( vMass.at(i)*vMass.at(i) / pow(Const::fMW, 2) );

	double MEGbranch = 3 * Const::fAem * std::norm(MEGamp) / (32 * Const::fPi);

	//bool MEG = MEGbranch < 4.2e-13;		//present
	bool MEG = MEGbranch < 5e-14;		//future

	return MEG;
}

//return true if satisfies unitarity by NSI constraints
bool InverseMatrix::NSI(std::vector<double> &vMass, Eigen::MatrixXcd &VA)
{
	Eigen::Matrix3d Unit;
	Eigen::Matrix3cd Kab;

	Unit <<	4.0e-3,	1.2e-4,	3.2e-3,
		1.2e-4,	1.6e-3,	2.1e-3,
		3.2e-3,	2.1e-3,	5.3e-3;

	for (unsigned int i = 4; i < nM(); ++i)
		for (unsigned int c = 0; c < Kab.cols(); ++c)
			for (unsigned int r = 0; r < c+1; ++r)
				Kab(r, c) += VA(r, i) * std::conj(VA(c, i));

	bool NSI = true;
	for (unsigned int c = 0; c < Kab.cols(); ++c)
		for (unsigned int r = 0; r < c+1; ++r)
			NSI *= Kab.cwiseAbs()(r, c) < Unit(r, c);

	return NSI;
}
