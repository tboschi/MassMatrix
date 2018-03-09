#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <complex>
#include <random>
#include <chrono>
#include <vector>
#include <algorithm>
#include <getopt.h>

#include "TTree.h"
#include "TFile.h"

#include "Tools.h"

int main(int argc, char** argv)
{
	std::ifstream InFile;
	std::ofstream OutFile;
	if (argc > 1)
		InFile.open(argv[1]);
	else
	{
		std::cout << "Give us an argument!";
		return 1;
	}
	if (argc > 2)
		OutFile.open(argv[2]);

	std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	std::string Line;
	std::stringstream ssL;

	unsigned int Grid = 100;
	unsigned int Hist_e[Grid][Grid] = {}, Hist_m[Grid][Grid] = {}, Hist_t[Grid][Grid] = {};

	double Mm, Ue, Um, Ut;
	double RM, Re, Rm, Rt;
	unsigned int Mm_i, Ue_i, Um_i, Ut_i;
	double lMmin = -4, lMmax = 0.7, lUmin = -16, lUmax = -6;
	double M_d = (lMmax - lMmin)/Grid, U_d = (lUmax - lUmin)/Grid;
	while (std::getline(InFile, Line))
	{
		ssL.str("");
		ssL.clear();

		ssL << Line;
		ssL >> Mm >> Ue >> Um >> Ut;
		Mm *= 1e-9;

		Mm_i = 2*Grid;
		Ue_i = 2*Grid;
		Um_i = 2*Grid;
		Ut_i = 2*Grid;

		//std::cout << "H0 " << M_d << "\t" << U_d << std::endl;
		for (unsigned int j = 0; j < Grid; ++j)
		{
			RM = (log10(Mm) - lMmin)/M_d - j;
			Re = (log10(Ue) - lUmin)/U_d - j;
			Rm = (log10(Um) - lUmin)/U_d - j;
			Rt = (log10(Ut) - lUmin)/U_d - j;
			//std::cout << "Ratio " << log10(Mm) << "\t" << RM << std::endl;
		       	//<< "\t" << Re << "\t" << Rm << "\t" << Rt << std::endl;

			if (RM > 0 && RM <= 1)
				Mm_i = j;

			if (Re > 0 && Re <= 1)
				Ue_i = j;
			if (Rm > 0 && Rm <= 1)
				Um_i = j;
			if (Rt > 0 && Rt <= 1)
				Ut_i = j;
		}

		//std::cout << "H1 " << Mm_i << "\t" << Ue_i << "\t" << Um_i << "\t" << Ut_i << std::endl;
		if (Mm_i < Grid && Ue_i < Grid && Um_i < Grid && Ut_i < Grid)
		{
			++Hist_e[Mm_i][Ue_i];
			++Hist_m[Mm_i][Um_i];
			++Hist_t[Mm_i][Ut_i];
		}
	}

	double Mass, Uu;
	unsigned int i = 0, j = 0;
	for (double logMass = lMmin; logMass < lMmax; logMass += M_d, ++i)
	{
		Mass = pow(10, logMass);
		unsigned int j = 0;
		for (double logUu = lUmin; logUu < lUmax; logUu += U_d, ++j)
		{
			Uu = pow(10, logUu);
			Out << Mass << "\t" << Uu << "\t"; 
			Out << Hist_e[i][j] << "\t"; 
			Out << Hist_m[i][j] << "\t"; 
			Out << Hist_t[i][j] << std::endl;
		}
	}

	return 0;
}
