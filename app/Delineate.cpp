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

	unsigned int Grid = 500;
	unsigned int Hist_e[Grid][Grid] = {}, Hist_m[Grid][Grid] = {}, Hist_t[Grid][Grid] = {};

	double Mm, Ue, Um, Ut;
	double RM, Re, Rm, Rt;
	unsigned int Mm_i, Ue_i, Um_i, Ut_i;
	double lMmin = -3, lMmax = 1, lUmin = -30, lUmax = 0;
	double M_d = (lMmax - lMmin)/Grid, U_d = (lUmax - lUmin)/Grid;
	while (std::getline(InFile, Line))
	{
		std::stringstream ssL(Line);

		double tmp;
		Mm = 0;
		Ue = 0;
		Um = 0;
		Ut = 0;

		bool Err = false;
		ssL >> Mm >> Ue >> Um >> Ut;
		/*
		std::cout << "In " << Line << std::endl;
		//dump all ssL
		for (unsigned int i = 0; !ssL.eof(); ++i)
		{
			double *pp;
			switch(i)
			{
				case 0:
					pp = &Mm;
					break;
				case 1:
					pp = &Ue;
					break;
				case 2:
					pp = &Um;
					break;
				case 3:
					pp = &Ut;
					break;
			}
			if (ssL >> tmp)
			{
				Err = false;
				*pp = tmp;
			}
			else if (!ssL.eof())
			{
				if (!Err)
					Err = true;
				else
					--i;

				*pp = sqrt(-1);
				ssL.clear();
				ssL.ignore(1);
			}
		}
		std::cout << "Out " << Mm << "\t" << Ue << "\t" << Um << "\t" << Ut << std::endl;
		//Mm *= 1e-9;
		//
		*/

		if (log10(Mm) < lMmin || log10(Mm) > lMmax)
			continue;

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
			//std::cout << "rest " << RM << "\t" << Re << "\t" << Rm << "\t" << Rt << std::endl;
			//std::cout << "Ratio " << log10(Mm) << "\t" << RM << std::endl;
		       	//<< "\t" << Re << "\t" << Rm << "\t" << Rt << std::endl;

			if (RM > -0.5 && RM < 0.5)
				Mm_i = j;

			if (Re > -0.5 && Re < 0.5)
				Ue_i = j;
			if (Rm > -0.5 && Rm < 0.5)
				Um_i = j;
			if (Rt > -0.5 && Rt < 0.5)
				Ut_i = j;
		}
		//std::cout << Mm << "\t" << Ue << "\t" << Um << "\t" << Ut << std::endl;
		//std::cout << "rest " << RM << "\t" << Re << "\t" << Rm << "\t" << Rt << std::endl;

		//std::cout << "H1 " << Mm_i << "\t" << Ue_i << "\t" << Um_i << "\t" << Ut_i << std::endl;
		if (Mm_i < Grid)
		{
			if (Ue_i < Grid)
				++Hist_e[Mm_i][Ue_i];
			if (Um_i < Grid)
				++Hist_m[Mm_i][Um_i];
			if (Ut_i < Grid)
				++Hist_t[Mm_i][Ut_i];
		}
	}

	double Mass, Uu;
	unsigned int i = 0, j = 0;
	for (double logMass = lMmin; logMass < lMmax && i < Grid; logMass += M_d, ++i)
	{
		Mass = pow(10, logMass);
		unsigned int j = 0;
		for (double logUu = lUmin; logUu < lUmax && j < Grid; logUu += U_d, ++j)
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
