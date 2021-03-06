#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <vector>
#include <cmath>
#include <getopt.h>
#include <cstdlib>

void Usage(char* argv0);
int main(int argc, char** argv)
{
	const struct option longopts[] = 
	{
		{"threshold", 	required_argument,	0, 't'},
		{"massdep", 	required_argument,	0, 'm'},
		{"start", 	required_argument,	0, 'A'},
		{"end", 	required_argument,	0, 'B'},
		{"input", 	no_argument,		0, 'i'},
		{"background", 	required_argument,	0, 'W'},
		{"output", 	required_argument,	0, 'o'},
		{"help", 	no_argument,	 	0, 'h'},
		{0,	0, 	0,	0},
	};

	int index; 
	int iarg = 0;
	opterr = 1;
	
	//Initialize variables
	std::ofstream OutFile;
	std::ifstream InFile;
	double Threshold = 2.44;
	bool UeFlag = false, UmFlag = false, UtFlag = false;
	
	while((iarg = getopt_long(argc,argv, "i:o:t:m:EMTh", longopts, &index)) != -1)
	{
		switch(iarg)
		{
			case 'i':
				InFile.open(optarg);
				break;
			case 'o':
				OutFile.open(optarg);
				break;
			case 't':
				Threshold = strtod(optarg, NULL);
				break;
			case 'E':
				UeFlag = true;
				break;
			case 'M':
				UmFlag = true;
				break;
			case 'T':
				UtFlag = true;
				break;
			case 'h':
				Usage(argv[0]);
				return 1;
			default:
				return 1;
				break;
		}
	}

	std::ostream &Out = (OutFile.is_open()) ? OutFile : std::cout;

	bool First = true;
	bool Reached = false;
	unsigned int rE, rM, rT;
	double rMass, rUu2, rEvt;
	double fMass, fUu2, fEvt;
	std::vector <double> vMass, vUu2, vEvt;

	std::string Line;
	std::stringstream ssL;
	while (std::getline(InFile, Line))
	{
		if (Line[0] == '#') continue;

		ssL.str("");
		ssL.clear();
		ssL << Line;

		ssL >> rMass >> rUu2 >> rE >> rM >> rT;
		vMass.push_back(rMass);
		vUu2.push_back(rUu2);
		if (UeFlag)
			vEvt.push_back(rE);
		else if (UmFlag)
			vEvt.push_back(rM);
		else if (UtFlag)
			vEvt.push_back(rT);
	}

	double Thr;
	unsigned int i;
	double Mass_ = -1.0;
	for (unsigned int j = 0; j < 2*vMass.size(); ++j)
	{
		if (j < vMass.size())
			i = j;
		else
			i = 2*vMass.size() - 1 - j;

		if (vMass.at(i) != Mass_ || j == vMass.size())
		//if (vMass.at(i) != Mass_)
		{
			Reached = false;
			Mass_ = vMass.at(i);
		}

		Thr = Threshold;

		if (!Reached && vEvt.at(i) > Thr && (i < vEvt.size() ? vEvt.at(i+1) > Thr : true))
		//if (!Reached && vEvt.at(i) > Thr)
		{
			if (First)
			{
				fMass = vMass.at(i);
				fUu2 = vUu2.at(i);
				fEvt = vEvt.at(i);
				First = false;
			}

			Out << vMass.at(i) << "\t" << vUu2.at(i) << "\t" << vEvt.at(i) << std::endl;
			Reached = true;
		}
	}

	Out << fMass << "\t" << fUu2 << "\t" << fEvt << std::endl;

	InFile.close();
	OutFile.close();

	return 0;
}

void Usage(char* argv0)
{
	std::cout << "Description" << std::endl;
	std::cout << "Usage : " << std::endl;
	std::cout << argv0 << " [OPTIONS]" << std::endl;
	std::cout <<"\n  -i,  --input" << std::endl;
	std::cout << "\t\tInput file, tabulated as Mass\tUu\tnEvt" << std::endl;
	std::cout <<"\n  -o,  --output" << std::endl;
	std::cout << "\t\tOutput file" << std::endl;
	std::cout <<"\n  -t,  --threshold" << std::endl;
	std::cout << "\t\tEvent threshold for signal. As mass dependent threshold, t is the y-intercept" << std::endl;
	std::cout <<"\n  -m,  --massdep" << std::endl;
	std::cout << "\t\tIf specified, a mass dependance is considered for threshold, [thr] = t + m * [Mass]" << std::endl;
	std::cout <<"\n  -h,  --help" << std::endl;
	std::cout << "\t\tPrint this message and exit" << std::endl;
}
