#include <iostream>
#include <ctime>
#include <cstdlib>
#include "Param.hpp"
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

Param::Param(int argc, char **argv):
	modfile(""), snvfile(""), hapfile(""), indfile(""), seed(0), chr("A"), recEvent("1"), mutEvent("1"), selEvent("1"), givenRecRate(false), givenMutRate(false), givenPop(true), givenGen(true)
{
	if (argc > 1 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help"))
	{
		help();
		exit(0);
	}
	
	if (argc < 3)
	{
		std::cerr << "Error: Need more arguments than provided, use -h/--help to get more help!" << std::endl;
		exit(1);
	}
	bool givenSeed = false;
	bool givenOut = false;
	
	for (int i = 0; i < argc; ++i)
	{
		std::string arg(argv[i]);
		if (arg == "-h" || arg == "--help")
		{
			help();
			exit(0);
		}
		else if (arg == "-in" || arg == "--inPrefix")
		{
			std::string inPrefix = std::string(argv[++i]);
			modfile = inPrefix+".mod";
			snvfile = inPrefix+".snv";
			hapfile = inPrefix+".hap";
			indfile = inPrefix+".ind";
		}
		else if (arg == "-mod" || arg == "--modfile")
		{
			modfile = std::string(argv[++i]);
		}
		else if (arg == "-snv" || arg == "--snvfile")
		{
			snvfile = std::string(argv[++i]);
		}
		else if (arg == "-hap" || arg == "--hapfile")
		{
			hapfile = std::string(argv[++i]);
		}
		else if (arg == "-ind" || arg == "--indfile")
		{
			indfile = std::string(argv[++i]);
		}
		else if (arg == "-p" || arg == "--poplist")
		{
			std::string temp = std::string(argv[++i]);
			boost::split(popList, temp, boost::is_any_of(","), boost::token_compress_on);
		}
		else if (arg == "-g" || arg == "--gen")
		{
			std::string temp = std::string(argv[++i]);
			std::vector<std::string> temp1;
			boost::split(temp1, temp, boost::is_any_of(","), boost::token_compress_on);
			popGen.resize(temp1.size());
			for (int j = 0; j < temp1.size(); j++)
			{
				popGen[j] = std::atoi(temp1.at(j).c_str());
			}
		}
		else if (arg == "-n" || arg == "--nInd")
		{
			std::string temp = std::string(argv[++i]);
			std::vector<std::string> temp1;
			boost::split(temp1, temp, boost::is_any_of(","), boost::token_compress_on);
			popNumber.resize(temp1.size());
			for (int j = 0; j < temp1.size(); j++)
			{
				popNumber[j] = std::atoi(temp1.at(j).c_str());
			}
		}
		else if (arg == "-s" || arg == "--seed")
		{
			seed = unsigned(atol(argv[++i]));
			givenSeed = true;
		}
		else if (arg == "-c" || arg == "--chr")
		{
			chr = std::string(argv[++i]);
		}
		else if (arg == "-rec" || arg == "--recRate")
		{
			givenRecRate = true;
			recRate = atof(argv[++i]);
		}
		else if (arg == "-mut" || arg == "--mutRate")
		{
			givenMutRate = true;
			mutRate = atof(argv[++i]);
		}
		else if (arg == "-out" || arg == "--outPrefix")
		{
			outPrefix = std::string(argv[++i]);
			givenOut = true;
		}
		else if (arg == "--no-rec")
		{
			recEvent = "0";	
		}	
		else if (arg == "--no-mut")
		{
			mutEvent = "0";
		}
		else if (arg == "--no-sel")
		{
			selEvent = "0";
		}
	}
	
	if (modfile.size() == 0)
	{
		std::cerr << "Error: Model description file must be specified!" << std::endl;
		exit(1);
	}
	if (snvfile.size() == 0)
	{
		std::cerr << "Error: SNV file must be specified!" << std::endl;
		exit(1);
	}
	if (hapfile.size() == 0)
	{
		std::cerr << "Error: Anchap file must be specified!" << std::endl;
		exit(1);
	}
	if (indfile.size() == 0)
	{
		std::cerr << "Error: Individual file must be specified!" << std::endl;
		exit(1);
	}
	if (popList.size() == 0)
	{
		givenPop = false;
		popList.resize(1);
		std::cout << "Warning: Output population is not specified and final generated population is output by default!" << std::endl;
	}
	if (popGen.size() == 0)
	{
		givenGen = false;
		popGen.resize(1);
		std::cout << "Warning: Output generation index is not specified and last generation is output by default!" << std::endl;
	}		
	if (popNumber.size() == 0)
	{
		popNumber.resize(1);
		popNumber[0] = 10;
		std::cout << "Warning: Output population size is not specified and population size is set as 10 by default!" << std::endl;
	}
	if (popList.size() != popNumber.size() || popList.size() != popGen.size())
	{
		std::cerr << "Error: The sampled number or generation index setting of output population is incorrect!" << std::endl;
		exit(1);
	}
	if (givenPop == false)
	{
		popList.clear();
	}
	if (givenGen == false)
	{
		popGen.clear();
	}
	for (int i = 0; i < popGen.size(); i++)
	{
		if (popGen.at(i) < 0)
		{
			std::cerr << "Error: The generation index \"" << popGen.at(i) << "\" is smaller than zero!" << std::endl;
			exit(1);
		}
	}
	for (int i = 0; i < popNumber.size(); i++)
	{
		if (popNumber.at(i) <= 0)
		{
			std::cerr << "Error: The number of individuals sampled \"" << popNumber.at(i) << "\" should be larger than zero!" << std::endl;
			exit(1);
		}
	}
	for (int i = 0; i < popList.size(); i++)
	{	
		if (count(popList.begin(), popList.end(), popList.at(i)) > 1)
		{
			std::vector<int> genIndex;
			for (int j = 0; j < popList.size(); j++)
			{
				if (popList.at(j) == popList.at(i))
				{
					genIndex.push_back(popGen.at(j));
				}
			}
	
			for (int j = 0; j < genIndex.size(); j++)
			{
				if (count(genIndex.begin(), genIndex.end(), genIndex.at(j)) > 1)
				{
					std::cerr << "Error: For output, population \"" << popList.at(i) << "\" and generation index \"" << genIndex.at(j) << "\" repeats!" << std::endl;
					exit(1);
				}
			}
		}
	}
	if (!givenSeed)
	{
		seed = unsigned(time(0));
	}
	srand(seed);
	if (!givenOut)
	{
		outPrefix = "out";
	}
	if (chr != "A" && chr != "X")
	{
		std::cerr << "Error: The chromosome setting can only be \"A\" or \"X\"!" << std::endl;
		exit(1);
	}
}

Param::~Param()
{
}

void Param::help()
{
	std::cout << "-------------------------------------------------------------------------" << std::endl;
	std::cout << kProgram << " " << kVersion << std::endl;
	std::cout << "Arguments: " << std::endl;
	std::cout << "\t-in/--inPrefix\t<string>\tThe prefix of four input files [required]" << std::endl;
	std::cout << "\t-mod/--modfile\t<string>\tModel description file [required]" << std::endl;	
	std::cout << "\t-ind/--indfile\t<string>\tIndividual information file [required]" << std::endl;
	std::cout << "\t-hap/--hapfile\t<string>\tAncestral haplotype data file [required]" << std::endl;
	std::cout << "\t-snv/--snvfile\t<string>\tSNV information file [required]" << std::endl;
	std::cout << std::endl;
	std::cout << "\t-p/--poplist\t<string>\tOutput population list, separated by comma [optional, default = the final generated population]" << std::endl;
	std::cout << "\t-g/--gen\t[integer]\tThe generation index of each output population, separated by comma [optional, default = the last generation in simulation]" << std::endl;
	std::cout << "\t-n/--nInd\t[integer]\tThe sampled number of each output population, separated by comma [optional, default = 10]" << std::endl;
	std::cout << std::endl;
	std::cout << "\t-s/--seed\t[integer]\tSeed of random generation [optional, default = time]" << std::endl;
	std::cout << "\t-c/--chr\t[string]\tType of chromosome to simulate: A(autosome), X(X-chromosome) [optional, default = A]"<< std::endl; 
	std::cout << std::endl;
	std::cout << "\t-rec/--recRate\t[double]\tThe uniform recombination rate for all loci (unit: Morgan per base pair) [optional]" << std::endl;
	std::cout << "\t-mut/--mutRate\t[double]\tThe uniform mutation rate for all loci (unit: per generation per site) [optional]" << std::endl;
	std::cout << "\t-out/--outPrefix\t<string>\tThe prefix of all output files [optional, default = out]" << std::endl;
	
	std::cout << std::endl << "Options:" << std::endl;
	std::cout << "\t--no-rec\tNo recombination events in the simulation [optional]" << std::endl;
	std::cout << "\t--no-mut\tNo mutation events in the simulation [optional]" << std::endl;
	std::cout << "\t--no-sel\tNo selection events in the simulation [optional]" << std::endl; 
	std::cout << std::endl;
	std::cout << "\t-h/--help\tPrint help message [optional]" << std::endl;
	std::cout << "-------------------------------------------------------------------------" << std::endl;	
}

std::string Param::getModfile() const
{
	return modfile;
}

std::string Param::getSnvfile() const
{
	return snvfile;
}

std::string Param::getHapfile() const
{
	return hapfile;
}

std::string Param::getIndfile() const
{
	return indfile;
}

std::vector<std::string> Param::getpopList() const
{
	return popList;
}

std::vector<int> Param::getpopGen() const 
{
	return popGen;
}

std::vector<int> Param::getpopNumber() const
{
	return popNumber;
}

unsigned long Param::getSeed() const
{
	return seed;
}

std::string Param::getChr() const
{
	return chr;
}

std::string Param::getRecEvent() const
{
	return  recEvent;
}

std::string Param::getMutEvent() const
{
	return mutEvent;
}

std::string Param::getSelEvent() const
{
	return selEvent;
}

double Param::getRecRate() const
{
	return recRate;
}

double Param::getMutRate() const
{
	return mutRate;
}

bool Param::getGivenRecRate() const
{
	return givenRecRate;
}

bool Param::getGivenMutRate() const
{
	return givenMutRate;
}

bool Param::getGivenPop() const
{
	return givenPop;
}

bool Param::getGivenGen() const
{
	return givenGen;
}

std::string Param::getoutPrefix() const
{
	return outPrefix;
}


