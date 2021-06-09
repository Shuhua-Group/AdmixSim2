/****************************************
 *	Project Name: AdmixSim 2
 *	File Name: AdmixSim2.cpp
 *	Author: Rui Zhang
 *	Date: July 7, 2020
 ****************************************/

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <algorithm>
#include <ctime>
#include <cmath>
#include <numeric>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include "Population.hpp"
#include "Param.hpp"
#include "Utils.hpp"

int main(int argc, char **argv)
{
	std::string modfile;
	std::string snvfile;
	std::string ancfile;
	std::string indfile;
	unsigned long seed;
	std::string chr;
	std::vector<std::string> popList;
	std::vector<int> popGen;
	std::vector<int> popNumber;
	std::string recEvent;
	std::string mutEvent;
	std::string selEvent;
	std::string outPrefix;
	bool givenRecRate;
	bool givenMutRate;	
	bool givenPop;
	bool givenGen;

	std::string outPop;
	int outGen; 
	
	Param *par;

	par = new Param(argc, argv);
	modfile = par->getModfile();
	snvfile = par->getSnvfile();
	ancfile = par->getHapfile();
	indfile = par->getIndfile();
	popList = par->getpopList();
	popGen = par->getpopGen();
	popNumber = par->getpopNumber();
	seed = par->getSeed();
	chr = par->getChr();
	recEvent = par->getRecEvent();
	mutEvent = par->getMutEvent();
	selEvent = par->getSelEvent();
	givenRecRate = par->getGivenRecRate();
	givenMutRate = par->getGivenMutRate();
	givenPop = par->getGivenPop();
	givenGen = par->getGivenGen();

	outPrefix = par->getoutPrefix();

	std::string SELOut = outPrefix+".sel";
	std::ofstream selout(SELOut.c_str());
	if (!selout.good()) 
	{
		std::cerr << "Error: Cannot open file \"" << SELOut << "\"!" << std::endl;
		exit(1);
	}
	selout << "Gen\tPop\tPosition(s)\tCondition\tFrequency\tMale_Frequency\tFemale_Frequency" << std::endl;

	std::string logfile = outPrefix+".log";
	std::ofstream logout(logfile.c_str());
	if (!logout.good())
	{
		std::cerr << "Error: Cannot open file \"" << logfile << "\"!" << std::endl;
		exit(1);
	}

	logout << "AdmixSim 2 " << std::endl;
	logout << "Arguments and Options:" << std::endl;
	logout << "  Modfile = " << modfile << std::endl;
	logout << "  SNVfile = " << snvfile << std::endl;
	logout << "  Hapfile = " << ancfile << std::endl;
	logout << "  Indfile = " << indfile << std::endl;
	if (givenPop == true && givenGen == true)
	{
		logout << "  Output population list and corresponding generation and population size: ";
		for (int i = 0; i < popList.size(); i++)
		{
			logout << popList.at(i) << " " << popGen.at(i) << "g " << popNumber.at(i) << "; ";	
		}
		logout << std::endl;	
	}
	
	time_t startTime = time(0);
	char *dstart = ctime(&startTime);

	/**********************
     *  reading snv file  *
     **********************/
	clock_t snvstart = clock();
	std::ifstream sfin(snvfile.c_str());
	if (!sfin.good())
	{
		std::cerr << "Error: Cannot open file \"" << snvfile << "\"!" << std::endl;
		return 1;
	}
	int phyPos = -1;					/* physical position */
	std::vector<int> phyPoss;
	double genPos = -1;				/* genetic position / accumulative recombination distribution(according to uniform recombination rate) */
	std::vector<double> genPoss;
	std::vector<double> psegenPoss;
	double mutPos;				/* accumulative mutation distribution (according to uniform mutation rate) */
	std::vector<double> mutPoss;
	std::vector<double> mutrate;
	std::map<int, int> phyPos_index;
	int Np = 0;
	std::string snvLine;	
	while (getline(sfin, snvLine))
	{
		std::istringstream iss(snvLine);
		iss >> phyPos;
		if (phyPos <= 0)
		{
			std::cerr << "Error: The physical position is less than zero or not given in SNV information file!" << std::endl;
			exit(1);
		}
		if (phyPoss.size() != 0 && phyPos <= phyPoss.back())
		{
			std::cerr << "Error: The position " << phyPos << " is not in strictly ascending order!" << std::endl;
			exit(1);
		}	
		phyPoss.push_back(phyPos);
		phyPos_index[phyPos] = Np;
		Np += 1;

		iss >> genPos;
		if (recEvent == "1")
		{
			if (givenRecRate == false)
			{
				if (genPos < 0)
				{
					std::cerr << "Error: The genetic distance " << std::setprecision(10) << genPos << " in snv file is incorrect or uniform recombination rate is not given!" << std::endl;
					exit(1);
				}
				if (genPoss.size() != 0 && genPos <= genPoss.back())
				{
					std::cerr << "Error: The genetic distance " << std::setprecision(10) << genPos << " in snv file is not in strictly ascending order!" << std::endl;
					exit(1);
				}
				genPoss.push_back(genPos);
			}
			else		/* calculating genetic distance according to uniform recombination rate */
			{
				if (genPoss.size() == 0)
				{
					genPoss.push_back(0);
				}
				else
				{
					double recRate = par->getRecRate();
					genPos = genPoss.back() + recRate*(phyPoss.back() - phyPoss.at(phyPoss.size()-2));	
					genPoss.push_back(genPos);
				}
			}
		}
		else
		{
			psegenPoss.push_back(genPos);
		}

		if (mutEvent == "1")
		{
			if (givenMutRate == false)	/* calculate mutation distance according to mutation rate in the third column of snv file */
			{
				double mut = -1;
				iss >> mut;
				mutrate.push_back(mut);
				if (mut < 0)
				{
					std::cerr << "Error: The mutation rate " << std::setprecision(10) << mut << " in snv file is incorrect or uniform mutation rate is not given!" << std::endl;
					exit(1);
				}
				if (mutPoss.size() == 0)
				{
					mutPoss.push_back(0);
				}
				else
				{
					mutPos = mutPoss.back() + mut*(phyPoss.back() - phyPoss.at(phyPoss.size()-2));
					mutPoss.push_back(mutPos);
				}
			}
			else						/* calculate mutation distance according to uniform mutation rate */
			{
				if (mutPoss.size() == 0)
				{
					mutPoss.push_back(0);
				}
				else
				{
					double mutRate = par->getMutRate();
					mutPos = mutPoss.back() + mutRate*(phyPoss.back() - phyPoss.at(phyPoss.size()-2));
					mutPoss.push_back(mutPos);
				}
			}
		}
		else
		{
			double mut = -1;
			iss >> mut;
			mutrate.push_back(mut);
		}
	}
	sfin.close();
	clock_t snvend = clock();


	/********************* 
 	 *  reading indfile  *
 	 *********************/
	clock_t indstart = clock();
	bool AncSex = 1;
	std::ifstream ifin(indfile.c_str());
	if (!ifin.good())
	{
		std::cerr << "Error: Cannot open file \"" << indfile << "\"!" << std::endl;
		return 1;
	}
	std::map<std::string, int> AncNe; 			/* ancestry name and corresponding Ne */
	std::map<std::string, std::string> IndAnc;	/* individual name and corresponding Anc name */ 
	std::map<int, std::string> IndexInd;		/* int start from 0 */
	std::map<int, int> IndexSex;				/* index int start from 0 */
	std::vector<std::string> InitialRef;
	std::string indLine;
	std::vector<std::string> indIDSet;	
	std::map<std::string, int> PopMaleSex;
	std::map<std::string, int> PopFemaleSex;
	int index = 0;
	while (getline(ifin, indLine))
	{
		std::istringstream iss(indLine);
		std::string indID;
		iss >> indID;
		if (indID == "")
		{
			std::cerr << "Error: Individual ID is not given in individual information file!" << std::endl;
			exit(1);
		}
		if (std::find(indIDSet.begin(), indIDSet.end(), indID) != indIDSet.end())
		{
			std::cerr << "Error: Individual ID \"" << indID << "\" Repeat!" << std::endl;
			exit(1);
		}
		else
		{
			indIDSet.push_back(indID);
		}

		std::string Pop;
		iss >> Pop;
		if (Pop == "")
		{
			std::cerr << "Error: Population label is not given in individual information file!" << std::endl;
			exit(1);
		}
		if (PopMaleSex.find(Pop) == PopMaleSex.end())
		{
			PopMaleSex[Pop] = 0;
		}
		if (PopFemaleSex.find(Pop) == PopFemaleSex.end())
		{
			PopFemaleSex[Pop] = 0;
		}

		if (AncNe.find(Pop) != AncNe.end())
		{
			AncNe.at(Pop) += 1;
		}
		else
		{
			AncNe[Pop] = 1;
		}
		
		IndAnc[indID] = Pop;
		IndexInd[index] = indID;
		
		if (std::find(InitialRef.begin(), InitialRef.end(), Pop) == InitialRef.end())
		{
			InitialRef.push_back(Pop);
		}
	
		int Sex = -1;
		iss >> Sex;
		if (Sex != 0 && Sex != 1 && Sex != 2)
		{
			std::cerr << "Error: The sex of \"" << indID << "\" is incorrect or not given in individual information file!" << std::endl;
			exit(1);
		}
		if (chr == "X")
		{
			if (Sex == 0)
			{
				std::cerr << "Error: The sex of \"" << indID << "\" should be specified for X chromosome simulation!" << std::endl;
				exit(1);
			}
		}	
		if (AncSex == 1 && Sex == 0)
		{
			AncSex = 0;
		}
		IndexSex[index] = Sex;
		if (Sex == 1)
		{
			PopMaleSex.at(Pop) += 1;
		}  
		else if (Sex == 2)
		{
			PopFemaleSex.at(Pop) += 1;
		}
		index += 1;
	}
	ifin.close();
	clock_t indend = clock();


	/*********************
 	 *  reading modfile  *
 	 *********************/
	clock_t modstart = clock(); 
	std::ifstream mfin(modfile.c_str());
	if (!mfin.good())
	{
		std::cerr << "Error: Cannot open the model file \"" << modfile << "\"!" << std::endl;
		exit(1);
	}

	std::string modLine;
	std::set<std::string> AllPopList;
	std::map<std::string, std::vector<int> > PopGenSE;	/* population name || generation start and end index */
	std::vector<std::string> AllModel;
	std::vector<int> ModelIndex;
	int newMod = 0;
	while (getline(mfin, modLine))
	{
		/* pure comments */
		if ((modLine.size() > 0 && modLine.at(0) == '#') || modLine.size() == 0)
		{
			continue;
		}
		/* remove comments */
		size_t compos = modLine.find('#');
		if (compos != std::string::npos)
		{
			modLine = modLine.substr(0, compos);
		}
		
		AllModel.push_back(modLine);
		if (modLine.at(0) == '*')
		{	
			ModelIndex.push_back(newMod);
			std::istringstream iss(modLine);
			std::string aa;

			iss >> aa;
			size_t GenSplit = aa.find("-");
			if (GenSplit == std::string::npos)
			{
				std::cerr << "Error: The generation setting is incorrect in \"" << modLine << "\"!" << std::endl;
				exit(1);
			}
			int GenStart = std::atoi(aa.substr(1, GenSplit-1).c_str());
			int GenEnd = std::atoi(aa.substr(GenSplit+1, aa.size()-GenSplit).c_str());	

			if (givenGen == false)
			{
				outGen = GenEnd;
			}
		
			iss >> aa;
			std::vector<std::string> temp;
			boost::split(temp, aa, boost::is_any_of(","), boost::token_compress_on);
			for (int i = 0; i < temp.size(); i++)
			{
				AllPopList.insert(temp.at(i));
				if (PopGenSE.find(temp.at(i)) == PopGenSE.end())
				{
					PopGenSE[temp.at(i)].push_back(GenStart);
					PopGenSE[temp.at(i)].push_back(GenEnd);
				}
				else
				{
					int start, end;
					min(GenStart, PopGenSE[temp.at(i)].front(), start);
					max(GenEnd, PopGenSE[temp.at(i)].back(), end);
					std::vector<int> tempIndex;
					tempIndex.resize(2);
					tempIndex[0] = start;
					tempIndex[1] = end;
					PopGenSE[temp.at(i)] = tempIndex;
				}
			}
			if (givenPop == false)
			{
				outPop = temp.back();	
			}
		}
		newMod += 1;
	}
	ModelIndex.push_back(AllModel.size());
	
	if (givenGen == false)
	{
		popGen.push_back(outGen);	
	}
	if (givenPop == false)
	{
		popList.push_back(outPop);
	}
	
	/* maximum between sampled generation and model generation */
	for (int i = 0; i < popList.size(); i++)
	{
		if (PopGenSE.find(popList.at(i)) == PopGenSE.end())
		{	
			std::cerr << "Error: The population " << popList.at(i) << " in the argument \" -p/--poplist\" doesn't exist!" << std::endl;
			exit(1); 
		}
		std::vector<int> tempIndex;
		tempIndex.resize(2);
		tempIndex[0] = PopGenSE.at(popList.at(i)).front();
		max(popGen.at(i), PopGenSE.at(popList.at(i)).back(), tempIndex[1]);
		PopGenSE[popList.at(i)] = tempIndex;
	}

	clock_t modend = clock();
	
	if (givenPop == false || givenGen == false)
	{
		logout << "  Output population list and corresponding generation and population size: ";
 		for (int i = 0; i < popList.size(); i++)	
		{
			logout << popList.at(i) << " " << popGen.at(i) << "g " << popNumber.at(i) << "; ";
		}
		logout << std::endl;	
	}

	logout << "  Seed = " << seed << std::endl;
	logout << "  Chr = " << chr << std::endl;

	if (recEvent == "0")
	{
		logout << "  Recombination Setting: No recombination events in the simulation" << std::endl;
	}
	if (givenRecRate == true)
	{
		logout << "  Recombination Setting: Uniform recombination rate = " << par->getRecRate() << std::endl;
	}
	if (givenRecRate == false && recEvent == "1")
	{
		logout << "  Recombinaion Setting: Using locus-specific genetic distance in snv file" << std::endl;
	}

	if (mutEvent == "0")
	{
		logout << "  Mutation Setting: No mutation events in the simulation" << std::endl;
	}
	if (givenMutRate == true)
	{
		logout << "  Mutation Setting: Uniform mutation rate = " << par->getMutRate() << std::endl;
	}
	if (givenMutRate == false && mutEvent == "1")
	{
		logout << "  Mutation Setting: Using locus-specific mutation rate in snv file" << std::endl;
	}
	
	if (selEvent == "0")
	{
		logout << "  Selection Setting: No selection events in the simulation" << std::endl;
	}
	if (outPrefix == "")
	{
		logout << "  Out prefix = out (by defaut)" << std::endl;
	}
	else
	{
		logout << "  Out prefix = " << outPrefix << std::endl;
	}

	logout << "==========================================================================================================" << std::endl;

	logout << std::endl << "Simulation start! ";
	logout << dstart << std::endl;
	

	logout << "Reading snvfile ";
	logout << "time: " << (double)(snvend-snvstart)/CLOCKS_PER_SEC << "s" << std::endl;
	
	logout << "Reading indfile ";
	logout << "time: " << (double)(indend-indstart)/CLOCKS_PER_SEC << "s" << std::endl;

	logout << "Reading modfile ";
	logout << "time: " << (double)(modend-modstart)/CLOCKS_PER_SEC << "s" << std::endl;

	/**********************************
 	 *  delete useless ancestry data  *
 	 **********************************/
	for (std::map<std::string, int>::iterator iter = AncNe.begin(); iter != AncNe.end();)
	{
		if (AllPopList.find(iter->first) == AllPopList.end())
		{
			PopMaleSex.erase(iter->first);
			PopFemaleSex.erase(iter->first);
			PopGenSE.erase(iter->first);
			AncNe.erase(iter++);
		}
		else
		{
			iter++;
		}
	}

	/**********************
 	 *  construct AllPop  *
 	 **********************/
	std::map<std::string, Population> AllPop;					/* AncName || population structure */
	for (std::map<std::string, int>::iterator iter = AncNe.begin(); iter != AncNe.end(); iter++)	/* input reference population */
	{
		Population pop;
		Generation gen;
		pop.Gens[PopGenSE.at(iter->first).front()-1] = gen;
		AllPop[iter->first] = pop;

		if (PopMaleSex.at(iter->first) == iter->second)
		{
			std::cout << "Warning: Individuals in reference population " << iter->first << " are all male!" << std::endl;
		}
		else if (PopFemaleSex.at(iter->first) == iter->second)
		{
			std::cout << "Warning: Individuals in reference population " << iter->first << " are all female!" << std::endl;
		}
	}

	/*******************************************************************
 	 *  initial input reference population first generation minus one  *
 	 *******************************************************************/
	for (std::map<std::string, int>::iterator iter = AncNe.begin(); iter != AncNe.end(); iter++)
	{
		std::vector<int> temp;
		temp.resize(2);
		temp[0] = PopGenSE.at(iter->first).front()-1;
		temp[1] = PopGenSE.at(iter->first).back();
		PopGenSE[iter->first] = temp;
		AllPop.at(iter->first).setIndex(temp);

		std::vector<int> tempNe;
		tempNe.resize(1);
		tempNe[0] = iter->second;
		AllPop.at(iter->first).addNe(PopGenSE.at(iter->first).front(), tempNe);

		std::map<std::string, std::vector<double> > tempAncP;
		std::vector<double> tempP;
		tempP.resize(1);
		tempP[0] = 1;
		tempAncP[iter->first] = tempP;
		AllPop.at(iter->first).addAnc(PopGenSE.at(iter->first).front(), tempAncP);
			
		std::vector<int> tempRealNe;
		tempRealNe.resize(2);
		tempRealNe[0] = PopMaleSex.at(iter->first);
		tempRealNe[1] = PopFemaleSex.at(iter->first);
		AllPop.at(iter->first).addRealNe(PopGenSE.at(iter->first).front(), tempRealNe);
	}

	/**************************
 	 *  setting evolve model  *
 	 **************************/
	std::vector<int> AllSpos; 	/* ascending order */
	for (int j = 0; j < ModelIndex.size()-1; j++)
	{
		bool newPop = false;
		std::vector<std::vector<int> > SphyPos;
		std::vector<std::string> UserSpos; 
		std::vector<std::vector<std::string> > SAllele;
		std::vector<std::string> AdditionRule;

		int FirstColNum = 0; 
		std::string FirstLine = AllModel.at(ModelIndex.at(j));

		std::istringstream iss(FirstLine);

		std::string aa;
		iss >> aa;	
		FirstColNum += 1;
		size_t GenSplit = aa.find("-");
		int GenStart = std::atoi(aa.substr(1, GenSplit-1).c_str());
		int GenEnd = std::atoi(aa.substr(GenSplit+1, aa.size()-GenSplit).c_str());

		if (ModelIndex.at(j+1)-ModelIndex.at(j)-1 == 0)
		{
			std::cerr << "Error: The model unit with the beginning \"" << FirstLine << "\" doesn't have related evolution information!" << std::endl;
			exit(1);
		}

		if ((GenEnd-GenStart+1) != (ModelIndex.at(j+1)-ModelIndex.at(j)-1))
		{
			std::cerr << "Error: The number of generation is not equal to the setting in \"" << FirstLine << "\"!" << std::endl;
			exit(1);
		}

		iss >> aa;
		FirstColNum += 1;
		std::vector<std::string> TempPop;
		boost::split(TempPop, aa, boost::is_any_of(","), boost::token_compress_on);

		for (int i = 0; i < TempPop.size(); i++)
		{
			if (AllPop.find(TempPop.at(i)) == AllPop.end())
			{
				if (i != TempPop.size()-1 || TempPop.size() == 1)
				{
					std::cerr << "Error: The population \"" << TempPop.at(i) << "\" doesn't exist or more than one admixed population is given in one model unit!" << std::endl;
					exit(1);
				}
				else if (TempPop.size() == 2)
				{
					std::cerr << "Error: There should be more than two reference population when generating a new population!" << std::endl;
					exit(1);
				}
				else if (i == TempPop.size()-1 && TempPop.size() != 1)
				{
					Population pop;
					pop.setIndex(PopGenSE.at(TempPop.at(i)));
					AllPop[TempPop.at(i)] = pop;
					newPop = true;
				}
			}
		}
		
		if (selEvent == "1")
		{
			while (iss >> aa)
			{
				FirstColNum += 1;
				std::vector<std::string> temp1;
				std::vector<int> sphyPos;
				std::vector<std::string> sAllele;
				std::string additionRule;

				boost::split(temp1, aa, boost::is_any_of(":"), boost::token_compress_on);
				if (temp1.size() == 2)
				{
					additionRule = "1";
				}
				else if (temp1.size() == 3)
				{
					if (temp1.at(2) != "1" && temp1.at(2) != "2")
					{
						std::cerr << "Error: The addition rule setting of \"" << aa << "\" is incorrect!" << std::endl;
						exit(1);
					}
					additionRule = temp1.at(2);
				}
				else
				{
					std::cerr << "Error: The selection setting \"" << aa << "\" is incorrect!" << std::endl;
					exit(1);
				}
				std::string phyPosSet = temp1.at(0);
				std::string alleleSet = temp1.at(1);
				std::vector<std::string> tempPos;
				boost::split(tempPos, phyPosSet, boost::is_any_of(","), boost::token_compress_on);

				for (int i = 0; i < tempPos.size(); i++)
				{
					std::vector<std::string> start_end;
					boost::split(start_end, tempPos.at(i), boost::is_any_of("-"), boost::token_compress_on);
					if (start_end.size() == 2)
					{
						int start = std::atoi(start_end.front().c_str());
						int end = std::atoi(start_end.back().c_str());
						if (std::binary_search(phyPoss.begin(), phyPoss.end(), start) == false)
						{
							std::cerr << "Error: The position " << start << " is not in the SNV information file!" << std::endl;
							exit(1);
						}
						if (std::binary_search(phyPoss.begin(), phyPoss.end(), end) == false)
						{
							std::cerr << "Error: The position " << end << " is not in the SNV information file!" << std::endl;
							exit(1);
						}
						int startIndex, endIndex;
						findIntPos(phyPoss, start, startIndex);
						findIntPos(phyPoss, end, endIndex);
						if (end == phyPoss.back())
						{
							endIndex -= 1;
						}
						for (int k = startIndex; k <= endIndex; k++)
						{
							sphyPos.push_back(phyPoss.at(k));
						}	
					}
					else
					{
						int pos = std::atoi(start_end.front().c_str());
						if (std::binary_search(phyPoss.begin(), phyPoss.end(), pos) == false)
						{
							std::cerr << "Error: The position " << pos << " is not in the SNV information file!" << std::endl;
							exit(1);
						}
						sphyPos.push_back(pos);
					}
				}
			 
				if (AllPopList.find(alleleSet) != AllPopList.end() && std::find(InitialRef.begin(), InitialRef.end(), alleleSet) == InitialRef.end())
				{
					std::cerr << "Error: The selected ancestry must be one of the initial reference population, cannot be " << alleleSet << "!" << std::endl;
					exit(1);
				}
				else if (std::find(InitialRef.begin(), InitialRef.end(), alleleSet) != InitialRef.end())
				{
					sAllele.push_back(alleleSet);
				}
				else
				{
					boost::split(sAllele, alleleSet, boost::is_any_of(","), boost::token_compress_on);
					for (int i = 0; i < sAllele.size(); i++)
					{
						std::string temp = sAllele.at(i);
						if (AllPopList.find(temp) != AllPopList.end())
						{
							std::cerr << "Error: When involves ancestry selection, the setting \"" << aa << "\" should be written in two columns." << std::endl;
							exit(1);
						}
						if (sphyPos.size() != temp.length())
						{
							std::cerr << "Error: The selection condition \"" << aa << "\" is incorrect!" << std::endl;
							exit(1);		
						}
						if (count(sAllele.begin(), sAllele.end(), temp) != 1)
						{
							std::cerr << "Error: The selected condition of \"" << aa << "\" replicates!" << std::endl;
							exit(1);
						}
					}	 
				}	
				for (int i = 0; i < sphyPos.size(); i++)
				{
					if (std::binary_search(AllSpos.begin(), AllSpos.end(), sphyPos.at(i)) == false)
					{
						AllSpos.insert(lower_bound(AllSpos.begin(), AllSpos.end(), sphyPos.at(i)), sphyPos.at(i));
					}
				}	
				SphyPos.push_back(sphyPos);
				UserSpos.push_back(phyPosSet);
				SAllele.push_back(sAllele);
				AdditionRule.push_back(additionRule);	
			}
		}		

		for (int k = ModelIndex.at(j)+1; k < ModelIndex.at(j+1); k++)
		{
			int GenIndex = GenStart+k-ModelIndex.at(j)-1;
			
			std::string ModelLine = AllModel.at(k);
			std::istringstream iss (ModelLine);

			int LaterColNum = 0;
			std::vector<std::vector<int> > ne;	/* population size */
			std::vector<std::string> tempNe;
			std::string popNe;
			iss >> popNe;
			LaterColNum += 1;
			boost::split(tempNe, popNe, boost::is_any_of(":"), boost::token_compress_on);
			if (tempNe.size() != 1 && tempNe.size() != 2)
			{
				std::cerr << "Error: The setting of population size \"" << popNe << "\" is incorrect!" << std::endl;
				exit(1);
			}

			if (chr == "A" && AncSex == 0 && tempNe.size() == 2)
			{
				std::cerr << "Error: Population size \"" << popNe << "\" cannot be sex-specified for autosome simulation and no sex setting in initial ancestry population!"<< std::endl;
				exit(1);
			}

			std::vector<std::string> tempNe1;
			boost::split(tempNe1, tempNe.front(), boost::is_any_of(","), boost::token_compress_on);
			if (tempNe1.size() != 1 && tempNe1.size() != TempPop.size())
			{
				std::cerr << "Error: The setting of population size \"" << popNe << "\" is incorrect!" << std::endl;
				exit(1);
			}

			if (tempNe1.size() == 1)
			{
				for (int i = 1; i < TempPop.size(); i++)
				{
					tempNe1.push_back(tempNe1.front());
				}
			}

			int PN;
			for (int i = 0; i < tempNe1.size(); i++)
			{
				PN = std::atoi(tempNe1.at(i).c_str());
				if (PN < 0)
				{
					std::cerr << "Error: Population size \"" << PN << "\" can not be negative!" << std::endl;
					exit(1);
				}
				std::vector<int> temp;
				temp.push_back(PN);
				ne.push_back(temp);
			}

			if (tempNe.size() == 2)
			{
				std::vector<std::string> tempNe2;
				boost::split(tempNe2, tempNe.back(), boost::is_any_of(","), boost::token_compress_on);
				if (tempNe2.size() != 1 && tempNe2.size() != TempPop.size())
				{
					std::cerr << "Error: The setting of population size \"" << popNe << "\" is incorrect!" << std::endl;
					exit(1);
				}
				
				if (tempNe2.size() == 1)
				{
					for (int i = 1; i < TempPop.size(); i++)
					{
						tempNe2.push_back(tempNe2.front());
 					}	
				}

				for (int i = 0; i < tempNe2.size(); i++)
				{
					PN = std::atoi(tempNe2.at(i).c_str());
					if (PN < 0)
					{
						std::cerr << "Error: Population size \"" << PN << "\" can not be negative!" << std::endl;
						exit(1);
					}
					ne.at(i).push_back(PN);
				}
			}

			for (int i = 0; i < TempPop.size(); i++)
			{
				if (AllPop.at(TempPop.at(i)).GNe.count(GenIndex) > 0 && AllPop.at(TempPop.at(i)).GNe.at(GenIndex) != ne.at(i))
				{
					std::cerr << "Error: The population size setting for population " << TempPop.at(i) << " in generation " << GenIndex << " conflicts!" << std::endl;
					exit(1);
				}
				else if (AllPop.at(TempPop.at(i)).GNe.count(GenIndex) == 0)
				{
					AllPop.at(TempPop.at(i)).addNe(GenIndex, ne.at(i));
				}
			}	

			/* admixture proportion */
			std::string propSet;
			iss >> propSet;
			LaterColNum += 1;
			std::vector<std::string> tempProp;

			int nAnc; 			
			if (newPop == false)		/* no admixed population */
			{
				nAnc = TempPop.size();
			}
			else						/* new admixed population */
			{
				nAnc = TempPop.size()-1;
		
				std::map<std::string, std::vector<double> > AncProp;
				boost::split(tempProp, propSet, boost::is_any_of(":"), boost::token_compress_on);
				if (tempProp.size() != 1 && tempProp.size() != 2)
				{
					std::cerr << "Error: The proportion setting \"" << propSet << "\" is incorrect!" << std::endl;
					exit(1);
				}				

				if (chr == "A" && AncSex == 0 && tempProp.size() == 2)
				{
					std::cerr << "Error: Mixture proportion \"" << propSet  << "\" cannot be sex-specified for autosome and no sex setting in initial ancestry population!" << std::endl;
					exit(1); 
				}			

				std::vector<std::string> tempProp1;
				boost::split(tempProp1, tempProp.front(), boost::is_any_of(","), boost::token_compress_on);
				if (tempProp1.size() != TempPop.size())
				{	
					std::cerr << "Error: The proportion setting \"" << propSet << "\" is incorrect!" << std::endl;
					exit(1); 
				}
				if (k == ModelIndex.at(j)+1 && std::atof(tempProp1.back().c_str()) != 0)
				{
					std::cerr << "Error: Proportion from admixed population in founder generation should be zero!" << std::endl;
					exit(1);
				}
				double propSum = 0;
				for (int i = 0; i < tempProp1.size(); i++)
				{
					std::vector<double> temp;
					temp.push_back(std::atof(tempProp1.at(i).c_str()));
					if (temp.at(0) < 0 || temp.at(0) > 1)
					{
						std::cerr << "Error: Admixture proportion \"" << temp.at(0) << "\" must be between 0 and 1!" << std::endl;
						exit(1);
					}
		
					propSum += temp.at(0);
					AncProp[TempPop.at(i)] = temp;
				}
					
				if (tempProp.size() == 2)
				{
					std::vector<std::string> tempProp2;
					boost::split(tempProp2, tempProp.back(), boost::is_any_of(","), boost::token_compress_on);
					if (tempProp1.size() != tempProp2.size())
					{
						std::cerr << "Error: Male proportion setting number is not equal to female! " << propSet << std::endl;
						exit(1);  
					}
					if (k == ModelIndex.at(j)+1 && std::atof(tempProp2.back().c_str()) != 0)
					{
						std::cerr << "Error: Proportion from admixed population in founder generation should be zero!" << std::endl;
						exit(1);
					} 
					for (int i = 0; i < tempProp2.size(); i++)
					{
						if (std::atof(tempProp2.at(i).c_str()) < 0 || std::atof(tempProp2.at(i).c_str()) > 1)
						{
							std::cerr << "Error: Admixture proportion \"" << std::atof(tempProp2.at(i).c_str()) << "\" must be between 0 and 1!" << std::endl;
							exit(1);
						}
						AncProp.at(TempPop.at(i)).push_back(std::atof(tempProp2.at(i).c_str()));
						propSum += std::atof(tempProp2.at(i).c_str());
					}
				}
				if (fabs(propSum-1.0) > 1e-12)
				{
					std::cerr << "Error: The admixture proportion \"" << propSet << "\" must sum to 1!"<< std::endl;
					exit(1);
				}
				AllPop.at(TempPop.back()).addAnc(GenIndex, AncProp);
			}

			for (int i = 0; i < nAnc; i++)
			{
				std::map<std::string, std::vector<double> > AncProp;
				std::vector<double> prop;
				prop.resize(1);
				prop[0] = 1;
				AncProp[TempPop.at(i)] = prop;
				if (AllPop.at(TempPop.at(i)).GAncP.count(GenIndex) != 0)
				{
					if (AncProp != AllPop.at(TempPop.at(i)).GAncP.at(GenIndex))
					{
						std::cerr << "Error: The source population and corresponding proportion setting for population " << TempPop.at(i) << " in generation " << GenIndex << " conflicts!" << std::endl;
						exit(1);
					}
				}
				else
				{
					AllPop.at(TempPop.at(i)).addAnc(GenIndex, AncProp);
				}
			}

			/* selection coefficients */
			if (SphyPos.size() != 0)
			{
				std::map<std::string, std::vector<std::vector<double> > > Pop_Degree;
				for (int i = 0; i < TempPop.size(); i++)	
				{	
					std::vector<std::vector<double> > tempDegrees;
					for (int mm = 0; mm < SphyPos.size(); mm++)
					{
						std::vector<double> tempdegree;
						tempDegrees.push_back(tempdegree);
					}
					Pop_Degree[TempPop.at(i)] = tempDegrees;								
				}											

				std::string Coefficient;
				int index = 0;
				while (iss >> Coefficient)
				{
					LaterColNum += 1;
					if (LaterColNum > FirstColNum)
					{
						std::cerr << "Error: The number of selection coefficients setting columns is not equal to the number of selected segment!" << std::endl;
						exit(1);
					}
					std::vector<std::string> tempCo;
					boost::split(tempCo, Coefficient, boost::is_any_of(":"), boost::token_compress_on);
					if (tempCo.size() != 1 && tempCo.size() != 2)
					{
						std::cerr << "Error: The setting of selection coefficients \"" << Coefficient << "\" is incorrect!" << std::endl;
						exit(1);
					}

					if (chr == "A" && AncSex == 0 && tempCo.size() == 2)
					{
						std::cerr << "Error: Selection coefficients \"" << Coefficient << "\" cannot be sex-specified for autosome and no sex setting in reference population!" << std::endl;
						exit(1);
					}
					
					std::vector<std::string> tempCo1;
					boost::split(tempCo1, tempCo.front(), boost::is_any_of(","), boost::token_compress_on);

					if (tempCo1.size() != 1 && tempCo1.size() != TempPop.size())
					{
						std::cerr << "Error: The setting of selection coefficients \"" << Coefficient << "\" is incorrect!" << std::endl;
						exit(1);
					}
					
					if (tempCo1.size() == 1)
					{
						for (int i = 1; i < TempPop.size(); i++)
						{
							tempCo1.push_back(tempCo1.front());
						}
					}
		
					double coefficient;
					for (int i = 0; i < tempCo1.size(); i++)
					{
						coefficient = std::atof(tempCo1.at(i).c_str());
						if (coefficient < -1)
						{
							std::cerr << "Error: The selection coefficient \"" << coefficient << "\" is smaller than -1!" << std::endl;
							exit(1); 
						}
						Pop_Degree.at(TempPop.at(i)).at(index).push_back(coefficient);
					}

					if (tempCo.size() == 2)
					{		
						std::vector<std::string> tempCo2;
						boost::split(tempCo2, tempCo.back(), boost::is_any_of(","), boost::token_compress_on);
						if (tempCo2.size() != 1 && tempCo2.size() != TempPop.size())
						{
							std::cerr << "Error: The setting of selection coefficients \"" << Coefficient << "\" is incorrect!" << std::endl;
							exit(1);
						}

						if (tempCo2.size() == 1)
						{
							for (int i = 1; i < TempPop.size(); i++)
							{
								tempCo2.push_back(tempCo2.front());
							}
						}
					
						for (int i = 0; i < tempCo2.size(); i++)
						{
							coefficient = std::atof(tempCo2.at(i).c_str());
							if (coefficient < -1)
							{
								std::cerr << "Error: The selection coefficient \"" << coefficient << "\" is smaller than -1!" << std::endl;
								exit(1);
							}
							Pop_Degree.at(TempPop.at(i)).at(index).push_back(coefficient);
						}
					}
					index += 1;
				}

				if (LaterColNum != FirstColNum)
				{
					std::cerr << "Error: The number of selection coefficients setting columns is not equal to the number of selected segment!" << std::endl;
					exit(1);
				}

				for (int i = 0; i < TempPop.size(); i++)
				{
					for (int m = 0; m < Pop_Degree.at(TempPop.at(i)).size(); m++)
					{
						if (Pop_Degree.at(TempPop.at(i)).at(m).front() != 0 || Pop_Degree.at(TempPop.at(i)).at(m).back() != 0)
						{
							AllPop.at(TempPop.at(i)).addSel(GenIndex, SphyPos.at(m), UserSpos.at(m), SAllele.at(m), AdditionRule.at(m), Pop_Degree.at(TempPop.at(i)).at(m));
						}	
					}
				}
			}
		}
	}


	/*********************
 	 *  reading hapfile  *
 	 *********************/
	logout << "Reading hapfile ";
	clock_t hapstart = clock();
	std::ifstream hfin(ancfile.c_str());
	if (!hfin.good())
	{
		std::cerr << "Error: Cannot open file \"" << ancfile << "\"!" << std::endl;
		return 1;
	}

	std::map<std::string, std::vector<std::string> > ancInds;	/* AncName || IndName */
	for (std::map<std::string, int>::iterator iter = AncNe.begin(); iter != AncNe.end(); iter++)
	{
		std::vector<std::string> IndNames;
		ancInds[iter->first] = IndNames;
	}

	for (int i = 0; i < IndexInd.size(); i++)
	{
		std::string hapLine1;
		std::string hapLine2;
		getline(hfin, hapLine1);
		getline(hfin, hapLine2);
	
		if (hapLine1 == "" || hapLine2 == "")
		{
			std::cerr << "Error: The number of individual in hap file is not equal to that in ind file!" << std::endl;
			exit(1); 
		}
		
		if (chr == "A" && (hapLine1.at(0) == '9' || hapLine2.at(0) == '9'))
		{
			std::cerr << "Error: The sequence of \"" << IndexInd.at(i) << "\" is incorrect for autosome simulation!" << std::endl;
			exit(1);
		}
	
		if (hapLine1.size() != phyPoss.size() || hapLine2.size() != phyPoss.size())
		{
			std::cerr << "Error: The haplotype length of \"" << IndexInd.at(i) << "\" is not equal to the number of snv!" << std::endl;
			exit(1);
		}
		
		int Sex = IndexSex.at(i);
		
		if (chr == "X" && ((Sex == 1 && hapLine1.at(0) != '9' && hapLine2.at(0) != '9') || (Sex == 2 && (hapLine1.at(0) == '9' || hapLine2.at(0) == '9'))))
		{
			std::cerr << "Error: The sex setting of \"" << IndexInd.at(i) << "\" is incorrect!" << std::endl;
			exit(1);
		}

		if (chr == "X" && Sex == 1 && hapLine1.at(0) == '9')
		{
			std::string temphap = hapLine2;
			hapLine2 = hapLine1;
			hapLine1 = temphap;	
		}

		if (AllPopList.find(IndAnc.at(IndexInd.at(i))) != AllPopList.end())
		{
			int hap1Name = 2*i;
			int hap2Name = 2*i+1;
		
			ancInds.at(IndAnc.at(IndexInd.at(i))).push_back(IndexInd.at(i));

			std::vector<Segment> segs1, segs2;
			segs1.resize(1);
			segs2.resize(1);
			segs1[0] = Segment(phyPoss.front(), phyPoss.back(), hap1Name);
			segs2[0] = Segment(phyPoss.front(), phyPoss.back(), hap2Name);
		
			std::vector<int> mutationpoints1;
			std::vector<int> mutationpoints2;		
			std::map<int, char> Sallele1;
			std::map<int, char> Sallele2;

			Haplotype hap1(segs1, mutationpoints1, Sallele1);
			Haplotype hap2(segs2, mutationpoints2, Sallele2);
			hap1.setSallele(hapLine1, phyPos_index, AllSpos);
			hap2.setSallele(hapLine2, phyPos_index, AllSpos);
		
			Individual ind(hap1, hap2, Sex);

			AllPop.at(IndAnc.at(IndexInd.at(i))).Gens.at(PopGenSE.at(IndAnc.at(IndexInd.at(i))).front()).inds.push_back(ind);
		}
	}
	hfin.close();
	clock_t hapend = clock();
	logout << "time: " << (double)(hapend-hapstart)/CLOCKS_PER_SEC << "s" << std::endl << std::endl;


	for (std::map<std::string, Population>::iterator iter = AllPop.begin(); iter != AllPop.end(); iter++)
	{
		(iter->second).complement(iter->first);
	}	

	for (std::map<std::string, Population>::iterator iter = AllPop.begin(); iter != AllPop.end(); iter++)
	{
		std::map<int, std::map<std::string, std::vector<double> > > TempGAncP = (iter->second).GAncP;
		for (int i = (iter->second).GIndices.front(); i <= (iter->second).GIndices.back(); i++)
		{
			for (std::map<std::string, std::vector<double> >::iterator iter1 = TempGAncP.at(i).begin(); iter1 != TempGAncP.at(i).end(); iter1++)
			{
				if ((iter1->second).front() != 0 && (iter1->second).back() != 0)
				{
					if ((i-1) > PopGenSE.at(iter1->first).back())
					{
						std::vector<int> tempIndex;
						tempIndex.resize(2);
						tempIndex[0] = PopGenSE.at(iter1->first).front();
						tempIndex[1] = i-1;
						PopGenSE[iter1->first] = tempIndex;	
					}
				}	
			} 
		}
	}

	std::map<int, std::vector<std::string> > IndexPopSet;  /* Generation index || population name set */
	for (std::map<std::string, Population>::iterator iter = AllPop.begin(); iter != AllPop.end(); iter++)
	{
		AllPop.at(iter->first).setIndex(PopGenSE[iter->first]);
		(iter->second).complement(iter->first);

		for (int i = (iter->second).GIndices.front(); i <= (iter->second).GIndices.back(); i++)	
		{
			if (IndexPopSet.find(i) == IndexPopSet.end())
			{
				std::vector<std::string> temp;
				temp.resize(1);
				temp[0] = iter->first;
				IndexPopSet[i] = temp;
			}
			else if (IndexPopSet.find(i) != IndexPopSet.end())
			{
				IndexPopSet.at(i).push_back(iter->first);
			}
		}

		for (std::map<int, std::vector<int> >::iterator iter1 = (iter->second).GNe.begin(); iter1 != (iter->second).GNe.end(); iter1++)
		{
			if ((iter1->second).size() == 2 && (iter1->second).front() == 0 && (iter1->second).back() != 0)
			{
				std::cout << "Warning: Individuals in population " << (iter->first) << " at generation " << (iter1->first) << " are all female!" << std::endl;
			}
			if ((iter1->second).size() == 2 && (iter1->second).back() == 0 && (iter1->second).front() != 0)
			{
				std::cout << "Warning: Individuals in population " << (iter->first) << " at generation " << (iter1->first) << " are all male!" << std::endl;
			}
		}
	}

	/******************************************************************
 	 *  check the output arguments -p/--poplist  -g/--gen  -n/--nInd  *
 	 ******************************************************************/
	for (int i = 0; i < popList.size(); i++)
	{
		if (AllPop.find(popList.at(i)) == AllPop.end())
		{
			std::cerr << "Error: The population " << popList.at(i) << " in the argument \"-p/--poplist\" doesn't exist!" << std::endl;
			exit(1); 
		}
		if (popGen.at(i) < AllPop.at(popList.at(i)).GIndices.front())
		{
			std::cerr << "Error: The population " << popList.at(i) << " doesn't have generation index " << popGen.at(i) << " in the argument \"-g/--gen\"!" << std::endl;
			exit(1);
		}
		int temp = accumulate(AllPop.at(popList.at(i)).GNe.at(popGen.at(i)).begin(), AllPop.at(popList.at(i)).GNe.at(popGen.at(i)).end(), 0);
		if (popNumber.at(i) > temp)
		{
			std::cerr << "Error: The output individuals of population " << popList.at(i) << " at generation " << popGen.at(i) << " in the argument \"-n/--nInd\" is larger than its population size!" << std::endl;
			exit(1);
		}
	}

	std::map<std::pair<std::string, int>, int> ToSample;
	for (int i = 0; i < popList.size(); i++)
	{
		ToSample[std::pair<std::string, int>(popList.at(i), popGen.at(i))] = popNumber.at(i);
	}

	/**********************************
     *  first generation setIndRange  *
     **********************************/
	for (std::map<std::string, int>::iterator iter = AncNe.begin(); iter != AncNe.end(); iter++)
	{
		AllPop.at(iter->first).Gens.at(PopGenSE.at(iter->first).front()).setIndRange();
	}	

	/****************************
	*  Initial ancestry sample  *
	*****************************/
	std::map<std::pair<std::string, int>, std::vector<Individual> > SampledInd;
	for (std::map<std::string, int>::iterator iter = AncNe.begin(); iter != AncNe.end(); iter++)
	{
		if (ToSample.count(std::pair<std::string, int>(iter->first, AllPop.at(iter->first).GIndices.front())) != 0)
		{
			std::vector<Individual> sample;
			AllPop.at(iter->first).Gens.at(AllPop.at(iter->first).GIndices.front()).sample(ToSample.at(std::pair<std::string, int>(iter->first, AllPop.at(iter->first).GIndices.front())), sample);
			SampledInd[std::pair<std::string, int>(iter->first, AllPop.at(iter->first).GIndices.front())] = sample;
		}
	}

	/********************
     *  evolve process  *
     ********************/
	clock_t eStart = clock();	//
	std::vector<int> IndexSet;
	for (std::map<int, std::vector<std::string> >::iterator iter = IndexPopSet.begin(); iter != IndexPopSet.end(); iter++)
	{
		IndexSet.push_back(iter->first);
	}
	std::sort(IndexSet.begin(), IndexSet.end());
	
	for (int i = IndexSet.front(); i <= IndexSet.back(); i++)
	{
		std::cout << "Generation " << i << std::endl;	//

		clock_t	genstart = clock(); 
		if (i != IndexSet.front())
		{
			for (int j = 0; j < IndexPopSet.at(i).size(); j++)
			{
				std::string popname = IndexPopSet.at(i).at(j);
				if ((AncNe.find(popname) != AncNe.end() && i != AllPop.at(popname).GIndices.front()) || AncNe.find(popname) == AncNe.end())
				{
					if (accumulate(AllPop.at(popname).GNe.at(i).begin(), AllPop.at(popname).GNe.at(i).end(), 0) != 0)
					{
				//		std::cout << popname << " evolve " << std::endl;	//
						std::vector<std::string> alleleString;
						AllPop.at(popname).evolve(popname, i, AllPop, chr, recEvent, phyPoss, genPoss, InitialRef, IndAnc, AncSex, alleleString, AllSpos, IndexInd);	
						for (int k = 0; k < alleleString.size(); k++)
						{
							selout << i << "\t" << popname << "\t" << alleleString.at(k) << std::endl;
						}
					}
					else
					{
						Generation nullGen;
						AllPop.at(popname).Gens[i] = nullGen;
						std::vector<int> nullRealNe;
						nullRealNe.resize(2);
						nullRealNe[0] = 0;
						nullRealNe[1] = 0;
						AllPop.at(popname).addRealNe(i, nullRealNe);	
					}
					/* Sample out */
					if (ToSample.count(std::pair<std::string, int>(popname, i)) != 0)
					{
						std::vector<Individual> sample;
 						AllPop.at(popname).Gens.at(i).sample(ToSample.at(std::pair<std::string, int>(popname, i)), sample);
						//std::vector<Individual> sample = AllPop.at(popname).Gens.at(i).inds; //
						SampledInd[std::pair<std::string, int>(popname, i)] = sample;
					}
				}
			}
		}
		
		/* add de novo mutation to populations in current generation*/
		clock_t mustart = clock();		
		if (mutEvent == "1" && i != IndexSet.back())
		{
			std::set<int> ExistPoss;
			if (i != IndexSet.front()) 	 /* addExistPoss */
			{
				for (int j = 0; j < IndexPopSet.at(i-1).size(); j++)
				{
					for (int k = 0; k < AllPop.at(IndexPopSet.at(i-1).at(j)).Gens.at(i-1).inds.size(); k++)
					{
						int m = 2;
						if (chr == "X" && AllPop.at(IndexPopSet.at(i-1).at(j)).Gens.at(i-1).inds.at(k).Sex == 1)
						{
							m = 1;
						}
						for (int n = 0; n < m; n++)
						{
							std::vector<int> Mutations;
							AllPop.at(IndexPopSet.at(i-1).at(j)).Gens.at(i-1).inds.at(k).getHapMut(n, Mutations);
							for (int nn = 0; nn < Mutations.size(); nn++)
							{
								ExistPoss.insert(Mutations.at(nn));
							}
						}
					}
				}
			}

			for (int j = 0; j < IndexPopSet.at(i).size(); j++)
			{
			//	std::cout << IndexPopSet.at(i).at(j) << " Mutation" << std::endl;	// 
				AllPop.at(IndexPopSet.at(i).at(j)).Gens.at(i).GenMuPoints(phyPoss, mutPoss, chr, ExistPoss);
			}
		}
		
		clock_t genend = clock();
	//	logout << "evolve time: " << (double)(mustart-genstart)/CLOCKS_PER_SEC << "s; " << std::endl;	//
	//	std::cout << "mutation time: " << (double)(genend-mustart)/CLOCKS_PER_SEC << "s; ";		//
	//	std::cout << std::endl;	//

		clock_t deleStart = clock();	//		
		/* delete useless data */
		if (i > IndexSet.front() && i != IndexSet.back())
		{
			for (std::map<std::string, Population>::iterator iter = AllPop.begin(); iter != AllPop.end();)
			{
				if (std::find((iter->second).GIndices.begin(), (iter->second).GIndices.end(), i-1) != (iter->second).GIndices.end())
				{
					(iter->second).Gens.erase(i-1);
					(iter->second).GNe.erase(i-1);
					(iter->second).GAncP.erase(i-1);
					(iter->second).GSphyPos.erase(i-1);
					(iter->second).GphyPos_allele.erase(i-1);
					(iter->second).GphyPos_Addition.erase(i-1);
					(iter->second).GSdegree.erase(i-1);
					(iter->second).GRealNe.erase(i-1);
						
					if ((iter->second).Gens.size() == 0)
					{
						AllPop.erase(iter++);
					}
					else
					{
						iter++;
					}
				}
				else
				{
					iter++;
				}
			}
		}		
		clock_t deleEnd = clock();	//
	//	logout << "delete time: " << (double)(deleEnd-deleStart)/CLOCKS_PER_SEC << "s; ";	//
	//	logout << "generation time: " << (double)(genend-genstart)/CLOCKS_PER_SEC << "s" << std::endl;	
	}
	clock_t eEnd = clock();
	logout << "Simulation time: " << (double)(eEnd-eStart)/CLOCKS_PER_SEC << "s" << std::endl << std::endl;

	AllPop.clear();

	/************************************
	 *  saving ancestral sequence data  *
	 ************************************/
	clock_t saveStart = clock();
	std::map<std::string, std::vector<std::string> > anchaps;	/* AncName || sequence */
	for (std::map<std::string, int>::iterator iter = AncNe.begin(); iter != AncNe.end(); iter++)
	{
		std::vector<std::string> haps;
		haps.reserve((iter->second)*2);
		anchaps[iter->first] = haps;
	}

	std::ifstream hfin1(ancfile.c_str());
	for (int i = 0; i < IndexInd.size(); i++)
	{
		std::string hapLine1;
		std::string hapLine2;
		getline(hfin1, hapLine1);
		getline(hfin1, hapLine2);

		if (AllPopList.find(IndAnc.at(IndexInd.at(i))) != AllPopList.end())
		{		
			anchaps.at(IndAnc.at(IndexInd.at(i))).push_back(hapLine1);
			anchaps.at(IndAnc.at(IndexInd.at(i))).push_back(hapLine2);
		}
	}
	clock_t saveEnd = clock();
	logout << "Saving ancestral sequence data time: " << (double)(saveEnd-saveStart)/CLOCKS_PER_SEC << "s" << std::endl << std::endl;

	/****************************
	 *  Update SNV information  *
	 ****************************/
	std::vector<int> FinalPoss;     	/* only mutation position of output population */
	std::vector<int> MutIndex;  		/* the index of mapped mutation points in phyPoss */
	std::map<int, std::vector<char> > MutAlleleSet;

	std::string SnvOut = outPrefix+".snv";
	std::ofstream snvout(SnvOut.c_str());
	if (!snvout.good())
	{
		std::cerr << "Error: Cannot open file \"" << SnvOut << "\"!" << std::endl;
		exit(1);
	}

	clock_t setStart = clock();	
	for (std::map<std::pair<std::string, int>, std::vector<Individual> >::iterator iter = SampledInd.begin(); iter != SampledInd.end(); iter++)
	{	
		for (int i = 0; i < (iter->second).size(); i++)
		{
			int k = 2;
			if (chr == "X" && (iter->second).at(i).Sex == 1)
			{
				k = 1;
			}
			for (int j = 0; j < k; j++)
			{
				std::vector<int> Mutations;
				(iter->second).at(i).getHapMut(j, Mutations);
				for (int m = 0; m < Mutations.size(); m++)
				{
					if (std::binary_search(FinalPoss.begin(), FinalPoss.end(), Mutations.at(m)) == false)
					{
						FinalPoss.insert(lower_bound(FinalPoss.begin(), FinalPoss.end(), Mutations.at(m)), Mutations.at(m));
					}
				}
			}
		}
	}

	MutIndex.resize(FinalPoss.size());
	for (int i = 0; i < FinalPoss.size(); i++)
	{		
		findIntPos(phyPoss, FinalPoss.at(i), MutIndex[i]);
	}

	/* setting the reference and alternative allele of each new mutation point */
	char dataType = 'A';
	if (anchaps.at(IndAnc.at(IndexInd.at(0))).at(0).at(0) == '0' || anchaps.at(IndAnc.at(IndexInd.at(0))).at(0).at(0) == '1')
	{
		dataType = 'B';
	}

	if (dataType == 'A')
	{
		std::string Alleles = "AGCT";
		for (int i = 0; i < FinalPoss.size(); i++)
		{
			std::vector<char> tempAlleleSet;
			tempAlleleSet.resize(2);
			int rand1 = rand()%4;
			int rand2 = rand()%4;
			while (rand2 == rand1)
			{
				rand2 = rand()%4;
			}
			tempAlleleSet[0] = Alleles.at(rand1);
			tempAlleleSet[1] = Alleles.at(rand2);
			MutAlleleSet[FinalPoss.at(i)] = tempAlleleSet;
		} 
	}
	else if (dataType == 'B')
	{
		std::vector<char> tempAlleleSet;
		tempAlleleSet.resize(2);
		tempAlleleSet[0] = '0';
		tempAlleleSet[1] = '1';
		for (int i = 0; i < FinalPoss.size(); i++)
		{
			MutAlleleSet[FinalPoss.at(i)] = tempAlleleSet;
		}
	}

	int start_index = 0;
	double glMutRate = 0;
	if (givenMutRate == true)	
	{
		glMutRate = par->getMutRate();
	}
	for (int i = 0; i < MutIndex.size(); i++)
	{	
		int end_index = MutIndex.at(i);
		for (int j = start_index; j < end_index; j++)
		{
			snvout << phyPoss.at(j) << "\t" << std::setprecision(10) << genPoss.at(j) << "\t";
			if (givenMutRate == false)
			{
				snvout << std::setprecision(10) << mutrate.at(j);
			}
			else
			{
				snvout << std::setprecision(10) << glMutRate;
			}
			snvout << "\t" << "F" << std::endl;
		}	
		snvout << FinalPoss.at(i) << "\t" << std::setprecision(10) << (genPoss.at(end_index)-genPoss.at(end_index-1))/(phyPoss.at(end_index)-phyPoss.at(end_index-1))*(FinalPoss.at(i)-phyPoss.at(end_index-1))+genPoss.at(end_index-1) << "\t";
		if (givenMutRate == false)
		{
			snvout << std::setprecision(10) << mutrate.at(end_index);
		}	
		else
		{
			snvout << std::setprecision(10) << glMutRate;
		}
		snvout << "\t" << "T," << MutAlleleSet.at(FinalPoss.at(i)).at(0) << "/" << MutAlleleSet.at(FinalPoss.at(i)).at(1) << std::endl;
		start_index = end_index;
	}

	for (int j = start_index; j < phyPoss.size(); j++)
	{
		snvout << phyPoss.at(j) << "\t";
		if (recEvent == "1")
		{
			snvout << std::setprecision(10) << genPoss.at(j) << "\t";
		}
		else
		{
			snvout << std::setprecision(10) << psegenPoss.at(j) << "\t";
		}
		
		if (mutEvent == "0" || givenMutRate == false)
		{
			snvout << std::setprecision(10) << mutrate.at(j);
		}
		else
		{
			snvout << std::setprecision(10) << glMutRate;
		}
		snvout << "\t" << "F" << std::endl;
	}

	snvout.close();
	clock_t setEnd = clock();
	logout << "Final add de novo mutation and output snvfile time: " << (double)(setEnd-setStart)/CLOCKS_PER_SEC << "s" << std::endl << std::endl;	

//	std::cout << std::endl;	//

	/************
	 *  output  *
	 ************/
	std::string HapOut = outPrefix+".hap"; 
	std::string SegOut = outPrefix+".seg";
	std::string IndOut = outPrefix+".ind"; 
	std::ofstream hapout(HapOut.c_str());
	std::ofstream segout(SegOut.c_str());
	std::ofstream indout(IndOut.c_str());
	
//	std::string MutOut = outPrefix+".mut";	//
//	std::ofstream mutout(MutOut.c_str());	//
//	std::string RecOut = outPrefix+".rec"; 	//
//	std::ofstream recout(RecOut.c_str());	//

	if (!hapout.good())
	{
		std::cerr << "Error: Cannot open file \"" << HapOut << "\"!" << std::endl;
		exit(1);
	}
	if (!segout.good())
	{
		std::cerr << "Error: Cannot open file \"" << SegOut << "\"!" << std::endl;
		exit(1);
	}
	if (!indout.good())
	{
		std::cerr << "Error: Cannot open file \"" << IndOut << "\"!" << std::endl;
		exit(1);
	}
	
	int kkk = 1;
	for (int ii = 0; ii < popList.size(); ii++)
	{
		clock_t outstart = clock();
		logout << "Population " << popList.at(ii) << " at generation " << popGen.at(ii) << " output ";
		
		for (int i = 0; i < popNumber.at(ii); i++)
		{
			indout << "Ind" << kkk << "\t" << popList.at(ii)+"_" << popGen.at(ii) << "\t" << SampledInd.at(std::pair<std::string, int>(popList.at(ii), popGen.at(ii))).at(i).Sex << std::endl;	
			int nhap = 2;
			if (chr == "X" && (SampledInd.at(std::pair<std::string, int>(popList.at(ii), popGen.at(ii))).at(i)).Sex == 1)
			{
				nhap = 1;
			}
			for (int j = 0; j < nhap; j++)
			{
				Haplotype hap;
				SampledInd.at(std::pair<std::string, int>(popList.at(ii), popGen.at(ii))).at(i).getHap(j, hap);

			//	mutout << hap.mutationpoints.size() << std::endl;	//
			//	recout << hap.segments.size()-1 << std::endl;	//

				segout << "Ind" << kkk << " Hap " << (j+1) << std::endl;

				std::string outStr = "";
				outStr.reserve(FinalPoss.size()+phyPoss.size());
			
				int SegSize = hap.segments.size();		
				for (int k = 0; k < SegSize; k++)
				{
					Segment seg = hap.segments.at(k);
					int label = seg.label;
					std::string anc = IndAnc.at(IndexInd.at(label/2));
					std::vector<std::string> tempInds = ancInds.at(anc);
					int index = 2*(std::distance(tempInds.begin(), std::find(tempInds.begin(), tempInds.end(), IndexInd.at(label/2))))+label%2;
					int left, right;
					findIntPos(phyPoss, seg.start, left);
					findIntPos(phyPoss, seg.end, right);
					outStr += anchaps.at(anc).at(index).substr(left, right-left);
				}

				if (FinalPoss.size() != 0)
				{
					std::vector<int> hapMutation = hap.mutationpoints;
					int start_index = 0;
					int muIndex = 0;
					for (int k = 0, len = MutIndex.size(); k < len; k++)
					{
						int end_index = MutIndex.at(k);
						hapout << outStr.substr(start_index, end_index-start_index);
						
						if (hapMutation.size() != 0 && hapMutation.at(muIndex) == FinalPoss.at(k))
						{
							hapout << MutAlleleSet.at(FinalPoss.at(k))[1];
							if (muIndex != hapMutation.size()-1)
							{
								muIndex++;
							}
						}
						else
						{
							hapout << MutAlleleSet.at(FinalPoss.at(k))[0];
						}
						start_index = end_index;
					}
					hapout << outStr.substr(start_index, outStr.size()-start_index) << std::endl;
				}
				else
				{
					hapout << outStr << std::endl;	
				}	

				if (recEvent == "1")
				{
					hap.smooth(IndAnc, IndexInd);	
					
					int SegSize = hap.segments.size();
					for (int k = 0; k < SegSize; k++)
					{
						Segment seg = hap.segments.at(k);

						int start_index, end_index;
						findIntPos(phyPoss, seg.start, start_index);
						findIntPos(phyPoss, seg.end, end_index);
						double start_gen = genPoss.front();
						double end_gen = genPoss.back();
						if (start_index != 0)
						{
							start_gen = (genPoss.at(start_index)-genPoss.at(start_index-1))/(phyPoss.at(start_index)-phyPoss.at(start_index-1))*(seg.start-phyPoss.at(start_index-1))+genPoss.at(start_index-1);
						}				
						if (end_index != phyPoss.size())
						{
							end_gen = (genPoss.at(end_index)-genPoss.at(end_index-1))/(phyPoss.at(end_index)-phyPoss.at(end_index-1))*(seg.end-phyPoss.at(end_index-1))+genPoss.at(end_index-1);
						}
	
						segout << seg.start << "\t" << std::setprecision(10) << start_gen << "\t" << seg.end << "\t" << std::setprecision(10) << end_gen << "\t" << IndAnc.at(IndexInd.at(seg.label/2)) << std::endl;		
					//	segout << seg.start << "\t" << std::setprecision(10) << start_gen << "\t" << seg.end << "\t" << std::setprecision(10) << end_gen << "\t" << IndAnc.at(IndexInd.at(seg.label/2)) << "\t" << IndexSex.at(seg.label/2) << std::endl;	//
					}
				}
				else
				{
					Segment seg = hap.segments.at(0);
					segout << seg.start << "\t" << std::setprecision(10) << psegenPoss.front() << "\t" << seg.end << "\t" << std::setprecision(10) << psegenPoss.back() << "\t" << IndAnc.at(IndexInd.at(seg.label/2)) << std::endl;

				}
			}
			if (nhap == 1)
			{
				std::string outStr(FinalPoss.size()+phyPoss.size(), '9');
				hapout << outStr << std::endl;
			}
			kkk += 1;
		}	
		clock_t outend = clock();
		logout << "time: " << (double)(outend-outstart)/CLOCKS_PER_SEC << "s" << std::endl;
	}
	hapout.close();
	segout.close();
	indout.close();
	logout << std::endl << "Simulation end! ";
	time_t endTime = time(0);
	char *dend = ctime(&endTime);
	logout << dend << std::endl;
	logout.close();
}


