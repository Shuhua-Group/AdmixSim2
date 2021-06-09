/****************************************
 *	Project Name: AdmixSim 2
 *	File Name: Generation.cpp
 *	Author: Rui Zhang
 *	Date: July 7, 2020
 ****************************************/

#include <cstdlib>
#include <algorithm>
#include <iostream>	
#include "Generation.hpp"
#include "Haplotype.hpp"
#include "Utils.hpp"

Generation::Generation()
{
}

Generation::~Generation()
{
}

void Generation::GenMuPoints(const std::vector<int> &phyPoss, const std::vector<double> &mutPoss, const std::string &chr, std::set<int> &ExistPoss)
{
	for (int i = 0; i < inds.size(); i++)
	{
		inds.at(i).haplotype1.HapMuPoints(getPoissonNumb(mutPoss.back()), phyPoss, mutPoss, ExistPoss);

		if (chr == "A" || (chr == "X" && inds.at(i).Sex == 2))
		{
			inds.at(i).haplotype2.HapMuPoints(getPoissonNumb(mutPoss.back()), phyPoss, mutPoss, ExistPoss);
		}
	}
}

/************************************************
 *  set individual sample range with selection  *  
 ************************************************/
void Generation::setIndRange(const std::string &chr, const std::vector<std::string> &InitialRef, const std::map<std::string, std::string> &IndAnc, const std::vector<std::vector<int> > &SphyPos, const std::vector<std::vector<std::string> > &Sallele, const std::vector<std::string> &Saddition, const std::vector<std::vector<double> > &Sdegree, const std::map<int, std::string> &IndexInd)
{
	IndRange.clear();
	int size = static_cast<int>(inds.size());
	std::vector<double> IndexSdegree;
	IndexSdegree.resize(size);
	if (size != 0)
	{
		IndRange.resize(size+1);
	}
	else
	{
		IndRange.resize(0);
	}

	double Sum = 0;
	for (int j = 0; j < size; j++)
	{
		int TempAdd = 0;
		double Indchance = 1;

		int max = 2;
		if (chr == "X" && inds.at(j).Sex == 1)
		{
			max = 1;
		}
		for (int k = 0; k < max; k++)
		{
			for (int p = 0; p < SphyPos.size(); p++)
			{
				std::string Anc;
				bool selectAnc = 0;
				if (Sallele.at(p).size() == 1 && std::find(InitialRef.begin(), InitialRef.end(), Sallele.at(p).front()) != InitialRef.end())
				{
					Anc = Sallele.at(p).front();
					selectAnc = 1;
				}

				for (int i = 0; i < Sallele.at(p).size(); i++)
				{
					int Nsame = 0;
					if (selectAnc == 1)
					{
						for (int n = 0; n < SphyPos.at(p).size(); n++)
						{
							int tempLabel;
							inds.at(j).getLabelName(SphyPos.at(p).at(n), tempLabel, k);
							if (IndAnc.at(IndexInd.at(tempLabel/2)) == Anc)
							{
								Nsame += 1;
							}
						}
					}
					else if (selectAnc == 0)
					{
						for (int n = 0; n < SphyPos.at(p).size(); n++)
						{
							if (inds.at(j).getSallele(SphyPos.at(p).at(n), k) == Sallele.at(p).at(i).at(n))
							{
								Nsame += 1;
							}
						}
					}

					if (Nsame == SphyPos.at(p).size())
					{
						if (k == 0)
						{
							TempAdd += 1;
							if (Sdegree.at(p).size() == 1 || inds.at(j).Sex == 1)
							{
								Indchance += Sdegree.at(p).front();
							}
							else if (Sdegree.at(p).size() == 2 && inds.at(j).Sex == 2)
							{
								Indchance += Sdegree.at(p).back();
							}
						}
						else if ((TempAdd == 1 && Saddition.at(p) == "1") || (TempAdd == 0))
						{
							if (Sdegree.at(p).size() == 1 || inds.at(j).Sex == 1)
							{
								Indchance += Sdegree.at(p).front();
							}
							else if (Sdegree.at(p).size() == 2 && inds.at(j).Sex == 2)
							{
								Indchance += Sdegree.at(p).back();
							}
						}	
					} 
				}
			}
		}
		if (Indchance < 0)
		{
			Indchance = 0;
		}
		IndexSdegree[j] = Indchance;
		Sum += Indchance;			
	}

	for (int i = 0; i < size; i++)
	{
		IndexSdegree.at(i) /= Sum;
	}
		
	double start = 0;
	if (size != 0)
	{
		IndRange[0] = start;
	}
	for (int j = 0; j < size; j++)
	{
		start += IndexSdegree.at(j);
		IndRange[j+1] = start;	
	}
}

/***************************************************
 *  set Individual sample range without selection  *
 ***************************************************/
void Generation::setIndRange() 
{
	IndRange.clear();
	int size = static_cast<int>(inds.size());
	if (size != 0)
	{
		IndRange.resize(size+1);
	}
	else
	{
		IndRange.resize(0);
	}
	double start = 0;
	double range = 1/(size*1.0);
	if (size != 0)
	{
		IndRange[0] = start;
	}
	for (int i = 0; i < size; i++)
	{
		start += range;
		IndRange[i+1] = start;
	}
}

/******************************
 *  sample nsamp individuals  *
 ******************************/ 
void Generation::sample(const int &nsamp, std::vector<Individual> &samp) const
{
	samp.clear();
	samp.resize(nsamp);
	
	std::vector<int> IndIndex; 
	IndIndex.resize(nsamp, inds.size()+2);	

	for (int i = 0; i < nsamp;)
	{
		double randfloat = rand()*1.0/RAND_MAX;
		while (randfloat == 1)
		{
			randfloat = rand()*1.0/RAND_MAX;
		}
		int index;
		findFloatPos(IndRange, randfloat, index);
		index -= 1;

		if (std::find(IndIndex.begin(), IndIndex.end(), index) == IndIndex.end())
		{
			Individual ind = inds.at(index);
			samp[i] = ind;
			IndIndex[i] = index;
			i++;
		}
	}
}

/***************************
 *	sample one individual  *
 ***************************/
void Generation::sample(int &index) const
{
	double randfloat = rand()*1.0/RAND_MAX;	
	while (randfloat == 1)
	{
		randfloat = rand()*1.0/RAND_MAX;
	}
	findFloatPos(IndRange, randfloat, index);
	index -=1;
}

