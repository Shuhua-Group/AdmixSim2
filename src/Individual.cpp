/****************************************
 *	Project Name: AdmixSim 2
 *	File Name: Individual.cpp
 *	Author: Rui Zhang
 *	Date: July 7, 2020
 ****************************************/

#include <cstdlib>
#include <algorithm>
#include <vector>
#include <cmath>
#include "Individual.hpp"
#include "Utils.hpp"

Individual::Individual()
{
}

Individual::Individual(const Haplotype &haplotype1, const Haplotype &haplotype2, const int &Sex):
	haplotype1(haplotype1), haplotype2(haplotype2), Sex(Sex)
{
}

Individual::~Individual()
{
}

int getPoissonNumb(const double &lambda)
{
	double length = exp(-lambda);
	double prob = 1.0;
	int numb = 0;
	do
	{
		numb++;
		prob *= rand()*1.0/RAND_MAX;	
	}while (prob > length);
	return numb - 1;
}

/******************************
 *  set recombination points  *
 ******************************/
void Individual::breakPoints(int *recompoints, const int &n, const std::vector<int> &phyPoss, const std::vector<double> &genPoss) const
{
	recompoints[0] = phyPoss.front();
	for (int i = 1; i < n-1; i++)
	{
		double recRandom = rand()*genPoss.back()/RAND_MAX;
		MapToPhyPos(recRandom, genPoss, phyPoss, recompoints[i]);
	}
	recompoints[n-1] = phyPoss.back();

	std::sort(recompoints, recompoints+n);
}

/******************************************
 *  get haplotype according to the index  *
 ******************************************/
void Individual::getHap(const int &index, Haplotype &hap) const
{
	if (index == 0)
	{
		hap = haplotype1;
	}
	else
	{
		hap = haplotype2;
	}
}

/****************************************************
 * 	get haplotype mutations according to the index  *
 ****************************************************/
void Individual::getHapMut(const int &index, std::vector<int> &Mutations) const
{
	if (index == 0)
	{
		Mutations = haplotype1.mutationpoints;
	}
	else
	{
		Mutations = haplotype2.mutationpoints;
	}
}

/*********************************************************
 * return segment label of specified physical position  *
 *********************************************************/
void Individual::getLabelName(const int &pos, int &Label, const int &index) const
{
	if (index == 0)
	{
		haplotype1.getLabelName(pos, Label);
	}
	else if (index == 1)
	{
		haplotype2.getLabelName(pos, Label);
	}
}

/**********************************
 *  return real allele situation  *
 **********************************/
char Individual::getSallele(const int &pos, const int &index) const
{
	char Allele;
	if (index == 0)
	{
		haplotype1.getSallele(pos, Allele);
	}
	else if (index == 1)
	{	
		haplotype2.getSallele(pos, Allele);
	}
	return Allele;
}

/*********************************
 *  return a haplotype randomly  *
 *********************************/
void Individual::recombine(const std::vector<int> &phyPoss, const std::vector<double> &genPoss, const std::vector<int> &AllSpos, Haplotype &hap)
{
	hap.segments.clear();
	hap.mutationpoints.clear();
	hap.Sallele.clear();

	hap.segments.reserve(haplotype1.segments.size()+haplotype2.segments.size());
	hap.mutationpoints.reserve(haplotype1.mutationpoints.size()+haplotype2.mutationpoints.size());

	int hapIndex = rand()%2;

	int numb = getPoissonNumb(genPoss.back())+2;
	int bps[numb];

	breakPoints(bps, numb, phyPoss, genPoss);

	if (hapIndex == 0)
	{
		int start = bps[0];
		for (int i = 1; i < numb; i++)
		{
			int end = bps[i];
			if (i % 2)		/* if i is not the multiple of 2 */
			{
				std::vector<Segment> tmps;
				haplotype1.extSegment(start, end, tmps);
				for (size_t j = 0; j < tmps.size(); j++)
				{
					hap.segments.push_back(tmps.at(j));
				}
				std::vector<int> tmpsmu;
				haplotype1.getMutations(start, end, tmpsmu);
				for (size_t k = 0; k < tmpsmu.size(); k++)
				{
					hap.mutationpoints.push_back(tmpsmu.at(k));
				}
				std::map<int, char> tmpS;
				haplotype1.getSallele(start, end, AllSpos, phyPoss, tmpS);
				for (std::map<int, char>::iterator iter = tmpS.begin(); iter != tmpS.end(); iter++)
				{
					hap.Sallele[iter->first] = iter->second;
				}
			}
			else			/* if i is the multiple of 2 */
			{
				std::vector<Segment> tmps;
				haplotype2.extSegment(start, end, tmps);
				for (size_t j = 0; j < tmps.size(); j++)
				{
					hap.segments.push_back(tmps.at(j));
				}			
				std::vector<int> tmpsmu;
				haplotype2.getMutations(start, end, tmpsmu);
				for (size_t k = 0; k < tmpsmu.size(); k++)
				{
					hap.mutationpoints.push_back(tmpsmu.at(k));
				}
				std::map<int, char> tmpS;
				haplotype2.getSallele(start, end, AllSpos, phyPoss, tmpS); 
				for (std::map<int, char>::iterator iter = tmpS.begin(); iter != tmpS.end(); iter++)
				{
					hap.Sallele[iter->first] = iter->second;
				}
			}
			start = end;
		}
	}
	else
	{
		int start = bps[0];
		for (int i = 1; i < numb; i++)
		{
			int end = bps[i];
			if (i % 2)		/* if i is not the multiple of 2 */
			{
				std::vector<Segment> tmps;
				haplotype2.extSegment(start, end, tmps);
				for (size_t j = 0; j < tmps.size(); j++)
				{
					hap.segments.push_back(tmps.at(j));
				}
				std::vector<int> tmpsmu;
				haplotype2.getMutations(start, end, tmpsmu);
				for (size_t k = 0; k < tmpsmu.size(); k++)
				{
					hap.mutationpoints.push_back(tmpsmu.at(k));
				}
				std::map<int, char> tmpS;
				haplotype2.getSallele(start, end, AllSpos, phyPoss, tmpS);
				for (std::map<int, char>::iterator iter = tmpS.begin(); iter != tmpS.end(); iter++)
				{
					hap.Sallele[iter->first] = iter->second;
				}
			}
			else			/* if i is the multiple of 2 */
			{
				std::vector<Segment> tmps;
				haplotype1.extSegment(start, end, tmps);
				for (size_t j = 0; j < tmps.size(); j++)
				{
					hap.segments.push_back(tmps.at(j));
				}
				std::vector<int> tmpsmu;
				haplotype1.getMutations(start, end, tmpsmu);
				for (size_t k = 0; k < tmpsmu.size(); k++)
				{
					hap.mutationpoints.push_back(tmpsmu.at(k));
				}
				std::map<int, char> tmpS;
				haplotype1.getSallele(start, end, AllSpos, phyPoss, tmpS);
				for (std::map<int, char>::iterator iter = tmpS.begin(); iter != tmpS.end(); iter++)
				{
					hap.Sallele[iter->first] = iter->second;
				}
			}
			start = end;
		}
	}

	std::vector<Segment>(hap.segments).swap(hap.segments);
	std::vector<int>(hap.mutationpoints).swap(hap.mutationpoints);	
	hap.updateBreaks();
}


