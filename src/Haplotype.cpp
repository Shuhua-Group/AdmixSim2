/****************************************
 *	Project Name: AdmixSim 2  
 *	File Name: Haplotype.cpp
 *	Author: Rui Zhang
 *	Date: July 7, 2020
 ****************************************/

#include <algorithm>
#include "Haplotype.hpp"
#include "Utils.hpp"

Haplotype::Haplotype()
{
}

Haplotype::Haplotype(const std::vector<Segment> &segments, const std::vector<int> &mutationpoints, const std::map<int, char> &Sallele):
		segments(segments), mutationpoints(mutationpoints), Sallele(Sallele)
{
	/* recompoints doesn't have start and last point */
	recompoints.clear();
	if (segments.size() > 1)
	{
		recompoints.resize(segments.size()-1);
	}
	else
	{
		recompoints.resize(0);
	}
	for (size_t i = 0; i < segments.size()-1; i++)
	{
		recompoints[i] = segments.at(i).end;
	}
}

Haplotype::~Haplotype()
{
}

int Haplotype::indexOf(const int &pos, int &index) const
{
	int left = 0;
	int right = static_cast<int>(recompoints.size());
	if (recompoints.size() == 0)
	{
		index = left;
		return 0;
	}
	else
	{
		if (pos <= recompoints.front())
		{
			index = left;
			return 0;
		}	
		if (pos > recompoints.back())
		{
			index = right;
			return 0;
		}

		int mid = (left + right + 1) / 2;
		while (left < right)
		{
			if (pos > recompoints.at(mid))
			{
				left = mid;
			}
			else	
			{
				right = mid - 1;
			}
			mid = (left + right + 1) / 2;
		}
		index = left + 1;
		return 0;
	}
}

/********************************************
 *  extract segments between start and end  *
 ********************************************/
void Haplotype::extSegment(const int &start, const int &end, std::vector<Segment> &extSegs) const
{
	extSegs.clear();
	extSegs.reserve(segments.size());
	if (start != end)
	{
		int startIndex, endIndex;
		indexOf(start, startIndex);
		indexOf(end, endIndex);
		/* if the start position is on the end of a segment, then move to next */
		if (start == segments.at(startIndex).end)
		{
			startIndex++;
		}
		/* if the start and end positions are on the same segment, then cut part of it */
		if (startIndex == endIndex)
		{
			extSegs.push_back(Segment(start, end, segments.at(startIndex).label));
		}
		else
		{
			extSegs.push_back(Segment(start, segments.at(startIndex).end, segments.at(startIndex).label));
			for (int i = startIndex+1; i < endIndex; i++)
			{
				extSegs.push_back(segments.at(i));
			}
			extSegs.push_back(Segment(segments.at(endIndex).start, end, segments.at(endIndex).label));
		}
	}
	std::vector<Segment>(extSegs).swap(extSegs);
}

/*************************************
 *  merge segments of some ancestry  *
 *************************************/
void Haplotype::smooth(const std::map<std::string, std::string> &IndAnc, const std::map<int, std::string> &IndexInd)
{
	std::vector<Segment> tmpSegments;
	tmpSegments.reserve(segments.size());
	Segment seg1 = segments.front();
	for (size_t i = 0; i < segments.size(); i++)
	{
		Segment seg2 = segments.at(i);
		if (IndAnc.at(IndexInd.at(seg1.label/2)) == IndAnc.at(IndexInd.at(seg2.label/2)))
		{
			seg1.setEnd(seg2.end);
		}
		else
		{
			tmpSegments.push_back(seg1);
			seg1 = seg2;
		}
	}	
	tmpSegments.push_back(seg1);
	segments = tmpSegments;
	updateBreaks();
}

/****************************************
 *  update recombination points vector  *
 ****************************************/
void Haplotype::updateBreaks()
{
	recompoints.clear();
	if (segments.size() > 1)
	{
		recompoints.resize(segments.size()-1);
	}
	else
	{
		recompoints.resize(0);
	}
	for (size_t i = 0; i < segments.size()-1; i++)
	{
		recompoints[i] = segments.at(i).end;
	}
}

/*********************************************************
 *	return segment label of specified physical position  *
 *********************************************************/
void Haplotype::getLabelName(const int &pos, int &Label) const
{
	int index;
	indexOf(pos, index);
	if (index < recompoints.size() && pos == recompoints.at(index))
	{
		index += 1;
	}
	Label = segments.at(index).label;
}

/***************************************************
 *  extract mutation points between start and end  *
 ***************************************************/
void Haplotype::getMutations(const int &start, const int &end, std::vector<int> &RangeMutations) const
{
	RangeMutations.clear();
	
	if (mutationpoints.size() != 0)
	{
		int startIndex,endIndex;
		findIntPos(mutationpoints, start, startIndex);
		findIntPos(mutationpoints, end, endIndex);
	
		RangeMutations.resize(endIndex-startIndex);

		for (int i = startIndex; i < endIndex; i++)
		{
			RangeMutations[i-startIndex] = mutationpoints.at(i);
		}
	}
}

/**************************************************************
 *  extract haplotype real allele state between start and end  *
 **************************************************************/
void Haplotype::getSallele(const int &start, const int &end, const std::vector<int> &AllSpos, const std::vector<int> &phyPoss, std::map<int, char> &RangeSallele) const
{
	RangeSallele.clear();

	if (AllSpos.size() != 0)
	{ 
		int startIndex, endIndex;
		findIntPos(AllSpos, start, startIndex);
		findIntPos(AllSpos, end, endIndex);

		if (AllSpos.size() == 1 && end == phyPoss.back() && AllSpos.front() == end)
		{
			endIndex += 1;
		}
	
		for (int i = startIndex; i < endIndex; i++)
		{
			RangeSallele[AllSpos.at(i)] = Sallele.at(AllSpos.at(i));
		}
	}
}

/**********************************
 *	return real allele situation  *
 **********************************/
void Haplotype::getSallele(const int &pos, char &Allele) const
{
	Allele = Sallele.at(pos);
}

/*******************************
 *  set real allele situation  *
 *******************************/
void Haplotype::setSallele(const std::string &hapseq, const std::map<int, int> &phyPos_index, const std::vector<int> &AllSpos)
{
	for (int i = 0; i < AllSpos.size(); i++)
	{
		Sallele[AllSpos.at(i)] = hapseq.at(phyPos_index.at(AllSpos.at(i)));	
	}
}

/******************************************
 *	add de novo mutations to a haplotype  *
 ******************************************/
void Haplotype::HapMuPoints(const int &n, const std::vector<int> &phyPoss, const std::vector<double> &mutPoss, std::set<int> &ExistPoss)
{
	mutationpoints.reserve(mutationpoints.size()+n);
	for (int i = 0; i < n;)
	{
		double mutRandom = rand()*mutPoss.back()/RAND_MAX;
		int newMut;
		MapToPhyPos(mutRandom, mutPoss, phyPoss, newMut);

		if (ExistPoss.find(newMut) == ExistPoss.end() && std::binary_search(phyPoss.begin(), phyPoss.end(), newMut) == false)
		{
			mutationpoints.insert(lower_bound(mutationpoints.begin(), mutationpoints.end(), newMut), newMut);
			ExistPoss.insert(newMut);
			i++;
		}
	}
}


