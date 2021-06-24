/****************************************
 *	Project Name: AdmixSim 2
 *	File Name: Haplotype.hpp
 *	Author: Rui Zhang
 *	Date: July 7, 2020
 ****************************************/

#ifndef HAPLOTYPE_HPP_
#define HAPLOTYPE_HPP_

#include <vector>
#include <map>
#include <set>
#include "Segment.hpp"


/**************************************************************
 * Each haplotype consists of a list of segments,			  *
 * a list recording the recombination points,				  *
 * a list recoding the de novo mutation points,				  *
 * and a map recording the physical positions under selection *
 * and corresponding states of allele						  *
 **************************************************************/

class Haplotype
{
public:
	std::vector<Segment> segments;
	std::vector<int> recompoints;       /* int is the mapped physical position */
	std::vector<int> mutationpoints;    /* int is the mapped physical position */
	std::map<int, char> Sallele;        /* int is the selected physical position; char is the real allele */
	Haplotype();
	Haplotype(const std::vector<Segment> &segments, const std::vector<int> &mutationpoints, const std::map<int, char> &Sallele);
	virtual ~Haplotype();
	int indexOf(const int &pos, int &index) const;
	void extSegment(const int &start, const int &end, std::vector<Segment> &extSegs) const;
	void smooth(const std::map<std::string, std::string> &IndAnc, const std::map<int, std::string> &IndexInd);
	void updateBreaks();
	void getLabelName(const int &pos, int &Label) const;

	void getMutations(const int &start, const int &end, std::vector<int> &RangeMutations) const;

	void getSallele(const int &start, const int &end, const std::vector<int> &AllSpos, const std::vector<int> &phyPoss, std::map<int, char> &RangeSallele) const;
	void getSallele(const int &pos, char &Allele) const;
	void setSallele(const std::string &hapseq, const std::map<int, int> &phyPos_index, const std::vector<int> &AllSpos);
	void HapMuPoints(const int &n, const std::vector<int> &phyPoss, const std::vector<double> &mutPoss, std::set<int> &ExistPoss);		
};

#endif /* HAPLOTYPE_HPP_ */
