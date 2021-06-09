/****************************************
 *	Project: AdmixSim 2
 *	File Name: Generation.hpp
 *	Author: Rui Zhang
 *	Date: July 7, 2020
 ****************************************/

#ifndef GENERATION_HPP_
#define GENERATION_HPP_

#include "Individual.hpp"

class Generation
{
private:
	std::vector<double> IndRange;

public:
	std::vector<Individual> inds;
	Generation();
	virtual ~Generation();

	void GenMuPoints(const std::vector<int> &phyPoss, const std::vector<double> &mutPoss, const std::string &chr, std::set<int> &ExistPoss); 
	void setIndRange(const std::string &chr, const std::vector<std::string> &InitialRef, const std::map<std::string, std::string> &IndAnc, const std::vector<std::vector<int> > &SphyPos, const std::vector<std::vector<std::string> > &Sallele, const std::vector<std::string> &Saddition, const std::vector<std::vector<double> > &Sdegree, const std::map<int, std::string> &IndexInd);
	void setIndRange();

	void sample(const int &nsamp, std::vector<Individual> &samp) const;
	void sample(int &index) const;
};

#endif /* GENERATION_HPP_ */

