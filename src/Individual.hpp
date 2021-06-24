/****************************************
 *	Project Name: AdmixSim 2
 *	File Name: Individual.hpp
 * 	Rui Zhang
 *	Date: July 7, 2020
 ****************************************/

#ifndef INDIVIDUAL_HPP_
#define INDIVIDUAL_HPP_

#include "Haplotype.hpp"

/**********************************************
 * Each Individual consists of two haplotype  *
 * and corresponding sex information          *
 **********************************************/

class Individual
{
private:
	void breakPoints(int *recompoints, const int &n, const std::vector<int> &phyPoss, const std::vector<double> &genPoss) const;
	
public:
	Haplotype haplotype1;
	Haplotype haplotype2;
	int Sex;
	Individual();
	Individual(const Haplotype &haplotype1, const Haplotype &haplotype2, const int &Sex);
	virtual ~Individual();
	void getHap(const int &index, Haplotype &hap) const;
	void getHapMut(const int &index, std::vector<int> &Mutations) const;
	void getLabelName(const int &pos, int &Label, const int &index) const;
	char getSallele(const int &pos, const int &index) const;
	void setSex(const int &sex);
	void recombine(const std::vector<int> &phyPoss, const std::vector<double> &genPoss, const std::vector<int> &AllSpos, Haplotype &hap);
};

int getPoissonNumb(const double &lambda);

#endif /* INDIVIDUAL_HPP_ */

