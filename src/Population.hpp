/****************************************
 *	Project Name: AdmixSim 2
 *	File Name: Population.hpp 
 *	Author: Rui Zhang
 *	Date: July 7, 2020
 ****************************************/

#ifndef POPULATION_HPP_
#define POPULATION_HPP_

#include "Generation.hpp"

/*************************************************************************
 * Each population consists of individuals spreading several generations *
 * specified in the demographic model									 *
 *************************************************************************/

class Population
{
public:
	std::map<int, Generation> Gens;

	std::vector<int> GIndices;
	std::map<int, std::vector<int> > GNe;									/* Generation index and corresponding Ne */
	std::map<int, std::map<std::string, std::vector<double> > > GAncP; 		/* Generation index and corresponding ancestry and proportions */	


	std::map<int, std::vector<std::vector<int> > > GSphyPos;					/* Generation index and selected physical position(s) allow replication */
	std::map<int, std::vector<std::string> > GUserSpos;						/* Generation index and user defined selected physical position */
	std::map<int, std::vector<std::vector<std::string> > > GphyPos_allele;	/* Generation Index and selected allele //allow replication */
	std::map<int, std::vector<std::string> > GphyPos_Addition;				/* Generation Index and addition rule (one+1,two+2/one,two+1) */
	std::map<int, std::vector<std::vector<double> > > GSdegree;				/* Generation Index and selection degree */

	std::map<int, std::vector<int> > GRealNe;								/* Generation index and real male/female size */

	Population();
	virtual ~Population();

	void setIndex(const std::vector<int> &Indices);
	void addNe(const int &Index, const std::vector<int> &Ne);
	void addRealNe(const int &Index, const std::vector<int> &RealNe);
	void addAnc(const int &Index, const std::map<std::string, std::vector<double> > &AncP);
	
	void addSel(const int &Index, const std::vector<int> &sphyPos, const std::string &userSpos, const std::vector<std::string> &allele, const std::string &Addition, const std::vector<double> &Degree);

	void complement(const std::string &Name);
	void GenerateInd(Individual &ind1, Individual &ind2, const std::string &chr, const std::string &recEvent, const std::vector<int> &phyPoss, const std::vector<double> &genPoss, const bool &AncSex, const std::vector<int> &AllSpos, Individual &admInd);
	void evolve(const std::string &name, const int &Index, std::map<std::string, Population> &AllPop, const std::string &chr, const std::string &recEvent, const std::vector<int> &phyPoss, const std::vector<double> &genPoss, const std::vector<std::string> &InitialRef, const std::map<std::string, std::string> &IndAnc, const bool &AncSex, std::vector<std::string> &alleleString, const std::vector<int> &AllSpos, const std::map<int, std::string> &IndexInd);
};

#endif /* POPULATION_HPP_ */

