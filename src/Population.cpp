/****************************************
 *	Project Name: AdmixSim 2
 *	File Name: Population.cpp
 *	Author: Rui Zhang
 *	Date: July 7, 2020
 ****************************************/

#include <iostream>
#include <algorithm>
#include <numeric>
#include "Population.hpp"
#include "Utils.hpp"

Population::Population()
{
}

Population::~Population()
{
}

void Population::setIndex(const std::vector<int> &Indices)
{
	GIndices.clear();
	GIndices.resize(Indices.back()-Indices.front()+1);
	for (int i = Indices.front(); i <= Indices.back(); i++)
	{
		GIndices[i-Indices.front()] = i;
	}
}

void Population::addNe(const int &Index, const std::vector<int> &Ne)
{
	GNe[Index] = Ne;
}

void Population::addRealNe(const int &Index, const std::vector<int> &RealNe)
{
	GRealNe[Index] = RealNe;
}

void Population::addAnc(const int &Index, const std::map<std::string, std::vector<double> > &AncP)
{
	GAncP[Index] = AncP;
}

void Population::addSel(const int &Index, const std::vector<int> &sphyPos, const std::string &userSpos, const std::vector<std::string> &allele, const std::string &Addition, const std::vector<double> &Degree)
{
	if (GSphyPos.find(Index) == GSphyPos.end())
	{
		std::vector<std::vector<int> > temp1;
		GSphyPos[Index] = temp1;

		std::vector<std::string> temp2;
		GUserSpos[Index] = temp2;	

		std::vector<std::vector<std::string> > temp3;
		GphyPos_allele[Index] = temp3;

		std::vector<std::string> temp4;
		GphyPos_Addition[Index] = temp4;

		std::vector<std::vector<double> > temp5;
		GSdegree[Index] = temp5;
	}

	GSphyPos.at(Index).push_back(sphyPos);
	GUserSpos.at(Index).push_back(userSpos);
	GphyPos_allele.at(Index).push_back(allele);
	GphyPos_Addition.at(Index).push_back(Addition);
	GSdegree.at(Index).push_back(Degree);
}

/********************************************************************************
 *  add unspecific Ne/admixture proportion/selection same with last generation  *
 ********************************************************************************/
void Population::complement(const std::string &Name)
{
	for (int i = GIndices.front(); i <= GIndices.back(); i++)
	{
		if (GNe.find(i) == GNe.end())
		{
			GNe[i] = GNe.at(i-1);
			GAncP[i] = GAncP.at(i-1);

			if (GSphyPos.find(i-1) != GSphyPos.end())
			{
				GSphyPos[i] = GSphyPos.at(i-1);
				GUserSpos[i] = GUserSpos.at(i-1);
				GphyPos_allele[i] = GphyPos_allele.at(i-1);
				GphyPos_Addition[i] = GphyPos_Addition.at(i-1);
				GSdegree[i] = GSdegree.at(i-1);
			}
		}	
	}
}

/***************************************
 *  Generate offspring of two parents  *
 ***************************************/
void Population::GenerateInd(Individual &ind1, Individual &ind2, const std::string &chr, const std::string &recEvent, const std::vector<int> &phyPoss, const std::vector<double> &genPoss, const bool &AncSex, const std::vector<int> &AllSpos, Individual &admInd)
{
	int sex1 = ind1.Sex;
	int sex2 = ind2.Sex;
	
	int GeneratedSex = 0;

	if (chr == "A")
	{
		if (recEvent == "0")
		{
			int hapIndex1 = rand()%2;
			if (hapIndex1 == 0)
			{
				admInd.haplotype1 = ind1.haplotype1;
			}
			else
			{
				admInd.haplotype1 = ind1.haplotype2;
			}

			int hapIndex2 = rand()%2;
			if (hapIndex2 == 0)
			{
				admInd.haplotype2 = ind2.haplotype1;
			}
			else
			{
				admInd.haplotype2 = ind2.haplotype2;	
			}
		}
		else
		{
			ind1.recombine(phyPoss, genPoss, AllSpos, admInd.haplotype1);
			ind2.recombine(phyPoss, genPoss, AllSpos, admInd.haplotype2);
		}
		
		if (AncSex == 1)
		{
			int temp = rand()%2;
			if (temp == 0)
			{
				GeneratedSex = 1;
			}
			else if (temp == 1)
			{
				GeneratedSex = 2;
			}
		}
		admInd.Sex = GeneratedSex;
	}
	else if (chr == "X")
	{
		int temp = rand()%2;
		if (temp == 0)
		{
			GeneratedSex = 1;
		}
		else if (temp == 1)
		{
			GeneratedSex = 2;
		}

		admInd.Sex = GeneratedSex;
		
		if (GeneratedSex == 1)		/* generate a boy */
		{
			if (sex1 == 1)
			{
				if (recEvent == "0")
				{
					int hapIndex = rand()%2;
					if (hapIndex == 0)
					{
						admInd.haplotype1 = ind2.haplotype1;
					}
					else
					{
						admInd.haplotype1 = ind2.haplotype2;
					}
				}
				else
				{
					ind2.recombine(phyPoss, genPoss, AllSpos, admInd.haplotype1);
				}
				admInd.haplotype2 = ind1.haplotype2;
			}
			else if (sex2 == 1)
			{
				if (recEvent == "0")
				{
					int hapIndex = rand()%2;
					if (hapIndex == 0)
					{
						admInd.haplotype1 = ind1.haplotype1;
					}
					else
					{
						admInd.haplotype1 = ind1.haplotype2;
					}
				}
				else
				{
					ind1.recombine(phyPoss, genPoss, AllSpos, admInd.haplotype1);
				}
				admInd.haplotype2 = ind2.haplotype2;
			}
		}
		else if (GeneratedSex == 2)	/* generate a girl */
		{
			if (sex1 == 1)
			{
				admInd.haplotype1 = ind1.haplotype1;
				if (recEvent == "0")
				{
					int hapIndex = rand()%2;
					if (hapIndex == 0)
					{
						admInd.haplotype2 = ind2.haplotype1;
					}
					else
					{
						admInd.haplotype2 = ind2.haplotype2;
					}
				}
				else
				{
					ind2.recombine(phyPoss, genPoss, AllSpos, admInd.haplotype2);
				}
			}
			else if (sex2 == 1)
			{
				if (recEvent == "0")
				{
					int hapIndex = rand()%2;
					if (hapIndex == 0)
					{
						admInd.haplotype1 = ind1.haplotype1;
					}
					else
					{
						admInd.haplotype1 = ind1.haplotype2;
					}
				}
				else
				{
					ind1.recombine(phyPoss, genPoss, AllSpos, admInd.haplotype1);
				}
				admInd.haplotype2 = ind2.haplotype1;
			}
		}
	}
}

void Population::evolve(const std::string &name, const int &Index, std::map<std::string, Population> &AllPop, const std::string &chr, const std::string &recEvent, const std::vector<int> &phyPoss, const std::vector<double> &genPoss, const std::vector<std::string> &InitialRef, const std::map<std::string, std::string> &IndAnc, const bool &AncSex, std::vector<std::string> &alleleString, const std::vector<int> &AllSpos, const std::map<int, std::string> &IndexInd)
{
	/* Set sampling range of each population according to model file */
	double start = 0;
	double middle;
	double end;
	std::map<std::string, std::vector<double> > range;	/* ancestry name || start,middle,end / start,end */
	int TotalAncMale = 0;
	int TotalAncFemale = 0;
	int TotalNum = 0;

	bool SexP = 0;
	if (AncSex == 1)
	{
		for (std::map<std::string, std::vector<double> >::iterator iter = GAncP.at(Index).begin(); iter != GAncP.at(Index).end(); iter++)
		{
			if ((iter->second).size() == 2)
			{
				SexP = 1;
				break;	
			}
			else
			{
				break;
			}
		}
	}

	/* male: 1/2*mi*pi/sum(mi*pi)   female: 1/2*mi*(1-pi)/sum(mi*(1-pi))   mi for admixture proportion  pi for male proportion */
	if (AncSex == 1 && SexP == 0)
	{
		double maleSum = 0;
		double femaleSum = 0;
		for (std::map<std::string, std::vector<double> >::iterator iter = GAncP.at(Index).begin(); iter != GAncP.at(Index).end(); iter++)
		{
			double tempsum = 0;
			for (int i = 0; i < (iter->second).size(); i++)
			{
				tempsum += (iter->second).at(i);
			}
			if (tempsum != 0)
			{
				int maleN = AllPop.at(iter->first).GRealNe.at(Index-1).front();
				int femaleN = AllPop.at(iter->first).GRealNe.at(Index-1).back();
				int sum = maleN + femaleN;
				maleSum += (iter->second).front() * maleN * 1.0 / sum;  
				femaleSum += (iter->second).front()* femaleN * 1.0 / sum;
			}
		}
		
		for (std::map<std::string, std::vector<double> >::iterator iter = GAncP.at(Index).begin(); iter != GAncP.at(Index).end(); iter++)
		{
			double tempsum = 0;
			for (int i = 0; i < (iter->second).size(); i++)
			{
				tempsum += (iter->second).at(i);
			}
			if (tempsum != 0)
			{
				double maleProp = 0;
				double femaleProp = 0;
				int maleN = AllPop.at(iter->first).GRealNe.at(Index-1).front();
				int femaleN = AllPop.at(iter->first).GRealNe.at(Index-1).back();
				int sum = maleN + femaleN;
				maleProp = 0.5 * (iter->second).front() * (maleN * 1.0 / sum) / maleSum;
				femaleProp = 0.5 * (iter->second).front() * (femaleN * 1.0 / sum) / femaleSum;
				std::vector<double> tempProp;
				tempProp.resize(2);
				tempProp[0] = maleProp;
				tempProp[1] = femaleProp;
				GAncP.at(Index).at(iter->first) = tempProp;
			}
		}
	}

	for (std::map<std::string, std::vector<double> >::iterator iter = GAncP.at(Index).begin(); iter != GAncP.at(Index).end(); iter++)
	{
		std::vector<double> start_end;
		if ((iter->second).size() == 2)
		{
			start_end.resize(3);
		}
		else
		{
			start_end.resize(2);
		}
		start_end[0] = start;

		if ((iter->second).size() == 2)
		{
			if ((iter->second).front() != 0 && AllPop.at(iter->first).GRealNe.at(Index-1).front() == 0)
			{
				std::cerr << "Error: Individuals of " << (iter->first) << " at generation " << (Index-1) << " are all female and cannot provide male for " << name << " at generation " << Index << "!" << std::endl;
				exit(1);
			}
			if ((iter->second).back() != 0 && AllPop.at(iter->first).GRealNe.at(Index-1).back() == 0)
			{
				std::cerr << "Error: Individuals of " << (iter->first) << " at generation " << (Index-1) << " are all male and cannot provide female for " << name << " at generation " << Index << "!" << std::endl;
				exit(1);
			}
			
			middle = start + (iter->second).front();
			start_end[1] = middle;
			start = middle;
		}
		end = start + (iter->second).back();
		start_end[start_end.size()-1] = end;

		range[iter->first] = start_end;
		start = end;

		if ((iter->second).front() != 0 || (iter->second).back() != 0)
		{
			TotalAncMale += AllPop.at(iter->first).GRealNe.at(Index-1).front();
			TotalAncFemale += AllPop.at(iter->first).GRealNe.at(Index-1).back();	
			TotalNum += accumulate(AllPop.at(iter->first).GNe.at(Index-1).begin(), AllPop.at(iter->first).GNe.at(Index-1).end(), 0);
		}

		if ((iter->second).front() != 0 && (iter->second).back() != 0 && accumulate(AllPop.at(iter->first).GNe.at(Index-1).begin(), AllPop.at(iter->first).GNe.at(Index-1).end(), 0) == 0)
		{
			std::cerr << "Error: The population size of " << (iter->first) << " at generation " << (Index-1) << " is zero and cannot provide individual for " << name << " at generation " << Index << "!" << std::endl;
			exit(1);
		}	
	}

	if (TotalAncMale == TotalNum)
	{
		std::cerr << "Error: The ancestry populations of " << name << " at generation " << Index << " are all male!" << std::endl;
		exit(1);
	}
	if (TotalAncFemale == TotalNum)
	{
		std::cerr << "Error: The ancestry populations of " << name << " at generation " << Index << " are all female!" << std::endl;
		exit(1);
	}

	Generation admixedGen;
	int GenMaleNum = 0;
	int GenFemaleNum = 0;
	
	int IndNum = accumulate(GNe.at(Index).begin(), GNe.at(Index).end(), 0);
	admixedGen.inds.resize(IndNum);

	int SetMaleNum = -1;
	int SetFemaleNum = -1;
	
	if (GNe.at(Index).size() == 2)
	{
		SetMaleNum = GNe.at(Index).front();
		SetFemaleNum = GNe.at(Index).back();
	}

	/* produce admixed population */
	for (int k = 0; k < IndNum; k++)
	{
	//	std::cout << "Ind " << k << std::endl;	//
		std::vector<std::string> PopLab;
		std::vector<int> ParentSex;
		for (int m = 0; m < 2; m++)
		{
			double rand1 = rand()*1.0/RAND_MAX;
			while (rand1 == 1)
			{
				rand1 = rand()*1.0/RAND_MAX;
			}
			std::string poplab;
			for (std::map<std::string, std::vector<double> >::iterator iter = range.begin(); iter != range.end(); iter++)
			{
				if (rand1 >= (iter->second).front() && rand1 < (iter->second).back())
				{
					poplab = iter->first;
				}
			}
			PopLab.push_back(poplab);
		
			if (range.at(poplab).size() == 3)
			{
				if (rand1 <= range.at(poplab).at(1))
				{
					ParentSex.push_back(1);
				}
				else if (rand1 > range.at(poplab).at(1))
				{
					ParentSex.push_back(2);
				}
			}
		}

		while (ParentSex.size() == 2 && ParentSex.front() == ParentSex.back())
		{
			PopLab.erase(PopLab.end()-1);
			ParentSex.erase(ParentSex.end()-1);
			double rand1 = rand()*1.0/RAND_MAX;
			while (rand1 == 1)
			{
				rand1 = rand()*1.0/RAND_MAX;
			}
			std::string poplab;
			for (std::map<std::string, std::vector<double> >::iterator iter = range.begin(); iter != range.end(); iter++)
			{
				if (rand1 >= (iter->second).front() && rand1 < (iter->second).back())
				{
					poplab = iter->first;
				}
			}
			PopLab.push_back(poplab);
			if (rand1 <= range.at(poplab).at(1))
			{
				ParentSex.push_back(1);
			}
			else if (rand1 > range.at(poplab).at(1))
			{
				ParentSex.push_back(2);
			}
		}

		int ind1Index, ind2Index;
		AllPop.at(PopLab.at(0)).Gens.at(Index-1).sample(ind1Index);
		AllPop.at(PopLab.at(1)).Gens.at(Index-1).sample(ind2Index);

		int sex1 = AllPop.at(PopLab.at(0)).Gens.at(Index-1).inds.at(ind1Index).Sex;
		int sex2 = AllPop.at(PopLab.at(1)).Gens.at(Index-1).inds.at(ind2Index).Sex;

		//ancestry no sex setting
		while (PopLab.at(0) == PopLab.at(1) && ind1Index == ind2Index)
		{
			AllPop.at(PopLab.at(0)).Gens.at(Index-1).sample(ind1Index);
			AllPop.at(PopLab.at(1)).Gens.at(Index-1).sample(ind2Index);
			sex1 = AllPop.at(PopLab.at(0)).Gens.at(Index-1).inds.at(ind1Index).Sex;
			sex2 = AllPop.at(PopLab.at(1)).Gens.at(Index-1).inds.at(ind2Index).Sex;	
		}

		while (ParentSex.size() == 2 && sex1 != ParentSex.front())
		{
			AllPop.at(PopLab.at(0)).Gens.at(Index-1).sample(ind1Index);
			sex1 = AllPop.at(PopLab.at(0)).Gens.at(Index-1).inds.at(ind1Index).Sex;
		}

		while (ParentSex.size() == 2 && sex2 != ParentSex.back())
		{
			AllPop.at(PopLab.at(1)).Gens.at(Index-1).sample(ind2Index);
			sex2 = AllPop.at(PopLab.at(1)).Gens.at(Index-1).inds.at(ind2Index).Sex;
		}

	//	std::cout << "poplab1 = " << PopLab.at(0) << " indIndex = " << ind1Index << std::endl;	//
	//	std::cout << "poplab2 = " << PopLab.at(1) << " indIndex = " << ind2Index << std::endl;	//

		Individual ind1 = AllPop.at(PopLab.at(0)).Gens.at(Index-1).inds.at(ind1Index);
		Individual ind2 = AllPop.at(PopLab.at(1)).Gens.at(Index-1).inds.at(ind2Index);
		Individual admixedInd;

		GenerateInd(ind1, ind2, chr, recEvent, phyPoss, genPoss, AncSex, AllSpos, admixedInd);

		if (admixedInd.Sex == 1)
		{
			GenMaleNum += 1;
		}
		else if (admixedInd.Sex == 2)
		{
			GenFemaleNum += 1;
		}
		
		while (SetMaleNum != -1 && (GenMaleNum > SetMaleNum || GenFemaleNum > SetFemaleNum))
		{ 
			if (GenMaleNum > SetMaleNum)
			{
				GenMaleNum -= 1;
			}
			else if (GenFemaleNum > SetFemaleNum)
			{
				GenFemaleNum -= 1;
			}
			
			GenerateInd(ind1, ind2, chr, recEvent, phyPoss, genPoss, AncSex, AllSpos, admixedInd);	
			if (admixedInd.Sex == 1)
			{
				GenMaleNum += 1;
			}
			else if (admixedInd.Sex == 2)
			{
				GenFemaleNum += 1;
			}
		}

		/* no sex-specified Ne setting     avoid the situation of all male or all female */		
		while (SetMaleNum == -1 && (GenMaleNum == IndNum || GenFemaleNum == IndNum))
		{
			if (GenMaleNum == IndNum)
			{
				GenMaleNum -= 1;
			}
			else if (GenFemaleNum == IndNum)
			{
				GenFemaleNum -= 1;
			}

			GenerateInd(ind1, ind2, chr, recEvent, phyPoss, genPoss, AncSex, AllSpos, admixedInd);
			if (admixedInd.Sex == 1)
			{
				GenMaleNum += 1;
			}
			else if (admixedInd.Sex == 2)
			{
				GenFemaleNum += 1;
			}
		}
		admixedGen.inds[k] = admixedInd;
	}			

	/* set IndRange for new generated Gen  */
	if (GSphyPos.size() != 0 && GSphyPos.find(Index) != GSphyPos.end())
	{
		admixedGen.setIndRange(chr, InitialRef, IndAnc, GSphyPos.at(Index), GphyPos_allele.at(Index), GphyPos_Addition.at(Index), GSdegree.at(Index), IndexInd);

		/* calculate allele frequency */
		for (int i = 0; i < GSphyPos.at(Index).size(); i++)
		{
			bool selectAnc = 0;
			std::string Anc;
			if (GphyPos_allele.at(Index).at(i).size() == 1 && std::find(InitialRef.begin(), InitialRef.end(), GphyPos_allele.at(Index).at(i).front()) != InitialRef.end())
			{
				Anc = GphyPos_allele.at(Index).at(i).front();
				selectAnc = 1;
			}

			for (int m = 0; m < GphyPos_allele.at(Index).at(i).size(); m++)
			{
				std::string tempaf;
				tempaf += GUserSpos.at(Index).at(i);
				tempaf += "\t";

				double frequency = 0;
				double male_f = 0;
				double female_f = 0;
				double TotalHap = 0;
				for (int j = 0; j < IndNum; j++)
				{
					int nhap = 2;
					if (chr == "X" && admixedGen.inds.at(j).Sex == 1)
					{
						nhap = 1;	
					}
					TotalHap += nhap;
					for (int k = 0; k < nhap; k++)
					{
						int Nsame = 0;
						if (selectAnc == 1)
						{
							for (int n = 0; n < GSphyPos.at(Index).at(i).size(); n++)
							{
								int tempLabel;
								admixedGen.inds.at(j).getLabelName(GSphyPos.at(Index).at(i).at(n), tempLabel, k);
								if (IndAnc.at(IndexInd.at(tempLabel/2)) == Anc)
								{
									Nsame += 1;
								}
							}
						}
						else if (selectAnc == 0)
						{
							for (int n = 0; n < GphyPos_allele.at(Index).at(i).at(m).size(); n++)
							{
								if (admixedGen.inds.at(j).getSallele(GSphyPos.at(Index).at(i).at(n), k) == GphyPos_allele.at(Index).at(i).at(m).at(n))
								{
									Nsame += 1;
								}
							}
						}

						if (Nsame == GSphyPos.at(Index).at(i).size())
						{
							frequency += 1;
							if (admixedGen.inds.at(j).Sex == 1)
							{
								male_f += 1;
							}
							else if (admixedGen.inds.at(j).Sex == 2)
							{
								female_f += 1;
							}
						}				
					}
				}
				tempaf += GphyPos_allele.at(Index).at(i).at(m);
				tempaf += "\t";
				frequency /= TotalHap;
				if (chr == "X")
				{
					male_f /= GenMaleNum;
					female_f /= (GenFemaleNum*2);
				}
				else if (chr == "A" && AncSex == 1)
				{
					male_f /= (GenMaleNum*2);
					female_f /= (GenFemaleNum*2);
				}

				tempaf += ToString(frequency);
				tempaf += "\t";
				if (AncSex == 1)
				{
					tempaf += ToString(male_f);
					tempaf += "\t";
					tempaf += ToString(female_f);	
				}
				else
				{
					tempaf += "-\t-";
				}
				alleleString.push_back(tempaf);
			}	
		}	
	}
	else
	{
		admixedGen.setIndRange();
	}

	AllPop.at(name).Gens[Index] = admixedGen;
	std::vector<int> tempRealNe;
	tempRealNe.resize(2);
	tempRealNe[0] = GenMaleNum;
	tempRealNe[1] = GenFemaleNum;
	AllPop.at(name).addRealNe(Index, tempRealNe);	
}

