/****************************************
 *	Project Name: AdmixSim 2
 *	File Name: Utils.hpp
 *	Author: Rui Zhang
 *	Date: July 7, 2020
 ****************************************/

#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <vector>
#include <string>

int findFloatPos(const std::vector<double> &poss, const double &pos, int &index);
void MapToPhyPos(const double &newRandom, const std::vector<double> &IniPoss, const std::vector<int> &phyPoss, int &newPhyPos);
int findIntPos(const std::vector<int> &poss, const int &pos, int &index); 
void min(const int&a, const int&b, int &c);
void max(const int&a, const int&b, int &c);
std::string ToString(const double &a); 

#endif /* UTILS_HPP_ */
