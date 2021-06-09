/****************************************
 *	Project Name: AdmimxSim 2
 *	File Name: Utils.cpp
 *	Author: Rui Zhang
 *	Date: July 7, 2020
 ****************************************/

#include "Utils.hpp"
#include <cmath>
#include <sstream>
#include <iomanip>

int findFloatPos(const std::vector<double> &poss, const double &pos, int &index)
{
	int left = 0;
	int right = static_cast<int>(poss.size());
	if (pos < poss.front())
	{
		index = left;
		return 0;
	}

	int mid = (left + right + 1)/2;
	while (left < right)
	{
		if (pos > poss.at(mid))
			left = mid;
		else
			right = mid - 1;
		mid = (left + right + 1)/2;
	}
	index = left + 1;
	return 0;
}

/****************************************************************************
 *  map generated new mutation or recombination point to physical distance  *
 ****************************************************************************/
void MapToPhyPos(const double &newRandom, const std::vector<double> &IniPoss, const std::vector<int> &phyPoss, int &newPhyPos)
{
	int k;
	findFloatPos(IniPoss, newRandom, k); 
	if (k == 0)
	{
		newPhyPos = phyPoss.front();
	}
	else
	{
		newPhyPos = phyPoss.at(k-1) + int((phyPoss.at(k) - phyPoss.at(k-1))*(newRandom - IniPoss.at(k-1))/(IniPoss.at(k) - IniPoss.at(k-1)));
	}
}

int findIntPos(const std::vector<int> &poss, const int &pos, int &index)
{
	int left = 0;
	int right = static_cast<int>(poss.size());
	if (pos <= poss.front())
	{
		index = left;
		return 0;
	}
	if (pos >= poss.back())
	{
		index = right;
		return 0;
	}
	int mid = (left + right + 1)/2;
	while (left < right)
	{
		if (pos > poss.at(mid))
			left = mid;
		else
			right = mid - 1;
		mid = (left + right + 1)/2;
	}
	index = left+1;	
	return 0;
}

void min(const int&a, const int &b, int &c)
{
	if (a <= b)
	{
		c = a;
	}
	else
	{
		c = b;
	}
}

void max(const int&a, const int &b, int &c)
{
	if (a >= b)
	{
		c = a;
	}
	else
	{
		c = b;
	}
}

std::string ToString(const double &a)
{
	std::stringstream stream;
	stream << a;
	return stream.str();
}

