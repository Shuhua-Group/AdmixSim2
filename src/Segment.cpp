/****************************************
 *	Project Name: AdmixSim2
 *	File Name: Segment.cpp
 *	Author: Rui Zhang
 *	Data: July 7, 2020
 ****************************************/

#include "Segment.hpp"

Segment::Segment():start(0), end(0), label(1)
{
}

Segment::Segment(const int &start, const int &end, const int &label):start(start), end(end), label(label)
{	
}

Segment::~Segment()
{
}

void Segment::setEnd(const int &end)
{
	this -> end = end;
}

