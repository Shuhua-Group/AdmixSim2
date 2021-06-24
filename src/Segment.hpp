/****************************************
 *	Project Name: AdmixSim 2
 *	File Name: Segment.hpp
 *	Author: Rui Zhang
 *	Date: July 7, 2020
 ****************************************/
 
#ifndef SEGMENT_HPP_
#define SEGMENT_HPP_

#include <string>

/************************************************************
 *	Each segment of haplotype is identified as 				*
 *	the start physical position, the end physical position, *
 *	and the label it originates from						*
 ************************************************************/

class Segment
{
public:
	int start;
	int end;
	int label;
	Segment();		
	Segment(const int &start, const int &end, const int &label);
	virtual ~Segment();	
	void setEnd(const int &end);
};

#endif /* SEGMENT_HPP_ */
