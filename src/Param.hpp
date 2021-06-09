#ifndef PARAM_HPP_
#define PARAM_HPP_

#include <string>
#include <vector>
#include <map>

const std::string kProgram = "AdmixSim";
const std::string kVersion = "2";

class Param
{
private:
	std::string modfile;
	std::string snvfile;
	std::string hapfile;
	std::string indfile;

	std::vector<std::string> popList;
	std::vector<int> popGen;
	std::vector<int> popNumber;
	unsigned long seed;
	std::string chr;

	std::string recEvent;
	std::string mutEvent;
	std::string selEvent;
	double recRate;
	double mutRate;
	std::string outPrefix;
	
	bool givenRecRate;
	bool givenMutRate;
	bool givenPop;
	bool givenGen;

	void help();

public:
	Param(int argc, char **argv);
	virtual ~Param();
	std::string getModfile() const;
	std::string getSnvfile() const;
	std::string getHapfile() const;
	std::string getIndfile() const;

	std::vector<std::string> getpopList() const;
	std::vector<int> getpopGen() const;
	std::vector<int> getpopNumber() const;
	unsigned long getSeed() const;
	std::string getChr() const;

	std::string getRecEvent() const;
	std::string getMutEvent() const;
	std::string getSelEvent() const;
	double getRecRate() const;
	double getMutRate() const;
	bool getGivenRecRate() const;
	bool getGivenMutRate() const;
	bool getGivenPop() const;
	bool getGivenGen() const;
	std::string getoutPrefix() const;
};

#endif /* PARAM_HPP_ */
