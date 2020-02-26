/*
 * RandomImmigrantSGMXGA.cpp
 *
 *  Created on: Aug 19, 2018
 *      Author: oakile
 */

#include <iostream>
#include <vector>
#include "RandomImmigrantSGMXGA.h"

namespace gmxga {
void RandomImmigrantSGMXGA::NoGenomeInitializer(GAGenome& genome) {

}
RandomImmigrantSGMXGA::RandomImmigrantSGMXGA(const GAGenome& gaGenome) : SimpleGMXGA(gaGenome) {
	period = 1;
	curStage = 0;
	numInitGen = 0;
	numImmigrateGen = 0;
	immigrationRate = 0;
	genomePrototype = NULL;
}

RandomImmigrantSGMXGA::~RandomImmigrantSGMXGA() {
	delete genomePrototype;
}

GAGenome* RandomImmigrantSGMXGA::getGenomePrototype() {
	return genomePrototype;
}

void RandomImmigrantSGMXGA::setGenomePrototype(GAGenome* genome) {
	this->genomePrototype = genome;
}

double RandomImmigrantSGMXGA::getImmigrationRate() const {
	return immigrationRate;
}

void RandomImmigrantSGMXGA::setImmigrationRate(double immigrationRate) {
	this->immigrationRate = immigrationRate;
}

int RandomImmigrantSGMXGA::getNumImmigrateGen() const {
	return numImmigrateGen;
}

void RandomImmigrantSGMXGA::setNumImmigrateGen(int numImmigrateGen) {
	this->numImmigrateGen = numImmigrateGen;
}

int RandomImmigrantSGMXGA::getNumInitGen() const {
	return numInitGen;
}

void RandomImmigrantSGMXGA::setNumInitGen(int numInitGen) {
	this->numInitGen = numInitGen;
}

int RandomImmigrantSGMXGA::getPeriod() const {
	return period;
}

void RandomImmigrantSGMXGA::setPeriod(int period) {
	this->period = period;
}

int RandomImmigrantSGMXGA::nGenerations(unsigned int n) {
	return -1;
}

void RandomImmigrantSGMXGA::initialize(unsigned int seed)
{
	//std::cout<<"STAGE: "<<curStage<<"/"<<period<<" ";
	if(curStage % period == 0 || immigrationRate==1.0) {
		for(int i=0;i<pop->size(); i++) {
			pop->individual(i).initializer(genomePrototype->initializer());
		}
		SimpleGMXGA::nGenerations(numInitGen);
	} else {
		//REPLACE WORST INDIVIDUALS
		int populationSize = pop->size();
		int numGenomesToReplace = immigrationRate*populationSize;

		pop->worst(0);
		std::vector<GA1DArrayAlleleGenome<int>*> worstGenomes;
		for(int i=0; i<numGenomesToReplace; i++) {
			worstGenomes.push_back((GA1DArrayAlleleGenome<int>*)&pop->individual(populationSize-1-i));
		}
		for(int i=0; i<numGenomesToReplace; i++) {
			GA1DArrayAlleleGenome<int>* oldGenome = worstGenomes.at(i);
			for(int j=0; j<oldGenome->length(); j++) {
				oldGenome->gene(j, GARandomInt(0, oldGenome->alleleset().size()-1));
			}
		}

		if(curStage % period != 0) {
			for(int i=0;i<pop->size(); i++) {
				pop->individual(i).initializer(NoGenomeInitializer);
			}
		}
		SimpleGMXGA::nGenerations(numImmigrateGen);
	}
	SimpleGMXGA::initialize(seed);
}

void RandomImmigrantSGMXGA::step() {
	GASimpleGA::step();
	/*
	//REPLACE WORST INDIVIDUALS
	int populationSize = pop->size();
	int numGenomesToReplace = immigrationRate*populationSize;

	pop->worst(0);
	std::vector<GA1DArrayAlleleGenome<int>*> worstGenomes;
	for(int i=0; i<numGenomesToReplace; i++) {
		worstGenomes.push_back((GA1DArrayAlleleGenome<int>*)&pop->individual(populationSize-1-i));
	}
	for(int i=0; i<numGenomesToReplace; i++) {
		GA1DArrayAlleleGenome<int>* oldGenome = worstGenomes.at(i);
		for(int j=0; j<oldGenome->length(); j++) {
			oldGenome->gene(j, GARandomInt(0, oldGenome->alleleset().size()-1));
		}
	}
	stats.update(*pop);
	*/
}
void RandomImmigrantSGMXGA::evolve(unsigned int seed){
	GASimpleGA::evolve(seed);
	curStage++;
}
} /* namespace gmxga */
