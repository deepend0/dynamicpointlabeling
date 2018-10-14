/*
 * RandomImmigrantSGMXGA.cpp
 *
 *  Created on: Aug 19, 2018
 *      Author: oakile
 */

#include <iostream>
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
	if(curStage % period == 0) {
		for(int i=0;i<pop->size(); i++) {
			pop->individual(i).initializer(genomePrototype->initializer());
		}
		SimpleGMXGA::nGenerations(numInitGen);
	} else {
		//REPLACE WORST INDIVIDUALS
		int populationSize = pop->size();
		int numGenomesToReplace = immigrationRate*populationSize;

		for(int i=0; i<numGenomesToReplace; i++) {
			GA1DArrayAlleleGenome<int>* oldGenome = (GA1DArrayAlleleGenome<int>*) &pop->worst(i);
			GA1DArrayAlleleGenome<int>* newGenome = new GA1DArrayAlleleGenome<int>(oldGenome->size(), oldGenome->alleleset(), oldGenome->evaluator(),
					NULL);
			newGenome->crossover(oldGenome->sexual());
			newGenome->initializer(genomePrototype->initializer());
			newGenome->initialize();
			newGenome->initializer(NoGenomeInitializer);
			pop->replace(newGenome, oldGenome);
		}

		if(curStage % period == 1) {
			for(int i=0;i<pop->size(); i++) {
				pop->individual(i).initializer(NoGenomeInitializer);
			}
		}
		SimpleGMXGA::nGenerations(numImmigrateGen);
	}
	SimpleGMXGA::initialize(seed);
}

void RandomImmigrantSGMXGA::evolve(unsigned int seed){
	GASimpleGA::evolve(seed);
	curStage++;
}
} /* namespace gmxga */
