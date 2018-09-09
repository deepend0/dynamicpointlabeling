/*
 * RandomImmigrantSGMXGA.cpp
 *
 *  Created on: Aug 19, 2018
 *      Author: oakile
 */

#include "RandomImmigrantSGMXGA.h"

namespace gmxga {

RandomImmigrantSGMXGA::RandomImmigrantSGMXGA(const GAGenome& gaGenome) : SimpleGMXGA(gaGenome) {
	period = 1;
	curStage = 0;
	numInitGen = 0;
	numImmigrateGen = 0;
	immigrationRate = 0;
}

RandomImmigrantSGMXGA::~RandomImmigrantSGMXGA() {
	// TODO Auto-generated destructor stub
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
	GARandomSeed(seed);

	if(curStage % period == 0) {
		pop->initialize();
		SimpleGMXGA::nGenerations(numInitGen);
	} else {
		//REPLACE WORST INDIVIDUALS
		int populationSize = pop->size();
		int numGenomesToReplace = immigrationRate*populationSize;
		GA1DArrayAlleleGenome<int>* genomesToReplace[numGenomesToReplace];
		for(int i=0; i<numGenomesToReplace; i++) {
			genomesToReplace[i] = (GA1DArrayAlleleGenome<int>*) &pop->worst(i);
		}
		for(int i=0; i<numGenomesToReplace; i++) {
			GA1DArrayAlleleGenome<int>* genomeToReplace = genomesToReplace[i];
			for(int j=0; j<genomeToReplace->size(); j++) {
				genomeToReplace->gene(j, genomeToReplace->alleleset(j).allele());
			}
		}
		SimpleGMXGA::nGenerations(numImmigrateGen);
	}
	pop->evaluate(gaTrue);

	stats.reset(*pop);

	if(!scross)
		GAErr(GA_LOC, className(), "initialize", gaErrNoSexualMating);
}

void RandomImmigrantSGMXGA::evolve(unsigned int seed){
	GASimpleGA::evolve(seed);
	curStage++;
}
} /* namespace gmxga */
