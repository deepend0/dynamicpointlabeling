/*
 * GMXSimpleGA.cpp
 *
 *  Created on: Aug 13, 2018
 *      Author: oakile
 */

#include "SimpleGMXGA.h"

namespace gmxga {

SimpleGMXGA::SimpleGMXGA(const GAGenome& gaGenome) : GASimpleGA(gaGenome) {
	gmx = NULL;
}

SimpleGMXGA::~SimpleGMXGA() {
	delete gmx;
}

GroupMaskingCrossover* SimpleGMXGA::getCrossover() {
	return this->gmx;
}
void SimpleGMXGA::setCrossover(GroupMaskingCrossover* gmx) {
	this->gmx = gmx;
	GASimpleGA::crossover(GMXDelegate(gmx));
}

GAGenome::SexualCrossoverFn SimpleGMXGA::crossover(GAGenome::SexualCrossoverFn f) {
	return NULL;
}

}
