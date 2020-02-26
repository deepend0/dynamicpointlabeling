/*
 * PresentFileSetDCGG.h
 *
 *  Created on: Nov 4, 2018
 *      Author: oakile
 */

#ifndef SIMULATION_CONFGRAPHGEN_PRESENTFILESETDCGG_H_
#define SIMULATION_CONFGRAPHGEN_PRESENTFILESETDCGG_H_

#include <string>

#include "ConflictGraph.h"
#include "ConflictGraphGenerator.h"

class PresentFileSetDCGG : public ConflictGraphGenerator {
public:
	struct Parameters {
		Parameters(std::string path);
		std::string path;
	};
	PresentFileSetDCGG(Parameters parameters);
	virtual ~PresentFileSetDCGG();
	virtual labelplacement::ConflictGraph* generate(int forwardStates);
private:
	std::string path;

	int currentState;
};

#endif /* SIMULATION_CONFGRAPHGEN_PRESENTFILESETDCGG_H_ */
