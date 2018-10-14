/*
 * GADPLPSimulation.h
 *
 *  Created on: Sep 29, 2018
 *      Author: oakile
 */

#ifndef _SIMULATION_GADPLPSIMULATION_H_
#define _SIMULATION_GADPLPSIMULATION_H_

#include "DynamicConflictGraphGenerator.h"
#include "GADPLP.h"
#include "Solution.h"

#include <vector>

struct GADPLPSimulationParameters {
	GADPLPSimulationParameters(int problemUpdatePeriod, int optimizationPeriod, int numberOfPeriods);
	int problemUpdatePeriod;
	int optimizationPeriod;
	int numberOfPeriods;
};
class GADPLPSimulation {
private:

public:
	struct SolutionContext {
		SolutionContext(labelplacement::ConflictGraph* sourcePrb, labelplacement::ConflictGraph* targetPrb, labelplacement::Solution* solution, int optimizationTime);
		~SolutionContext();
		labelplacement::ConflictGraph* sourcePrb = NULL;
		labelplacement::ConflictGraph* targetPrb = NULL;
		labelplacement::Solution* solution = NULL;
		int optimizationTime = 0;
	};
	GADPLPSimulation(DynamicConflictGraphGenerator* confgraphgen, GADPLP::GADPLP* gadplp, GADPLPSimulationParameters parameters);
	GADPLPSimulation(GADPLP::GADPLPParameters& gadplpParameters, DynamicConflictGraphGeneratorParameters& cggParameters, GADPLPSimulationParameters& simulationParameters);
	virtual ~GADPLPSimulation();
	std::vector<SolutionContext>* runSimulation();
private:
	DynamicConflictGraphGenerator* confgraphgen;
	GADPLP::GADPLP* gadplp;
	int problemUpdatePeriod;
	int optimizationPeriod;
	int numberOfPeriods;
};


#endif /* _SIMULATION_GADPLPSIMULATION_H_ */
