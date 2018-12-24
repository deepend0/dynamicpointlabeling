/*
 * GADPLPSimulation.h
 *
 *  Created on: Sep 29, 2018
 *      Author: oakile
 */

#ifndef _SIMULATION_GADPLPSIMULATION_H_
#define _SIMULATION_GADPLPSIMULATION_H_

#include "GADPLP.h"
#include "Solution.h"

#include <functional>
#include <vector>
#include "ConflictGraphGenerator.h"

struct GADPLPSimulationParameters {
	GADPLPSimulationParameters(int mode, int problemUpdatePeriod, int optimizationPeriod, int numberOfPeriods);
	int mode;
	int problemUpdatePeriod;
	int optimizationPeriod;
	int numberOfPeriods;
};
class GADPLPSimulation {
private:

public:
	struct GADPLPSimulationPeriodNotification {
		GADPLPSimulationPeriodNotification();
		~GADPLPSimulationPeriodNotification();
		labelplacement::ConflictGraph* conflictGraph = NULL;
		int conflictGraphChange = -1;
		std::vector<labelplacement::Solution*>* solutions = NULL;
		int periodStartTime = 0;
		int periodEndTime = 0;
		bool optimizationOccured=false;
		std::vector<int>* solutionTimes = NULL;
	};

	GADPLPSimulation(ConflictGraphGenerator* confgraphgen, GADPLP::GADPLP* gadplp, GADPLPSimulationParameters parameters);
	virtual ~GADPLPSimulation();
	void runSimulation(std::function<void(GADPLPSimulationPeriodNotification*)>& simulationObserver);
private:
	ConflictGraphGenerator* confgraphgen;
	GADPLP::GADPLP* gadplp;
	int mode;
	int problemUpdatePeriod;
	int optimizationPeriod;
	int numberOfPeriods;
};


#endif /* _SIMULATION_GADPLPSIMULATION_H_ */
