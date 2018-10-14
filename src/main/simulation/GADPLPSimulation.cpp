/*
 * GADPLPSimulation.cpp
 *
 *  Created on: Sep 29, 2018
 *      Author: oakile
 */

#include "GADPLPSimulation.h"

#include <iostream>
#include <cstdlib>

GADPLPSimulationParameters::GADPLPSimulationParameters(int problemUpdatePeriod, int optimizationPeriod, int numberOfPeriods) {
	this->problemUpdatePeriod = problemUpdatePeriod;
	this->optimizationPeriod = optimizationPeriod;
	this->numberOfPeriods = numberOfPeriods;
}
GADPLPSimulation::SolutionContext::SolutionContext(labelplacement::ConflictGraph* sourcePrb, labelplacement::ConflictGraph* targetPrb, labelplacement::Solution* solution, int optimizationTime) {
	this->sourcePrb = sourcePrb;
	this->targetPrb = targetPrb;
	this->solution = solution;
	this->optimizationTime = optimizationTime;
}
GADPLPSimulation::SolutionContext::~SolutionContext() {
}
GADPLPSimulation::GADPLPSimulation(DynamicConflictGraphGenerator* confgraphgen, GADPLP::GADPLP* gadplp, GADPLPSimulationParameters parameters) {
	this->confgraphgen = confgraphgen;
	this->gadplp = gadplp;
	this->problemUpdatePeriod = parameters.problemUpdatePeriod;
	this->optimizationPeriod = parameters.optimizationPeriod;
	this->numberOfPeriods = parameters.numberOfPeriods;
}
GADPLPSimulation::GADPLPSimulation(GADPLP::GADPLPParameters& gadplpParameters, DynamicConflictGraphGeneratorParameters& cggParameters, GADPLPSimulationParameters& simulationParameters) {
	this->gadplp = new GADPLP::GADPLP();
	this->gadplp->init(gadplpParameters);
	this->confgraphgen = new DynamicConflictGraphGenerator(cggParameters);
	this->problemUpdatePeriod = simulationParameters.problemUpdatePeriod;
	this->optimizationPeriod = simulationParameters.optimizationPeriod;
	this->numberOfPeriods = simulationParameters.numberOfPeriods;
}

GADPLPSimulation::~GADPLPSimulation() {
	delete gadplp;
	delete confgraphgen;
}

std::vector<GADPLPSimulation::SolutionContext>* GADPLPSimulation::runSimulation() {
	std::vector<SolutionContext>* result = new std::vector<SolutionContext>();
	int lastProblemUpdatePeriod = -1;

	for(int i=0; i<numberOfPeriods; i++) {
		int sourceUpdatePeriod =  (i*optimizationPeriod) / problemUpdatePeriod;
		std::cout<<"GADPLPSimulation::runSimulation\t"<<"Current Optimization Period"<<i<<" Source Update Period: "<<sourceUpdatePeriod<<std::endl;
		labelplacement::ConflictGraph* lastConfGraph = confgraphgen->generate(sourceUpdatePeriod - lastProblemUpdatePeriod);
		lastProblemUpdatePeriod = sourceUpdatePeriod;

		clock_t startG = clock();
		labelplacement::Solution& solution = gadplp->optimize(*lastConfGraph);
		clock_t endG = clock();
		int optimizationTime = 1000 * (endG - startG) / CLOCKS_PER_SEC;
		int targetUpdatePeriod = (i*optimizationPeriod + optimizationTime) / problemUpdatePeriod;
		std::cout<<"GADPLPSimulation::runSimulation"<<"\tOptimization Time:"<<optimizationTime<<"\tTarget Update Period: "<<targetUpdatePeriod<<std::endl;
		labelplacement::ConflictGraph* currentConflictGraph = confgraphgen->generate(targetUpdatePeriod-lastProblemUpdatePeriod);
		result->push_back(SolutionContext(lastConfGraph, currentConflictGraph, &solution, optimizationTime));
		lastProblemUpdatePeriod = targetUpdatePeriod;
	}
	return result;
}
