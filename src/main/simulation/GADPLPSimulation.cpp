/*
 * GADPLPSimulation.cpp
 *
 *  Created on: Sep 29, 2018
 *      Author: oakile
 */

#include "GADPLPSimulation.h"

#include <iostream>
#include <cstdlib>
#include <cmath>

GADPLPSimulationParameters::GADPLPSimulationParameters(int mode, int problemUpdatePeriod, int optimizationPeriod, int numberOfPeriods, int numberOfOptimizations) {
	this->mode = mode;
	this->problemUpdatePeriod = problemUpdatePeriod;
	this->optimizationPeriod = optimizationPeriod;
	this->numberOfPeriods = numberOfPeriods;
	this->numberOfOptimizations = numberOfOptimizations;
}
GADPLPSimulation::GADPLPSimulationPeriodNotification::GADPLPSimulationPeriodNotification() {
}
GADPLPSimulation::GADPLPSimulationPeriodNotification::~GADPLPSimulationPeriodNotification() {
	delete solutions;
	delete solutionTimes;
}
GADPLPSimulation::GADPLPSimulation(ConflictGraphGenerator* confgraphgen, std::function<Optimizer> optimizer, GADPLPSimulationParameters parameters) {
	this->confgraphgen = confgraphgen;
	this->optimizer = optimizer;
	this->mode = parameters.mode;
	this->problemUpdatePeriod = parameters.problemUpdatePeriod;
	this->optimizationPeriod = parameters.optimizationPeriod;
	this->numberOfPeriods = parameters.numberOfPeriods;
	this->numberOfOptimizations = parameters.numberOfOptimizations;
}


GADPLPSimulation::~GADPLPSimulation() {
	delete confgraphgen;
}

/*std::vector<GADPLPSimulation::SolutionContext>* GADPLPSimulation::runSimulationPeriodic() {
	std::vector<SolutionContext>* result = new std::vector<SolutionContext>();
	int lastProblemUpdatePeriod = -1;

	for(int i=0; i<numberOfPeriods; i++) {
		int sourceUpdatePeriod =  (i*optimizationPeriod) / problemUpdatePeriod;
		labelplacement::ConflictGraph* lastConfGraph = confgraphgen->generate(sourceUpdatePeriod - lastProblemUpdatePeriod);
		lastProblemUpdatePeriod = sourceUpdatePeriod;
		std::cout<<"GADPLPSimulation::runSimulation\t"<<"Optimization Period: "<<i<<" Source Update Period: "<<sourceUpdatePeriod;

		clock_t startG = clock();
		labelplacement::Solution& solution = gadplp->optimize(*lastConfGraph);
		clock_t endG = clock();
		int optimizationTime = 1000 * (endG - startG) / CLOCKS_PER_SEC;
		int targetUpdatePeriod = (i*optimizationPeriod + optimizationTime) / problemUpdatePeriod;
		labelplacement::ConflictGraph* currentConflictGraph = confgraphgen->generate(targetUpdatePeriod-lastProblemUpdatePeriod);
		result->push_back(SolutionContext(lastConfGraph, currentConflictGraph, &solution, optimizationTime));
		lastProblemUpdatePeriod = targetUpdatePeriod;

		std::cout<<" Optimization Time: "<<optimizationTime<<"\tTarget Update Period: "<<targetUpdatePeriod;
		int diff = DynamicConflictGraphGenerator::compareConflictGraphs(lastConfGraph->getConflictGraphOfPositions()->getVertexNumber(), lastConfGraph->getConflictGraphOfPositions()->getAdjacencyList(), currentConflictGraph->getConflictGraphOfPositions()->getAdjacencyList());
		std::cout<<" Diff: "<<diff<<std::endl;

		if(optimizationTime>optimizationPeriod) {
			std::cout<<"GADPLPSimulation::runSimulation::ERROR\tOptimization time is greater than optimization period."<<std::endl;
			for(std::vector<GADPLPSimulation::SolutionContext>::iterator it = result->begin(); it!=result->end(); it++) {
				GADPLPSimulation::SolutionContext solutionCtx = *it;
				delete solutionCtx.sourcePrb;
				delete solutionCtx.targetPrb;
				delete solutionCtx.solution;
			}
			delete result;
			return NULL;
		}
	}
	return result;
}*/

void GADPLPSimulation::runSimulation(std::function<void(GADPLPSimulationPeriodNotification*)>& simulationObserver, labelplacement::Solution* initialSolution) {
	//Last occured solution in timeline
	labelplacement::Solution* lastSolution = NULL;
	int lastSolutionTime=-1;
	//Future solution
	int nextOptimizationTime=0;
	int currentNumberOfOptimizations=0;
	std::vector<labelplacement::Solution*> forwardSolutions;
	std::vector<int> forwardSolutionTimes;

	for(int i=0; i<numberOfPeriods && (currentNumberOfOptimizations<numberOfOptimizations || nextOptimizationTime>=i*problemUpdatePeriod); i++) {
		GADPLPSimulationPeriodNotification* periodNotification = new GADPLPSimulationPeriodNotification();
		periodNotification->conflictGraph = confgraphgen->generate(1);
		periodNotification->solutions = new std::vector<labelplacement::Solution*>();
		periodNotification->solutionTimes = new std::vector<int>();
		//Period start and end time range inclusive
		periodNotification->periodStartTime = i*problemUpdatePeriod;
		periodNotification->periodEndTime = (i+1)*problemUpdatePeriod-1;

		/*if(i==0) {
			lastSolution = initialSolution;
			lastSolutionTime = 0;
			initialSolution->setConflictGraph(periodNotification->conflictGraph);
			forwardSolutions.push_back(lastSolution);
			forwardSolutionTimes.push_back(lastSolutionTime);
		}*/
		while(nextOptimizationTime>=periodNotification->periodStartTime && nextOptimizationTime<=periodNotification->periodEndTime
				&& currentNumberOfOptimizations<numberOfOptimizations){
			//Optimization time must be greater than 0.
			//std::cout<<"Next Opt Time: "<<nextOptimizationTime<<std::endl;
			clock_t startG = clock();
			labelplacement::Solution* forwardSolution = &optimizer(*periodNotification->conflictGraph);
			clock_t endG = clock();
			int optimizationTime = ceil((double)1000 * ((double)endG - (double)startG) / (double)CLOCKS_PER_SEC);
			int nextSolutionTime = nextOptimizationTime+optimizationTime;
			forwardSolution->setConflictGraph(periodNotification->conflictGraph);
			periodNotification->optimizationOccured=true;
			forwardSolutions.push_back(forwardSolution);
			forwardSolutionTimes.push_back(nextSolutionTime);

			currentNumberOfOptimizations++;
			if(mode==0) {
				nextOptimizationTime+=optimizationTime;
			} else if(mode==1) {
				nextOptimizationTime+=optimizationPeriod;
			}
		}

		if(forwardSolutionTimes.size()==0 || forwardSolutionTimes[0]!=periodNotification->periodStartTime){
			periodNotification->solutions->push_back(lastSolution);
			periodNotification->solutionTimes->push_back(lastSolutionTime);
		}

		while(forwardSolutionTimes.size()>0 && forwardSolutionTimes[0]>=periodNotification->periodStartTime && forwardSolutionTimes[0]<=periodNotification->periodEndTime){
			int forwardSolutionTime=forwardSolutionTimes[0];
			labelplacement::Solution* forwardSolution = forwardSolutions[0];

			//std::cout<<"Fwd Sol Time: "<<forwardSolutionTime<<std::endl;
			periodNotification->solutions->push_back(forwardSolution);
			periodNotification->solutionTimes->push_back(forwardSolutionTime);

			lastSolution = forwardSolution;
			lastSolutionTime = nextOptimizationTime;

			forwardSolutionTimes.erase(forwardSolutionTimes.begin());
			forwardSolutions.erase(forwardSolutions.begin());
		}
		simulationObserver(periodNotification);
	}
}
