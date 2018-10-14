/*
 * Test.cpp
 *
 *  Created on: Aug 27, 2018
 *      Author: oakile
 */

#include <cstdlib>
#include <string>
#include <vector>

#include "ConflictGraphLoader.h"
#include "GADPLP.h"
#include "GADPLPSimulation.h"
#include "Solution.h"

void runGADPLPWithSimulation(GADPLP::GADPLPParameters& gadplpParameters, DynamicConflictGraphGeneratorParameters& cggParameters, GADPLPSimulationParameters& simulationParameters);
void runGADPLPFromFile(std::string problemPath, int numInstances, int* instanceNumbers, GADPLP::GADPLPParameters gaplpParameters);
void testDynamicConflictGraphGenerator();

int main(int argc, char** argv) {
	if(argc>0){
		std::string testType=std::string(argv[1]);
		if(testType=="fromfile" || testType=="withsimulation" || testType=="immigranttest") {
			std::cout<<"Params";
			for(int i=1; i<argc; i++){
				std::cout<<" "<<argv[i];
			}
			std::cout<<std::endl;
			double coefPop=atof(argv[2]);
			double crossoverRate = atof(argv[3]);
			double mutationRate = atof(argv[4]);
			double coefInitGen=atof(argv[5]);
			double coefImmGen=atof(argv[6]);
			double immRate=atof(argv[7]);
			double groupProportion=atof(argv[8]);
			double groupProportionMargin=atof(argv[9]);
			int numPoints = atoi(argv[10]);
			int numPositionsPerPoint = atoi(argv[11]);
			if(testType=="fromfile"){
				std::string problemPath = std::string(argv[12]);
				int numInstances = atoi(argv[13]);

				GADPLP::GADPLPParameters gadplpParameters(numPoints, numPositionsPerPoint, coefPop*numPoints, crossoverRate, mutationRate, 0, true, numInstances,
						coefInitGen*numPoints, coefImmGen*numPoints, immRate, 0, groupProportion, groupProportionMargin, 0, 0, 0);
				int instanceNumbers [numInstances];
				for(int j=0; j<numInstances; j++) {
					instanceNumbers[j]=j+1;
				}
				runGADPLPFromFile(problemPath, numInstances, instanceNumbers, gadplpParameters);
			} else if(testType=="withsimulation") {
				int areaWidth = atoi(argv[12]);
				int areaHeight = atoi(argv[13]);
				int labelWidth = atoi(argv[14]);
				int labelHeight = atoi(argv[15]);
				int problemUpdatePeriod = atoi(argv[16]);
				int optimizationPeriod = atoi(argv[17]);
				int numberOfPeriods = atoi(argv[18]);

				GADPLP::GADPLPParameters gadplpParameters(numPoints, numPositionsPerPoint, coefPop*numPoints, crossoverRate, mutationRate, 0, true, optimizationPeriod,
						coefInitGen*numPoints, coefImmGen*numPoints, immRate, 0, groupProportion, groupProportionMargin, 0, 0, 0);
				DynamicConflictGraphGeneratorParameters cggParameters(areaWidth, areaHeight, numPoints, numPositionsPerPoint, labelWidth, labelHeight);
				GADPLPSimulationParameters simulationParameters(problemUpdatePeriod, optimizationPeriod, numberOfPeriods);
				runGADPLPWithSimulation(gadplpParameters, cggParameters, simulationParameters);
			} else if(testType=="immigranttest") {
				std::string problemPath = std::string(argv[12]);
				int numInstances = atoi(argv[13]);

				GADPLP::GADPLPParameters gadplpParameters(numPoints, numPositionsPerPoint, coefPop*numPoints, crossoverRate, mutationRate, 0, true, numInstances,
						coefInitGen*numPoints, coefImmGen*numPoints, immRate, 0, groupProportion, groupProportionMargin, 0, 0, 0);
				int instanceNumbers [numInstances];
				for(int j=0; j<numInstances; j++) {
					instanceNumbers[j]=1;
				}
				runGADPLPFromFile(problemPath, numInstances, instanceNumbers, gadplpParameters);
			} else {
				std::cout<<"Unidentified type of test!"<<std::endl;
			}
		} else if(testType=="cggtest") {
			testDynamicConflictGraphGenerator();
		}
	} else {
		std::cout<<"Enter test type"<<std::endl;
	}
}

void runGADPLPWithSimulation(GADPLP::GADPLPParameters& gadplpParameters, DynamicConflictGraphGeneratorParameters& cggParameters, GADPLPSimulationParameters& simulationParameters) {
	GADPLPSimulation gadplpSimulation(gadplpParameters, cggParameters, simulationParameters);
	std::vector<GADPLPSimulation::SolutionContext>* solutions = gadplpSimulation.runSimulation();
	for(std::vector<GADPLPSimulation::SolutionContext>::iterator it = solutions->begin(); it!=solutions->end(); it++) {
		GADPLPSimulation::SolutionContext solutionCtx = *it;
		int conflictsForSourceProblem = solutionCtx.solution->getConflictSize();
		solutionCtx.solution->setConflictGraph(solutionCtx.targetPrb);
		int conflictsForTargetProblem = solutionCtx.solution->getConflictSize();
		std::cout<<"OT: "<<solutionCtx.optimizationTime<<" CSP:"<<conflictsForSourceProblem<<" CTP:"<<conflictsForTargetProblem<<std::endl;
		delete solutionCtx.sourcePrb;
		delete solutionCtx.targetPrb;
		delete solutionCtx.solution;
	}
	delete solutions;
}

void runGADPLPFromFile(std::string problemPath, int numInstances, int* instanceNumbers, GADPLP::GADPLPParameters gaplpParameters){
	labelplacement::ConflictGraphLoader cgl;
	GADPLP::GADPLP gadplp;
	gadplp.init(gaplpParameters);
	for(int i=0; i<numInstances; i++) {
		int instanceNumber = instanceNumbers[i];
		std::string filename = (instanceNumber < 10 ? "0" : "") + std::to_string(instanceNumber) + ".dat";
		std::cout << "I: " << problemPath + filename;
		std::string instanceFullPath = problemPath + filename;
		labelplacement::ConflictGraph* cg = cgl.load(instanceFullPath);
		clock_t startG = clock();
		labelplacement::Solution& solution = gadplp.optimize(*cg);
		clock_t endG = clock();
		int gaTimeCost = 1000 * (endG - startG) / CLOCKS_PER_SEC;
		std::cout << "\tC: " << solution.getConflictSize();
		std::cout << "\tT: " << gaTimeCost << std::endl;
		delete &solution;
		delete cg;
	}
}

void testDynamicConflictGraphGenerator() {
	DynamicConflictGraphGeneratorParameters params(800, 600, 100, 4, 20, 5);
	DynamicConflictGraphGenerator cgg(params);
	for(int i=0;i<10;i++) {
		std::cout<<i<<std::endl;
		cgg.generate(1);
	}
}
