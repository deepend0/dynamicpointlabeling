/*
 * Test.cpp
 *
 *  Created on: Aug 27, 2018
 *      Author: oakile
 */

#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include "ConflictGraphLoader.h"
#include "ConflictGraphGenerator.h"
#include "PresentFileSetDCGG.h"
#include "RandomDynamicConflictGraphGenerator.h"
#include "Solution.h"
#include "GADPLP.h"
#include "GADPLPSimulation.h"
#include "FHPLPOptimizerAdaptor.h"
#include "fhplp/FastHeuristicPLP.hpp"

void runGADPLPSimulation(GADPLP::GADPLPParameters& gadplpParameters, RandomDynamicConflictGraphGeneratorParameters& cggParameters, GADPLPSimulationParameters& simulationParameters);
void runGADPLPSimulation(GADPLP::GADPLPParameters& gadplpParameters, PresentFileSetDCGG::Parameters& cggParameters, GADPLPSimulationParameters& simulationParameters);
void runDPLPSimulation(std::function<GADPLPSimulation::Optimizer> optimizer, ConflictGraphGenerator* cgg, GADPLPSimulationParameters& simulationParameters, labelplacement::Solution* initialSolution);
void runGADPLPFromFile(std::string problemPath, int numInstances, int* instanceNumbers, GADPLP::GADPLPParameters gaplpParameters);
void runFHPLPSimulation(fhplp::FastHeuristicPLPParameters& fhplpParameters, PresentFileSetDCGG::Parameters& cggParameters, GADPLPSimulationParameters& simulationParameters);
void generateConflictGraphToFile(RandomDynamicConflictGraphGeneratorParameters& cggParameters, std::string outputPath, int numStates);
void testDynamicConflictGraphGenerator();

struct ConflictSizeSimulationObserver {
	~ConflictSizeSimulationObserver() {
		if(lastSolution!=NULL){
			delete lastSolution->getConflictGraph();
			delete lastSolution;
			if(lastConflictGraphFree) {
				delete lastConflictGraph;
			}
		}
	}
	labelplacement::Solution* lastSolution = NULL;
	labelplacement::ConflictGraph* lastConflictGraph = NULL;
	bool lastConflictGraphFree = false;
	double averageConflictSize=0;
	double averageConflictSizeToFirstTarget=0;
	double totalTime=0;
	int currentPeriod=0;
	double averageConflictSizeToSource=0;
	double averageConflictGraphChange=0;
	std::vector<int> conflictGraphSizes;
	std::vector<int> conflictGraphChanges;
	int numSolutions=0;
	std::vector<int> solutionTimes;
	void operator()(GADPLPSimulation::GADPLPSimulationPeriodNotification* periodNotification) {
		currentPeriod++;
		labelplacement::ConflictGraph* cg = periodNotification->conflictGraph;
		if(lastConflictGraph!=NULL) {
			periodNotification->conflictGraphChange = labelplacement::ConflictGraph::compareConflictGraphs(lastConflictGraph->getConflictGraphOfPositions()->getVertexNumber(),
					lastConflictGraph->getConflictGraphOfPositions()->getAdjacencyList(),
					cg->getConflictGraphOfPositions()->getAdjacencyList());
		}

		std::cout<<"CONFLICTS ALL: "<<periodNotification->conflictGraph->getConflictGraphOfPositions()->getEdgeNumber()<<"\tCHG: "<<periodNotification->conflictGraphChange<<std::endl;
		if(periodNotification->conflictGraphChange!=-1) {
			averageConflictGraphChange=(averageConflictGraphChange*(currentPeriod-2) + (double)periodNotification->conflictGraphChange/(double)periodNotification->conflictGraph->getConflictGraphOfPositions()->getEdgeNumber())/(currentPeriod-1);
			conflictGraphSizes.push_back(periodNotification->conflictGraph->getConflictGraphOfPositions()->getEdgeNumber());
			conflictGraphChanges.push_back(periodNotification->conflictGraphChange);
		}
		std::vector<int*>* solutionIntervals = new std::vector<int*>();
		int* lastSolutionInterval = NULL;

		for(int i=0; i<periodNotification->solutions->size(); i++) {

			labelplacement::Solution* solution = periodNotification->solutions->at(i);
			int solutionTime = periodNotification->solutionTimes->at(i);
			if(i!=0){
				lastSolutionInterval[1]=solutionTime-1;
			}
			lastSolutionInterval = new int[4];
			solutionIntervals->push_back(lastSolutionInterval);
			if(i==0) {
				lastSolutionInterval[0]=periodNotification->periodStartTime;
			} else {
				lastSolutionInterval[0]=solutionTime;
			}
			if(solution!=NULL) {
				lastSolutionInterval[2]=solution->getConflictSize(cg);
				lastSolutionInterval[3]=solution->getConflictSize();
			} else {
				lastSolutionInterval[2]=-1;
				lastSolutionInterval[3]=-1;
			}

			if(lastSolution!=solution) {
				if(lastSolution!=NULL){
					if(lastSolution->getConflictGraph()!=solution->getConflictGraph()){
						delete lastSolution->getConflictGraph();
					}
					delete lastSolution;
				}
				lastSolution = solution;
				solutionTimes.push_back(solutionTime);

				averageConflictSizeToSource=(averageConflictSizeToSource*numSolutions+(double)lastSolutionInterval[3])/(numSolutions+1);
				averageConflictSizeToFirstTarget=(averageConflictSizeToFirstTarget*numSolutions+(double)lastSolutionInterval[2])/(numSolutions+1);
				numSolutions++;
			}
		}
		lastSolutionInterval[1]=periodNotification->periodEndTime;

		double weightedConflictSizeOfPeriod=0;
		double firstIntervalStart=-1;
		double periodTime = 0;
		for(std::vector<int*>::iterator it=solutionIntervals->begin(); it!=solutionIntervals->end(); it++){
			int* solutionInterval = *it;
			double intervalStart = (double)solutionInterval[0]/1000.0;
			double intervalEnd = (double)solutionInterval[1]/1000.0;
			std::cout<<intervalStart<<"\t"<<intervalEnd<<"\t"<<solutionInterval[3]<<"\t"<<solutionInterval[2]<<std::endl;
			if(solutionInterval[2]!=-1){
				if(firstIntervalStart==-1){
					firstIntervalStart=intervalStart;
				}
				int intervalTime = intervalEnd-intervalStart+1.0;
				periodTime += intervalTime;
				weightedConflictSizeOfPeriod+=intervalTime*(double)solutionInterval[2];
			}
		}
		int* solutionInterval = *(solutionIntervals->end()-1);
		double lastIntervalEnd = (double)solutionInterval[1]/1000.0;

		if(periodTime!=0) {
			averageConflictSize = (averageConflictSize*totalTime+weightedConflictSizeOfPeriod)/(totalTime+periodTime);
			totalTime+=periodTime;
		}

		if(lastConflictGraphFree) {
			delete lastConflictGraph;
		}
		lastConflictGraph = cg;
		if(!periodNotification->optimizationOccured) {
			lastConflictGraphFree=true;
		} else {
			lastConflictGraphFree=false;
		}
		delete periodNotification;
	}
};

int main(int argc, char** argv) {
	if(argc>0){
		std::string testType=std::string(argv[1]);
		if(testType=="fromfile" || testType=="withsimulation" || testType=="withsimulationfromfile" || testType=="immigranttest") {
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
				int mode = atoi(argv[16]);
				int problemUpdatePeriod = atoi(argv[17]);
				int optimizationPeriod = atoi(argv[18]);
				int numberOfPeriods = atoi(argv[19]);
				int numberOfOptimizations = atoi(argv[20]);

				GADPLP::GADPLPParameters gadplpParameters(numPoints, numPositionsPerPoint, coefPop*numPoints, crossoverRate, mutationRate, 0, true, optimizationPeriod,
						coefInitGen*numPoints, coefImmGen*numPoints, immRate, 0, groupProportion, groupProportionMargin, 0, 0, 0);
				RandomDynamicConflictGraphGeneratorParameters cggParameters(areaWidth, areaHeight, numPoints, numPositionsPerPoint, labelWidth, labelHeight);
				GADPLPSimulationParameters simulationParameters(mode, problemUpdatePeriod, optimizationPeriod, numberOfPeriods, numberOfOptimizations);
				runGADPLPSimulation(gadplpParameters, cggParameters, simulationParameters);
			} else if(testType=="withsimulationfromfilecoef") {
				std::string problemPath = std::string(argv[12]);
				int mode = atoi(argv[13]);
				int problemUpdatePeriod = atoi(argv[14]);
				int optimizationPeriod = atoi(argv[15]);
				int numberOfPeriods = atoi(argv[16]);
				int numberOfOptimizations = atoi(argv[17]);

				GADPLP::GADPLPParameters gadplpParameters(numPoints, numPositionsPerPoint, coefPop*numPoints, crossoverRate, mutationRate, 0, true, optimizationPeriod,
						coefInitGen*numPoints, coefImmGen*numPoints, immRate, 0, groupProportion, groupProportionMargin, 0, 0, 0);
				PresentFileSetDCGG::Parameters cggParameters(problemPath);
				GADPLPSimulationParameters simulationParameters(mode, problemUpdatePeriod, optimizationPeriod, numberOfPeriods, numberOfOptimizations);
				runGADPLPSimulation(gadplpParameters, cggParameters, simulationParameters);
			} else if(testType=="withsimulationfromfile") {
				std::string problemPath = std::string(argv[12]);
				int mode = atoi(argv[13]);
				int problemUpdatePeriod = atoi(argv[14]);
				int optimizationPeriod = atoi(argv[15]);
				int numberOfPeriods = atoi(argv[16]);
				int numberOfOptimizations = atoi(argv[17]);

				GADPLP::GADPLPParameters gadplpParameters(numPoints, numPositionsPerPoint, (int)coefPop, crossoverRate, mutationRate, 0, true, optimizationPeriod,
						(int)coefInitGen, (int)coefImmGen, immRate, 0, groupProportion, groupProportionMargin, 0, 0, 0);
				PresentFileSetDCGG::Parameters cggParameters(problemPath);
				GADPLPSimulationParameters simulationParameters(mode, problemUpdatePeriod, optimizationPeriod, numberOfPeriods, numberOfOptimizations);
				runGADPLPSimulation(gadplpParameters, cggParameters, simulationParameters);
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
		} else if(testType=="fhplpsimulation") {
			int numPoints = atoi(argv[2]);
			int numPositionsPerPoint = atoi(argv[3]);
			int greedyMethod = atoi(argv[4]);
			int greedySubmethod = atoi(argv[5]);
			double a = atof(argv[6]);
			double b = atof(argv[7]);
			int hc = atoi(argv[8]);
			int localSearch = atoi(argv[9]);
			int mode = atoi(argv[10]);
			int problemUpdatePeriod = atoi(argv[11]);
			int optimizationPeriod = atoi(argv[12]);
			int numberOfPeriods = atoi(argv[13]);
			int numberOfOptimizations = atoi(argv[14]);
			std::string problemPath = std::string(argv[15]);
			fhplp::FastHeuristicPLPParameters fhplpParameters(numPoints, numPositionsPerPoint, greedyMethod, greedySubmethod, a, b, hc, localSearch);
			PresentFileSetDCGG::Parameters cggParameters(problemPath);
			GADPLPSimulationParameters simulationParameters(mode, problemUpdatePeriod, optimizationPeriod, numberOfPeriods, numberOfOptimizations);
			runFHPLPSimulation(fhplpParameters, cggParameters, simulationParameters);
		} else if(testType=="cggtest") {
			testDynamicConflictGraphGenerator();
		} else if(testType=="tofile") {
			int areaWidth = atoi(argv[2]);
			int areaHeight = atoi(argv[3]);
			int numberOfPoints = atoi(argv[4]);
			int labelWidth = atoi(argv[5]);
			int labelHeight = atoi(argv[6]);
			int positionsPerPoint = atoi(argv[7]);
			int numStates = atoi(argv[8]);
			std::string outputPath(argv[9]);
			RandomDynamicConflictGraphGeneratorParameters cggParameters(areaWidth, areaHeight, numberOfPoints, positionsPerPoint, labelWidth, labelHeight);
			if(argc == 16) {
				cggParameters.mpParamsIndicator = true;
				cggParameters.meanSpd = atof(argv[10]);
				cggParameters.stdSpd = atof(argv[11]);
				cggParameters.maxSpd = atof(argv[12]);
				cggParameters.meanDir = atof(argv[13]);
				cggParameters.stdDir = atof(argv[14]);
				cggParameters.maxDir = atof(argv[15]);
			}
			generateConflictGraphToFile(cggParameters, outputPath, numStates);
		}
	} else {
		std::cout<<"Enter test type"<<std::endl;
	}
}

void runGADPLPSimulation(GADPLP::GADPLPParameters& gadplpParameters, RandomDynamicConflictGraphGeneratorParameters& cggParameters, GADPLPSimulationParameters& simulationParameters) {
	GADPLP::GADPLP* gadplp = new GADPLP::GADPLP();
	gadplp->init(gadplpParameters);
	labelplacement::Solution* initialSolution = new labelplacement::Solution(NULL);
	int* labelPlacements = new int[gadplpParameters.numPoints];
	for(int i=0; i<gadplpParameters.numPoints; i++) {
		labelPlacements[i]=i%gadplpParameters.numPositionsPerPoint;
	}
	initialSolution->setLabelPlacements(labelPlacements);
	ConflictGraphGenerator* confgraphgen = new RandomDynamicConflictGraphGenerator(cggParameters);
	std::function<GADPLPSimulation::Optimizer> optimizer = std::bind(&GADPLP::GADPLP::optimize, gadplp, std::placeholders::_1);
	runDPLPSimulation(optimizer, confgraphgen, simulationParameters, initialSolution);
	delete gadplp;
}

void runGADPLPSimulation(GADPLP::GADPLPParameters& gadplpParameters, PresentFileSetDCGG::Parameters& cggParameters, GADPLPSimulationParameters& simulationParameters) {
	GADPLP::GADPLP* gadplp = new GADPLP::GADPLP();
	gadplp->init(gadplpParameters);
	labelplacement::Solution* initialSolution = new labelplacement::Solution(NULL);
	int* labelPlacements = new int[gadplpParameters.numPoints];
	for(int i=0; i<gadplpParameters.numPoints; i++) {
		labelPlacements[i]=i%gadplpParameters.numPositionsPerPoint;
	}
	initialSolution->setLabelPlacements(labelPlacements);
	ConflictGraphGenerator* confgraphgen = new PresentFileSetDCGG(cggParameters);
	std::function<GADPLPSimulation::Optimizer> optimizer = std::bind(&GADPLP::GADPLP::optimize, gadplp, std::placeholders::_1);
	runDPLPSimulation(optimizer, confgraphgen, simulationParameters, initialSolution);
	delete gadplp;
}

void runFHPLPSimulation(fhplp::FastHeuristicPLPParameters& fhplpParameters, PresentFileSetDCGG::Parameters& cggParameters, GADPLPSimulationParameters& simulationParameters) {
	fhplp::FastHeuristicPLP fhplp;
	fhplp.init(fhplpParameters);
	fhplp::FHPLPOptimizerAdaptor fhplpOptimizerAdaptor(&fhplp);
	labelplacement::Solution* initialSolution = new labelplacement::Solution(NULL);
	int* labelPlacements = new int[fhplpParameters.numPoints];
	for(int i=0; i<fhplpParameters.numPoints; i++) {
		labelPlacements[i]=i%fhplpParameters.numPositionsPerPoint;
	}
	initialSolution->setLabelPlacements(labelPlacements);
	ConflictGraphGenerator* confgraphgen = new PresentFileSetDCGG(cggParameters);
	std::function<GADPLPSimulation::Optimizer> optimizer = std::bind(&fhplp::FHPLPOptimizerAdaptor::optimize, fhplpOptimizerAdaptor, std::placeholders::_1);
	runDPLPSimulation(optimizer, confgraphgen, simulationParameters, initialSolution);
}

void runDPLPSimulation(std::function<GADPLPSimulation::Optimizer> optimizer, ConflictGraphGenerator* cgg, GADPLPSimulationParameters& simulationParameters, labelplacement::Solution* initialSolution) {
	GADPLPSimulation gadplpSimulation(cgg, optimizer, simulationParameters);
	std::function<void(GADPLPSimulation::GADPLPSimulationPeriodNotification*)> conflictSizeSimulationObserverFn = ConflictSizeSimulationObserver();
	ConflictSizeSimulationObserver* conflictSizeSimulationObserver = conflictSizeSimulationObserverFn.target<ConflictSizeSimulationObserver>();

	std::cout<<"START\tEND\tCONF SRC\tCONF TGT"<<std::endl;
	gadplpSimulation.runSimulation(conflictSizeSimulationObserverFn, initialSolution);
	std::cout<<"SOL TMS";
	double totalSolTime=0;
	int numSolutions=conflictSizeSimulationObserver->solutionTimes.size();
	for(int i=0; i<numSolutions;i++){
		double solTime = (double)conflictSizeSimulationObserver->solutionTimes.at(i)/1000.0;
		std::cout<<"\t"<<solTime;
		totalSolTime+=solTime;
	}
	std::cout<<std::endl;
	double lastSolTime = (double)conflictSizeSimulationObserver->solutionTimes.at(numSolutions-1)/1000.0;
	double avgSoltTime = lastSolTime/(double)numSolutions;
	double confGraphChangeSum=0;
	int confGraphChangeCount=conflictSizeSimulationObserver->conflictGraphChanges.size();
	for(int i=0; i<confGraphChangeCount; i++) {
		confGraphChangeSum+=(double)conflictSizeSimulationObserver->conflictGraphChanges[i]/(double)conflictSizeSimulationObserver->conflictGraphSizes[i];
	}
	double avgConfGraphChange = confGraphChangeSum/(double)confGraphChangeCount;
	double meanDistanceSquaredSum=0;
	for(int i=0; i<confGraphChangeCount; i++) {
		meanDistanceSquaredSum+=pow((double)conflictSizeSimulationObserver->conflictGraphChanges[i]/(double)conflictSizeSimulationObserver->conflictGraphSizes[i]-avgConfGraphChange,2);
	}
	double confGraphChangeStdDev = sqrt(meanDistanceSquaredSum/(double)confGraphChangeCount);
	std::cout<<"AVG CONF\tAVG CONF SRC\tAVG CONF TARGET\tAVG SOL TIME\tAVG CHG\tST DEV CHG"<<std::endl;
	std::cout<<"OVERALL\t"<<conflictSizeSimulationObserver->averageConflictSize<<"\t"<<conflictSizeSimulationObserver->averageConflictSizeToSource<<"\t"
			<<conflictSizeSimulationObserver->averageConflictSizeToFirstTarget<<"\t"<<avgSoltTime<<"\t"<<avgConfGraphChange<<"\t"<<confGraphChangeStdDev<<std::endl;
}

void runGADPLPFromFile(std::string problemPath, int numInstances, int* instanceNumbers, GADPLP::GADPLPParameters gaplpParameters){
	labelplacement::ConflictGraphLoader cgl;
	GADPLP::GADPLP gadplp;
	gadplp.init(gaplpParameters);
	for(int i=0; i<numInstances; i++) {
		int instanceNumber = instanceNumbers[i];
		std::string filename = "d" + std::to_string(gaplpParameters.numPoints) + "_" + (instanceNumber < 10 ? "0" : "") + std::to_string(instanceNumber) + ".dat";
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
	RandomDynamicConflictGraphGeneratorParameters params(800, 600, 100, 4, 20, 5);
	RandomDynamicConflictGraphGenerator cgg(params);
	for(int i=0;i<10;i++) {
		std::cout<<i<<std::endl;
		cgg.generate(1);
	}
}

void generateConflictGraphToFile(RandomDynamicConflictGraphGeneratorParameters& cggParameters, std::string outputPath, int numStates)
{
	RandomDynamicConflictGraphGenerator cgg(cggParameters);
	labelplacement::ConflictGraph* lastCG = NULL;
	std::vector<int> targetConflicts;
	std::vector<int> diffs;
	for(int i=0; i<numStates; i++) {
		labelplacement::ConflictGraph* cg = cgg.generate(1);
		std::string filePath = outputPath + "/" + (i+1<10?"0":"") + std::to_string(i+1) + std::string(".dat");
		labelplacement::ConflictGraph::writeConflictGraphIntoFile(cg->getConflictGraphOfPositions()->getAdjacencyList(), cggParameters.numberOfPoints, cggParameters.positionsPerPoint, filePath);
		if(lastCG!=NULL) {
			int diff = labelplacement::ConflictGraph::compareConflictGraphs(lastCG->getConflictGraphOfPositions()->getVertexNumber(),
					lastCG->getConflictGraphOfPositions()->getAdjacencyList(),
					cg->getConflictGraphOfPositions()->getAdjacencyList());
			targetConflicts.push_back(cg->getConflictGraphOfPositions()->getEdgeNumber());
			diffs.push_back(diff);
			delete lastCG;
		}
		lastCG=cg;
	}
	delete lastCG;

	std::ofstream ofsParams(outputPath+"/params.txt", std::ios_base::out);
	ofsParams<<"Area Height:\t"<<cggParameters.areaHeight<<std::endl;
	ofsParams<<"Area Width:\t"<<cggParameters.areaWidth<<std::endl;
	ofsParams<<"Label Height:\t"<<cggParameters.labelHeight<<std::endl;
	ofsParams<<"Label Width:\t"<<cggParameters.labelWidth<<std::endl;
	ofsParams<<"Number of Points:\t"<<cggParameters.numberOfPoints<<std::endl;
	ofsParams<<"Number of Positions Per Point:\t"<<cggParameters.positionsPerPoint<<std::endl;
	if(cggParameters.mpParamsIndicator) {
		ofsParams<<"Mean Speed:\t"<<cggParameters.meanSpd<<std::endl;
		ofsParams<<"Std Dev Speed:\t"<<cggParameters.stdSpd<<std::endl;
		ofsParams<<"Max Speed:\t"<<cggParameters.maxSpd<<std::endl;
		ofsParams<<"Mean Direction:\t"<<cggParameters.meanDir<<std::endl;
		ofsParams<<"Std Dev Direction:\t"<<cggParameters.stdDir<<std::endl;
		ofsParams<<"Max Direction:\t"<<cggParameters.maxDir<<std::endl;
	}
	ofsParams.close();

	std::ofstream ofsDiff(outputPath+"/diffs.txt", std::ios_base::out);
	for(std::vector<int>::iterator it = diffs.begin(), it2 = targetConflicts.begin(); it!=diffs.end(); it++, it2++) {
		ofsDiff<<*it2<<"\t"<<*it<<std::endl;
	}
	ofsDiff.close();

}
