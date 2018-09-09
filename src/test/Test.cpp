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
#include "Solution.h"

//Problem Sets
std::string problemSetsPath = "/home/oakile/Workspace/DynamicPointLabeling/test/input_data/";
const int numProblemsInSet6 = 5;
std::string pathsProblemSet6 [numProblemsInSet6] = {
	"dynamic_6/d100/",
	"dynamic_6/d250/",
	"dynamic_6/d500/",
	"dynamic_6/d750/",
	"dynamic_6/d1000/"
};

int numInstancesProblemSet6 = 25;
int numPointsProblemSet6 [numProblemsInSet6] = {100, 250, 500, 750, 1000};
int numPosPerPointProblemSet6 = 4;

void runGADPLPForInstanceSet(std::string basePath, std::string instanceSetPath, int numInstances, int* instanceNumbers, GADPLP::GADPLPParameters gaplpParameters);
void runGADPLPForProblemsOfDifferentSize(double coefPop, double coefInitGen, double coefImmGen, double immRate, double groupProportion, double groupProportionMargin);
void testImmigrationIsEffective(double coefPop, double coefInitGen, double coefImmGen, double immRate, double groupProportion, double groupProportionMargin);

int main(int argc, char** argv) {
	//GA Parameters
	/*double coefPop=1;
	double coefInitGen=0.33;
	double coefImmGen=0.15;
	double immRate=0.20;
	double groupProportion=0.5;
	double groupProportionMargin=0.2;*/
	//runGADPLPForProblemsOfDifferentSize(coefPop, coefInitGen, coefImmGen, immRate, groupProportion, groupProportionMargin);
	if(argc==8){
		int testType=atoi(argv[1]);
		if(testType==1 || testType==2) {
			double coefPop=atof(argv[2]);
			double coefInitGen=atof(argv[3]);
			double coefImmGen=atof(argv[4]);
			double immRate=atof(argv[5]);
			double groupProportion=atof(argv[6]);
			double groupProportionMargin=atof(argv[7]);
			std::cout<<"Params "<<coefPop<<" "<<coefInitGen<<" "<<coefImmGen<<" "<<immRate<<" "<<groupProportion<<" "<<groupProportionMargin<<std::endl;
			if(testType==1){
				runGADPLPForProblemsOfDifferentSize(coefPop, coefInitGen, coefImmGen, immRate, groupProportion, groupProportionMargin);
			} else if(testType==2) {
				testImmigrationIsEffective(coefPop, coefInitGen, coefImmGen, immRate, groupProportion, groupProportionMargin);
			} else {
				std::cout<<"Unidentified type of test!"<<std::endl;
			}
		}
	} else {
		std::cout<<"Wrong parameters!"<<std::endl;
	}
}

void runGADPLPForProblemsOfDifferentSize(double coefPop, double coefInitGen, double coefImmGen, double immRate, double groupProportion, double groupProportionMargin) {
	for(int i=0; i<numProblemsInSet6; i++)
	{
		GADPLP::GADPLPParameters parameters(numPointsProblemSet6[i], numPosPerPointProblemSet6, coefPop*numPointsProblemSet6[i], 0, true, numInstancesProblemSet6,
				coefInitGen*numPointsProblemSet6[i], coefImmGen*numPointsProblemSet6[i], immRate,
			0, groupProportion, groupProportionMargin,
			0, 0, 0);
		int instanceNumbers [numInstancesProblemSet6];
		for(int j=0; j<numInstancesProblemSet6; j++) {
			instanceNumbers[j]=j+1;
		}
		runGADPLPForInstanceSet(problemSetsPath, pathsProblemSet6[i], numInstancesProblemSet6, instanceNumbers, parameters);
	}
}

void runGADPLPForInstanceSet(std::string basePath, std::string instanceSetPath, int numInstances, int* instanceNumbers, GADPLP::GADPLPParameters gaplpParameters){
	labelplacement::ConflictGraphLoader cgl;
	GADPLP::GADPLP gadplp;
	gadplp.init(gaplpParameters);
	for(int i=0; i<numInstances; i++) {
		int instanceNumber = instanceNumbers[i];
		std::string filename = (instanceNumber < 10 ? "0" : "") + std::to_string(instanceNumber) + ".dat";
		std::cout << "I: " << instanceSetPath + filename;
		std::string instanceFullPath = basePath + instanceSetPath + filename;
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

void testImmigrationIsEffective(double coefPop, double coefInitGen, double coefImmGen, double immRate, double groupProportion, double groupProportionMargin) {

	for(int i=0; i<numProblemsInSet6; i++)
	{
		GADPLP::GADPLPParameters parameters(numPointsProblemSet6[i], numPosPerPointProblemSet6, coefPop*numPointsProblemSet6[i], 0, true, numInstancesProblemSet6,
				coefInitGen*numPointsProblemSet6[i], coefImmGen*numPointsProblemSet6[i], immRate,
			0, groupProportion, groupProportionMargin,
			0, 0, 0);
		int instanceNumbers [numInstancesProblemSet6];
		for(int j=0; j<10; j++) {
			instanceNumbers[j]=1;
		}
		runGADPLPForInstanceSet(problemSetsPath, pathsProblemSet6[i], 10, instanceNumbers, parameters);
	}
}
