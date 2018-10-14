/*
 * GADPLP.cpp
 *
 *  Created on: Aug 20, 2018
 *      Author: oakile
 */

#include "GADPLP.h"

namespace GADPLP {

GADPLPParameters::GADPLPParameters(){
	this->numPoints = 0;
	this->numPositionsPerPoint = 0;
	this->populationSize = 0;
	this->crossoverRate = 0;
	this->mutationRate = 0;
	this->selector = 0;
	this->elitist = 0;
	this->period = 0;
	this->numInitGen = 0;
	this->numImmigrateGen = 0;
	this->immigrationRate = 0;
	this->crossoverType = 0;
	this->groupProportion = 0;
	this->groupProportionMargin = 0;
	this->individualProportion = 0;
	this->individualProportionMargin = 0;
	this->individualScoreLimit = 0;
}

GADPLPParameters::GADPLPParameters(int numPoints, int numPositionsPerPoint, int populationSize, double crossoverRate, double mutationRate, int selector,
		bool elitist, int period, int numInitGen, int numImmigrateGen, double immigrationRate,
		int crossoverType, double groupProportion, double groupProportionMargin,
		double individualProportion, double individualProportionMargin, double individualScoreLimit) {
	this->numPoints = numPoints;
	this->numPositionsPerPoint = numPositionsPerPoint;
	this->populationSize = populationSize;
	this->crossoverRate = crossoverRate;
	this->mutationRate = mutationRate;
	this->selector = selector;
	this->elitist = elitist;
	this->period = period;
	this->numInitGen = numInitGen;
	this->numImmigrateGen = numImmigrateGen;
	this->immigrationRate = immigrationRate;
	this->crossoverType = crossoverType;
	this->groupProportion = groupProportion;
	this->groupProportionMargin = groupProportionMargin;
	this->individualProportion = individualProportion;
	this->individualProportionMargin = individualProportionMargin;
	this->individualScoreLimit = individualScoreLimit;
}

GADPLP::GADPLP(){
	ga = NULL;
	objective.conflictGraph = new labelplacement::ConflictGraph*();

}

GADPLP::~GADPLP() {
	delete objective.conflictGraph;
}

GADPLP::Objective::Objective() {
	conflictGraph = NULL;
}
float GADPLP::Objective::operator()(GAGenome & g) {
	GA1DArrayAlleleGenome<int> & genome = (GA1DArrayAlleleGenome<int> &) g;
	float objectiveValue;
	float potentialConflictSize = 0;
	float conflictSize = 0;
	int pointCount =
			(*conflictGraph)->getConflictGraphOfPoints()->getVertexNumber();
	int positionNumber = (*conflictGraph)->getPositionNumberPerPoint();
	for (int i = 0; i < genome.length(); i++) {
		int pointNumber = i;
		int position1 = genome.gene(i);
		int posIx1 = pointNumber * positionNumber + position1;

		std::vector<int>* conflictingPoints =
				(*conflictGraph)->getConflictGraphOfPoints()->getAdjacencyList()
						+ pointNumber;
		for (std::vector<int>::iterator it = conflictingPoints->begin();
				it != conflictingPoints->end(); it++) {
			int conflictingPointNo = *it;
			int conflictingPointIx = conflictingPointNo;
			int position2 = genome.gene(conflictingPointIx);
			int posIx2 = conflictingPointNo * positionNumber + position2;
			if ((*conflictGraph)->getConflictGraphOfPositions()->getAdjacencyMatrix()[posIx1][posIx2]) {
				conflictSize++;
				break;
			}
		}

		potentialConflictSize +=
				(*conflictGraph)->getConflictGraphOfPositions()->getAdjacencyList()[posIx1].size();
	}
	objectiveValue = conflictSize /*+ 2.0*(1-(double)ga->generation()/(double)ga->nGenerations())*potentialConflictSize*/;
	return conflictSize;
}

GADPLP::AlleleEvaluator::AlleleEvaluator() {
	conflictGraph = NULL;
}

float GADPLP::AlleleEvaluator::operator()(const GA1DArrayGenome<int> & g, int i) {
	int conflictSize = 0;
	int position1 = g.gene(i);
	int positionNumber = conflictGraph->getPositionNumberPerPoint();
	int posIx1 = i * positionNumber + position1;

	std::vector<int>* conflictingPoints =
			conflictGraph->getConflictGraphOfPoints()->getAdjacencyList() + i;
	for (std::vector<int>::iterator it = conflictingPoints->begin();
			it != conflictingPoints->end(); it++) {
		int conflictingPointIx = *it;
		int position2 = g.gene(conflictingPointIx);
		int posIx2 = conflictingPointIx * positionNumber + position2;
		if (conflictGraph->getConflictGraphOfPositions()->getAdjacencyMatrix()[posIx1][posIx2]) {
			conflictSize++;
			break;
		}
	}
	return conflictSize;
}

void GADPLP::init(GADPLPParameters& parameters) {
	//INITIALIZE PROBLEM
	int alleleSetVals[parameters.numPositionsPerPoint];
	for (int i = 0; i < parameters.numPositionsPerPoint; i++)
		alleleSetVals[i] = i;
	GAAlleleSet<int> alleleSet(parameters.numPositionsPerPoint, alleleSetVals);
	GA1DArrayAlleleGenome<int>* genome = new GA1DArrayAlleleGenome<int>(parameters.numPoints, alleleSet, objective, NULL);
	genome->initializer(GA1DArrayAlleleGenome<int>::UniformInitializer);
	ga = new gmxga::RandomImmigrantSGMXGA(*genome);
	ga->setGenomePrototype(genome);
	ga->minimaxi(GAGeneticAlgorithm::MINIMIZE);

	this->parameters = parameters;
	//SET GA PARAMETERS
	ga->terminator(GAGeneticAlgorithm::TerminateUponGeneration);

	if (parameters.populationSize % 2 == 1)
		parameters.populationSize++;
	ga->populationSize(parameters.populationSize);
	if(parameters.selector == 0) {
		ga->selector(GATournamentSelector());
	}
	ga->elitist(parameters.elitist?GABoolean::gaTrue:GABoolean::gaFalse);
	ga->setPeriod(parameters.period);
	ga->setNumInitGen(parameters.numInitGen);
	ga->setNumImmigrateGen(parameters.numImmigrateGen);
	ga->setImmigrationRate(parameters.immigrationRate);
	if(parameters.crossoverType==0) {
		UniformRankedGroupMaskingCrossover* crossover = new UniformRankedGroupMaskingCrossover();
		crossover->setParams(parameters.groupProportion, parameters.groupProportionMargin,
				parameters.individualProportion, parameters.individualProportionMargin, parameters.individualScoreLimit);
		crossover->setGA(ga);
		std::function<UniformRankedGroupMaskingCrossover::AlleleEvaluator>* evaluator = new std::function<UniformRankedGroupMaskingCrossover::AlleleEvaluator>(AlleleEvaluator());
		crossover->setEvaluator(evaluator);
		crossover->setGrouper(BoundedDiameterMinCutClusterGrouper(1));
		ga->setCrossover(crossover);
	}
	ga->set(gaNpCrossover, parameters.crossoverRate);
	ga->set(gaNpMutation, parameters.mutationRate);
}

void GADPLP::updateProblem(labelplacement::ConflictGraph& conflictGraph) {
	/*for(int i=0; i<ga->population().size(); i++) {
		Objective* objectiveOfIndividual = ga->population().individual(i).evaluator().target<Objective>();
		objectiveOfIndividual->conflictGraph = conflictGraph;
	}*/

	*objective.conflictGraph=&conflictGraph;
	if(parameters.crossoverType==0) {
		ga->getCrossover()->initialize(conflictGraph.getConflictGraphOfPoints());
		GMXDelegate* gmxDelegate = ga->sexual().target<GMXDelegate>();
		UniformRankedGroupMaskingCrossover* uniformRankedGroupMaskingCrossover=(UniformRankedGroupMaskingCrossover*)gmxDelegate->gmx;
		((AlleleEvaluator*)uniformRankedGroupMaskingCrossover->getEvaluator())->conflictGraph = &conflictGraph;
	}
}

labelplacement::Solution& GADPLP::buildSolution(labelplacement::ConflictGraph& conflictGraph, const GA1DArrayAlleleGenome<unsigned int> & individual) {
	labelplacement::Solution* solution = new labelplacement::Solution(&conflictGraph);
	labelplacement::Solution& solutionRef = *solution;
	int* labelPlacements = new int[individual.length()];
	solutionRef.setLabelPlacements(labelPlacements);
	for (int i = 0; i < individual.length(); i++) {
		labelPlacements[i] = individual.gene(i);
	}
	return solutionRef;
}

labelplacement::Solution& GADPLP::optimize(labelplacement::ConflictGraph& conflictGraph) {
	updateProblem(conflictGraph);
	ga->evolve();
	const GA1DArrayAlleleGenome<unsigned int> & bestIndividual =
				(GA1DArrayAlleleGenome<unsigned int> &) ga->statistics().bestIndividual();
	return buildSolution(conflictGraph, bestIndividual);
}

} /* namespace GADPLP */
