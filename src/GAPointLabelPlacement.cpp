/*
 * PointLabelingGA.cpp
 *
 *  Created on: Dec 18, 2017
 *      Author: oakile
 */
#include "GAPointLabelPlacement.h"

#include <cstdlib>
#include <iostream>
#include <vector>

#include <ga/GASimpleGA.h>
#include <ga/GAPopulation.h>

using namespace std;

namespace labelplacement {

GAPointLabelPlacement::GAPointLabelPlacement() {
	conflictGraph = __null;
}

GAPointLabelPlacement::~GAPointLabelPlacement() {
	// TODO Auto-generated destructor stub
}

ConflictGraph* GAPointLabelPlacement::conflictGraph = NULL;
bool* GAPointLabelPlacement::mask = NULL;
int* GAPointLabelPlacement::groupedMask = NULL;
int GAPointLabelPlacement::numberOfMaskGroupes = 0;
int GAPointLabelPlacement::lastGeneration = -1;
GASimpleGA* GAPointLabelPlacement::ga = NULL;

ConflictGraph* GAPointLabelPlacement::getConflictGraph() {
	return conflictGraph;
}
void GAPointLabelPlacement::setConflictGraph(ConflictGraph* conflictGraph) {
	GAPointLabelPlacement::conflictGraph = conflictGraph;
}
bool* GAPointLabelPlacement::getMask() {
	return GAPointLabelPlacement::mask;
}
int* GAPointLabelPlacement::getGroupedMask() {
	return GAPointLabelPlacement::groupedMask;
}
int GAPointLabelPlacement::getNumberOfMaskGroupes() {
	return GAPointLabelPlacement::numberOfMaskGroupes;
}
float GAPointLabelPlacement::objective(GAGenome & g) {
	GA1DArrayAlleleGenome<int> & genome = (GA1DArrayAlleleGenome<int> &) g;
	float objectiveValue;
	float potentialConflictSize = 0;
	float conflictSize = 0;
	for(int i=0; i<genome.length(); i++)
	{
			int position1 = genome.gene(i);
			int posIx1 = i*conflictGraph->getPositionNumber()+position1;

			vector<int>* conflictingPoints = conflictGraph->getConflictingPoints()+i;
			for(vector<int>::iterator it=conflictingPoints->begin(); it!=conflictingPoints->end(); it++)
			{
				int conflictingPointIx = *it;
				int position2 = genome.gene(conflictingPointIx);
				int posIx2 = conflictingPointIx*conflictGraph->getPositionNumber()+position2;
				if(conflictGraph->getAdjacencyMatrix()[posIx1][posIx2])
				{
					conflictSize++;
					break;
				}
			}

			potentialConflictSize += conflictGraph->getConflictingPositions()[posIx1].size();
	}
	objectiveValue = conflictSize + 0.25*potentialConflictSize;
	return conflictSize;
}
void GAPointLabelPlacement::generateMask(int depth, double effectiveCoverage) {
	generateGroupedMask(depth, effectiveCoverage);
	if(mask != NULL)
		delete [] mask;
	mask = new bool[conflictGraph->getPointNumber()];
	for(int i=0; i<conflictGraph->getPointNumber();i++) {
		mask[i] = (groupedMask[i]==0?false:true);
	}
}
void GAPointLabelPlacement::generateGroupedMask(int depth, double effectiveCoverage) {
	const int pointNumber = conflictGraph->getPointNumber();
	const int effectiveMaskSize = pointNumber*effectiveCoverage;
	if(groupedMask != NULL)
			delete [] groupedMask;
	groupedMask = new int[pointNumber];
	for(int i=0; i<pointNumber; i++) {
		groupedMask[i] = 0;
	}
	if(depth>-1) {
		vector<int> pointsNonEffectivelyMasked;
		for(int i=0; i<pointNumber; i++) {
			pointsNonEffectivelyMasked.push_back(i);
		}
		vector<int> pointsEffectivelyMasked;
		int lastDepthCurStart = -1;
		int lastDepthCurEnd = -1;
		int groupNumber = 0;
		for(int curEffectiveMaskSize = 0; curEffectiveMaskSize<effectiveMaskSize;) {
			groupNumber++;
			int pointToMaskIx = GARandomInt(0,pointsNonEffectivelyMasked.size()-1);
			int pointToMask = pointsNonEffectivelyMasked[pointToMaskIx];
			pointsNonEffectivelyMasked.erase(pointsNonEffectivelyMasked.begin()+pointToMaskIx);
			pointsEffectivelyMasked.push_back(pointToMask);
			groupedMask[pointToMask]=groupNumber;
			curEffectiveMaskSize++;
			lastDepthCurStart=lastDepthCurEnd=pointsEffectivelyMasked.size()-1;
			for(int curDepth=0; curDepth<depth && curEffectiveMaskSize<effectiveMaskSize; curDepth++) {
				for(;lastDepthCurStart<=lastDepthCurEnd && curEffectiveMaskSize<effectiveMaskSize; lastDepthCurStart++) {
					int pointToFindRivals = pointsEffectivelyMasked[lastDepthCurStart];
					vector<int> conflictingPositions = conflictGraph->getConflictingPositions()[pointToFindRivals];
					for(vector<int>::iterator it=conflictingPositions.begin(); it!=conflictingPositions.end() && curEffectiveMaskSize<effectiveMaskSize; it++) {
						int conflictingPosition = *it;
						int rivalPoint = conflictingPosition/conflictGraph->getPositionNumber();
						bool existInPointsEffectivelyMasked = false;
						for(vector<int>::iterator it2=pointsEffectivelyMasked.begin(); it2!=pointsEffectivelyMasked.end(); it2++) {
							if(*it2 == rivalPoint) {
								existInPointsEffectivelyMasked = true;
								break;
							}
						}
						if(!existInPointsEffectivelyMasked) {
							pointsEffectivelyMasked.push_back(rivalPoint);
							groupedMask[rivalPoint]=groupNumber;
							curEffectiveMaskSize++;
						}
					}
				}
				lastDepthCurEnd = pointsEffectivelyMasked.size()-1;
			}
		}
		numberOfMaskGroupes = groupNumber;
	}
}
int GAPointLabelPlacement::uniformCrossoverWithMasking(const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2)
{
	double maskCoverage = 0.3 + 0.5 * (double)ga->generation()/(double)ga->nGenerations();

	if(lastGeneration<ga->generation()) {
		generateMask(2, maskCoverage);
		lastGeneration=ga->generation();
	}
	/*for(int i=0; i<conflictGraph->getPointNumber();i++) {
		cout<<mask[i]?1:0;
	}
	cout<<endl;*/
	const GA1DArrayGenome<int> &mom=DYN_CAST(const GA1DArrayGenome<int> &, p1);
	const GA1DArrayGenome<int> &dad=DYN_CAST(const GA1DArrayGenome<int> &, p2);

	GA1DArrayGenome<int> &sis=DYN_CAST(GA1DArrayGenome<int> &, *c1);
	GA1DArrayGenome<int> &bro=DYN_CAST(GA1DArrayGenome<int> &, *c2);

	for(int i=mom.length()-1; i>=0; i--){
		if(mask[i]){
			  sis.gene(i, mom.gene(i));
			  bro.gene(i, dad.gene(i));
		} else {
			if(GARandomBit()){
			  sis.gene(i, mom.gene(i));
			  bro.gene(i, dad.gene(i));
			} else {
			  sis.gene(i, dad.gene(i));
			  bro.gene(i, mom.gene(i));
			}
		}
	}
	return 2;
}
int GAPointLabelPlacement::uniformCrossoverWithGroupedMasking(const GAGenome& p1, const GAGenome& p2,
			      GAGenome* c1, GAGenome* c2) {
	int maskGroupParents[numberOfMaskGroupes];
	for(int i=0;i<numberOfMaskGroupes; i++) {
		maskGroupParents[i] = GARandomBit();
	}
	double maskCoverage = 0.1 + 0.8 * (double)ga->generation()/(double)ga->nGenerations();

	if(lastGeneration<ga->generation()) {
		generateGroupedMask(1, 0.7);
		lastGeneration=ga->generation();
	}

	const GA1DArrayGenome<int> &mom=DYN_CAST(const GA1DArrayGenome<int> &, p1);
	const GA1DArrayGenome<int> &dad=DYN_CAST(const GA1DArrayGenome<int> &, p2);

	GA1DArrayGenome<int> &sis=DYN_CAST(GA1DArrayGenome<int> &, *c1);
	GA1DArrayGenome<int> &bro=DYN_CAST(GA1DArrayGenome<int> &, *c2);

	for(int i=mom.length()-1; i>=0; i--){
		if(groupedMask[i]!=0){
			int maskGroupParent = maskGroupParents[groupedMask[i]-1];
			if(maskGroupParent==1) {
				sis.gene(i, mom.gene(i));
				bro.gene(i, dad.gene(i));
			}
			else {
				  sis.gene(i, dad.gene(i));
				  bro.gene(i, mom.gene(i));
			}
		} else {
			if(GARandomBit()){
			  sis.gene(i, mom.gene(i));
			  bro.gene(i, dad.gene(i));
			} else {
			  sis.gene(i, dad.gene(i));
			  bro.gene(i, mom.gene(i));
			}
		}
	}
	return 2;
}
Solution& GAPointLabelPlacement::optimize(ConflictGraph& conflictGraph)
{
	this->conflictGraph = &conflictGraph;
	int* alleleSet = new int[conflictGraph.getPositionNumber()];
	for(int i=0;i<conflictGraph.getPositionNumber(); i++)
		alleleSet[i] = i;
	GAAlleleSet<int> positionAlleleSet(conflictGraph.getPositionNumber(), alleleSet);
	GA1DArrayAlleleGenome<int> genome(conflictGraph.getPointNumber(), positionAlleleSet, objective, NULL);
	genome.initializer(GA1DArrayAlleleGenome<int>::UniformInitializer);
	//genome.crossover(GA1DArrayAlleleGenome<int>::UniformCrossover);
	genome.crossover(GAPointLabelPlacement::uniformCrossoverWithGroupedMasking);
	ga=new GASimpleGA(genome);
	ga->selector(GARankSelector());
	GAParameterList params;
	GASimpleGA::registerDefaultParameters(params);
	params.set(gaNscoreFilename, "bog.dat");
	params.set(gaNflushFrequency, 50);
	params.set(gaNpMutation, 0.001);
	params.set(gaNpCrossover, 0.90);

	int pointNumber = conflictGraph.getPointNumber();
	int populationSize = pointNumber/10;
	if(populationSize%2==1)
		populationSize++;
	ga->populationSize(populationSize);
	ga->minimaxi(GAGeneticAlgorithm::MINIMIZE);
	ga->nGenerations(pointNumber*3);
	ga->elitist(GABoolean::gaTrue);
	ga->terminator(GAGeneticAlgorithm::TerminateUponGeneration);
	ga->evolve();

	const GA1DArrayAlleleGenome<unsigned int> & bestIndividual = (GA1DArrayAlleleGenome<unsigned int> &) ga->statistics().bestIndividual();
	Solution* solution = new Solution(this->conflictGraph);
	Solution& solutionRef = *solution;
	int* labelPlacements = new int[conflictGraph.getPointNumber()];
	solution->setLabelPlacements(labelPlacements);
	for(int i=0; i<bestIndividual.length(); i++)
	{
		labelPlacements[i] = bestIndividual.gene(i);
	}
	lastGeneration=-1;
	delete alleleSet;
	delete ga;
	return solutionRef;
}

} /* namespace labelplacement */
