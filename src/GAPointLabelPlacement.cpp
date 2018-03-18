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

#include "BoundDiameterMinCutClustering.h"

using namespace std;
using namespace graph;
using namespace clustering;

namespace labelplacement {

GAPointLabelPlacement::GAPointLabelPlacement() {
	conflictGraph = __null;
}

GAPointLabelPlacement::~GAPointLabelPlacement() {
}

ConflictGraph* GAPointLabelPlacement::conflictGraph = NULL;
bool* GAPointLabelPlacement::mask = NULL;
double GAPointLabelPlacement::initialSelectedGroupPortion = 0.1;
double GAPointLabelPlacement::addedSelectedGroupPortion = 0.8;
int GAPointLabelPlacement::alleleNumber = 0;
int* GAPointLabelPlacement::alleleGroups = NULL;
int GAPointLabelPlacement::numberOfGroups = 0;
vector<int>* GAPointLabelPlacement::groups = NULL;
int GAPointLabelPlacement::lastGeneration = -1;
vector<int>** GAPointLabelPlacement::selectedGroups = NULL;
int* GAPointLabelPlacement::selectedGroupNumbers = NULL;
int GAPointLabelPlacement::numberOfSelectedGroups = 0;
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
vector<int>* GAPointLabelPlacement::getGroups() {
	return groups;
}
int GAPointLabelPlacement::getNumberOfGroups() {
	return GAPointLabelPlacement::numberOfGroups;
}
int GAPointLabelPlacement::getNumberOfSelectedGroups() {
	return GAPointLabelPlacement::numberOfSelectedGroups;
}
void GAPointLabelPlacement::initialize(ConflictGraph& conflictGraph) {
	//cout<<"GAPointLabelPlacement::initialize>start"<<endl;
	GAPointLabelPlacement::conflictGraph = &conflictGraph;
	alleleNumber = conflictGraph.getConflictGraphOfPoints()->getVertexNumber();
	if(alleleGroups!=NULL) {
		delete alleleGroups;
	}
	alleleGroups = new int[alleleNumber];
	if(groups!=NULL) {
		delete [] groups;
	}
	BoundDiameterMinCutClustering clustering;
	//cout<<"GAPointLabelPlacement::initialize>>clustering start"<<endl;
	vector<Graph*>* clusters = clustering.cluster(conflictGraph.getConflictGraphOfPoints(),2);
	numberOfGroups = clusters->size();
	cout<<"GAPointLabelPlacement::initialize>>number of groups: "<<numberOfGroups<<endl;
	groups = new vector<int>[numberOfGroups];
	for(int i=0; i<numberOfGroups; i++) {
		Graph* cluster = clusters->at(i);
		groups[i] = vector<int>(*cluster->getVertices());
		delete cluster;
	}
	delete clusters;
	//cout<<"GAPointLabelPlacement::initialize>end"<<endl;
}
int* GAPointLabelPlacement::getAlleleGroups() {
	//cout<<"GAPointLabelPlacement::getAlleleGroups>start"<<endl;
	if(lastGeneration<ga->generation()) {
		double currentSelectedGroupPortion = initialSelectedGroupPortion + addedSelectedGroupPortion * (double)ga->generation()/(double)ga->nGenerations();
		numberOfSelectedGroups = numberOfGroups * currentSelectedGroupPortion;
		//cout<<"GAPointLabelPlacement::getAlleleGroups>num selected groups: "<<numberOfSelectedGroups<<endl;
	}
	selectedGroups = new vector<int>*[numberOfSelectedGroups];
	selectedGroupNumbers = new int[numberOfSelectedGroups];
	vector<int> groupNumbersToSelect;
	for(int i=0;i<numberOfGroups;i++){
		groupNumbersToSelect.push_back(i);
	}
	//cout<<"GAPointLabelPlacement::getAlleleGroups>determine random group numbers"<<endl;
	for(int i=0; i<numberOfSelectedGroups; i++) {
		int groupNumberToSelectIndex = GARandomInt(0, groupNumbersToSelect.size()-1);
		int groupNumberToSelect = groupNumbersToSelect[groupNumberToSelectIndex];
		selectedGroupNumbers[i] = groupNumberToSelect;
		groupNumbersToSelect.erase(groupNumbersToSelect.begin()+groupNumberToSelectIndex);
	}

	for(int i=0; i<alleleNumber; i++) {
		alleleGroups[i] = 0;
	}
	//cout<<"GAPointLabelPlacement::getAlleleGroups>assign selected groups"<<endl;
	for(int i=0; i<numberOfSelectedGroups; i++) {
		int selectedGroupNumber = selectedGroupNumbers[i];
		selectedGroups[i] = &getGroups()[selectedGroupNumber];
		vector<int>* selectedGroup = selectedGroups[i];
		for(int j=0; j<selectedGroup->size(); j++) {
			int alleleNumber =  selectedGroup->at(j);
			alleleGroups[alleleNumber] = i;
		}
	}
	delete [] selectedGroups;
	delete [] selectedGroupNumbers;
	selectedGroups = NULL;
	//cout<<"GAPointLabelPlacement::getAlleleGroups>end"<<endl;
	lastGeneration=ga->generation();
	return alleleGroups;
}

float GAPointLabelPlacement::objective(GAGenome & g) {
	GA1DArrayAlleleGenome<int> & genome = (GA1DArrayAlleleGenome<int> &) g;
	float objectiveValue;
	float potentialConflictSize = 0;
	float conflictSize = 0;
	int pointNumber = conflictGraph->getConflictGraphOfPoints()->getVertexNumber();
	int positionNumber = conflictGraph->getPositionNumber();
	for(int i=0; i<genome.length(); i++)
	{
			int position1 = genome.gene(i);
			int posIx1 = i*positionNumber+position1;

			vector<int>* conflictingPoints = conflictGraph->getConflictGraphOfPoints()->getAdjacencyList()+i;
			for(vector<int>::iterator it=conflictingPoints->begin(); it!=conflictingPoints->end(); it++)
			{
				int conflictingPointIx = *it;
				int position2 = genome.gene(conflictingPointIx);
				int posIx2 = conflictingPointIx*positionNumber+position2;
				if(conflictGraph->getConflictGraphOfPositions()->getAdjacencyMatrix()[posIx1][posIx2])
				{
					conflictSize++;
					break;
				}
			}

			potentialConflictSize += conflictGraph->getConflictGraphOfPositions()->getAdjacencyList()[posIx1].size();
	}
	objectiveValue = conflictSize + 0.25*potentialConflictSize;
	return conflictSize;
}
void GAPointLabelPlacement::generateMask(int depth, double effectiveCoverage) {
	generateGroupedMask(depth, effectiveCoverage, -1);
	mask = new bool[alleleNumber];
	for(int i=0; i<alleleNumber;i++) {
		mask[i] = (alleleGroups[i]==0?false:true);
	}
}
void GAPointLabelPlacement::generateGroupedMask(int depth, double effectiveCoverage, int branchingLimit) {
	const int pointNumber = alleleNumber;
	const int effectiveMaskSize = pointNumber*effectiveCoverage;
	alleleGroups = new int[pointNumber];
	for(int i=0; i<pointNumber; i++) {
		alleleGroups[i] = 0;
	}
	if(depth>-1) {
		vector<int> pointsToBeEffectivelyMaskedLimited;
		vector<int> pointsToBeEffectivelyMaskedOverLimit;
		for(int i=0; i<pointNumber; i++)
		{
			if(branchingLimit == 0)
			{
				pointsToBeEffectivelyMaskedLimited.push_back(i);
			}
			else
			{
				int pointConflictSize = conflictGraph->getConflictGraphOfPoints()->getAdjacencyList()[i].size();
				if(pointConflictSize > branchingLimit)
				{
					pointsToBeEffectivelyMaskedOverLimit.push_back(i);
				}
				else
				{
					pointsToBeEffectivelyMaskedLimited.push_back(i);
				}
			}
		}
		vector<int> pointsEffectivelyMasked;
		bool pointsEffectivelyMaskedIndicator[pointNumber];
		for(int i=0;i<pointNumber; i++)
		{
			pointsEffectivelyMaskedIndicator[i] = false;
		}
		int lastDepthCurStart = -1;
		int lastDepthCurEnd = -1;
		int groupNumber = 0;
		for(int curEffectiveMaskSize = 0; curEffectiveMaskSize<effectiveMaskSize;)
		{
			groupNumber++;
			vector<int>* pointsToBeEffectivelyMasked;
			if(!pointsToBeEffectivelyMaskedLimited.empty())
			{
				pointsToBeEffectivelyMasked = &pointsToBeEffectivelyMaskedLimited;
			}
			else
			{
				pointsToBeEffectivelyMasked = &pointsToBeEffectivelyMaskedOverLimit;
			}
			int pointToMaskIx = GARandomInt(0,pointsToBeEffectivelyMasked->size()-1);
			int pointToMask = pointsToBeEffectivelyMasked->at(pointToMaskIx);

			pointsToBeEffectivelyMasked->erase(pointsToBeEffectivelyMasked->begin()+pointToMaskIx);
			pointsEffectivelyMasked.push_back(pointToMask);
			pointsEffectivelyMaskedIndicator[pointToMask] = true;
			alleleGroups[pointToMask]=groupNumber;
			curEffectiveMaskSize++;

			lastDepthCurStart=lastDepthCurEnd=pointsEffectivelyMasked.size()-1;
			for(int curDepth=0; curDepth<depth && curEffectiveMaskSize<effectiveMaskSize; curDepth++) {
				for(;lastDepthCurStart<=lastDepthCurEnd && curEffectiveMaskSize<effectiveMaskSize; lastDepthCurStart++) {
					int pointToFindRivals = pointsEffectivelyMasked[lastDepthCurStart];
					vector<int> conflictingPoints = conflictGraph->getConflictGraphOfPoints()->getAdjacencyList()[pointToFindRivals];
					if(conflictingPoints.size()>branchingLimit && !pointsToBeEffectivelyMaskedLimited.empty())
					{
						continue;
					}
					for(vector<int>::iterator it=conflictingPoints.begin(); it!=conflictingPoints.end() && curEffectiveMaskSize<effectiveMaskSize; it++) {
						int conflictingPoint = *it;
						if(!pointsEffectivelyMaskedIndicator[conflictingPoint]) {
							for(vector<int>::iterator it2=pointsToBeEffectivelyMasked->begin(); it2!=pointsToBeEffectivelyMasked->end(); it2++)
							{
								if(*it2==conflictingPoint)
								{
									pointsToBeEffectivelyMasked->erase(it2);
									break;
								}
							}
							pointsEffectivelyMasked.push_back(conflictingPoint);
							pointsEffectivelyMaskedIndicator[conflictingPoint] = true;
							alleleGroups[conflictingPoint]=groupNumber;
							curEffectiveMaskSize++;
						}
					}
				}
				lastDepthCurEnd = pointsEffectivelyMasked.size()-1;
			}
		}
		numberOfGroups = groupNumber;
	}
}
int GAPointLabelPlacement::uniformCrossoverWithMasking(const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2)
{
	double maskCoverage = 0.3 + 0.7 * (double)ga->generation()/(double)ga->nGenerations();

	if(lastGeneration<ga->generation()) {
		if(mask != NULL)
			delete [] mask;
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
int GAPointLabelPlacement::uniformCrossoverWithRandomPartialGroupedMasking(const GAGenome& p1, const GAGenome& p2,
			      GAGenome* c1, GAGenome* c2) {
	int* alleleGroups = getAlleleGroups();
	int alleleGroupParents[numberOfSelectedGroups];
	for(int i=0;i<numberOfSelectedGroups; i++) {
		alleleGroupParents[i] = GARandomBit();
	}

	const GA1DArrayGenome<int> &mom=DYN_CAST(const GA1DArrayGenome<int> &, p1);
	const GA1DArrayGenome<int> &dad=DYN_CAST(const GA1DArrayGenome<int> &, p2);

	GA1DArrayGenome<int> &sis=DYN_CAST(GA1DArrayGenome<int> &, *c1);
	GA1DArrayGenome<int> &bro=DYN_CAST(GA1DArrayGenome<int> &, *c2);

	for(int i=mom.length()-1; i>=0; i--){
		if(alleleGroups[i]!=0){
			int maskGroupParent = alleleGroupParents[alleleGroups[i]];
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
	initialize(conflictGraph);
	int pointNumber = conflictGraph.getConflictGraphOfPoints()->getVertexNumber();
	int positionsPerPointNumber = conflictGraph.getPositionNumber();

	int* alleleSet = new int[positionsPerPointNumber];
	for(int i=0;i<positionsPerPointNumber; i++)
		alleleSet[i] = i;
	GAAlleleSet<int> positionAlleleSet(positionsPerPointNumber, alleleSet);
	GA1DArrayAlleleGenome<int> genome(pointNumber, positionAlleleSet, objective, NULL);
	genome.initializer(GA1DArrayAlleleGenome<int>::UniformInitializer);
	//genome.crossover(GA1DArrayAlleleGenome<int>::UniformCrossover);
	genome.crossover(GAPointLabelPlacement::uniformCrossoverWithRandomPartialGroupedMasking);
	ga=new GASimpleGA(genome);
	ga->selector(GARankSelector());
	GAParameterList params;
	GASimpleGA::registerDefaultParameters(params);
	params.set(gaNscoreFilename, "bog.dat");
	params.set(gaNflushFrequency, 50);
	params.set(gaNpMutation, 0.001);
	params.set(gaNpCrossover, 0.90);

	int populationSize = pointNumber/10;
	if(populationSize%2==1)
		populationSize++;
	ga->populationSize(populationSize);
	ga->minimaxi(GAGeneticAlgorithm::MINIMIZE);
	ga->nGenerations(15*pointNumber);
	ga->elitist(GABoolean::gaTrue);
	ga->terminator(GAGeneticAlgorithm::TerminateUponGeneration);
	ga->evolve();

	const GA1DArrayAlleleGenome<unsigned int> & bestIndividual = (GA1DArrayAlleleGenome<unsigned int> &) ga->statistics().bestIndividual();
	Solution* solution = new Solution(this->conflictGraph);
	Solution& solutionRef = *solution;
	int* labelPlacements = new int[pointNumber];
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
