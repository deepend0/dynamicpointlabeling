/*
 * PointLabelingGA.cpp
 *
 *  Created on: Dec 18, 2017
 *      Author: oakile
 */
#include "GAPointLabelPlacement.h"

#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <vector>

#include <ga/GAPopulation.h>
#include <ga/GASimpleGA.h>
#include <ga/GAStatistics.h>

#include "BoundDiameterMinCutClustering.h"

using namespace std;
using namespace graph;
using namespace clustering;

namespace labelplacement {

struct ComparePairIntDouble {
	bool operator()(pair<int, double> p1, pair<int, double> p2) {
		return p1.second < p2.second;
	}
};
struct ComparePairIntDoubleMax {
	bool operator()(pair<int, double> p1, pair<int, double> p2) {
		return p1.second > p2.second;
	}
};
struct CompareVerticesAccToObjectivePart {
	const GA1DArrayGenome<int>* parent;
	bool operator()(int a1, int a2) {
		return GAPointLabelPlacement::scoreOfAnAllele(*parent,a1) < GAPointLabelPlacement::scoreOfAnAllele(*parent, a2);
	}
};
struct CompareVertexGroupsByAvgIncidence {
	Graph* graph;
	bool operator()(vector<int> g1, vector<int> g2) {
		double g1Score = 0;
		for(int i=0; i<g1.size(); i++) {
			int v = g1[i];
			g1Score += graph->getAdjacencyList()[v].size();
		}
		g1Score /= (double)g1.size();
		double g2Score = 0;
		for(int i=0; i<g2.size(); i++) {
			int v = g2[i];
			g2Score += graph->getAdjacencyList()[v].size();
		}
		g2Score /= (double)g2.size();
		return g1Score < g2Score;
	}
};
GAPointLabelPlacement::GAPointLabelPlacement() {
	conflictGraph = __null;
}

GAPointLabelPlacement::~GAPointLabelPlacement() {
}
void NoInitializer(GAGenome & c){
}

ConflictGraph* GAPointLabelPlacement::conflictGraph = NULL;
bool* GAPointLabelPlacement::mask = NULL;
double GAPointLabelPlacement::initialSelectedGroupPortion = 0.00;
double GAPointLabelPlacement::addedSelectedGroupPortion = 0.00;
double GAPointLabelPlacement::selectedIndividualPortion = 0.30;
double GAPointLabelPlacement::addedSelectedIndividualPortion = 0.0;
int GAPointLabelPlacement::alleleNumber = 0;
int* GAPointLabelPlacement::alleleGroups = NULL;
int* GAPointLabelPlacement::allelePointNumbers = NULL;
unordered_map<int, int>* GAPointLabelPlacement::alleleNumberIndexMap = NULL;
int GAPointLabelPlacement::numberOfGroups = 0;
vector<int>* GAPointLabelPlacement::groups = NULL;
int GAPointLabelPlacement::lastGeneration = -1;
int GAPointLabelPlacement::numberOfSelectedGroups = 0;
GASimpleGA* GAPointLabelPlacement::ga = NULL;
GAPopulation* GAPointLabelPlacement::lastPopulation = NULL;

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
void GAPointLabelPlacement::initialize(ConflictGraph& conflictGraph,
		bool sequenceGroups, int clusteringType) {
	//cout<<"GAPointLabelPlacement::initialize>start"<<endl;
	alleleNumber = conflictGraph.getConflictGraphOfPoints()->getVertexNumber();
	if (alleleGroups != NULL) {
		delete alleleGroups;
	}
	alleleGroups = new int[alleleNumber];
	if (groups != NULL) {
		delete[] groups;
	}
	if (clusteringType == 0) {
		BoundDiameterMinCutClustering clustering;
		//cout<<"GAPointLabelPlacement::initialize>>clustering start"<<endl;
		vector<Graph*>* clusters = clustering.cluster(
				conflictGraph.getConflictGraphOfPoints(), 2);
		numberOfGroups = clusters->size();
		//cout<<"GAPointLabelPlacement::initialize>>number of groups: "<<numberOfGroups<<endl;
		groups = new vector<int> [numberOfGroups];
		for (int i = 0; i < numberOfGroups; i++) {
			Graph* cluster = clusters->at(i);
			groups[i] = vector<int>(*cluster->getVertices());
			delete cluster;
		}
		delete clusters;
	} else if (clusteringType == 1) {
		vector<vector<int>*>* clusters = clusterPivotal(2, 10);
		numberOfGroups = clusters->size();
		groups = new vector<int> [numberOfGroups];
		for (int i = 0; i < numberOfGroups; i++) {
			groups[i] = vector<int>(*clusters->at(i));
			delete clusters->at(i);
		}
		delete clusters;
	} else if(clusteringType == -1 ) {
		numberOfGroups = 1;
		groups = new vector<int> [numberOfGroups];
		for(int i=0;i<alleleNumber; i++) {
			groups[0].push_back(i);
		}
	}

	CompareVertexGroupsByAvgIncidence compare;
	compare.graph = conflictGraph.getConflictGraphOfPoints();
	sort(groups, groups+numberOfGroups, compare);
//	cout<<"Groups ranking: ";
//	for(int i=0;i<numberOfGroups; i++) {
//		vector<int>* g1 = groups + i;
//		double g1Score = 0;
//		for(int i=0; i<g1->size(); i++) {
//			int v = g1->at(i);
//			g1Score += conflictGraph.getConflictGraphOfPoints()->getAdjacencyList()[v].size();
//		}
//		g1Score /= (double)g1->size();
//		cout<<g1Score<<" ";
//	}
//	cout<<endl;
	if (sequenceGroups) {
		allelePointNumbers = new int[alleleNumber];
		alleleNumberIndexMap = new unordered_map<int, int>();
		int ix = 0;
		for (int i = 0; i < numberOfGroups; i++) {
			vector<int>* group = groups + i;
			for (vector<int>::iterator it = group->begin(); it != group->end();
					it++) {
				int pointNumber = *it;
				allelePointNumbers[ix] = pointNumber;
				alleleNumberIndexMap->insert(pair<int, int>(pointNumber, ix));
				ix++;
			}
		}
	}
	//cout<<"GAPointLabelPlacement::initialize>end"<<endl;
}
int* GAPointLabelPlacement::getRandomPartialAlleleGroups() {
	//cout<<"GAPointLabelPlacement::getRandomPartialAlleleGroups>start"<<endl;
	if (lastGeneration < ga->generation()) {
		double currentSelectedGroupPortion = initialSelectedGroupPortion
				+ addedSelectedGroupPortion * (double) ga->generation()
						/ (double) ga->nGenerations();
		numberOfSelectedGroups = numberOfGroups * currentSelectedGroupPortion;
		//cout<<"GAPointLabelPlacement::getRandomPartialAlleleGroups>num selected groups: "<<numberOfSelectedGroups<<endl;
	}
	vector<int>** selectedGroups = new vector<int>*[numberOfSelectedGroups];
	int* selectedGroupNumbers = new int[numberOfSelectedGroups];
	vector<int> groupNumbersToSelect;
	for (int i = 0; i < numberOfGroups; i++) {
		groupNumbersToSelect.push_back(i);
	}
	//cout<<"GAPointLabelPlacement::getRandomPartialAlleleGroups>determine random group numbers"<<endl;
	for (int i = 0; i < numberOfSelectedGroups; i++) {
		int groupNumberToSelectIndex = GARandomInt(0,
				groupNumbersToSelect.size() - 1);
		int groupNumberToSelect = groupNumbersToSelect[groupNumberToSelectIndex];
		selectedGroupNumbers[i] = groupNumberToSelect;
		groupNumbersToSelect.erase(
				groupNumbersToSelect.begin() + groupNumberToSelectIndex);
	}

	for (int i = 0; i < alleleNumber; i++) {
		alleleGroups[i] = 0;
	}
	//cout<<"GAPointLabelPlacement::getRandomPartialAlleleGroups>assign selected groups"<<endl;
	for (int i = 0; i < numberOfSelectedGroups; i++) {
		int selectedGroupNumber = selectedGroupNumbers[i];
		selectedGroups[i] = &groups[selectedGroupNumber];
		vector<int>* selectedGroup = selectedGroups[i];
		for (int j = 0; j < selectedGroup->size(); j++) {
			int alleleNumber = selectedGroup->at(j);
			alleleGroups[alleleNumber] = i;
		}
	}
	delete[] selectedGroups;
	delete[] selectedGroupNumbers;
	selectedGroups = NULL;
	//cout<<"GAPointLabelPlacement::getRandomPartialAlleleGroups>end"<<endl;
	lastGeneration = ga->generation();
	return alleleGroups;
}

int* GAPointLabelPlacement::getProblemRankedPartialAlleleGroups() {
	//cout<<"GAPointLabelPlacement::getProblemRankedPartialAlleleGroups>start"<<endl;
	if (lastGeneration < ga->generation()) {
		double currentSelectedGroupPortion = initialSelectedGroupPortion
				+ addedSelectedGroupPortion * (double) ga->generation()
						/ (double) ga->nGenerations();
		numberOfSelectedGroups = numberOfGroups * currentSelectedGroupPortion;
	}
	for (int i = 0; i < alleleNumber; i++) {
		alleleGroups[i] = -1;
	}
	//cout<<"GAPointLabelPlacement::getProblemRankedPartialAlleleGroups>assign selected groups"<<endl;
	for (int i = 0; i < numberOfSelectedGroups; i++) {
		vector<int>* selectedGroup = groups + i;
		for (int j = 0; j < selectedGroup->size(); j++) {
			int alleleNumber = selectedGroup->at(j);
			alleleGroups[alleleNumber] = i;
		}
	}
	//cout<<"Allele Groups: ";
	//for(int i=0; i<alleleNumber; i++) {
	//	cout<<alleleGroups[i]<<" ";
	//}
	//cout<<endl;
	//cout<<"GAPointLabelPlacement::getProblemRankedPartialAlleleGroups>end"<<endl;
	lastGeneration = ga->generation();
	return alleleGroups;
}
int** GAPointLabelPlacement::getSolutionRankedPartialAlleleGroups(const GA1DArrayGenome<int>& p1,
		const GA1DArrayGenome<int>& p2) {
	//cout<<"GAPointLabelPlacement::getRankedPartialAlleleGroups>start"<<endl;
	if (lastGeneration < ga->generation()) {
		double selectedGroupPortion = initialSelectedGroupPortion
				+ addedSelectedGroupPortion * (double) ga->generation()
						/ (double) ga->nGenerations();
		numberOfSelectedGroups = numberOfGroups * selectedGroupPortion;
		//cout<<"GAPointLabelPlacement::getRankedPartialAlleleGroups>num selected groups: "<<numberOfSelectedGroups<<endl;
	}

	int **alleleGroups = new int*[2];
	alleleGroups[0] = new int[alleleNumber];
	alleleGroups[1] = new int[alleleNumber];

	for (int i = 0; i < alleleNumber; i++) {
		alleleGroups[0][i] = -1;
		alleleGroups[1][i] = -1;
	}

	//cout<<"GAPointLabelPlacement::getRankedPartialAlleleGroups>assign selected groups"<<endl;
	vector<pair<int, double>> groupScoresOfP1;
	vector<pair<int, double>> groupScoresOfP2;
	for (int i = 0; i < numberOfGroups; i++) {
		vector<int>* group = groups + i;
		double groupScoreP1 = 0;
		double groupScoreP2 = 0;
		for (int j = 0; j < group->size(); j++) {
			int groupVertex = group->at(j);
			groupScoreP1 += scoreOfAnAllele(p1, groupVertex);
			groupScoreP2 += scoreOfAnAllele(p2, groupVertex);
		}
		groupScoreP1 /= (double) group->size();
		groupScoreP2 /= (double) group->size();

		groupScoresOfP1.push_back(pair<int, double>(i, groupScoreP1));
		groupScoresOfP2.push_back(pair<int, double>(i, groupScoreP2));
	}
	double selectedGroupPortion = (double) numberOfSelectedGroups
			/ (double) numberOfGroups;
	if (selectedGroupPortion < 0.5) {
		partial_sort(groupScoresOfP1.begin(),
				groupScoresOfP1.begin() + numberOfSelectedGroups,
				groupScoresOfP1.end(), ComparePairIntDouble());
		partial_sort(groupScoresOfP2.begin(),
				groupScoresOfP2.begin() + numberOfSelectedGroups,
				groupScoresOfP2.end(), ComparePairIntDouble());
		vector<pair<int, double>>::iterator itP1 = groupScoresOfP1.begin();
		vector<pair<int, double>>::iterator itP2 = groupScoresOfP2.begin();
		for (int i = 0; i < numberOfSelectedGroups; i++) {
			int selectedGroupNumberP1 = (*itP1++).first;
			int selectedGroupNumberP2 = (*itP2++).first;
			vector<int>* selectedGroup1 = groups + selectedGroupNumberP1;
			vector<int>* selectedGroup2 = groups + selectedGroupNumberP2;

			for (int j = 0; j < selectedGroup1->size(); j++) {
				int alleleNumber = selectedGroup1->at(j);
				alleleGroups[0][alleleNumber] = selectedGroupNumberP1;
			}

			for (int j = 0; j < selectedGroup2->size(); j++) {
				int alleleNumber = selectedGroup2->at(j);
				alleleGroups[1][alleleNumber] = selectedGroupNumberP2;
			}
		}
	} else {
		int numberOfNonSelectedGroups = numberOfGroups - numberOfSelectedGroups;
		partial_sort(groupScoresOfP1.begin(),
				groupScoresOfP1.begin() + numberOfNonSelectedGroups,
				groupScoresOfP1.end(), ComparePairIntDoubleMax());
		partial_sort(groupScoresOfP2.begin(),
				groupScoresOfP2.begin() + numberOfNonSelectedGroups,
				groupScoresOfP2.end(), ComparePairIntDoubleMax());

		vector<pair<int, double>>::iterator itP1 = groupScoresOfP1.end() - 1;
		vector<pair<int, double>>::iterator itP2 = groupScoresOfP2.end() - 1;
		for (int i = 0; i < numberOfSelectedGroups; i++) {
			int selectedGroupNumberP1 = (*itP1--).first;
			int selectedGroupNumberP2 = (*itP2--).first;
			vector<int>* selectedGroup1 = groups + selectedGroupNumberP1;
			vector<int>* selectedGroup2 = groups + selectedGroupNumberP2;

			for (int j = 0; j < selectedGroup1->size(); j++) {
				int alleleNumber = selectedGroup1->at(j);
				alleleGroups[0][alleleNumber] = selectedGroupNumberP1;
			}

			for (int j = 0; j < selectedGroup2->size(); j++) {
				int alleleNumber = selectedGroup2->at(j);
				alleleGroups[1][alleleNumber] = selectedGroupNumberP2;
			}
		}
	}

//	cout<<"Allele Groups 1: ";
//	for(int i=0; i<alleleNumber; i++) {
//		cout<<alleleGroups[0][i]<<" ";
//	}
//	cout<<endl;
//	cout<<"Allele Groups 2: ";
//	for(int i=0; i<alleleNumber; i++) {
//		cout<<alleleGroups[1][i]<<" ";
//	}
//	cout<<endl;
	//cout<<"GAPointLabelPlacement::getRankedPartialAlleleGroups>end"<<endl;
	lastGeneration = ga->generation();
	return alleleGroups;
}
int** GAPointLabelPlacement::getSolutionRankedPartialAlleleGroupsAndIndividuals(const GA1DArrayGenome<int>& p1,
		const GA1DArrayGenome<int>& p2) {
	//cout<<"GAPointLabelPlacement::getRankedPartialAlleleGroups>start"<<endl;
	if (lastGeneration < ga->generation()) {
		double selectedGroupPortion = initialSelectedGroupPortion
				+ addedSelectedGroupPortion * (double) ga->generation()
						/ (double) ga->nGenerations();
		numberOfSelectedGroups = numberOfGroups * selectedGroupPortion;
		//cout<<"GAPointLabelPlacement::getRankedPartialAlleleGroups>num selected groups: "<<numberOfSelectedGroups<<endl;
	}

	int **alleleGroups = new int*[2];
	alleleGroups[0] = new int[alleleNumber];
	alleleGroups[1] = new int[alleleNumber];

	for (int i = 0; i < alleleNumber; i++) {
		alleleGroups[0][i] = -1;
		alleleGroups[1][i] = -1;
	}

	//cout<<"GAPointLabelPlacement::getRankedPartialAlleleGroups>assign selected groups"<<endl;
	vector<pair<int, double>> groupScoresOfP1;
	vector<pair<int, double>> groupScoresOfP2;
	for (int i = 0; i < numberOfGroups; i++) {
		vector<int>* group = groups + i;
		double groupScoreP1 = 0;
		double groupScoreP2 = 0;
		for (int j = 0; j < group->size(); j++) {
			int groupVertex = group->at(j);
			groupScoreP1 += scoreOfAnAllele(p1, groupVertex);
			groupScoreP2 += scoreOfAnAllele(p2, groupVertex);
		}
		groupScoreP1 /= (double) group->size();
		groupScoreP2 /= (double) group->size();

		groupScoresOfP1.push_back(pair<int, double>(i, groupScoreP1));
		groupScoresOfP2.push_back(pair<int, double>(i, groupScoreP2));
	}
	vector<int> nonGroupedAllelesP1;
	vector<int> nonGroupedAllelesP2;

	/*vector<pair<int, double>> selectedGroupScoresOfP1;
	vector<pair<int, double>> selectedGroupScoresOfP2;
	for(vector<pair<int, double>>::iterator it = groupScoresOfP1.begin(); it!= groupScoresOfP1.end(); it++) {
		pair<int, double> groupScore = *it;
		if(groupScore.second > 0.001*(double)alleleNumber) {
			vector<int>* nonSelectedGroup = groups + groupScore.first;
			for (int j = 0; j < nonSelectedGroup->size(); j++) {
				int alleleNumber = nonSelectedGroup->at(j);
				nonGroupedAllelesP1.push_back(alleleNumber);
			}
		} else {
			selectedGroupScoresOfP1.push_back(groupScore);
		}
	}
	for (vector<pair<int, double>>::iterator it = groupScoresOfP1.begin(); it!= groupScoresOfP1.end(); it++) {
		int selectedGroupNumberP1 = it->first;
		vector<int>* selectedGroup1 = groups + selectedGroupNumberP1;

		for (int j = 0; j < selectedGroup1->size(); j++) {
			int alleleNumber = selectedGroup1->at(j);
			alleleGroups[0][alleleNumber] = selectedGroupNumberP1;
		}
	}
	for(vector<pair<int, double>>::iterator it = groupScoresOfP2.begin(); it!= groupScoresOfP2.end(); it++) {
		pair<int, double> groupScore = *it;
		if(groupScore.second > 0.001*(double)alleleNumber) {
			vector<int>* nonSelectedGroup = groups + groupScore.first;
			for (int j = 0; j < nonSelectedGroup->size(); j++) {
				int alleleNumber = nonSelectedGroup->at(j);
				nonGroupedAllelesP2.push_back(alleleNumber);
			}
		}
	}
	for (vector<pair<int, double>>::iterator it = groupScoresOfP2.begin(); it!= groupScoresOfP2.end(); it++) {
		int selectedGroupNumberP2 = it->first;
		vector<int>* selectedGroup2 = groups + selectedGroupNumberP2;

		for (int j = 0; j < selectedGroup2->size(); j++) {
			int alleleNumber = selectedGroup2->at(j);
			alleleGroups[1][alleleNumber] = selectedGroupNumberP2;
		}
	}*/
	double selectedGroupPortion = (double) numberOfSelectedGroups
			/ (double) numberOfGroups;
	if (selectedGroupPortion < 0.5) {
		partial_sort(groupScoresOfP1.begin(),
				groupScoresOfP1.begin() + numberOfSelectedGroups,
				groupScoresOfP1.end(), ComparePairIntDouble());
		partial_sort(groupScoresOfP2.begin(),
				groupScoresOfP2.begin() + numberOfSelectedGroups,
				groupScoresOfP2.end(), ComparePairIntDouble());
		vector<pair<int, double>>::iterator itP1 = groupScoresOfP1.begin();
		vector<pair<int, double>>::iterator itP2 = groupScoresOfP2.begin();
		for (int i = 0; i < numberOfSelectedGroups; i++) {
			int selectedGroupNumberP1 = (*itP1++).first;
			int selectedGroupNumberP2 = (*itP2++).first;
			vector<int>* selectedGroup1 = groups + selectedGroupNumberP1;
			vector<int>* selectedGroup2 = groups + selectedGroupNumberP2;

			for (int j = 0; j < selectedGroup1->size(); j++) {
				int alleleNumber = selectedGroup1->at(j);
				alleleGroups[0][alleleNumber] = selectedGroupNumberP1;
			}

			for (int j = 0; j < selectedGroup2->size(); j++) {
				int alleleNumber = selectedGroup2->at(j);
				alleleGroups[1][alleleNumber] = selectedGroupNumberP2;
			}
		}
		for(int i = numberOfSelectedGroups; i < numberOfGroups; i++) {
			int nonselectedGroupNumberP1 = (*itP1++).first;
			int nonselectedGroupNumberP2 = (*itP2++).first;
			vector<int>* nonselectedGroup1 = groups + nonselectedGroupNumberP1;
			vector<int>* nonselectedGroup2 = groups + nonselectedGroupNumberP2;
			for (int j = 0; j < nonselectedGroup1->size(); j++) {
				int alleleNumber = nonselectedGroup1->at(j);
				nonGroupedAllelesP1.push_back(alleleNumber);
			}

			for (int j = 0; j < nonselectedGroup2->size(); j++) {
				int alleleNumber = nonselectedGroup2->at(j);
				nonGroupedAllelesP2.push_back(alleleNumber);
			}
		}
	} else {
		int numberOfNonSelectedGroups = numberOfGroups - numberOfSelectedGroups;
		partial_sort(groupScoresOfP1.begin(),
				groupScoresOfP1.begin() + numberOfNonSelectedGroups,
				groupScoresOfP1.end(), ComparePairIntDoubleMax());
		partial_sort(groupScoresOfP2.begin(),
				groupScoresOfP2.begin() + numberOfNonSelectedGroups,
				groupScoresOfP2.end(), ComparePairIntDoubleMax());

		vector<pair<int, double>>::iterator itP1 = groupScoresOfP1.end() - 1;
		vector<pair<int, double>>::iterator itP2 = groupScoresOfP2.end() - 1;
		for (int i = 0; i < numberOfSelectedGroups; i++) {
			int selectedGroupNumberP1 = (*itP1--).first;
			int selectedGroupNumberP2 = (*itP2--).first;
			vector<int>* selectedGroup1 = groups + selectedGroupNumberP1;
			vector<int>* selectedGroup2 = groups + selectedGroupNumberP2;

			for (int j = 0; j < selectedGroup1->size(); j++) {
				int alleleNumber = selectedGroup1->at(j);
				alleleGroups[0][alleleNumber] = selectedGroupNumberP1;
			}

			for (int j = 0; j < selectedGroup2->size(); j++) {
				int alleleNumber = selectedGroup2->at(j);
				alleleGroups[1][alleleNumber] = selectedGroupNumberP2;
			}
		}

		for(int i = numberOfSelectedGroups; i < numberOfGroups; i++) {
			int nonselectedGroupNumberP1 = (*itP1--).first;
			int nonselectedGroupNumberP2 = (*itP2--).first;
			vector<int>* nonselectedGroup1 = groups + nonselectedGroupNumberP1;
			vector<int>* nonselectedGroup2 = groups + nonselectedGroupNumberP2;
			for (int j = 0; j < nonselectedGroup1->size(); j++) {
				int alleleNumber = nonselectedGroup1->at(j);
				nonGroupedAllelesP1.push_back(alleleNumber);
			}

			for (int j = 0; j < nonselectedGroup2->size(); j++) {
				int alleleNumber = nonselectedGroup2->at(j);
				nonGroupedAllelesP2.push_back(alleleNumber);
			}
		}
	}


	/*double numberOfSelectedAlleles = alleleNumber * (selectedIndividualPortion
			+ addedSelectedIndividualPortion * (double) ga->generation()
					/ (double) ga->nGenerations());

	if(numberOfSelectedAlleles > nonGroupedAllelesP1.size()) {
		numberOfSelectedAlleles = nonGroupedAllelesP1.size();
	}

	if(numberOfSelectedAlleles > nonGroupedAllelesP2.size()) {
		numberOfSelectedAlleles = nonGroupedAllelesP2.size();
	}*/

	int individualNumber = numberOfGroups;
	/*CompareVerticesAccToObjectivePart compareVerticesAccToObjectivePart;
	compareVerticesAccToObjectivePart.parent = &p1;
	partial_sort(nonGroupedAllelesP1.begin(), nonGroupedAllelesP1.begin() + numberOfSelectedAlleles, nonGroupedAllelesP1.end(), compareVerticesAccToObjectivePart);
	compareVerticesAccToObjectivePart.parent = &p2;
	partial_sort(nonGroupedAllelesP2.begin(), nonGroupedAllelesP2.begin() + numberOfSelectedAlleles, nonGroupedAllelesP2.end(), compareVerticesAccToObjectivePart);

	for(int i=0; i<numberOfSelectedAlleles; i++) {
		int alleleNumber = nonGroupedAllelesP1[i];
		alleleGroups[0][alleleNumber] = individualNumber;
		alleleNumber = nonGroupedAllelesP2[i];
		alleleGroups[1][alleleNumber] = individualNumber++;
	}*/
	vector<int> nonGroupedAllelesSelectedP1;
	vector<int> nonGroupedAllelesSelectedP2;
	double scoreLimit = alleleNumber * 0.000;
	for(vector<int>::iterator it = nonGroupedAllelesP1.begin(); it!=nonGroupedAllelesP1.end(); it++) {
		int allele = *it;
		float alleleScore = scoreOfAnAllele(p1, allele);
		if(alleleScore <= scoreLimit) {
			nonGroupedAllelesSelectedP1.push_back(allele);
		}
	}
	for(vector<int>::iterator it = nonGroupedAllelesSelectedP1.begin(); it!=nonGroupedAllelesSelectedP1.end(); it++) {
		int alleleNumber = *it;
		alleleGroups[0][alleleNumber] = individualNumber++;
	}
	for(vector<int>::iterator it = nonGroupedAllelesP2.begin(); it!=nonGroupedAllelesP2.end(); it++) {
		int allele = *it;
		float alleleScore = scoreOfAnAllele(p2, allele);
		if(alleleScore <= scoreLimit) {
			nonGroupedAllelesSelectedP2.push_back(allele);
		}
	}
	for(vector<int>::iterator it = nonGroupedAllelesSelectedP2.begin(); it!=nonGroupedAllelesSelectedP2.end(); it++) {
		int alleleNumber = *it;
		alleleGroups[1][alleleNumber] = individualNumber++;
	}
	//cout<<"GAPointLabelPlacement::getRankedPartialAlleleGroups>end"<<endl;
	lastGeneration = ga->generation();
	return alleleGroups;
}

float GAPointLabelPlacement::objective(GAGenome & g) {
	GA1DArrayAlleleGenome<int> & genome = (GA1DArrayAlleleGenome<int> &) g;
	float objectiveValue;
	float potentialConflictSize = 0;
	float conflictSize = 0;
	int pointCount =
			conflictGraph->getConflictGraphOfPoints()->getVertexNumber();
	int positionNumber = conflictGraph->getPositionNumber();
	for (int i = 0; i < genome.length(); i++) {
		int pointNumber =
				allelePointNumbers == NULL ? i : allelePointNumbers[i];
		int position1 = genome.gene(i);
		int posIx1 = pointNumber * positionNumber + position1;

		vector<int>* conflictingPoints =
				conflictGraph->getConflictGraphOfPoints()->getAdjacencyList()
						+ pointNumber;
		for (vector<int>::iterator it = conflictingPoints->begin();
				it != conflictingPoints->end(); it++) {
			int conflictingPointNo = *it;
			int conflictingPointIx;
			if (allelePointNumbers != NULL) {
				conflictingPointIx = alleleNumberIndexMap->find(
						conflictingPointNo)->second;
			} else {
				conflictingPointIx = conflictingPointNo;
			}
			int position2 = genome.gene(conflictingPointIx);
			int posIx2 = conflictingPointNo * positionNumber + position2;
			if (conflictGraph->getConflictGraphOfPositions()->getAdjacencyMatrix()[posIx1][posIx2]) {
				conflictSize++;
				break;
			}
		}

		potentialConflictSize +=
				conflictGraph->getConflictGraphOfPositions()->getAdjacencyList()[posIx1].size();
	}
	objectiveValue = conflictSize /*+ 2.0*(1-(double)ga->generation()/(double)ga->nGenerations())*potentialConflictSize*/;
	return conflictSize;
}
float GAPointLabelPlacement::scoreOfAnAllele(const GA1DArrayGenome<int> & g,
		int i) {
	int conflictSize = 0;
	int position1 = g.gene(i);
	int positionNumber = conflictGraph->getPositionNumber();
	int posIx1 = i * positionNumber + position1;

	vector<int>* conflictingPoints =
			conflictGraph->getConflictGraphOfPoints()->getAdjacencyList() + i;
	for (vector<int>::iterator it = conflictingPoints->begin();
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
vector<vector<int>*>* GAPointLabelPlacement::clusterPivotal(int depth,
		int maxBranchingFactor) {
	const int pointNumber = alleleNumber;
	vector<vector<int>*>* clusters = new vector<vector<int>*>();

	if (depth > -1) {
		vector<int> pointsToBeEffectivelyMaskedLimited;
		vector<int> pointsToBeEffectivelyMaskedOverLimit;
		for (int i = 0; i < pointNumber; i++) {
			if (maxBranchingFactor == 0) {
				pointsToBeEffectivelyMaskedLimited.push_back(i);
			} else {
				int pointConflictSize =
						conflictGraph->getConflictGraphOfPoints()->getAdjacencyList()[i].size();
				if (pointConflictSize > maxBranchingFactor) {
					pointsToBeEffectivelyMaskedOverLimit.push_back(i);
				} else {
					pointsToBeEffectivelyMaskedLimited.push_back(i);
				}
			}
		}
		vector<int> pointsEffectivelyMasked;
		bool pointsEffectivelyMaskedIndicator[pointNumber];
		for (int i = 0; i < pointNumber; i++) {
			pointsEffectivelyMaskedIndicator[i] = false;
		}
		int lastDepthCurStart = -1;
		int lastDepthCurEnd = -1;
		int groupNumber = 0;
		for (int curEffectiveMaskSize = 0; curEffectiveMaskSize < pointNumber;
				) {
			vector<int>* cluster = new vector<int>();
			clusters->push_back(cluster);
			groupNumber++;
			vector<int>* pointsToBeEffectivelyMasked;
			if (!pointsToBeEffectivelyMaskedLimited.empty()) {
				pointsToBeEffectivelyMasked =
						&pointsToBeEffectivelyMaskedLimited;
			} else {
				pointsToBeEffectivelyMasked =
						&pointsToBeEffectivelyMaskedOverLimit;
			}
			int pointToMaskIx = GARandomInt(0,
					pointsToBeEffectivelyMasked->size() - 1);
			int pointToMask = pointsToBeEffectivelyMasked->at(pointToMaskIx);

			pointsToBeEffectivelyMasked->erase(
					pointsToBeEffectivelyMasked->begin() + pointToMaskIx);
			pointsEffectivelyMasked.push_back(pointToMask);
			pointsEffectivelyMaskedIndicator[pointToMask] = true;
			cluster->push_back(pointToMask);
			curEffectiveMaskSize++;

			lastDepthCurStart = lastDepthCurEnd = pointsEffectivelyMasked.size()
					- 1;
			for (int curDepth = 0;
					curDepth < depth && curEffectiveMaskSize < pointNumber;
					curDepth++) {
				for (;
						lastDepthCurStart <= lastDepthCurEnd
								&& curEffectiveMaskSize < pointNumber;
						lastDepthCurStart++) {
					int pointToFindRivals =
							pointsEffectivelyMasked[lastDepthCurStart];
					vector<int> conflictingPoints =
							conflictGraph->getConflictGraphOfPoints()->getAdjacencyList()[pointToFindRivals];
					if (conflictingPoints.size() > maxBranchingFactor
							&& !pointsToBeEffectivelyMaskedLimited.empty()) {
						continue;
					}
					for (vector<int>::iterator it = conflictingPoints.begin();
							it != conflictingPoints.end()
									&& curEffectiveMaskSize < pointNumber;
							it++) {
						int conflictingPoint = *it;
						if (!pointsEffectivelyMaskedIndicator[conflictingPoint]) {
							for (vector<int>::iterator it2 =
									pointsToBeEffectivelyMasked->begin();
									it2 != pointsToBeEffectivelyMasked->end();
									it2++) {
								if (*it2 == conflictingPoint) {
									pointsToBeEffectivelyMasked->erase(it2);
									break;
								}
							}
							pointsEffectivelyMasked.push_back(conflictingPoint);
							pointsEffectivelyMaskedIndicator[conflictingPoint] =
									true;
							cluster->push_back(conflictingPoint);
							curEffectiveMaskSize++;
						}
					}
				}
				lastDepthCurEnd = pointsEffectivelyMasked.size() - 1;
			}
		}
	}
	return clusters;
}

int GAPointLabelPlacement::uniformCrossoverWithRandomPartialGroupedMasking(
		const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {
	int* alleleGroups = getRandomPartialAlleleGroups();
	int alleleGroupParents[numberOfSelectedGroups];
	for (int i = 0; i < numberOfSelectedGroups; i++) {
		alleleGroupParents[i] = GARandomBit();
	}

	const GA1DArrayGenome<int> &mom = DYN_CAST(const GA1DArrayGenome<int> &,
			p1);
	const GA1DArrayGenome<int> &dad = DYN_CAST(const GA1DArrayGenome<int> &,
			p2);

	GA1DArrayGenome<int> &sis = DYN_CAST(GA1DArrayGenome<int> &, *c1);
	GA1DArrayGenome<int> &bro = DYN_CAST(GA1DArrayGenome<int> &, *c2);

	for (int i = mom.length() - 1; i >= 0; i--) {
		if (alleleGroups[i] != 0) {
			int maskGroupParent = alleleGroupParents[alleleGroups[i]];
			if (maskGroupParent == 1) {
				sis.gene(i, mom.gene(i));
				bro.gene(i, dad.gene(i));
			} else {
				sis.gene(i, dad.gene(i));
				bro.gene(i, mom.gene(i));
			}
		} else {
			if (GARandomBit()) {
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
int GAPointLabelPlacement::uniformCrossoverWithProblemRankedPartialGroupedMasking(
		const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {
	int* alleleGroups = getProblemRankedPartialAlleleGroups();
	int alleleGroupParents[numberOfSelectedGroups];
	for (int i = 0; i < numberOfSelectedGroups; i++) {
		alleleGroupParents[i] = GARandomBit();
	}

	const GA1DArrayGenome<int> &mom = DYN_CAST(const GA1DArrayGenome<int> &,
			p1);
	const GA1DArrayGenome<int> &dad = DYN_CAST(const GA1DArrayGenome<int> &,
			p2);

	GA1DArrayGenome<int> &sis = DYN_CAST(GA1DArrayGenome<int> &, *c1);
	GA1DArrayGenome<int> &bro = DYN_CAST(GA1DArrayGenome<int> &, *c2);

	for (int i = mom.length() - 1; i >= 0; i--) {
		if (alleleGroups[i] != -1) {
			int maskGroupParent = alleleGroupParents[alleleGroups[i]];
			if (maskGroupParent == 1) {
				sis.gene(i, mom.gene(i));
				bro.gene(i, dad.gene(i));
			} else {
				sis.gene(i, dad.gene(i));
				bro.gene(i, mom.gene(i));
			}
		} else {
			if (GARandomBit()) {
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

int GAPointLabelPlacement::uniformCrossoverWithSolutionRankedPartialGroupedMasking(
		const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {

	const GA1DArrayGenome<int> &mom = DYN_CAST(const GA1DArrayGenome<int> &,
			p1);
	const GA1DArrayGenome<int> &dad = DYN_CAST(const GA1DArrayGenome<int> &,
			p2);

	GA1DArrayGenome<int> &sis = DYN_CAST(GA1DArrayGenome<int> &, *c1);
	GA1DArrayGenome<int> &bro = DYN_CAST(GA1DArrayGenome<int> &, *c2);

	int **alleleGroups = getSolutionRankedPartialAlleleGroups(mom, dad);

	for (int i = mom.length() - 1; i >= 0; i--) {
		if (alleleGroups[0][i] != -1 && alleleGroups[1][i] == -1) {
			sis.gene(i, mom.gene(i));
			bro.gene(i, mom.gene(i));
		} else if (alleleGroups[0][i] == -1 && alleleGroups[1][i] != -1) {
			sis.gene(i, dad.gene(i));
			bro.gene(i, dad.gene(i));
		} else if (alleleGroups[0][i] != -1 && alleleGroups[1][i] != -1) {
			sis.gene(i, mom.gene(i));
			bro.gene(i, dad.gene(i));
		} else {
			if (GARandomBit()) {
				sis.gene(i, mom.gene(i));
				bro.gene(i, dad.gene(i));
			} else {
				sis.gene(i, dad.gene(i));
				bro.gene(i, mom.gene(i));
			}
		}
	}
	delete alleleGroups[0];
	delete alleleGroups[1];
}
int GAPointLabelPlacement::uniformCrossoverWithSolutionRankedPartialGroupedAndIndividualMasking(
		const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {

	const GA1DArrayGenome<int> &mom = DYN_CAST(const GA1DArrayGenome<int> &,
			p1);
	const GA1DArrayGenome<int> &dad = DYN_CAST(const GA1DArrayGenome<int> &,
			p2);

	GA1DArrayGenome<int> &sis = DYN_CAST(GA1DArrayGenome<int> &, *c1);
	GA1DArrayGenome<int> &bro = DYN_CAST(GA1DArrayGenome<int> &, *c2);

	int **alleleGroups = getSolutionRankedPartialAlleleGroupsAndIndividuals(mom, dad);

	for (int i = mom.length() - 1; i >= 0; i--) {
		if (alleleGroups[0][i] != -1 && alleleGroups[1][i] == -1) {
			sis.gene(i, mom.gene(i));
			bro.gene(i, mom.gene(i));
		} else if (alleleGroups[0][i] == -1 && alleleGroups[1][i] != -1) {
			sis.gene(i, dad.gene(i));
			bro.gene(i, dad.gene(i));
		} else if (alleleGroups[0][i] != -1 && alleleGroups[1][i] != -1) {
			sis.gene(i, mom.gene(i));
			bro.gene(i, dad.gene(i));
		} else {
			if (GARandomBit()) {
				sis.gene(i, mom.gene(i));
				bro.gene(i, dad.gene(i));
			} else {
				sis.gene(i, dad.gene(i));
				bro.gene(i, mom.gene(i));
			}
		}
	}
	delete alleleGroups[0];
	delete alleleGroups[1];
}
Solution& GAPointLabelPlacement::optimize(ConflictGraph& conflictGraph) {
	GAPointLabelPlacement::conflictGraph = &conflictGraph;
	initialize(conflictGraph, false, -1);
	int pointNumber =
			conflictGraph.getConflictGraphOfPoints()->getVertexNumber();
	int positionsPerPointNumber = conflictGraph.getPositionNumber();

	int* alleleSet = new int[positionsPerPointNumber];
	for (int i = 0; i < positionsPerPointNumber; i++)
		alleleSet[i] = i;
	GAAlleleSet<int> positionAlleleSet(positionsPerPointNumber, alleleSet);
	GA1DArrayAlleleGenome<int> genome(pointNumber, positionAlleleSet, objective,
	NULL);
	//genome.crossover(GA1DArrayAlleleGenome<int>::UniformCrossover);
	genome.crossover(
			GAPointLabelPlacement::uniformCrossoverWithSolutionRankedPartialGroupedAndIndividualMasking);
	//genome.crossover(GA1DArrayAlleleGenome<int>::OnePointCrossover);
	ga = new GASimpleGA(genome);
	ga->selector(GATournamentSelector());
	GAParameterList params;
	GASimpleGA::registerDefaultParameters(params);
	params.set(gaNscoreFilename, "bog.dat");
	params.set(gaNflushFrequency, 50);
	params.set(gaNpMutation, 0.001);
	params.set(gaNpCrossover, 0.90);

	int populationSize = pointNumber / 10;
	if (populationSize % 2 == 1)
		populationSize++;
	ga->populationSize(populationSize);
	ga->minimaxi(GAGeneticAlgorithm::MINIMIZE);
	ga->elitist(GABoolean::gaTrue);
	ga->terminator(GAGeneticAlgorithm::TerminateUponGeneration);
	ga->selectScores(GAStatistics::AllScores);
	ga->scoreFilename(string("output3/scores_"+fileName+".txt").c_str());
	ga->flushFrequency(pointNumber);
	if(lastPopulation != NULL) {
		//cout<<"\tAve: "<<lastPopulation->ave()<<"\tDiv: "<<lastPopulation->div()<<"\tDev: "<<lastPopulation->dev();
		//REPLACE WORST INDIVIDUALS
		for(int i=0; i<0.00*pointNumber; i++) {
			int randomIndividual = GARandomInt(0, pointNumber-1);
			GA1DArrayAlleleGenome<int>* oldGenome = (GA1DArrayAlleleGenome<int>*) &lastPopulation->worst(i);
			GA1DArrayAlleleGenome<int>* newGenome = new GA1DArrayAlleleGenome<int>(pointNumber, positionAlleleSet, objective,
					NULL);
			newGenome->crossover(
					GAPointLabelPlacement::uniformCrossoverWithSolutionRankedPartialGroupedAndIndividualMasking);
			for(int j=0; j<newGenome->size(); j++) {
				newGenome->gene(j, newGenome->alleleset(j).allele());
			}
			lastPopulation->replace(newGenome, oldGenome);
			//delete oldGenome;
		}
		for(int i=0; i<lastPopulation->size(); i++) {
			lastPopulation->individual(i).initializer(NoInitializer);
		}
		lastPopulation->evaluate(GABoolean(true));
		ga->population(*lastPopulation);
		ga->nGenerations(0.033*pointNumber);
		delete lastPopulation;
	}
	else
	{
		ga->nGenerations(0.33*pointNumber);
		genome.initializer(GA1DArrayAlleleGenome<int>::UniformInitializer);
	}
	ga->evolve();
	lastPopulation = new GAPopulation(ga->population());
	const GA1DArrayAlleleGenome<unsigned int> & bestIndividual =
			(GA1DArrayAlleleGenome<unsigned int> &) ga->statistics().bestIndividual();
	ga->statistics().write(string("output3/stats_"+fileName+".txt").c_str());
	//ga->flushScores();
	Solution* solution = new Solution(this->conflictGraph);
	Solution& solutionRef = *solution;
	int* labelPlacements = new int[pointNumber];
	solution->setLabelPlacements(labelPlacements);
	for (int i = 0; i < bestIndividual.length(); i++) {
		int pointNo;
		if (allelePointNumbers != NULL) {
			pointNo = allelePointNumbers[i];
		} else {
			pointNo = i;
		}
		labelPlacements[pointNo] = bestIndividual.gene(i);
	}
	lastGeneration = -1;
	delete alleleSet;
	delete ga;
	if (allelePointNumbers != NULL) {
		delete allelePointNumbers;
		delete alleleNumberIndexMap;
	}
	return solutionRef;
}

} /* namespace labelplacement */
