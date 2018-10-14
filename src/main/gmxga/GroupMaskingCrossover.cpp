/*
 * GroupMaskingCrossover.cpp
 *
 *  Created on: Aug 6, 2018
 *      Author: oakile
 */

#include "GroupMaskingCrossover.h"

#include <algorithm>

using namespace std;

GroupMaskingCrossover::GroupMaskingCrossover() {
	groups = NULL;
}

GroupMaskingCrossover::~GroupMaskingCrossover() {
	for(vector<vector<int>*>::iterator it = groups->begin(); it!=groups->end(); it++) {
		delete *it;
	}
	delete groups;
}

void GroupMaskingCrossover::initialize(graph::Graph* graph) {
}

void GroupMaskingCrossover::setGrouper(std::function<AlleleGrouper> grouper) {
	this->grouper = grouper;
}

GMXDelegate::GMXDelegate(GroupMaskingCrossover* gmx) {
	this->gmx=gmx;
}

int GMXDelegate::operator()(const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {
	return (*gmx)(p1, p2, c1, c2);
}

UniformRankedGroupMaskingCrossover::UniformRankedGroupMaskingCrossover() {
	groupProportion = 0;
	groupProportionMargin = 0;
	individualProportion = 0;
	individualProportionMargin = 0;
	individualScoreLimit=0;
	groups = NULL;
	alleleNumber = 0;
	numberOfGroups = 0;
	curNumberOfGroups = 0;
	ga = NULL;
	evaluator = NULL;
}

UniformRankedGroupMaskingCrossover::~UniformRankedGroupMaskingCrossover() {
	delete evaluator;
}

void UniformRankedGroupMaskingCrossover::setParams(double groupProportion, double groupProportionMargin,
		double individualProportion, double individualProportionMargin, double individualScoreLimit) {
	this->groupProportion = groupProportion;
	this->groupProportionMargin = groupProportionMargin;
	this->individualProportion = individualProportion;
	this->individualProportionMargin = individualProportionMargin;
	this->individualScoreLimit=individualScoreLimit;
}

std::function<UniformRankedGroupMaskingCrossover::AlleleEvaluator>* UniformRankedGroupMaskingCrossover::getEvaluator() {
	return this->evaluator;
}

void UniformRankedGroupMaskingCrossover::setEvaluator(std::function<AlleleEvaluator>* evaluator) {
	this->evaluator = evaluator;
}

void UniformRankedGroupMaskingCrossover::setGA(GASimpleGA* ga) {
	this->ga = ga;
}

void UniformRankedGroupMaskingCrossover::initialize(graph::Graph* graph) {
	if(groups!=NULL){
		for(vector<vector<int>*>::iterator it = groups->begin(); it!=groups->end(); it++) {
			delete *it;
		}
		delete groups;
	}
	if(groupProportion!=0 && groupProportionMargin!=0) {
		groups = grouper(graph);
		numberOfGroups = groups->size();
	}
	GroupMaskingCrossover::initialize(graph);
	alleleNumber = graph->getVertices()->size();
}

int UniformRankedGroupMaskingCrossover::operator()(const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {
	curNumberOfGroups = (groupProportion + (double)ga->generation() / (double) ga->nGenerations()
			* groupProportionMargin) * numberOfGroups;

	const GA1DArrayGenome<int> &mom = DYN_CAST(const GA1DArrayGenome<int> &,
				p1);
	const GA1DArrayGenome<int> &dad = DYN_CAST(const GA1DArrayGenome<int> &,
			p2);

	GA1DArrayGenome<int> &sis = DYN_CAST(GA1DArrayGenome<int> &, *c1);
	GA1DArrayGenome<int> &bro = DYN_CAST(GA1DArrayGenome<int> &, *c2);

	int* alleleGroups[2];
	alleleGroups[0] = selectGroupsOfIndividual(mom);
	alleleGroups[1] = selectGroupsOfIndividual(dad);

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

	return -1;
}

int* UniformRankedGroupMaskingCrossover::selectGroupsOfIndividual(const GA1DArrayGenome<int>& individual) {
	int* alleleGroups = new int[alleleNumber];

	for (int i = 0; i < alleleNumber; i++) {
		alleleGroups[i] = -1;
	}

	vector<int> nonGroupedAlleles;

	if(groups != NULL) {
		vector<pair<int, double>> groupScores;
		for (int i = 0; i < numberOfGroups; i++) {
			vector<int>* group = groups->at(i);
			double groupScore = 0;
			for (int j = 0; j < group->size(); j++) {
				int groupVertex = group->at(j);
				groupScore += (*evaluator)(individual, groupVertex);
			}
			groupScore /= (double) group->size();

			groupScores.push_back(pair<int, double>(i, groupScore));
		}

		double selectedGroupPortion = (double) curNumberOfGroups
				/ (double) numberOfGroups;

		if (selectedGroupPortion < 0.5) {
			partial_sort(groupScores.begin(),
					groupScores.begin() + curNumberOfGroups,
					groupScores.end(), ComparePairIntDoubleMin());
			vector<pair<int, double>>::iterator itGroupScore = groupScores.begin();
			for (int i = 0; i < curNumberOfGroups; i++) {
				int selectedGroupNumber = (*itGroupScore++).first;
				vector<int>* selectedGroup = groups->at(selectedGroupNumber);

				for (int j = 0; j < selectedGroup->size(); j++) {
					int alleleNumber = selectedGroup->at(j);
						alleleGroups[alleleNumber] = selectedGroupNumber;
				}
			}
			for(int i = curNumberOfGroups; i < numberOfGroups; i++) {
				int nonselectedGroupNumber = (*itGroupScore++).first;
				vector<int>* nonselectedGroup = groups->at(nonselectedGroupNumber);
				for (int j = 0; j < nonselectedGroup->size(); j++) {
					int alleleNumber = nonselectedGroup->at(j);
					nonGroupedAlleles.push_back(alleleNumber);
				}
			}
		} else {
			int numberOfNonSelectedGroups = numberOfGroups - curNumberOfGroups;
			partial_sort(groupScores.begin(),
					groupScores.begin() + numberOfNonSelectedGroups,
					groupScores.end(), ComparePairIntDoubleMax());

			vector<pair<int, double>>::iterator itGroupScore = groupScores.end() - 1;
			for (int i = 0; i < curNumberOfGroups; i++) {
				int selectedGroupNumber = (*itGroupScore--).first;
				vector<int>* selectedGroup = groups->at(selectedGroupNumber);

				for (int j = 0; j < selectedGroup->size(); j++) {
					int alleleNumber = selectedGroup->at(j);
					alleleGroups[alleleNumber] = selectedGroupNumber;
				}
			}

			for(int i = curNumberOfGroups; i < numberOfGroups; i++) {
				int nonselectedGroupNumber = (*itGroupScore--).first;
				vector<int>* nonselectedGroup = groups->at(nonselectedGroupNumber);
				for (int j = 0; j < nonselectedGroup->size(); j++) {
					int alleleNumber = nonselectedGroup->at(j);
					nonGroupedAlleles.push_back(alleleNumber);
				}
			}
		}
	} else {
		for(int i=0;i<alleleNumber;i++) {
			nonGroupedAlleles.push_back(i);
		}
	}
	//Individual Selection Step
	int individualNumber = curNumberOfGroups;
	vector<int> nonGroupedAllelesSelected;
	for(vector<int>::iterator it = nonGroupedAlleles.begin(); it!=nonGroupedAlleles.end(); it++) {
		int allele = *it;
		float alleleScore = (*evaluator)(individual, allele);
		if(alleleScore <= individualScoreLimit) {
			nonGroupedAllelesSelected.push_back(allele);
		}
	}

	for(vector<int>::iterator it = nonGroupedAllelesSelected.begin(); it!=nonGroupedAllelesSelected.end(); it++) {
		int alleleNumber = *it;
		alleleGroups[alleleNumber] = individualNumber++;
	}
	return alleleGroups;
}
