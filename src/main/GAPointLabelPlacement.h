/*
 * GAPointLabelPlacement.h
 *
 *  Created on: Dec 18, 2017
 *      Author: oakile
 */

#ifndef _POINTLABELINGGA_H_
#define _POINTLABELINGGA_H_

#include <string>
#include <unordered_map>

#include <ga/GA1DArrayGenome.h>
#include <ga/GASimpleGA.h>

#include "ConflictGraph.h"
#include "Solution.h"

namespace labelplacement {

class GAPointLabelPlacement {
public:
	GAPointLabelPlacement();
	virtual ~GAPointLabelPlacement();
	static ConflictGraph* getConflictGraph();
	static void setConflictGraph(ConflictGraph* conflictGraph);
	static bool* getMask();
	static std::vector<int>* getGroups();
	static int getNumberOfGroups();
	static int getNumberOfSelectedGroups();
	static void initialize(ConflictGraph& conflictGraph, bool sequenceGroups,
			int clusteringType);
	static int* getRandomPartialAlleleGroups();
	static int** getSolutionRankedPartialAlleleGroups(const GA1DArrayGenome<int>& p1,
			const GA1DArrayGenome<int>& p2);
	static int** getSolutionRankedPartialAlleleGroupsAndIndividuals(const GA1DArrayGenome<int>& p1,
			const GA1DArrayGenome<int>& p2);
	static int* getProblemRankedPartialAlleleGroups();
	static float objective(GAGenome & g);
	static float scoreOfAnAllele(const GA1DArrayGenome<int> & g, int i);
	static std::vector<std::vector<int>*>* clusterPivotal(int depth,
			int maxBranchingFactor);
	static int uniformCrossoverWithMasking(const GAGenome&, const GAGenome&,
			GAGenome*, GAGenome*);
	static int uniformCrossoverWithRandomPartialGroupedMasking(const GAGenome&,
			const GAGenome&, GAGenome*, GAGenome*);
	static int uniformCrossoverWithProblemRankedPartialGroupedMasking(const GAGenome&,
			const GAGenome&, GAGenome*, GAGenome*);
	static int uniformCrossoverWithSolutionRankedPartialGroupedMasking(const GAGenome&,
			const GAGenome&, GAGenome*, GAGenome*);
	static int uniformCrossoverWithSolutionRankedPartialGroupedAndIndividualMasking(const GAGenome&,
			const GAGenome&, GAGenome*, GAGenome*);
	Solution& optimize(ConflictGraph& conflictGraph);
	std::string fileName;

private:
	static ConflictGraph* conflictGraph;
	static bool* mask;
	static double initialSelectedGroupPortion;
	static double addedSelectedGroupPortion;
	static double selectedIndividualPortion;
	static double addedSelectedIndividualPortion;
	static int numberOfGroups;
	static int numberOfSelectedGroups;
	static std::vector<int>* groups;
	static int* alleleGroups;
	static int* allelePointNumbers;
	static std::unordered_map<int, int>* alleleNumberIndexMap;
	static int alleleNumber;
	static int lastGeneration;
	static GASimpleGA* ga;
};

} /* namespace labelplacement */

#endif /* _POINTLABELINGGA_H_ */
