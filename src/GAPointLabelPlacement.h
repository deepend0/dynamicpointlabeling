/*
 * GAPointLabelPlacement.h
 *
 *  Created on: Dec 18, 2017
 *      Author: oakile
 */

#ifndef _POINTLABELINGGA_H_
#define _POINTLABELINGGA_H_

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
	static void initialize(ConflictGraph& conflictGraph);
	static int* getAlleleGroups();
	static float objective(GAGenome & g);
	static void generateMask(int depth, double coverage);
	static void generateGroupedMask(int depth, double coverage, int branchingLimit);
	static int uniformCrossoverWithMasking(const GAGenome&, const GAGenome&,
				      GAGenome*, GAGenome*);
	static int uniformCrossoverWithRandomPartialGroupedMasking(const GAGenome&, const GAGenome&,
				      GAGenome*, GAGenome*);
	Solution& optimize(ConflictGraph& conflictGraph);

private:
	static ConflictGraph* conflictGraph;
	static bool* mask;
	static double initialSelectedGroupPortion;
	static double addedSelectedGroupPortion;
	static int numberOfGroups;
	static int numberOfSelectedGroups;
	static std::vector<int>* groups;
	static std::vector<int>** selectedGroups;
	static int* selectedGroupNumbers;
	static int* alleleGroups;
	static int alleleNumber;
	static int lastGeneration;
	static GASimpleGA* ga;
};

} /* namespace labelplacement */

#endif /* _POINTLABELINGGA_H_ */
