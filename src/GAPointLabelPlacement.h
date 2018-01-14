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
	static int* getGroupedMask();
	static int getNumberOfMaskGroupes();
	static float objective(GAGenome & g);
	static void generateMask(int depth, double coverage);
	static void generateGroupedMask(int depth, double coverage);
	static int uniformCrossoverWithMasking(const GAGenome&, const GAGenome&,
				      GAGenome*, GAGenome*);
	static int uniformCrossoverWithGroupedMasking(const GAGenome&, const GAGenome&,
				      GAGenome*, GAGenome*);
	Solution& optimize(ConflictGraph& conflictGraph);

private:
	static ConflictGraph* conflictGraph;
	static bool* mask;
	static int numberOfMaskGroupes;
	static int* groupedMask;
	static int lastGeneration;
	static GASimpleGA* ga;
};

} /* namespace labelplacement */

#endif /* _POINTLABELINGGA_H_ */
