/*
 * Solution.h
 *
 *  Created on: Dec 18, 2017
 *      Author: oakile
 */

#ifndef SRC_SOLUTION_H_
#define SRC_SOLUTION_H_

#include "ConflictGraph.h"

namespace labelplacement {

class Solution {
public:
	Solution(ConflictGraph* conflictGraph);
	virtual ~Solution();
	ConflictGraph* getConflictGraph();
	void setConflictGraph(ConflictGraph* conflictGraph);
	int* getLabelPlacements();
	void setLabelPlacements(int* labelPlacements);
	int getConflictSize();
	int getConflictSize(ConflictGraph* conflictGraph);
private:
	ConflictGraph* conflictGraph;
	int* labelPlacements;
};

} /* namespace labelplacement */

#endif /* SRC_SOLUTION_H_ */
