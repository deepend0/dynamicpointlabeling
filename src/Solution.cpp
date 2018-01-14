/*
 * Solution.cpp
 *
 *  Created on: Dec 18, 2017
 *      Author: oakile
 */

#include "Solution.h"

namespace labelplacement {

Solution::Solution(ConflictGraph* conflictGraph) {
	this->conflictGraph = conflictGraph;
	this->labelPlacements = __null;
}

Solution::~Solution() {
	delete [] labelPlacements;
}

int* Solution::getLabelPlacements() {
	return this->labelPlacements;
}
void Solution::setLabelPlacements(int* labelPlacements) {
	this->labelPlacements = labelPlacements;
}

int Solution::getConflictSize()
{
	int conflictSize = 0;
	for(int i=0; i<conflictGraph->getPointNumber(); i++)
	{
		int position1 = labelPlacements[i];
		int posIx1 = i*conflictGraph->getPositionNumber()+position1;
		for(int j=0; j<conflictGraph->getPointNumber(); j++)
		{
			int position2 = labelPlacements[j];
			int posIx2 = j*conflictGraph->getPositionNumber()+position2;
			if(conflictGraph->getAdjacencyMatrix()[posIx1][posIx2])
			{
				conflictSize++;
				break;
			}
		}
	}
	return conflictSize;
}
} /* namespace labelplacement */
