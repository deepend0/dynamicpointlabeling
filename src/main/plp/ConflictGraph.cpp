/*
 * ConflictGraph.cpp
 *
 *  Created on: Dec 17, 2017
 *      Author: oakile
 */

#include "ConflictGraph.h"

namespace labelplacement {

ConflictGraph::ConflictGraph() {
	this->positionNumber = 0;
	this->conflictGraphOfPoints = __null;
	this->conflictGraphOfPositions = __null;
}

ConflictGraph::~ConflictGraph() {
	delete conflictGraphOfPoints;
	delete conflictGraphOfPositions;
}

Graph* ConflictGraph::getConflictGraphOfPoints() const {
	return conflictGraphOfPoints;
}

void ConflictGraph::setConflictGraphOfPoints(Graph* conflictGraphOfPoints) {
	this->conflictGraphOfPoints = conflictGraphOfPoints;
}

Graph* ConflictGraph::getConflictGraphOfPositions() const {
	return conflictGraphOfPositions;
}

void ConflictGraph::setConflictGraphOfPositions(
		Graph* conflictGraphOfPositions) {
	this->conflictGraphOfPositions = conflictGraphOfPositions;
}

int ConflictGraph::getPositionNumber() const {
	return this->positionNumber;
}
void ConflictGraph::setPositionNumber(int positionNumber) {
	this->positionNumber = positionNumber;
}
} /* namespace labelplacement */
