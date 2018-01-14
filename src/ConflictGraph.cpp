/*
 * ConflictGraph.cpp
 *
 *  Created on: Dec 17, 2017
 *      Author: oakile
 */

#include "ConflictGraph.h"

namespace labelplacement {

ConflictGraph::ConflictGraph() {

	this->pointNumber=0;
	this->positionNumber=0;
	this->adjacencyMatrix=__null;
	this->conflictingPositions=__null;
	this->conflictingPoints=__null;
}

ConflictGraph::~ConflictGraph() {
	for(int i=0;i<pointNumber*positionNumber;i++)
		delete [] adjacencyMatrix[i];
	delete [] adjacencyMatrix;

	delete [] conflictingPositions;
	delete [] conflictingPoints;
}

bool** ConflictGraph::getAdjacencyMatrix() const {
	return adjacencyMatrix;
}

void ConflictGraph::setAdjacencyMatrix(bool** adjacencyMatrix) {
	this->adjacencyMatrix = adjacencyMatrix;
}

int ConflictGraph::getPointNumber() const {
	return pointNumber;
}

void ConflictGraph::setPointNumber(int pointNumber) {
	this->pointNumber = pointNumber;
}

int ConflictGraph::getPositionNumber() const {
	return positionNumber;
}

void ConflictGraph::setPositionNumber(int positionNumber) {
	this->positionNumber = positionNumber;
}

std::vector<int>* ConflictGraph::getConflictingPositions() const {
	return this->conflictingPositions;
}

void ConflictGraph::setConflictingPositions(std::vector<int>* conflictingPositions) {
	this->conflictingPositions = conflictingPositions;
}

std::vector<int>* ConflictGraph::getConflictingPoints() const {
	return this->conflictingPoints;
}

void ConflictGraph::setConflictingPoints(std::vector<int>* conflictingPoints) {
	this->conflictingPoints = conflictingPoints;
}

} /* namespace labelplacement */
