/*
 * ConflictGraph.cpp
 *
 *  Created on: Dec 17, 2017
 *      Author: oakile
 */

#include "ConflictGraph.h"

namespace labelplacement {

ConflictGraph::ConflictGraph() {
	this->positionNumberPerPoint = 0;
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

int ConflictGraph::getPositionNumberPerPoint() const {
	return this->positionNumberPerPoint;
}
void ConflictGraph::setPositionNumberPerPoint(int positionNumberPerPoint) {
	this->positionNumberPerPoint = positionNumberPerPoint;
}

void ConflictGraph::initializeOf(int positionNumberPerPoint, int sizeOfAdjacencyListOfPositions, std::vector<int>* adjacencyListOfPositions) {
	this->positionNumberPerPoint = positionNumberPerPoint;

	int numberOfPositions = sizeOfAdjacencyListOfPositions;
	int numberOfPoints = numberOfPositions / positionNumberPerPoint;

	bool** adjacencyMatrixOfPositions = new bool*[numberOfPositions];
	for(int i=0; i<numberOfPositions; i++) {
		adjacencyMatrixOfPositions[i] = new bool[numberOfPositions];
		for(int j=0; j< numberOfPositions; j++) {
			adjacencyMatrixOfPositions[i][j] = false;
		}
	}

	conflictGraphOfPositions = new Graph();
	conflictGraphOfPositions->setVertexNumber(numberOfPositions);
	conflictGraphOfPositions->setAdjacencyList(adjacencyListOfPositions);
	conflictGraphOfPositions->setAdjacencyMatrix(adjacencyMatrixOfPositions);

	std::vector<int>* adjacencyListOfPoints = new std::vector<int>[numberOfPoints];
	bool** adjacencyMatrixOfPoints = new bool*[numberOfPoints];
	for(int i=0; i<numberOfPoints; i++) {
		adjacencyMatrixOfPoints[i] = new bool[numberOfPoints];
		for(int j=0; j< numberOfPoints; j++) {
			adjacencyMatrixOfPoints[i][j] = false;
		}
	}

	conflictGraphOfPoints = new Graph();
	conflictGraphOfPoints->setVertexNumber(numberOfPoints);
	conflictGraphOfPoints->setAdjacencyList(adjacencyListOfPoints);
	conflictGraphOfPoints->setAdjacencyMatrix(adjacencyMatrixOfPoints);

	int edgeNumberOfPositions = 0;
	int edgeNumberOfPoints = 0;
	for(int i=0; i<numberOfPositions; i++) {
		int curPointIx = i / positionNumberPerPoint;
		std::vector<int>* adjacencyListOfCurrentPosition = adjacencyListOfPositions+i;
		for(std::vector<int>::iterator it=adjacencyListOfCurrentPosition->begin(); it!=adjacencyListOfCurrentPosition->end(); it++){
			int adjacentPositionIx = *it;
			adjacencyMatrixOfPositions[i][adjacentPositionIx]=true;
			edgeNumberOfPositions++;
			int adjacentPointIx = adjacentPositionIx / positionNumberPerPoint;
			if(!adjacencyMatrixOfPoints[curPointIx][adjacentPointIx]) {
				adjacencyMatrixOfPoints[curPointIx][adjacentPointIx] = true;
				adjacencyListOfPoints[curPointIx].push_back(adjacentPointIx);
				edgeNumberOfPoints++;
			}
		}
	}
	edgeNumberOfPositions /= 2;
	edgeNumberOfPoints /= 2;
	conflictGraphOfPositions->setEdgeNumber(edgeNumberOfPositions);
	conflictGraphOfPoints->setEdgeNumber(edgeNumberOfPoints);

	std::vector<int>* verticesOfConflictGraphOfPositions = new std::vector<int>();
	for(int i=0;i<numberOfPositions;i++) {
		verticesOfConflictGraphOfPositions->push_back(i);
	}
	std::vector<int>* verticesOfConflictGraphOfPoints = new std::vector<int>();
	for(int i=0;i<numberOfPoints;i++) {
		verticesOfConflictGraphOfPoints->push_back(i);
	}
	conflictGraphOfPoints->setVertices(verticesOfConflictGraphOfPositions);
	conflictGraphOfPoints->setVertices(verticesOfConflictGraphOfPoints);
}
} /* namespace labelplacement */
