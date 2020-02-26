/*
 * ConflictGraph.cpp
 *
 *  Created on: Dec 17, 2017
 *      Author: oakile
 */

#include "ConflictGraph.h"

#include <fstream>
#include <sstream>

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

int ConflictGraph::compareConflictGraphs(int size, std::vector<int>* cg1, std::vector<int>* cg2) {
	int diff=0;
	int conflictSize1=0;
	int conflictSize2=0;
	for(int i=0; i<size; i++) {
		std::vector<int> adjList1 = std::vector<int>(*(cg1 + i));
		std::vector<int> adjList2 = std::vector<int>(*(cg2 + i));
		conflictSize1+=adjList1.size();
		conflictSize2+=adjList2.size();
		for(int j=0; j<adjList1.size();) {
			bool isFound = false;
			int k=0;
			for(; k < adjList2.size(); k++) {
				if(adjList1[j] == adjList2[k]) {
					isFound = true;
					break;
				}
			}
			if(isFound) {
				adjList1.erase(adjList1.begin()+j);
				adjList2.erase(adjList2.begin()+k);
			}
			else
			{
				j++;
			}
		}
		diff+= adjList2.size();
	}
	diff/=2;
	conflictSize1/=2;
	conflictSize2/=2;
	int diffConflictSize = (int)abs(conflictSize1-conflictSize2);
	//std::cout<<conflictSize1<<" "<<conflictSize2<<" "<<diffConflictSize<<endl;
	return diff;
}

void ConflictGraph::writeConflictGraphIntoFile(std::vector<int>* adjacencyListOfPositions, int numPoints, int positionsPerPoint, std::string filename) {
	int numPositions = numPoints * positionsPerPoint;
	std::ofstream ofs(filename, std::ios_base::out);
	ofs<<std::endl;
	ofs<<numPoints<<std::endl;
	ofs<<positionsPerPoint<<std::endl;
	for(int i=0; i<numPositions; i++) {
		std::vector<int>* adjacencyListOfPosition = adjacencyListOfPositions + i;
		ofs<<adjacencyListOfPosition->size()<<std::endl;
		for(std::vector<int>::iterator it=adjacencyListOfPosition->begin(); it!=adjacencyListOfPosition->end(); it++) {
			ofs<<*it+1;
			if(it+1!=adjacencyListOfPosition->end()) {
				ofs<<" ";
			}
		}
		ofs<<std::endl;
	}
	ofs.close();
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

void ConflictGraph::initializeOf(int positionNumberPerPoint, int numberOfPositions, std::vector<int>* adjacencyListOfPositions) {
	this->positionNumberPerPoint = positionNumberPerPoint;

	bool** adjacencyMatrixOfPositions = new bool*[numberOfPositions];
	for(int i=0; i<numberOfPositions; i++) {
		adjacencyMatrixOfPositions[i] = new bool[numberOfPositions];
		for(int j=0; j<numberOfPositions; j++) {
			adjacencyMatrixOfPositions[i][j] = false;
		}
	}

	conflictGraphOfPositions = new Graph();
	conflictGraphOfPositions->setVertexNumber(numberOfPositions);
	conflictGraphOfPositions->setAdjacencyList(adjacencyListOfPositions);
	conflictGraphOfPositions->setAdjacencyMatrix(adjacencyMatrixOfPositions);

	int numberOfPoints = numberOfPositions / positionNumberPerPoint;

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
