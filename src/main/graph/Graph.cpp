/*
 * Graph.cpp
 *
 *  Created on: Feb 1, 2018
 *      Author: oakile
 */

#include "Graph.h"
#include <iostream>
namespace graph {

Graph::Graph() {
	vertexNumber = 0;
	edgeNumber = 0;
	vertices = __null;
	adjacencyList = __null;
	adjacencyMatrix = __null;
}

Graph::~Graph() {
	delete vertices;
	delete[] adjacencyList;
	for (int i = 0; i < vertexNumber; i++) {
		delete[] adjacencyMatrix[i];
	}
	delete[] adjacencyMatrix;
}

int Graph::getVertexNumber() const {
	return this->vertexNumber;
}

void Graph::setVertexNumber(int vertexNumber) {
	this->vertexNumber = vertexNumber;
}

int Graph::getEdgeNumber() const {
	return this->edgeNumber;
}

void Graph::setEdgeNumber(int edgeNumber) {
	this->edgeNumber = edgeNumber;
}

std::vector<int>* Graph::getVertices() const {
	return vertices;
}

void Graph::setVertices(std::vector<int>* vertices) {
	this->vertices = vertices;
}

bool** Graph::getAdjacencyMatrix() const {
	return adjacencyMatrix;
}

void Graph::setAdjacencyMatrix(bool** adjacencyMatrix) {
	this->adjacencyMatrix = adjacencyMatrix;
}

std::vector<int>* Graph::getAdjacencyList() const {
	return this->adjacencyList;
}

void Graph::setAdjacencyList(std::vector<int>* adjacencyList) {
	this->adjacencyList = adjacencyList;
}

void Graph::generateAdjacencyMatrixFromAdjacencyList() {
	this->adjacencyMatrix = new bool*[vertexNumber];
	for (int i = 0; i < vertexNumber; i++) {
		adjacencyMatrix[i] = new bool[vertexNumber];
	}
	for (int i = 0; i < vertexNumber; i++) {
		for (int j = 0; j < vertexNumber; j++) {
			adjacencyMatrix[i][j] = false;
		}
	}
	for (int i = 0; i < vertexNumber; i++) {
		for (std::vector<int>::iterator it = adjacencyList[i].begin();
				it != adjacencyList[i].end(); it++) {
			adjacencyMatrix[i][*it] = true;
			edgeNumber++;
		}
	}
	edgeNumber /= 2;
}
} /* namespace labelplacement */

