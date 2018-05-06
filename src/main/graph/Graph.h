/*
 * Graph.h
 *
 *  Created on: Feb 1, 2018
 *      Author: oakile
 */

#ifndef _GRAPH_H_
#define _GRAPH_H_

#include <vector>

namespace graph {

class Graph {
public:
	Graph();
	virtual ~Graph();

	int getVertexNumber() const;

	void setVertexNumber(int vertexNumber);

	int getEdgeNumber() const;

	void setEdgeNumber(int edgeNumber);

	std::vector<int>* getVertices() const;

	void setVertices(std::vector<int>* vertices);

	bool** getAdjacencyMatrix() const;

	void setAdjacencyMatrix(bool** adjacencyMatrix);

	std::vector<int>* getAdjacencyList() const;

	void setAdjacencyList(std::vector<int>* adjacencyList);

	void generateAdjacencyMatrixFromAdjacencyList();

private:
	int vertexNumber;
	int edgeNumber;
	std::vector<int>* vertices;
	bool** adjacencyMatrix;
	std::vector<int>* adjacencyList;
};

} /* namespace labelplacement */

#endif /* _GRAPH_H_ */
